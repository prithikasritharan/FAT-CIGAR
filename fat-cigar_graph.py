import os
import sys
import argparse
import re
import time
import json
from collections import namedtuple
from itertools import groupby
try:
    import pysam
except:
    import subprocess
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pysam'])
    import pysam




def build_table(f, size):
    '''
    Computes the read name as a string hash for each line in the JSON file and adds it to a hash table and stores the line position in the file. If a collision occurs, all file positions for that key are stored in a list through chaining using the Node tuple. Returns the indexed hash table.
    '''
    start = time.time()
    #initialise the hash table
    table = [ None ] * size
    while True:
        #stores the position of f in the file 
        pos = f.tell()
        try:
            data = json.loads(f.readline())
        except ValueError:
            pass
            break
        readname = re.sub(r'^(\S+)\s+\S+$', r'\1', data[u'name'])
        if not readname: 
            break
        #computes hash by dividing the line by the size of the table 
        i = hash(readname) % size
        #stores hash key at index if empty, if not chains the hash key
        if table[i] is None:
            table[i] = pos
        else:
            table[i] = Node(pos, table[i])
    end = time.time()
    print("Built index: ",end - start," secs")
    return table



def create_cigar(line):
    '''
    Reads in the JSON data for a read (found through searching the hash table for the file position) and converts the individual path edits into the corresponding CIGAR string equivalent. The disjoined list of CIGAR strings are passed to the join_cigarstring function to join the list into a string. Returns the CIGAR2 string.
    ''' 
    data = json.loads(line)

    #Intialise list for appending the CIGARs based on the path mapping edits
    cigar2 = []

    max_node = len(data[u'path'][u'mapping'])
    
    #cycles through the path mapping edit information for each read and uses the to_length and from_length to recreate the CIGAR2 string
    for i in range(0, max_node):
        #gets the number of edits against each node in the path
        edits = len(data[u'path'][u'mapping'][i][u'edit'])

        for j in range(0, edits):
            if (("from_length" in data[u'path'][u'mapping'][i][u'edit'][j]) & ("to_length" in data[u'path'][u'mapping'][i][u'edit'][j])):
                #checks for matches by ensuring a sequence isn't specified  or whether the sequence is a N in mapping edit 
                if (("sequence" not in data[u'path'][u'mapping'][i][u'edit'][j]) or (data[u'path'][u'mapping'][i][u'edit'][j][u'sequence'] == 'N')):
                    match = str(data[u'path'][u'mapping'][i][u'edit'][j][u'from_length']) + '='
                    cigar2.append(match)
                else:
                    mismatch = str(data[u'path'][u'mapping'][i][u'edit'][j][u'from_length']) + 'X'
                    cigar2.append(mismatch)


            #evaluates soft-clips and insertions by checking if only to to_length is present in the mapping edit
            elif ("from_length" not in data[u'path'][u'mapping'][i][u'edit'][j]):
                #if insertion occurs at beginning or end of read sequence, the sequence is soft-clipped
                if (i == 0 and j == 0) or (i == max_node-1 and j == edits-1):
                    soft_clip = str(data[u'path'][u'mapping'][i][u'edit'][j][u'to_length']) + 'S'
                    cigar2.append(soft_clip)
                else:
                    insertion = str(data[u'path'][u'mapping'][i][u'edit'][j][u'to_length']) + 'I'
                    cigar2.append(insertion)  

            #evaluates deletions based on the absence of to_length in the mapping edit
            elif ("to_length" not in data[u'path'][u'mapping'][i][u'edit'][j]):
                deletion = str(data[u'path'][u'mapping'][i][u'edit'][j][u'from_length']) + 'D'
                cigar2.append(deletion) 


    #calls the create_cigar function to join the values in the cigar2 list to create the final CIGAR2 string 
    final_cigar = join_cigarstring(list(reversed(cigar2))) if 'is_reverse' in data[u'refpos'][0] else join_cigarstring(cigar2) 
    return final_cigar 



def join_cigarstring(cig_list):
    '''
    Joins the matches, mismatches, indels and soft clips from a list of disjoined CIGAR strings. Returns the final CIGAR2 string. 
    '''
    #initialize new list for the joined CIGAR strings
    new_cig2 = []

    #count for the total number of matched sequences
    match_sum = 0

    #Cycles through the list containing the individual CIGAR strings and joins them together in the new_cig2 list
    for index, i in enumerate(cig_list):
        num = re.findall(r'\d+', i)
        match_type = ''.join(re.findall(r'\D', i))
        match_sum += int(''.join(num))
        #if the index is smaller than the list length and the next value is also a match, moves onto the next iteration
        if (index < (len(cig_list)-1) and re.search(match_type, cig_list[index+1])):
            continue
        #if the index is the last position in the list or the next value in the list is not a match, the sum total of the matches are added to the new_cig2 list before resetting the counter
        if (index == (len(cig_list)-1) or not re.search(match_type, cig_list[index+1])):
            new_match = str(match_sum) + match_type  
            new_cig2.append(new_match)
            #resets the counter for the number of matches
            match_sum = 0
    
    #joins up the values in the list to create the final CIGAR2 string
    fin_cig2 = str(''.join(new_cig2))
    return fin_cig2



def norm_cigar(old_cigar_string):
    '''
    Replaces the surjected CIGAR string with the non-surjected normal CIGAR string when the XG tag is specified i.e. not CIGAR2 string. 
    '''
    #replaces any specified matches or mismatches with M
    for ch in ['=','X']:
        if ch in old_cigar_string:
            old_cigar_string = old_cigar_string.replace(ch,"M")
    num = re.findall(r'\d+\D', old_cigar_string)
    if len(num) > 1:
        old_cigar_string = join_cigarstring(num)
    #passed the replaced string to the join_cigarstring function to return the joint CIGAR string
    return old_cigar_string 



def search(rname, readseq, table, f):
    '''
    Computes the hash for the read name from the BAM file and looks it up in the hash table. If the key is found, reads the JSON data at that position from the file and compares the readname and sequence to verify you really have a match. If there are multiple positions, each one is checked until you find a match or none. The JSON data at that position is passed to the create_cigar function. Returns the CIGAR2 string.
    '''
    i = hash(rname) % len(table)
    entry = table[i]
    while entry is not None:
        pos = entry.pos if isinstance(entry, Node) else entry
        f.seek(pos)
        line = f.readline() 
        data = json.loads(line)
        if (data[u'name'] == rname and data[u'sequence'] == readseq):
            if ("path" in data):
                new_cigar = create_cigar(line)
                aln_score = data[u'score'] if "score" in data else 0 
                return new_cigar, aln_score
            else:
                return "*", 0
        entry = entry.next if isinstance(entry, Node) else None




def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])



def main():
    try:
        bam = pysam.AlignmentFile(args.bam_in, "rb")
    except FileNotFoundError:
        print('The specified BAM file does not exist')
 
    out_bam = pysam.AlignmentFile(args.bam_out, "wb", template=bam)
    
    if not bam.has_index():
        #checks if BAM file has index, if not, the file is sorted and indexed (both files written out to same directory as input
        sorted_name = args.bam_in.replace(".bam","_sorted.bam") 
        pysam.sort("-o", sorted_name, args.bam_in)
        pysam.index(sorted_name)
        bam = pysam.AlignmentFile(sorted_name, "rb")

    #size of hash table
    SIZE = 2**24 

    #need to read in from command line
    with open(args.json_file,'r') as f:
        #calls build_table function to build the hash table for the file based on the given size
        table = build_table(f, SIZE)

        #reads through the BAM file and stores the read name and sequence for each read. This is passed to the search function which searches through the hash table to find the corresponding read in the JSON file. Once a match is found, the JSON data is used to create the CIGAR2 string for the read alignments against a variation graph based on the information within the .path.mapping[].edit arrays. The original CIGAR string is modified to the CIGAR2 string and the read is written to the output BAM file.
        for read in bam.fetch():
            rname = read.query_name
            readseq = reverse_complement(read.query_sequence) if read.is_reverse else read.query_sequence
            cigar_string, score = search(rname, readseq, table, f)
            if args.xg_tag: 
                read.set_tag("XG", cigar_string)
                read.cigarstring = norm_cigar(cigar_string)
            else:
                read.cigarstring = cigar_string
            read.set_tag('AS', score)
            out_bam.write(read)

    bam.close()
    out_bam.close()



if __name__ == '__main__':
    #Instantiate the parser
    parser = argparse.ArgumentParser(description='Produces the FAT-CIGAR string for BAM files containing read alignments against a variation graph')

    #Required positional argument
    parser.add_argument('json_file', type=str,
                    help='A required argument stating the input json file (created from the alignment GAM file using vg view)')
    parser.add_argument('bam_in', type=str,
                    help='A required argument stating the input surjected BAM file (created from the alignment GAM file using vg surject)')
    parser.add_argument('bam_out', type=str,
                    help='A required argument stating the output BAM file for the CIGAR2 string')

    #Optional arguments
    parser.add_argument('-xg','--xg_tag', default=False, action='store_true',
                    help='Writes the modified CIGAR string as the \'XG\' tag (overwrites the surjected CIGAR string as default)')

    #Parse arguments
    args = parser.parse_args()

    #creates a hash table using the namedtuple function
    Node = namedtuple('Node', ['pos', 'next'])

    main()

