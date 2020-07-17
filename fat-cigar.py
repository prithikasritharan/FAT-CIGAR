import os
import sys
import argparse
import re
import time
import json
import pip
from collections import namedtuple
try:
    from itertools import groupby, imap, takewhile
except ImportError:
    imap=map
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
    Reads in the JSON data for a read (found through searching the hash table for the file position) and converts the individual path edits into the corresponding CIGAR string equivalent. The disjoined list of CIGAR strings are passed to the join_cigarstring function to join the list into a string. Returns the FAT-CIGAR string.
    ''' 
    data = json.loads(line)

    #Intialise list for appending the CIGARs based on the path mapping edits
    cigar2 = []

    max_node = len(data[u'path'][u'mapping'])
    
    #cycles through the path mapping edit information for each read and uses the to_length and from_length to recreate the FAT-CIGAR string
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


    #calls the create_cigar function to join the values in the FAT-CIGAR list to create the final FAT-CIGAR string 
    final_cigar = join_cigarstring(list(reversed(cigar2))) if 'is_reverse' in data[u'refpos'][0] else join_cigarstring(cigar2) 
    return final_cigar 



def global_score_calc(new_cig, n_base):
    '''
    Takes in the modified CIGAR string and calculates the global alignment score using the same scoring system utilised by BWA
    '''
    #splits up the CIGAR string into a list storing each operation as a separate value
    cig_split = ["".join(x) for _, x in groupby(new_cig, key=str.isdigit)]

    score = 0
    points = 0
    item = iter(cig_split)
    for i in item:
        j = next(item)
        if j == "S" or j == "H":
            points = int(i)*0  
        if j == "X":
            points = int(i)*-4
        if j == "=":
            points = int(i)*1
        if j == "I" or j == "D":
            points = -6 + (int(i)*-1)
        score += int(points)
    #adds a minus one penalty from the total matches based on the number of N's in the read seq
    score = score - (n_base*2)
    return score



def join_cigarstring(cig_list):
    '''
    Joins the matches, mismatches, indels and soft clips from a list of disjoined CIGAR strings. Returns the final FAT-CIGAR string. 
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
    
    #joins up the values in the list to create the final FAT-CIGAR string
    fin_cig2 = str(''.join(new_cig2))
    return fin_cig2


def join_cigarstring_linear(split_cigar):
    '''
    Joins the matches, mismatches, indels and soft clips from a list of disjoined CIGAR strings. Returns the final FAT-CIGAR string. 
    '''
    #initialize new list for the joined CIGAR strings
    joint_cigar = []
    #count for the total number of matched sequences
    match_sum = 0
 
    #cycles through the list containing the individual CIGAR strings and joins them together in the new_cig2 list
    for index, i in enumerate(split_cigar):
        if i.isdigit():
            match_sum += int(i)
            match_type = split_cigar[index+1]
            #if the index is smaller than the list length and the next value is also a match, moves onto the next iteration
            if (index < (len(split_cigar)-2) and match_type == split_cigar[index+3]):
                continue
            #if the index is the last position in the list or the next value in the list is not a match, the sum total of the matches are added to the new_cig2 list before resetting the counter
            if (index == (len(split_cigar)-2) or not match_type == split_cigar[index+3]):
                new_match = str(match_sum) + match_type  
                joint_cigar.append(new_match)
                #resets the counter for the number of matches
                match_sum = 0

    #joins up the values in the list to create the final CIGAR2 string
    return str(''.join(joint_cigar))



def make_cs_tag(final_cig, read_seq, ref_seq):
    '''
    Builds short form of cs tag using the reference and query sequences based of the modified FAT-CIGAR string
    '''
    ref_seq = ref_seq.lower()
    read_seq = read_seq.lower()

    #splits up the CIGAR string into a list storing each operation as a separate value
    cig_split = ["".join(x) for _, x in groupby(final_cig, key=str.isdigit)]

    #pointer for reference base
    ref_pointer = 0
    #pointer for base in read sequence
    read_pointer = 0

    cs_tag = []
    item = iter(cig_split)

    #loops through the modified CIGAR string and converts it to cs tag format by finding the corresponding reference/read sequence
    for i in item:
        j = next(item)
        if j == "=":
            ref_pointer += int(i)
            read_pointer += int(i)
            tag = ":" + i
            cs_tag.append(tag)
        elif j == "X":
            for m in range(0,int(i)):
                ref_pointer += 1
                read_pointer += 1
                tag = "*" + ref_seq[ref_pointer-1] + read_seq[read_pointer-1]
                cs_tag.append(tag)
        elif j == "I":
            tag = "+" + read_seq[read_pointer:read_pointer+int(i)]
            cs_tag.append(tag) 
            read_pointer += int(i)
        elif j == "D":
            tag = "-" + ref_seq[ref_pointer:ref_pointer+int(i)]
            cs_tag.append(tag) 
            ref_pointer += int(i)
        elif j == "S":
            read_pointer += int(i)
        else:
            continue 
    return "".join(cs_tag)



def mismatch(match):
    '''
    Takes mismatched bases from the mdtag and converts it into the CIGAR string format (X) 
    '''
    diff_match = match.group()
    return str(len(diff_match)) + 'X'



def my_deletion(match):
    '''
    Finds deletions from the mdtag and converts them into CIGAR string format (D) by replacing it with the sequence length
    '''
    match = match.group()
    match = re.sub('\^', '', match)
    return str(len(match)) + 'D'



def my_insertion(ins, cig, new_cig):
    '''
    Takes the position of insertion sequences from the CIGAR string and uses it to identify the insertion position in the md_tag.
    Breaks up the match sequence containing the masked insertion and calculates the new matches on either side of the insertion sequence.    
    '''
    divider = re.compile("\d+\D")
    #obtains a list containing each CIGAR string types as separate values 
    cig_split = divider.findall(cig)

    insertion_seq=[]
    for i, item in enumerate(cig_split):
        if (re.search(ins, item)):
            inserted = strip_list(cig_split[:i])
            #calculates the position before the insertion sequence  
            pos = sum_list(inserted)
  
            #split up the md tag into list
            new_mod = divider.findall(new_cig)

            #calls strip_list to remove any non-numbers and converts tag to a list of integer values
            mod_strip = strip_list(new_mod)

            total = 0
            j = 0
            while total <= pos:
                total += mod_strip[j]
                j += 1

            end_match = total - pos
            ins_pos = mod_strip[j-1]
            start_match = ins_pos - end_match
            end_match = re.sub(r'(\d+)', r'\1=', str(end_match))
            start_match = re.sub(r'(\d+)', r'\1=', str(start_match))
            insertion_seq = [item,end_match]
            new_mod[j-1] = start_match
            new_mod = ''.join(new_mod[:j] + insertion_seq + new_mod[j:]) 
            new_cig = re.sub(r'\d+(?=[ACGT])','', new_mod)
    return new_cig



def norm_cigar(old_cigar_string):
    '''
    Replaces the surjected CIGAR string with the non-surjected normal CIGAR string when the XG tag is specified i.e. not FAT-CIGAR string. 
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



def norm_cigar_linear(old_cigar_string):
    '''
    Replaces the surjected CIGAR string with the non-surjected normal CIGAR string when the XG tag is specified i.e. not CIGAR2 string. 
    '''
    #replaces any specified matches or mismatches with M
    for ch in ['=','X']:
        if ch in old_cigar_string:
            old_cigar_string = old_cigar_string.replace(ch,"M")

    num = re.findall(r'\d+\D', old_cigar_string)
    if len(num) > 1:
        old_cigar_string = join_cigarstring_linear(["".join(x) for _, x in groupby(old_cigar_string, key=str.isdigit)])

    #passed the replaced string to the join_cigarstring_linear function to return the joint CIGAR string
    return old_cigar_string 



def search(rname, readseq, table, f):
    '''
    Computes the hash for the read name from the BAM file and looks it up in the hash table. If the key is found, reads the JSON data at that position from the file and compares the readname and sequence to verify you really have a match. If there are multiple positions, each one is checked until you find a match or none. The JSON data at that position is passed to the create_cigar function. Returns the FAT-CIGAR string.
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



def strip_list(a):
    '''
    Strips a list to remove any non-digits
    '''
    stripped = list(imap(int, [re.sub('[\D]', '', i) for i in a]))
    return stripped


def sum_list(l):
    '''
    Sums up all the values within a list
    '''
    sum = 0
    for x in l:
        sum += x
    return int(sum)



def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])



def main():
    try:
        bam = pysam.AlignmentFile(args.input_bam, "rb")
    except FileNotFoundError:
        print('The specified BAM file does not exist')
 
    out_bam = pysam.AlignmentFile(args.output_bam, "wb", template=bam)
    
    if not bam.has_index():
        #checks if BAM file has index, if not, the file is sorted and indexed (both files written out to same directory as input
        sorted_name = args.input_bam.replace(".bam","_sorted.bam") 
        pysam.sort("-o", sorted_name, str(args.input_bam))
        pysam.index(sorted_name)
        bam = pysam.AlignmentFile(sorted_name, "rb")


    #VG toolkit
    if (args.mode == "graph_vg"):
        #size of hash table
        SIZE = 2**24 

        #need to read in from command line
        with open(args.json_file,'r') as f:
            #calls build_table function to build the hash table for the file based on the given size
            table = build_table(f, SIZE)

            #reads through the BAM file and stores the read name and sequence for each read. This is passed to the search function which searches through the hash table to find the corresponding read in the JSON file. Once a match is found, the JSON data is used to create the FAT-CIGAR string for the read alignments against a variation graph based on the information within the .path.mapping[].edit arrays. The original CIGAR string is modified to the FAT-CIGAR string and the read is written to the output BAM file.
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



    #SevenBridges Graph Genome toolkit
    if (args.mode == "graph_sb"):
        for read in bam.fetch():
            if not read.is_unmapped:
                try:
                    read.get_tag("XG")   
                except KeyError:
                    read.set_tag("XG", read.cigarstring.replace("M","="))
                read.cigarstring = norm_cigar_linear(read.get_tag("XG"))
            out_bam.write(read)



    #BWA
    if (args.mode == "linear"):
        for read in bam.fetch():
            if not read.is_unmapped:
                cigar = read.cigarstring
                try:
                    md_tag = read.get_tag("MD")
                except KeyError:
                    print("Tag 'MD' is missing for reads within the BAM file. The 'MD' tag is required for all reads.\nProgram Terminating")
                    os._exit(1)

                #removes the tag name from the tag 
                md_tag = md_tag.replace("MD:Z:","")
                #adds M for match next to any digits 
                md_mod = re.sub(r'(\d+)', r'\1=', md_tag)
                #replaces deletion sequences with the length of seq + D rather than the seq itself
                md_del = (re.sub(r'\^([ATGC]+)', my_deletion, md_mod)).strip()
                #removes any 0M from the md tag
                md_nomatch = re.sub(r'^0=|(?<=\D)0=', '', md_del)
                #replaces any mismatched bases with the length of the base + X
                new_cigar = (re.sub(r'([ATGC]+)', mismatch, md_nomatch)).strip()

                #adds both the soft and hard clippings from the CIGAR strings to the MD tag 
                start = re.compile("^\d+[SH]")
                #checks for clippings at the start of the CIGAR string  
                if (re.search(start, cigar)):
                    start_clip = start.search(cigar)
                    start_match = start_clip.group()
                    new_cigar = start_match + new_cigar

                #checks for clippings at the end of the CIGAR string
                end = re.compile("\d+[SH]$")
                if (re.search(end, cigar)):
                    end_clip = end.search(cigar)
                    end_match = end_clip.group()
                    new_cigar += end_match

                #calculates the position and adds any insertions into the md tag  
                insertion = re.compile("\d+I")
                if (re.search(insertion, cigar)):
                    new_cigar = my_insertion(insertion, cigar, new_cigar)
             
                #Swap N's in read sequence to matches
                n_seq = []

                if 'N' in read.query_sequence:
                    n_seq = [i for i, base in enumerate(read.query_sequence) if base == 'N']
                    for pos in n_seq:
                        #strips the CIGAR string into a list
                        cig_split = ["".join(x) for _, x in groupby(re.sub(r'\d+D', '', new_cigar), key=str.isdigit)]
                        total = 0
                        for index, i in enumerate(cig_split): 
                            if i.isdigit():
                                if total < pos:
                                    total += int(i)
                                    continue
                                #checks for a single mismatch at position
                                if cig_split[index] == "1" and cig_split[index+1] == "X":
                                    cig_split[index+1] = "="
                                    #pass to fuction that joins cigar together
                                    new_cigar = join_cigarstring_linear(cig_split)  
                                    break
                                #if there is more than one mismatch, calculates the exact position of the N base to convert to match
                                if cig_split[index] != "1" and cig_split[index+1] == "X":
                                    mismatch_split = int(i) * ["1X"]
                                    mismatch_split[(pos - total)] = "1=" 
                                    cig_split[index:index+2] = ["".join(x) for _, x in groupby("".join(mismatch_split), key=str.isdigit)]
                                    new_cigar = join_cigarstring_linear(cig_split)
                                    break
                            else:
                                continue
                            break         

                #assigns the new CIGAR string containing the mismatches to the read
                if args.xg_tag: 
                    read.set_tag("XG", new_cigar)
                else:
                    read.cigarstring = new_cigar

                #changes the alignment score to the final global alignment score
                if args.global_score: read.set_tag("AS", global_score_calc(new_cigar, len(n_seq)))
 
                #sets the short form of the cs tag if option selected by user  
                if args.cs_tag: read.set_tag("cs", make_cs_tag(new_cigar, read.query_sequence, read.get_reference_sequence()))

            out_bam.write(read)

    bam.close()
    out_bam.close()



if __name__ == '__main__':
    #Instantiate the parser
    parser = argparse.ArgumentParser(description='Produces the FAT-CIGAR string for BAM files containing read alignments against a variation graph')

    subparsers = parser.add_subparsers(help='commands', dest='commands')

    #Parser for working with vg
    vg_parser = subparsers.add_parser('graph_vg', help='Generates the FAT-CIGAR string for output from the vg toolkit')
    vg_parser.add_argument('json_file', type=str,
                    help='A required argument stating the input json file (created from the alignment GAM file using vg view)')
    vg_parser.add_argument('input_bam', type=str,
                    help='A required argument stating the input surjected BAM file (created from the alignment GAM file using vg surject)')
    vg_parser.add_argument('output_bam', type=str,
                    help='A required argument stating the output BAM file for the FAT-CIGAR string')
    vg_parser.add_argument('-xg','--xg_tag', default=False, action='store_true',
                    help='Writes the FAT-CIGAR string as the \'XG\' tag (overwrites the surjected CIGAR string as default)')    


    #Parser for working with SevenBridges graph genome toolkit
    sb_parser = subparsers.add_parser('graph_sb', help='Generates the non-surjected CIGAR string for output from the SevenBridges Graph Genome toolkit using the \'XG\' tag information')
    sb_parser.add_argument('input_bam', type=str,
                    help='A required argument stating the input BAM file')
    sb_parser.add_argument('output_bam', type=str,
                    help='A required argument stating the output BAM file to which the non-surjected CIGAR string will be written (overwrites the original CIGAR string)')


    #Parser for working with BWA
    linear_parser = subparsers.add_parser('linear', help='Generates the FAT-CIGAR string for output from BWA')
    linear_parser.add_argument('input_bam', type=str,
                    help='A required argument stating the input BAM file (The \'MD\' tag must be present for all reads within BAM file)')
    linear_parser.add_argument('output_bam', type=str,
                    help='A required argument stating the output BAM file to which the FAT-CIGAR string will be written')
    linear_parser.add_argument('-xg','--xg_tag', default=False, action='store_true',
                    help='Writes the FAT-CIGAR string as the \'XG\' tag (overwrites the original CIGAR string as default)')
    linear_parser.add_argument('-g','--global_score', default=False, action='store_true',
                    help='Calculates the global alignment score from the FAT-CIGAR string')
    linear_parser.add_argument('-cs','--cs_tag', default=False, action='store_true',
                    help='Adds the short form of the \'cs\' tag (also used to represent differences between query and reference sequence by Minimap2)')


    #Parse arguments
    args = parser.parse_args()

    #creates a hash table using the namedtuple function
    Node = namedtuple('Node', ['pos', 'next'])

    main()

