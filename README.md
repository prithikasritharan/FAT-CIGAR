# FAT-CIGAR

The FAT-CIGAR python script allows you to obtain exact CIGAR-like string representations of sequence read alignments against both linear and graph-based reference genomes, without the masking of any bases. We call these exact, extended representations FAT-CIGAR strings. Obtaining these strings is useful if you want to compare precisely how well a set of sequence reads map to both a linear and a graph reference. Scores are also given for each precise mapping, using the same penalty scheme. The input to the FAT-CIGAR script is a BAM file (see below for which read mappers provide suitable input) and, in the case of working with alignments from the vg toolkit, an additional JSON file. All output is again encoded within a BAM file, allowing for flexible downstream analysis.   
&nbsp;

&nbsp;


## Works With
The FAT-CIGAR script is able to work with alignments generated from the following three read mappers:
* BWA
* vg
* SevenBridges Graph Genome toolkit
&nbsp;

&nbsp;


## Dependencies
**Samtools**

For alignments against a reference genome with BWA, the MD tag must also be present for all the reads within the BAM file by running the samtools calmd function prior to running the script as follows:
```python
samtools calmd aln.bam ref.fasta
```
&nbsp;


## Usage

### BWA

For alignments produced by the BWA mapper, the script must be run using the `linear` command. It takes an input BAM file with an MD tag present and outputs a BAM file containing the FAT-CIGAR string. As the script is reliant on BAM files that contain this MD tag, it must be used after running the samtools routine calmd as described in the Dependencies section above. 

The script should be run with the following commands:
```python
python fat-cigar.py linear -h -xg -g -cs input_bam output_bam
```

**Default arguments:**

When the script is run with just the required `input_bam` (input BAM file) and `output_bam` (output BAM file) arguments, without any optional arguments, the original CIGAR string present in the input BAM file will be replaced with the FAT-CIGAR string in the output BAM file.

**Options:**

The `xg` option  preserves the original CIGAR string from the input BAM file while writing out the FAT-CIGAR string as the XG tag, thereby producing an output BAM file with both CIGAR and FAT-CIGAR strings. 

The `g` option can be used to obtain the global alignment scores of the read mapping. As the alignment scores output by BWA are calculated using the local alignment during the Smith-Waterman extension phase instead of the final banded global alignment, it may not be reflective of the actual alignment. The true global alignment scores will be calculated from the FAT-CIGAR string, replacing the AS tag in the output BAM file. 

The `cs` option is used to obtain the short-form of the CS string, a further alternative alignment representation as used by the Minimap2 mapper. 

&nbsp;

### SevenBridges Graph Genome Toolkit

Graphical read mappers such as vg and the SevenBridges Graph Genome toolkit are capable of producing both surjected and non-surjected alignments of a read to a reference structure. A surjected alignment refers to an alignment against a graph-based reference genome that is represented as if it were aligned against the linear reference genome. Conversely, a non-surjected alignment refers to an alignment against the graph itself. Alignments output from the Graph Genome toolkit contain the surjected CIGAR string and the non-surjected FAT-CIGAR string as the XG tag. The `graph_sb` command should be used in order to generate the non-surjected CIGAR string from the FAT-CIGAR string. 

The script should be run as follows:
```python
python fat-cigar.py graph_sb input_bam output_bam
```

The script requires both the `input_bam` (input BAM file) and `output_bam` (output BAM file) arguments which will replace the surjected CIGAR string with the non-surjected CIGAR string in the output BAM file. The FAT-CIGAR string may also be missing for reads that match exactly against the reference, therefore, the script also checks whether the XG tag is present and if not, writes the FAT-CIGAR string to the XG tag. 

&nbsp;  

### vg Toolkit
The `graph_sb` command should be used to produce FAT-CIGAR strings for alignments of sequence reads against variation graphs from the vg toolkit (https://github.com/vgteam/vg). BAM files produced by vg contain the surjected CIGAR string, with FAT-CIGAR output containing non-surjected FAT-CIGAR strings.

The script should be run as follows:     
```python
python fat-cigar.py -h -xg json_file bam_in bam_out
```


**Default arguments:**

At a minimum, the script requires as input the JSON file containing the alignment information, the surjected BAM file and an output BAM file. The default arguments replace an input (surjected) CIGAR string with an output (non-surjected) FAT-CIGAR string in the output BAM file. The alignment scores will also be written to the AS tag in the output BAM file as default.

**Options:**

The `xg` option can be specified to preserve the (surjected) CIGAR string and to write the (non-surjected) FAT-CIGAR string as the XG tag. 

**Obtaining the vg input files:**

The JSON file can be obtained from the vg GAM file using:
The script requires the JSON file containing the alignment information, the surjected BAM file and an output BAM file to be provided as input to run. The `xg` option can be specified to write the FAT-CIGAR string as the XG tag. The alignment scores will also be written to the AS tag in the output BAM file as default. 

The JSON file can be obtained from the GAM file using:
```python
vg view -aj alignments.gam > output.json
```

The surjected BAM file can be obtained using the XG index of the graph and the GAM file:   
```python
vg surject -x index.xg -b alignments.gam > output.bam
```

It is preferred that the GAM file is sorted using `vg gamsort` prior to obtaining the JSON and BAM files. If there are >5,000,000 sequence reads, it is recommended that the GAM file is split up into multiple smaller files using `vg chunk` before converting into JSON and BAM format in order to speed up runtime. The script should then be run on each individual file and the output BAM files can be merged into a single file with the `samtools merge` function. 
     
