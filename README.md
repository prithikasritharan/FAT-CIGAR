# FAT-CIGAR

The FAT-CIGAR python scripts allows you to obtain an exact CIGAR string representation, also referred to as the FAT-CIGAR string, of the alignment of sequence reads against both linear and graph-based reference genomes without the masking of any bases.

## Usage

The linear FAT-CIGAR script (fat-cigar_linear.py) is mainly utilised to obtain FAT-CIGAR string for alignment information within BAM files output by both the BWA mapper and the SevenBridges Graph Genome toolkit, used to create variation graphs. The script requires Samtools to be installed in order to run. For alignments against the traditional linear reference genome with BWA, the MD tag must also be present for all the reads within the BAM file by running the `samtools calmd` function prior to running the script. 

The script should be run with the following commands:
```python
python fat-cigar_linear.py -h -xg -c -g -cs bam_file out_file
```
When the script is run as default with just the required `bam_file` (input BAM file) and `out_file` (output BAM file) arguments, without any optional arguments, the original CIGAR string will be overwritten with the FAT-CIGAR string in the output BAM file. The `xg` option can be used to write out the FAT-CIGAR string as the XG tag in the output BAM file in order to prevent the original CIGAR string being overwritten. The `g` option can be used to obtain the global alignment scores. As the alignment scores output by BWA are calculated using the local alignment during the Smith-Waterman extension phase instead of the final banded global alignment, it may not be reflective of the actual alignment. The true global alignment scores will be calculated from the FAT-CIGAR string, replacing the AS tag in the output BAM file. The short-form of the CS tag, which is another form of alignment representation used by the Minimap2 mapper, can also be obtained using the `cs` option. 

Alignments output from the SevenBridges Graph Genome toolkit contain the surjected CIGAR string and the non-surject FAT-CIGAR string as the XG tag. A surjected alignment refers to alignments against a graph-based reference genome that is represented as if aligned against the linear reference genome whilst a non-surjected alignment refers to alignments against the graph itself. The `c` option allows you to obtain the non-surjected CIGAR string from the FAT-CIGAR string to be written to the output BAM file.  
  

The graph FAT-CIGAR script (fat-cigar_graph.py) is used specifically to produce the FAT-CIGAR string for alignments against variation graphs from vg as the BAM files produced by vg (https://github.com/vgteam/vg) contains the surjected CIGAR string. 

The script is run as follows:     
```python
python fat-cigar_graph.py -h -xg json_file bam_in bam_out
```
The script requires the JSON file containing the alignment information, the surjected BAM file and an output BAM file to be provided as input to run. The `xg` option can be specified to write the FAT-CIGAR string as the XG tag. The alignment scores will also be written to the AS tag in the output BAM file as default. 

The JSON file can be obtained from the GAM file using:
```vg view -aj alignments.gam > output.json```

The surjected BAM file can be obtained using the XG index of the graph and the GAM file:   
```vg surject -x index.xg -b alignments.gam > output.bam```

It is preferred that the GAM file is sorted using `vg gamsort` prior to obtaining the JSON and BAM files. If there are >5,000,000 sequence reads, it is recommended that the GAM file is split up into multiple smaller files using `vg chunk` before converting into JSON and BAM format in order to speed up runtime. The script should then be run on each individual file and the output BAM files can be merged into a single file with the `samtools merge` function. 
     
