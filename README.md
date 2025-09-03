# ECCFP - EccDNA Caller based on Consecutive Full Pass
A software package for identifying eccDNAs in the ONT sequencing data of RCA-amplified eccDNA
## Introduction
ECCFP utilizes all consecutive full passes to determine accurate eccDNA positions and generate consensus sequences, based on rolling circle amplification and nanopore sequencing.  
[![ECCFP: BioRxiv](https://img.shields.io/badge/DOI-10.1101/2025.05.13.653627-blue)](https://doi.org/10.1101/2025.05.13.653627)  
|![ECCFP](./images/ECCFP%20Figure%201.png)|
|:---------------------------------------:|
#### Reference
ECCFP is optimized from the [Flec](https://github.com/icebert/eccDNA_RCA_nanopore.git) software.
- [eccDNA_RCA_nanopore](https://github.com/icebert/eccDNA_RCA_nanopore.git)
- [minimap2](https://github.com/lh3/minimap2/)
#### Dependency (python packages)
- [numpy](https://numpy.org/)
- [pandas](https://pandas.pydata.org/)
- [pyfaidx](https://pypi.org/project/pyfaidx/)
- [pyfastx](https://pypi.org/project/pyfastx/)
- [Biopython](https://biopython.org)
## Installation
```
git clone https://github.com/WSG-Lab/ECCFP.git
cd ECCFP
conda env create --file environment.yml -y
conda activate ECCFP
pip install -e .

# Verify installation and display help message
eccfp --help
```
## Usage
#### simple usage
```
eccfp --fastq input.fastq --paf mapping.paf --reference ref.fa -o output
```
#### Parameter
```
usage: eccfp [options]

options:
  -h, --help            show this help message and exit

  Generate candidate eccDNA from rolling circle reads:

  --fastq FASTQ         input reads in fastq format
  --paf PAF             input alignments in PAF format
  --reference REFERENCE
                        reference genome sequences in fasta format
  --output OUTPUT, -o OUTPUT
                        output folder, the default is the current folder
  --maxOffset MAXOFFSET
                        maximum offset of start/end positions between two sub-reads to be considered as mapping to the same location, default is 20
  --minMapQual MINMAPQUAL
                        minimum mapping quality of sub-reads, default is 30

  Get accurate eccDNA location from candidate eccDNA:

  --fluctuate FLUCTUATE
                        maximum offset between the candidate eccDNA and the candidate eccDNA in a group, default is 20bp
  --nf NF               minimum number of fullpass to support candidate eccDNA, default is 2
  --cov COV             minimum fullpass used for merging a set of candidate eccDNA, default is 1
  --nc NC               minimum number of candidate eccDNA required to merge overlapped candidate eccDNAs, default is 2

  Generate consensus sequences and variants:

  --minDP MINDP         minimum depth to call variants, default is 4
  --minAF MINAF         minimum alternative allele frequency to call variants, default is 0.75

```
#### Example
```
cd example
minimap2 -cx map-ont ref.fa example.fastq --secondary=no -t 8 -o mapping.paf
eccfp --fastq example.fastq --paf mapping.paf --reference ref.fa -o output
```
## Output
Five files are generated following the completion of the pipeline: unit.txt and candidate_consolidated.txt(intermediate files), final_eccDNA.txt, consensus_sequence.fasta, and variant.txt (result files).   
|file|details|
|----|-------|  
|unit.txt|full pass alignment details for all candidate eccDNAs detected in reads|
|candidate_consolidated.txt|the consolidating steps used to derive accurate eccDNAs from the candidate eccDNAs|
|final_eccDNA.txt|accurate eccDNA information|
|consensus_sequence.fasta|the consensus sequences of accurate eccDNAs|
|variants.txt|variant profiles specific to these accurate eccDNAs|

###### Example final_eccDNA.txt file
|eccDNApos|Nfullpass|Nfragments|Nreads|refLength|seqLength|
|---------|---------|----------|------|---------|---------|
|chr16:757924-758274(+)|4|1|1|351|350|
|chr1:21828265-21828439(+)\|chr7:132167689-132167880(+)|24|2|2|367|367|

||description|
|---------|---------------|
|eccDNApos|eccDNA position|
|Nfullpass|Number of consecutive full pass for this eccDNA covered by all reads|
|Nfragments|Number of fragment that form this eccDNA|
|Nreads|Number of reads identified for the eccDNA|
|refLength|The length of reference genome that this eccDNA |
|seqLength|The length of consensus sequence that this eccDNA |

###### variants.txt file
||description|
|-------|-------|
|col1|chromsome|
|col2|position in the reference genome|
|col3|reference base|
|col4|variant|
|col5|supportive coverage depth|
|col6|total coverage depth|
|col7|type|
|col8|eccDNApos|

## Citation
> Zhang, T. et al. ECCFP: a consecutive full pass based bioinformatic analysis for eccDNA identification using Nanopore sequencing data. bioRxiv 2025.05.13.653627 (2025) doi:10.1101/2025.05.13.653627.
