# Orthomatic 

*Ortholog recovery using curated ortholog datasets and reciprocal best-match criterion.*

This pipeline uses a set of protein queries in a representative taxon (extracted from 
a user-provided set of orthologous peptides) to BLAST-search user-provided protein fasta files
(e.g., assembled & translated transcriptomes). 

Matches are then reciprocally BLASTed to the representative taxon's genome. Sequences 
satisfying the criterion of mutual-best-match are then pooled  and aligned by MAFFT.

-----

## Requirements & Install
The orthomatic script can be run directly on Unix machines, with the following dependencies in the user's $PATH.

### Dependencies (must be in $PATH )

 - NCBI BLAST 2.2.31+ (specifically: blastp, makeblastdb)		
 - GNU parallel																								
 - CD-HIT (cdhit), if -c invoked	
 - MAFFT (mafft)	
 - selectSeqs.pl 	(provided here, written by Naoki Takebayashi)
 - seqConverter.pl (for nexus input option, provided here, by Olaf Bininda-Emonds)

[Source for selectSeqs.pl](http://raven.wrrb.uaf.edu/~ntakebay/teaching/programming/perl-scripts/perl-scripts.html)
[Source fro seqConverter.pl](https://uol.de/systematik-evolutionsbiologie/programme)																
### Accessory scripts (included): 		
These two scripts must be in the working directory   
								
- pullOGsfromBlast.py											
- parse_recipBLAST.py		

### Prepare scripts 

> chmod +x orthomatic.sh
> chmod +x parse_recipBLAST.py
> chmod +x pullOGsfromBlast.py

----

## Inputs
- directory containing a set of orthologous peptides, in either separate fasta files or in NEXUS format
	- see example in OGs.zip at https://datadryad.org/stash/dataset/doi:10.5061/dryad.k6tq2
- reference taxon with directory containing fasta file of genome protein models 
	- the reference taxon must be represented in the ortholog dataset
- directory containing protein fastas in which to search orthologs
- file extensions: fasta files should end in ".fa" and nexus in ".nex"
- fasta headers: alphanumeric and _ characters only (no spaces)
	- See fasta example: http://datadryad.org/bitstream/handle/10255/dryad.98320/OGs.zip?sequence=1

----

## Usage

Example: 

> orthomatic.sh -T 8 -i OGs -r Nematostella -t taxa_blastdatabase -e 1e-20 -g ref_genome_blastdatabase


### Parameters

>	-h 		help
>	
>	-D 		option to delete temporary files 
>	
>	-c  	option to run CDHIT on FASTAs specified by -t or -m
>	
>	-T  	number of threads
>	
>	-i  	directory with either: 
>			one fasta per gene; each fasta has a sequence from Reference Species
>			-OR-
>			single nexus file with concatenated data; gene partitions specified in charset block
>				
>	-r 		Reference Species handle; no spaces, name must match fasta headerline in OGDIR/ files
>	
>	-t 		directory with peptide fasta files for each taxon to search, including any reference species
>				
>	-m     add new taxon dataset to existing results (put FASTA in TAXDB directory)
>               **Requires directory of previous alignments called OGfastasALN**
>                 
>	-g 		directory with peptide fasta genome for reference species. File must start with REFSP handle. 
>	-e 		e-value cut-off for blast

---

## Output

Fasta files for the recovered orthologs will be found the directories **OGfastasets** (unaligned) and **OGfastasALN** (MAFFT-aligned)

