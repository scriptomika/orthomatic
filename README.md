# Orthomatic 

*Ortholog recovery using curated ortholog datasets and reciprocal best-match criterion.*

This pipeline uses a set of protein (or nucleotide) queries in a representative taxon (extracted from 
a user-provided set of orthologous peptides) to BLAST-search user-provided protein (or nucleotide) fasta files
(e.g., assembled & translated transcriptomes, metagenome assemblies, etc). 

Matches are then reciprocally BLASTed to the representative taxon's genome. Sequences 
satisfying the criterion of mutual-best-match are then pooled  and aligned by MAFFT.

-----

## Requirements & Install
Once the local scripts are executable (chmod commands below), the orthomatic script can be run directly on Unix machines, with the following dependencies in the user's $PATH.

### Dependencies (must be in $PATH )

 - NCBI BLAST 2.2.31+ (specifically: blastp, makeblastdb)		
 - GNU parallel																								
 - CD-HIT (cd-hit), if -c invoked	
 - MAFFT (mafft)	
 - selectSeqs.pl 	(provided here, written by Naoki Takebayashi)
 - seqConverter.pl (in order to use nexus input option, provided here, by Olaf Bininda-Emonds)

[Source for selectSeqs.pl](http://raven.wrrb.uaf.edu/~ntakebay/teaching/programming/perl-scripts/perl-scripts.html)

[Source fro seqConverter.pl](https://uol.de/systematik-evolutionsbiologie/programme)																
### Accessory scripts (included): 		
These two scripts must be in the $PATH  
								
- pullOGsfromBlast.py											
- parse_recipBLAST.py		

### Prepare scripts 
```
$ chmod +x parse_recipBLAST.py
$ chmod +x pullOGsfromBlast.py
$ chmod +x orthomatic.sh
$ mv parse_recipBLAST.py pullOGsfromBlast.py orthomatic.sh ~/path/to/scripts/ #specify directory
$ export PATH=$PATH:~/path/to/scripts
```


----

## Inputs
1. directory containing a set of orthologous peptides or nueleotides, in either separate fasta files or in NEXUS format
	- see example in OGs.zip at https://datadryad.org/stash/dataset/doi:10.5061/dryad.k6tq2
2. reference taxon string
  	- example: 'Amphimedon' 
3. directory containing fastas to be searched
	- including: fasta of complete protein/nucleotide models (ideally a genome) from reference taxon
	- the reference taxon string must be found in the sequence headers
 	- example: >Amphimedon_seq1_c0_gi_1	 
   
- file extensions: fasta files should end in ".fa" and nexus in ".nex"
- fasta headers: alphanumeric and _ characters only (no spaces)
	- See fasta example: http://datadryad.org/bitstream/handle/10255/dryad.98320/OGs.zip?sequence=1

### Example: Renaming fasta headers based on filename
Example: rename headers on Transdecoder output
N.B. example file extension is 'transdecoder.fasta'.
```
#adds sample name and leading zeros to contig fastas
$ cd raw_fastas/ 
$ ls
Agelas_conifera.transdecoder.fasta
$ head -1 *fasta
>Agcon_Gene.1::g.1::m.1  ORF type:5prime_partial len:194 (+),score=50.48 Agcon_00001:3-584(+)

$ for fa in *transdecoder.fasta ; do f=${fa%%.*}; f1=$(basename $f); echo $f1;
$	awk '/^>/{printf ">" "%06d\n", ++i; next}{print}' < $fa |\ #replaces headers with numbers 0-999999
$	awk -v new=">${f1}_" '{sub(/>/, new)}1'   > ${f1}.fa 2>/dev/null; #appends file handle to numbered headers
$ done

$ ls
Agelas_conifera.transdecoder.fasta	Agelas_conifera.fa

$ head -1 *.fa
>Agelas_conifera_000001

$mv Agelas_conifera.fa ../newdata_fastas
```

### Example: Setting up Reference species
Example reference species: *Amphimedon queenslandica*

This species should have sequences represented in the curated gene sets supplied to Orthomatic AND a fasta file (preferably a genome) in the fasta directory.

REFSP string will be "Amphimedon", as that string is found in 1) the  gene set seq header, in 2) the species' fasta file name and in 3) its fasta sequence headers

```
$ grep "Amphimedon" ortholog_dir/gene1.fa
>Amphimedon_gene_xyz

# (the Reference species string can occur anywhere in the sequence ID: >genexyz_Amphimedon_acoi3776j3 would also work)

# fasta directory includes only 1 file begininning with "Amphimedon"
$ ls newdata_fastas/
Amphimedon_queenslandica.genome.fa
Amph_compressa.fa #other Amphimedons are labelled to avoid a match
Amph_viridis.fa

$ head -1 newdata_fastas/Amphimedon_queenslandica.genome.fa
>Amphimedon_queenslandica_00001 

$ orthomatic.sh -T 24 -r Amphimedon -t newdata_fastas -i ortholog_dir -e 1e-50
```

----

## Usage

Example: 
```
 orthomatic.sh -T 8 -i ortholog_dir -s prot -r Amphimedon -g Aque_genome -t newdata_fastas -e 1e-20 
```

### Parameters
```
	-h 		help	
	-c  	option to run CDHIT on FASTAs specified by -t or -m	
	-T  	number of threads	
	-s      sequence format (prot or nucl)	
	-i  	directory with either: 
			one fasta per gene; each fasta has a sequence from Reference Species
			-OR-
			single nexus file with concatenated data; gene partitions specified in charset block				
	-r 		Reference Species handle; no spaces, name must match fasta headerline in OGDIR/ files	
	-t 		directory with peptide fasta files for each taxon to search, including any reference species
	-m     option to add new taxon dataset to existing results (put FASTA in TAXDB directory)
              **Requires directory of previous alignments called orthomatic_alignments**               
	-g 		directory with peptide fasta genome for reference species. File must start with REFSP handle. 
	-e 		e-value cut-off for blast
	-k 		option to keep temporary files 	(for troubleshooting, discard is default)

```
---

## Output

MAFFT-Aligned fasta files for the recovered orthologs will be found the directory **orthomatic_alignments**.

