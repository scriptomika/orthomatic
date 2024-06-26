#!/bin/bash


# Usage info
show_help() {
cat << EOF
Example: ${0##*/} -T 8 -i OGs -r Nematostella -t taxa_blastdatabase -e 1e-20 -g ref_genome_blastdatabase

Usage: ${0##*/} [-h] [-D] [-c] [-m] [-T num] [-i OG_dir] [-r reference_species] [-t taxa_fastas_dir] [-g ref_genome_fasta_dir] [-e num]

	-h <help>	display this help and exit
	
	-D DEL		option to delete temporary files 
	
	-c CDHIT	option to run CDHIT on FASTAs specified by -t or -m
	
	-T MAXCPU	number of threads
	
	-i OGDIR	directory with either: 
				one fasta per gene; each fasta has a sequence from Reference Species
				-OR-
				single nexus file with concatenated data; gene partitions specified in charset block
				
	-r REFSP	Reference Species handle; no spaces, name must match fasta headerline in OGDIR/ files
	
	-t TAXDB	directory with CDS (coding nucleotide) fasta files for each taxon to search, including any reference species
				here they have been cdhit
				
	-m MORET    add new taxon dataset to existing results (put FASTA in TAXDB directory)
                **Requires directory of previous alignments called OGfastasALN**
                 
	-g REFDB	directory with CDS fasta genome for reference species. File must start with REFSP handle. 
	-e EVAL		e-value cut-off for blast

No spaces in any FASTA headers. 
Fasta file extension uses ".fa"
Nexus file extension uses ".nex"

See fasta example: http://datadryad.org/bitstream/handle/10255/dryad.98320/OGs.zip?sequence=1

#################################################################
#
#  Orthomatic BLASTs fastas (e.g., assembled & translated 
#  transcriptomes) using queries from designated reference taxon 
#  for a provided set of orthologous CDS.  Matches are 
#  reciprocally BLASTed to the reference taxon\'s genome. Sequences 
#  satisfying the criterion of mutual-best-match are then pooled 
#  and aligned by MAFFT.
#
#################################################################
# 
# Dependencies (must be in $PATH): 								
#   NCBI BLAST 2.2.31+ (specifically: blastn, makeblastdb)		
#   GNU parallel																								
#   CD-HIT (cdhit), if -c invoked																									
# 																	
# Dependencies (included): 		
#   
#   selectSeqs.pl 								
#   pullOGsfromBlast.py											
#   parse_recipBLAST.py		
#   seqConverter.pl (for nexus input option)
#									
#################################################################

EOF
}

# Initialize our own variables:
OGDIR=""
REFSP=""
TAXDB=""
REFDB=""
maxCPU=""
EVAL=""
CDHIT=""
DEL=""
MORET=""
DELY=""

OPTIND=1
Tflag=false; iflag=false; rflag=false; tflag=false; gflag=false; eflag=false
Terr="Missing: -T <number of threads>"
ierr="Missing: -i <directory with query gene fastas>"
rerr="Missing: -r <reference_species_name>"
terr="Missing: -t <directory with taxon fastas>"
gerr="Missing: -g <directory with reference genome fasta(s)>"
eerr="Missing: -e <e-value, e.g., 1e-20>"
merr="Missing: -m <name of new taxa FASTA file>"

while getopts hT:i:r:t:g:e:m:cD opt; do
    case $opt in
        h)  show_help; exit 0;;
        T)  maxCPU=$OPTARG; if [ ! $OPTARG ]; then echo $Terr; echo "boo"; show_help; exit; fi;;
        i)  OGDIR=$OPTARG; if [ ! -d $OPTARG ]; then echo $ierr; show_help; exit; fi;;
        r)  REFSP=$OPTARG; if [ ! $OPTARG ]; then echo $rerr; show_help; exit; fi;;
        t)  TAXDB=$OPTARG; if [ ! -d $OPTARG ]; then echo $terr; show_help; exit; fi;;
        g)  REFDB=$OPTARG; if [ ! -d $OPTARG ]; then echo $gerr; show_help; exit; fi;;
        e) 	EVAL=$OPTARG;  if [ ! $OPTARG ]; then echo $eerr; show_help; exit; fi;;
        m)  MORET=$OPTARG ;;
        c)  CDHIT=$OPTARG; echo "Input taxa fastas will be consolidated with CD-HIT" ;;
        D)  DEL=$OPTARG; echo "Clean-up selected: Temporary files will be deleted"; DELY="yes" ;;
        \?) echo "Unknown option: -$OPTARG" >&2; show_help; exit 1;;
        :) echo "Missing option argument for -$OPTARG" >&2; show_help; exit 1;;
        *) echo "Unimplemented option: -$OPTARG" >&2;show_help; exit 1;;
    esac
done
shift "$((OPTIND-1))" # Shift off the options and optional --.

if [[ -z $OGDIR || -z $TAXDB || -z $REFSP || -z $EVAL || -z $maxCPU || -z $REFDB ]]; then show_help; exit 1
elif ! $Tflag && [[ -d $1 ]]; then echo $Terr; exit 1
elif ! $iflag && [[ -d $1 ]]; then echo $ierr; exit 1
elif ! $rflag && [[ -d $1 ]]; then echo $rerr; exit 1
elif ! $tflag && [[ -d $1 ]]; then echo $terr; exit 1
elif ! $gflag && [[ -d $1 ]]; then echo $gerr; exit 1
elif ! $eflag && [[ -d $1 ]]; then echo $eerr; exit 1
fi

if [ $MORET ] && [ ! "$(ls ./${TAXADB}/${MORET}  2>/dev/null)" ]; then echo $merr; show_help; exit; fi;
merr2="Directory of previous alignments not found: ./OGfastasALN\ Try running without -m " 
if [ $MORET ] && [ ! -d ./OGfastasALN  ]; then echo $merr2 ; exit; fi;

echo "" "Running...." ""
echo "Reference species: $REFSP"
echo "Orignal gene datasets in: $OGDIR"
echo "New taxa datasets to search in: $TAXDB"
echo "Reference species genome in: $REFDB"
echo "Number of processors alloted: $maxCPU"


mkdir ./Seed_OGs
mkdir ./blast1_results
mkdir ./blast2_results
mkdir ./hit1_fasta
mkdir ./ref_OGs
mkdir ./OGfastasets
if [ ! $MORET ]; then mkdir ./OGfastasALN 2>/dev/null ; fi

###     Stage 1 : Pull OrthoGroup sequences for the Reference species 
###     and use them to blast the reference genome

#Grab REFSP seqs from each the Orthology Group, renumbers, tags species and write to file: REFSP_OGs_renamed.fa
if [ -s ./Seed_OGs/${REFSP}_OGs_renamed.fa ] #if file exists and is > 0 bytes
	then echo "OGs from $REFSP ready"
else
		if [ "$(ls -A $OGDIR/*.nex 2>/dev/null)" ]
		then
		#convert single nexus to multiple fasta
		orig="$(ls *.nex)"
		./seqConverter.pl -s -in -of -d*.nex;  rm ${orig%.nex}.fasta; rename 's/fasta/fa/' *.fasta
		fi
	for i in $(ls $OGDIR/*)
	do
	  echo "searching  $i ..."
	  selectSeqs.pl -m $REFSP  $i >> ${REFSP}_OGs.fa
	done
	echo "renaming FASTA headers" 			#Number each $REFSP OG seq so that there are unique identifiers
	#puts leading zeroes in numbering scheme
	awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < ${REFSP}_OGs.fa > ${REFSP}_OGs_renamed.fa
	rm *_OGs.fa
	#perl -p -i -e 's/>/>Nematostella|/g' ${REFSP}_OGs_renamed.fa
	#mv *renamed.fa ./Seed_OGs 
	awk -v new=">${REFSP}\|" '{sub(/>/, new)}1' ${REFSP}_OGs_renamed.fa > tmp 2>/dev/null
	mv tmp ./Seed_OGs/${REFSP}_OGs_renamed.fa
	rm ${REFSP}_OGs_renamed.fa
fi

#BLAST REFSP OGs against whole REFSP genome to get consistent gene names
if [ "$(ls -A ./ref_OGs/*ref_blastout 2>/dev/null)" ]  #(if not empty:)
then echo "$REFSP OGs BLAST against own genome: complete"
else
	for fasta in ./${REFDB}/${REFSP}*.fa
	
	do
	echo $fasta
	  if [ ! -s ${fasta}.nhr ]; then
		  echo "Formating for BLAST  $fasta ..." #Format ref_genome CDS datasets for BLAST
		  makeblastdb -in $fasta -parse_seqids -dbtype nucl
	  fi
	  echo "BLASTING  ./Seed_OGs/${REFSP}_OGs_renamed.fa against own genome:  $fasta ..."
	  blastn -query ./Seed_OGs/${REFSP}_OGs_renamed.fa -db $fasta -out ${fasta%.}_ref_blastout -evalue $EVAL -outfmt 6 # -max_target_seqs 1

	  #Remove redundant accessions
	  cat ${fasta%.}_ref_blastout| cut -f1 |uniq|  while read line ; do grep -m 1 $line ${fasta%.}_ref_blastout >> ${fasta%.}_ref_blastout_2; done 

 	  #Pull fasta seqs for each hit
      cut -f2 ${fasta%.}_ref_blastout_2 > file
 	  selectSeqs.pl -f file $fasta > ref_cds_genome_seq.fa
 	  rm file 
	  mv ./${REFDB}/*ref_blastout* ./ref_OGs
	  mv  ref_cds_genome_seq.fa ./hit1_fasta/${REFSP}.fa
	done
fi

###     Stage 2 : Blast(1) target species using OrthoGroup sequences for the Reference species; 
###     perform CDHIT on targets prior to BLAST to minimize paralogy/isoform noise during reciprocity step later

#CDHIT cds datasets before BLAST query
if [ ! $CDHIT ]; then echo "Skipping CDHIT step...";
else
	echo "Running CD-HIT..."
	if [ $MORET ]; then rename 's/\.fa/\.fas/' ./${TAXDB}/${MORET}
	cdhit -i ./${TAXDB}/${MORET%.fa}.fas -o ./${TAXDB}/${MORET%.fa}_cdht98.fa -c 0.98 -n 5
	NEWTX=./${TAXDB}/${MORET%.fa}_cdht98.fa
	else
	rename 's/\.fa/\.fas/' ./${TAXDB}/*.fa # change orig ext. to fas so it's ignored by later calls
	parallel --jobs $maxCPU 'cdhit -i {} -o {.}_cdht98.fa -c 0.98 -n 5' ::: ./${TAXDB}/*.fas	
	fi
fi
if [ ! $CDHIT ] && [ $MORET ]; then NEWTX=./${TAXDB}/${MORET}; fi

#Format nucleotide datasets for BLAST
if [ $MORET ]; then makeblastdb -in NEWTX -parse_seqids -dbtype nucl; fi

if [ "$(ls -A ./${TAXDB}/*phr  2>/dev/null)" ]; 
then echo "BLAST databases for new taxa found"
else
	echo "Formatting BLAST databases for new taxa..."
	parallel --jobs $maxCPU 'makeblastdb -in {} -parse_seqids -dbtype nucl' ::: ./${TAXDB}/*.fa
fi

#BLAST $REFSP OGs against each formated CDS dataset
if [ $MORET ]; then blastn -db NEWTX -query ./hit1_fasta/"${REFSP}".fa -out ${NEWTX}_blastout1 -evalue $EVAL -outfmt 6 # -max_target_seqs 1
	mv ./${TAXDB}/*_blastout1 ./blast1_results; 
fi

if [ "$(ls -A ./blast1_results)" ]; 
then echo "BLAST1 complete"
else
	echo "BLAST1: Querying $REFSP against new taxa..."
#	parallel --jobs $maxCPU 'blastn -db {} -query' ./hit1_fasta/"${REFSP}".fa '-out {.}_blastout1 -evalue '"$EVAL"' -outfmt 6 -max_target_seqs 1' ::: ./${TAXDB}/*.fa
	parallel --jobs $maxCPU 'blastn -db {} -query' ./hit1_fasta/"${REFSP}".fa '-out {.}_blastout1 -evalue '"$EVAL"' -outfmt 6 ' ::: ./${TAXDB}/*.fa
	mv ./${TAXDB}/*_blastout1 ./blast1_results
fi


#Pull fasta seqs for each hit to use in 2nd round of blasts
function pullseqs () {
	cat $1| cut -f1 |uniq|  while read line ; do grep -m 1 $line $hitfile >> ${1}_2; done 
	cut -f2 ${1}_2 |perl -p -e 's/\>//g' > temp
	local filenamefull=${1##*/}
	local filename="${filenamefull%_blastout1}.fa"
	echo "Pulling FASTA seqs from  $filename ..."
	selectSeqs.pl -f temp ./${TAXDB}/$filename >> ./hit1_fasta/${filename}
	rm temp
	}

if [ $MORET ]; then pullseqs ./blast1_results/${NEWTX}_blastout1; 
	else
	for hitfile in ./blast1_results/*blastout1
	do
	pullseqs $hitfile
	done
fi
echo "Blast1 hits pulled"


###    Stage 3:  Check reciprocity in Blast matches and make FASTAS:
###    -BLAST(2, reciprocal) query best hits from BLAST2 back to Reference species genome
###    -create fasta for each OrthoGroup containing seqs that satisfy BestReciprocalHit criterion

#BLAST REFSP OGs against best hits from each formated nucleotide dataset's 
echo "Running reciprocal BLASTs..."; echo ""
refgen=$(echo ./"${REFDB}"/${REFSP}*.fa)
if [ $MORET ]; then blastn -query ./hit1_fasta/${NEWTX}.fa -evalue $EVAL -db $refgen -out {1.}_blastout2 -outfmt 6 # -max_target_seqs 1; 
else 
#parallel --jobs $maxCPU 'blastp -query {1} -db {2} -out {1.}_blastout2 -evalue '"$EVAL"' -outfmt 6 -max_target_seqs 1' ::: ./hit1_fasta/*.fa ::: $refgen
parallel --jobs $maxCPU 'blastn -query {1} -db {2} -out {1.}_blastout2 -evalue '"$EVAL"' -outfmt 6' ::: ./hit1_fasta/*.fa ::: $refgen
fi
mv ./hit1_fasta/*_blastout2 ./blast2_results

function parsepull {
	#Remove redundant accessions
   cat $1| cut -f1 |uniq|  while read line ; do grep -m 1 $line $1 >> ${1}_2; done
   local basepath=${1%_blastout2}
   local base=${basepath##*/}
   # determine which seqfile to grab from using taxon in filename
   if [ -s ./${TAXDB}/${base}.fa ]; then seqfile=./${TAXDB}/${base}.fa
     ./parse_recipBLAST.py ${1}_2 ./blast1_results/${base}_blastout1_2 1 2 12 high rbh
     tail -n+2 rbh > rbh2
     ./pullOGsfromBlast.py $seqfile rbh2 2 1; rm rbh*
   elif [ -s ./${REFDB}/${base}.fa ]; then seqfile=./${REFDB}/${base}.fa
     ./pullOGsfromBlast.py $seqfile ./ref_OGs/${base}.fa_ref_blastout_2 2 2
   fi
  echo "Pulling Best Reciprocal Blast Hits from " $seqfile
}

if [ $MORET ]; then parsepull ./blast2_results/${NEWTX%.fa}_blastout2
else
for blastout in ./blast2_results/*blastout2
do
  parsepull $blastout
done; 
fi

#if multiple anchors
#look in first anchor's fasta-set
# grep deflines for all anchors
# search next anchor's fasta-set for co-occuring deflines: if found keep fasta, if not discard.
#
#only align kept fasta

##     Stage 4: align each OG
if [ $MORET ]; then parallel --jobs $maxCPU 'mafft --auto {} > {.}.aln' ::: *fa
else
mv *.fa ./OGfastasets/
parallel --jobs $maxCPU 'mafft --auto {} > {.}.aln' ::: ./OGfastasets/*fa
fi
mv ./OGfastasets/*aln ./OGfastasALN/
echo "Mafft alignment complete"

##restore original taxa fastanames [only applies if CD-HIT flagged]
if [ $CDHIT ]; then
	rename 's/\.fas/\.fa/' ./${TAXDB}/*.fa
fi

## Tidy up [optional]
if [ $DELY ]; then 
rm ./${TAXDB}/*.{phr,pin,pog,psd,psi,psq}
rm ./${REFDB}/*.{phr,pin,pog,psd,psi,psq}
rm -rf ./Seed_OGs
rm -rf ./blast1_results
rm -rf ./blast2_results
rm -rf ./hit1_fasta
rm -rf ./ref_OGs
rm -rf ./OGfastasets
rm -rf ./hit1_fasta
fi
