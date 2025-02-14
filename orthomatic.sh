#!/bin/bash
set -e

#version4: 
#omitting REFDB

# future dev: 
#   accept file containing ref species names instead of arg1; 
#   multiple anchor species.   
#add cd-hit-est

# Usage info
show_help() {
cat << EOF
Example: ${0##*/} -T 8 -i ortholog_dir/ -s prot -r Amphimedon -t newdata_fasta_dir/ -e 1e-20 

Usage: ${0##*/} [-h] [-c] [-m] [-T num] [-s seqtype] [-i OG_dir] [-r reference_species] [-t taxa_fastas_dir] [-e num] [-k]

	-h <help>	display this help and exit
	
	-c CDHIT	option to run CDHIT on FASTAs specified by -t or -m
	
	-T MAXCPU	number of threads
	
	-s SEQTYPE  prot or nucl (protein or nucleotide input)

	-i OGDIR	directory with either: 
				one fasta per gene (.fa)
				-OR-
				single nexus file with concatenated data; gene partitions specified in charset block (.nex)
				
				Reference Species must be included, with REFSP used in sequence id.
				
	-r REFSP	Reference Species handle; no spaces, name must be found in sequence ids of OGDIR/ files
				species should have good genomic data that will be used in reciprocal BLASTs with all other taxa
	
	-t TAXDB	directory with peptide fasta files for each taxon to search
				reference species file must include REFSP in file name
				NO OTHER FILES should include REFSP in file name
				all fasta must end with .fa
				
				
	-m MORET    add new taxon dataset to existing results (put FASTA in TAXDB directory)
                **Requires directory of previous alignments called orthomatic_alignments**
                 
	-e EVAL		e-value cut-off for blast
	
	-k KEEP		option to keep intermediate files (for troubleshooting, discarded by default) 

No spaces in any FASTA headers. 
Fasta file extension uses ".fa"
Nexus file extension uses ".nex"

See fasta example: http://datadryad.org/bitstream/handle/10255/dryad.98320/OGs.zip?sequence=1

#################################################################
#
#  Orthomatic BLASTs fastas (e.g., assembled & translated 
#  transcriptomes) using queries from designated reference taxon 
#  for a provided set of orthologous peptides.  Matches are 
#  reciprocally BLASTed to the reference taxon\'s genome. Sequences 
#  satisfying the criterion of mutual-best-match are then pooled 
#  and aligned by MAFFT.
#
#################################################################
# 
# Dependencies (must be in $PATH): 								
#   NCBI BLAST 2.2.31+ (specifically: blastp, makeblastdb)		
#   GNU parallel																
#   CD-HIT (cd-hit), if -c invoked														
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
maxCPU=""
EVAL=""
CDHIT=""
MORET=""
KEEP="yes"
SEQTYPE=""

OPTIND=1
Tflag=false; iflag=false; rflag=false; tflag=false;eflag=false; sflag=false;
Terr="Missing: -T <number of threads>"
serr="Missing: -s <prot or nucl>"
ierr="Missing: -i <directory with query gene fastas>"
rerr="Missing: -r <reference_species_name>"
terr="Missing: -t <directory with taxon fastas>"
eerr="Missing: -e <e-value, e.g., 1e-20>"
merr="Missing: -m <name of new taxa FASTA file>"

while getopts hT:i:s:r:t:e:m:cD opt; do
    case $opt in
        h)  show_help; exit 0;;
        T)  maxCPU=$OPTARG; if [ ! $OPTARG ]; then echo $Terr; echo "boo"; show_help; exit; fi;;
	    s)  SEQTYPE=$OPTARG; if [ ! $OPTARG ]; then echo $serr; show_help; exit; fi;; 
        i)  OGDIR=$OPTARG; if [ ! -d $OPTARG ]; then echo $ierr; show_help; exit; fi;;
        r)  REFSP=$OPTARG; if [ ! $OPTARG ]; then echo $rerr; show_help; exit; fi;;
        t)  TAXDB=$OPTARG; if [ ! -d $OPTARG ]; then echo $terr; show_help; exit; fi;;
        e)  EVAL=$OPTARG;  if [ ! $OPTARG ]; then echo $eerr; show_help; exit; fi;;
        m)  MORET=$OPTARG ;;
        c)  CDHIT=$OPTARG; echo "Input taxa fastas will be consolidated with CD-HIT" ;;
        k)  KEEP=$OPTARG; echo "Clean-up supressed: Temporary files will be retained"; KEEP="" ;;
        \?) echo "Unknown option: -$OPTARG" >&2; show_help; exit 1;;
        :) echo "Missing option argument for -$OPTARG" >&2; show_help; exit 1;;
        *) echo "Unimplemented option: -$OPTARG" >&2;show_help; exit 1;;
    esac
done
shift "$((OPTIND-1))" # Shift off the options and optional --.

if [[ -z $OGDIR || -z $TAXDB || -z $REFSP || -z $EVAL || -z $maxCPU ]]; then show_help; exit 1
elif ! $Tflag && [[ -d $1 ]]; then echo $Terr; exit 1
elif ! $iflag && [[ -d $1 ]]; then echo $ierr; exit 1
elif ! $rflag && [[ -d $1 ]]; then echo $rerr; exit 1
elif ! $tflag && [[ -d $1 ]]; then echo $terr; exit 1
elif ! $eflag && [[ -d $1 ]]; then echo $eerr; exit 1
elif ! $sflag && [[ -d $1 ]]; then echo $serr; exit 1
fi

if [ $MORET ] && [ ! "$(ls ./${TAXADB}/${MORET}*  2>/dev/null)" ]; then echo $merr; show_help; exit; fi;
merr2="Directory of previous alignments not found: ./orthomatic_alignments\ Try running without -m " 
if [ $MORET ] && [ ! -d ./orthomatic_alignments  ]; then echo $merr2 ; exit; fi;

echo "" "Running...." ""
echo "Reference species: $REFSP"
echo "Orignal gene datasets in: $OGDIR"
echo "New taxa datasets to search in: $TAXDB"
echo "Sequence format is $SEQTYPE"
echo "Number of processors alloted: $maxCPU"



###     Stage 1 : Pull OrthoGroup sequences for the Reference species 
###     and use them to blast the reference genome

#Grab REFSP seqs from each the Orthology Group, renumbers, tags species and write to file: REFSP_OGs_renamed.fa
if [ -s ./Seed_OGs/${REFSP}_OGs_renamed.fa ] #if file exists and is > 0 bytes
	then echo "OGs from $REFSP ready"; echo ""
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
        mkdir ./Seed_OGs/
        #perl -p -i -e 's/>/>Nematostella|/g' ${REFSP}_OGs_renamed.fa
	#mv *renamed.fa ./Seed_OGs 
	awk -v new=">${REFSP}\|" '{sub(/>/, new)}1' ${REFSP}_OGs_renamed.fa > tmp 2>/dev/null
	mv tmp ./Seed_OGs/${REFSP}_OGs_renamed.fa
	rm ${REFSP}_OGs_renamed.fa
    perl -pi -e 's/-//g' ./Seed_OGs/${REFSP}_OGs_renamed.fa #remove any gaps if using OGs from alignments
fi

#BLAST REFSP OGs against whole REFSP genome to get consistent gene names
if [ "$(ls -A ./ref_OGs/*ref_blastout 2>/dev/null)" ]  #(if not empty:)
then echo "$REFSP OGs BLAST against own genome: complete" ; echo ""
else
mkdir ./ref_OGs
mkdir ./hit1_fasta
	for fasta in ./${TAXDB}/${REFSP}*.fa
	
	do
	if [[ $SEQTYPE == "prot" ]]; then
	  if [ ! -s ${fasta}.phr ] ; then
		  echo "Formating for BLAST  $fasta ..." #Format ref_genome peptide datasets for BLAST
		  makeblastdb -in $fasta -parse_seqids -dbtype prot  #TO checkif exists first
	  fi
	  echo "BLASTING  ./Seed_OGs/${REFSP}_OGs_renamed.fa against own genome:  $fasta ..."
	  blastp -query ./Seed_OGs/${REFSP}_OGs_renamed.fa -db $fasta -out ${fasta%.}_ref_blastout -evalue $EVAL -outfmt 6 -max_target_seqs 10
    else
    if [ ! -s ${fasta}.nhr ]; then
		  echo "Formating for BLAST  $fasta ..." #Format ref_genome CDS datasets for BLAST
		  makeblastdb -in $fasta -parse_seqids -dbtype nucl
	  fi
	  echo "BLASTING  ./Seed_OGs/${REFSP}_OGs_renamed.fa against own genome:  $fasta ..."
	  blastn -query ./Seed_OGs/${REFSP}_OGs_renamed.fa -db $fasta -out ${fasta%.}_ref_blastout -evalue $EVAL -outfmt 6 # -max_target_seqs 1

    fi
    
	#Remove redundant accessions and retain top match
	  cat ${fasta%.}_ref_blastout |sort -k1,1 -k12,12nr > sorted
	cat sorted| cut -f1 |uniq|  while read line ; do grep -m 1 $line sorted >> ${fasta%.}_ref_blastout_2; done 
	rm sorted


 	  #Pull fasta seqs for each hit
         cut -f2 ${fasta%.}_ref_blastout_2 > file
 	  selectSeqs.pl -f file $fasta > ref_pep_genome_seq.fa
 	  rm file 
	  mv ./${TAXDB}/*ref_blastout* ./ref_OGs
	  mv  ref_pep_genome_seq.fa ./hit1_fasta/${REFSP}.fa
	done
fi

###     Stage 2 : Blast(1) target species using OrthoGroup sequences for the Reference species; 
###     perform CDHIT on targets prior to BLAST to minimize paralogy/isoform noise during reciprocity step later

#CDHIT peptide datasets before BLAST query
if [ ! $CDHIT ]; then echo "Skipping CDHIT step...";
else
	echo "Running CD-HIT..."
	if [ $MORET ]; then rename 's/\.fa/\.fas/' ./${TAXDB}/${MORET}
	cdhit -i ./${TAXDB}/${MORET%.fa}.fas -o ./${TAXDB}/${MORET%.fa}_cdht98.fa -c 0.98 -n 5
	NEWTX=./${TAXDB}/${MORET%.fa}_cdht98.fa
	else
	rename 's/\.fa/\.fas/' ./${TAXDB}/*.fa # change orig ext. to fas so it's ignored by later calls
	parallel --jobs $maxCPU 'cd-hit -i {} -o {.}_cdht98.fa -c 0.98 -n 5' ::: ./${TAXDB}/*.fas	
	wait
	fi
fi
if [ ! $CDHIT ] && [ $MORET ]; then NEWTX=./${TAXDB}/${MORET}; fi

#Format peptide datasets for BLAST

if [ $MORET ]; then makeblastdb -in NEWTX -parse_seqids -dbtype $SEQTYPE; fi

	echo "Formatting BLAST databases for new taxa..."
	parallel --jobs $maxCPU 'makeblastdb -in {} -parse_seqids -dbtype ' $SEQTYPE ::: ./${TAXDB}/*.fa
	wait

#BLAST $REFSP OGs against each formated peptide dataset
if [ $MORET ] && [[ $SEQTYPE == "prot" ]]; then blastp -db NEWTX -query ./hit1_fasta/"${REFSP}".fa -out ${NEWTX}_blastout1 -evalue $EVAL -outfmt 6 -max_target_seqs 10;fi
if [ $MORET ] && [[ $SEQTYPE == "nucl" ]]; then blastn -db NEWTX -query ./hit1_fasta/"${REFSP}".fa -out ${NEWTX}_blastout1 -evalue $EVAL -outfmt 6 ;fi


if [ "$(ls -A ./blast1_results  2>/dev/null))" ]; 
	then echo "BLAST1 complete"
else
        mkdir ./blast1_results

	echo "BLAST1: Querying $REFSP against new taxa..."
	if [[ $SEQTYPE == "prot" ]]; then
	parallel --jobs $maxCPU 'blastp -db {} -query' ./hit1_fasta/"${REFSP}".fa '-out {.}_blastout1 -evalue '"$EVAL"' -outfmt 6 -max_target_seqs 10' ::: ./${TAXDB}/*.fa
	fi
	if [[ $SEQTYPE == "nucl" ]]; then 
	parallel --jobs $maxCPU 'blastn -db {} -query' ./hit1_fasta/"${REFSP}".fa '-out {.}_blastout1 -evalue '"$EVAL"' -outfmt 6 ' ::: ./${TAXDB}/*.fa
	wait
	fi

mv ./${TAXDB}/*_blastout1 ./blast1_results
fi

#Pull fasta seqs for each hit to use in 2nd round of blasts
function pullseqs () {
#	cat $1| cut -f1 |uniq|  while read line ; do grep -m 1 $line $hitfile >> ${1}_2; done 
        cat $1 |sort -k1,1 -k12,12nr > ${1}.sorted
	cat ${1}.sorted| cut -f1 |uniq|  while read line ; do grep -m 1 $line ${1}.sorted >> ${1}_2; done 
	rm ${1}.sorted
	cut -f2 ${1}_2 |perl -p -e 's/\>//g' > ${1}.temp
	local filenamefull=${1##*/}
	local filename="${filenamefull%_blastout1}.fa"
#echo $filename
#echo ${2}/${filename}
	echo "Pulling FASTA seqs from  $filename ..."
	selectSeqs.pl -f ${1}.temp ./${2}/${filename} >> ./hit1_fasta/${filename}
	rm ${1}.temp
	}
export -f pullseqs

if [ "$(ls -A ./hit1_fasta 2>/dev/null)" ] ; then
     echo "Sequences ready for Blast1 hits in ./hit1_fasta"
else
	if [ $MORET ]; then pullseqs ./blast1_results/${NEWTX}_blastout1 $TAXDB; 
	else
	#for hitfile in ./blast1_results/*blastout1
#	do
	#pullseqs $hitfile
#	done
	parallel --jobs $maxCPU pullseqs  ::: ./blast1_results/*blastout1 ::: $TAXDB
    wait
   fi
fi
echo "Sequences pulled for Blast1 hits into ./hit1_fasta"


###    Stage 3:  Check reciprocity in Blast matches and make FASTAS:
###    -BLAST(2, reciprocal) query best hits from BLAST2 back to Reference species genome
###    -create fasta for each OrthoGroup containing seqs that satisfy BestReciprocalHit criterion

#BLAST REFSP OGs against best hits from each formated peptide dataset's 
if [ "$(ls -A ./blast2_results  2>/dev/null))" ]; 
        then echo "BLAST2 complete"; echo ""
else
echo "Running reciprocal BLASTs..."; echo ""
refgen=$(echo ./"${TAXDB}"/${REFSP}*.fa)
if [[ $SEQTYPE == "prot" ]]; then
	if [ $MORET ] ;  then blastp -query ./hit1_fasta/${NEWTX}.fa -evalue $EVAL -db $refgen -out {1.}_blastout2 -outfmt 6 -max_target_seqs 10; 
	else 
	parallel --jobs $maxCPU 'blastp -query {1} -db {2} -out {1.}_blastout2 -evalue '"$EVAL"' -outfmt 6 -max_target_seqs 10' ::: ./hit1_fasta/*.fa ::: $refgen
	wait
	fi
else

	if [ $MORET ]; then blastn -query ./hit1_fasta/${NEWTX}.fa -evalue $EVAL -db $refgen -out {1.}_blastout2 -outfmt 6 # -max_target_seqs 1; 
	else 
	parallel --jobs $maxCPU 'blastn -query {1} -db {2} -out {1.}_blastout2 -evalue '"$EVAL"' -outfmt 6' ::: ./hit1_fasta/*.fa ::: $refgen
	wait
	fi
fi
mkdir ./blast2_results
mv ./hit1_fasta/*_blastout2 ./blast2_results
fi

function parsepull {
	#Remove redundant accessions
   #cat $1| cut -f1 |uniq|  while read line ; do grep -m 1 $line $1 >> ${1}_2; done
   cat $1 |sort -k1,1 -k12,12nr > ${1}.sorted
	cat ${1}.sorted| cut -f1 |uniq|  while read line ; do grep -m 1 $line ${1}.sorted >> ${1}_2; done 
	rm ${1}.sorted   

local basepath=${1%_blastout2}
   local base=${basepath##*/}
   # determine which seqfile to grab from using taxon in filename
   # ref_OG becomes file prefix
   # blast hit sequences are pulled and appended to the same ref_OG fasta
   if [ -s ./${2}/${base}.fa ]; then seqfile=./${2}/${base}.fa
     ./parse_recipBLAST.py ${1}_2 ./blast1_results/${base}_blastout1_2 1 2 12 high ${1}.rbh
     tail -n+2 ${1}.rbh > ${1}.rbh2
     ./pullOGsfromBlast.py $seqfile ${1}.rbh2 2 1; rm ${1}.rbh*
 #  elif [ -s ./${3}/${base}.fa ]; then seqfile=./${3}/${base}.fa
 #    ./pullOGsfromBlast.py $seqfile ./ref_OGs/${base}.fa_ref_blastout_2 2 2
   fi
  echo "Pulling Best Reciprocal Blast Hits from " $seqfile
}
export -f parsepull

##     Stage 4: establish which blast hits are repriprocated, pull seqs to indiv OG fasta, align each OG
if [ -d OGfastasets ] && [ "$(ls -A ./OGfastasets  2>/dev/null))" ]; 
        then echo "Reciprocated hits grouped and aligned by OG"; echo ""
else
	echo "Checking reciprocity of best hits and grabbing those sequences..."; echo ""
	if [ $MORET ]; then 
		parsepull ./blast2_results/${NEWTX%.fa}_blastout2 $TAXDB #$REFDB

	else
	#for blastout in ./blast2_results/*blastout2
	## parsepull $blastout
	#done; 
	parallel --jobs $maxCPU parsepull  ::: ./blast2_results/*blastout2 ::: $TAXDB #::: $REFDB
	wait
	fi

    echo "Running MAFFT on each OG"; echo ""
    mkdir OGfastasets
    mv *.fa ./OGfastasets/
    parallel --jobs $maxCPU 'mafft --auto {} > {.}.aln' ::: ./OGfastasets/*fa
    wait
fi
mkdir orthomatic_alignments
mv ./OGfastasets/*aln ./orthomatic_alignments/
echo "Mafft alignment complete"

##restore original taxa fastanames [only applies if CD-HIT flagged]
if [ $CDHIT ]; then
	rename 's/\.fas/\.fa/' ./${TAXDB}/*.fa
fi

## Tidy up [optional]
if [ $KEEP ]; then 
rm ./${TAXDB}/*.{phr,pin,pog,psd,psi,psq}
rm -rf ./Seed_OGs
rm -rf ./blast1_results
rm -rf ./blast2_results
rm -rf ./hit1_fasta
rm -rf ./ref_OGs
rm -rf ./OGfastasets
rm -rf ./hit1_fasta
fi
