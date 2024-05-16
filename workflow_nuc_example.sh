#!/bin/bash

#SBATCH --partition=shared
#SBATCH -J orthoSetUP
#SBATCH --output setup_blast.log
#SBATCH --exclude=node117,node118
#SBATCH --open-mode=append

module load linuxbrew/colsa

#prep for orthomatic
mkdir ref_genome
mkdir gene_sets
mkdir samples

refsp_cds="path/to/Aoroides_cds.fasta"
protdata="/mnt/lz01/plachetzki/shared/metazoa2023/Alien_Indexed_maxine/Amphimedon_queenslandica.faa"

##REF_GENOME
#cp genome to ref_genome dir
cp $refsp_cds ref_genome/

# create duplicate fasta, cleaning up headers if needed
sed 's/ .*$//' ref_genome/${refsp_cds} > ref_genome/${refsp_cds%.fasta}.fa
#fasta file suffix is .fa

#example format of refsp_cds
# >g43249.t1
# atgagtgttttgacaaatgtgttaatgaaactacctgttttccctatatgtgaggaggagtgtgtttggt
# tgagctgtgctgaggttggtcagatgagccagcagccgtcctcctcgatgtgtacgagtgttgttcggag
# cagacgtagtgctgccgagtgcattactggggtgctgcgcttgaagggagcagaagcattgcccggctgg
# ctcggtagcagcaagaagacaattctggtggcccaaccgttgtaccagttgccggccactggtcttggaa


# GENE SETS
# 1: identify matches between focal species CDS and protein data set

blastx -query ref_genome/${refsp_cds%.fasta}.fa  -subject $protdata -evalue 1e-50 -outfmt 6 -num_threads 24 -out blast1
tblastn -query $protdata -subject ref_genome/${refsp_cds%.fasta}.fa -evalue 1e-50 -outfmt 6 -num_threads 24 -out blast2

./parse_recip_blast.py blast1 blast2 1 2 12 high rbh
tail -n+2 rbh |cut -f1 |sort|uniq >rbh2

selectSeqs.pl -f rbh2 ref_genome/${refsp_cds%.fasta}.fa> refsp_hits.fas


mv refsp_hits.fas gene_sets/
cd gene_sets/
#splits up multi-seq fasta to many single-seq fastas
cat gene_sets/refsp_hits.fas | awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")}
        print $0 > filename
}'

#Adds species to header line
#replace Aoroides with your ref species
for fa in refsp_hits/*.fa ; do sed 's/>g.*$/>Aoroides/' $fa > tmp; mv tmp $fa ; done


##SAMPLES
#mv cds fastas (.fa) to samples/


#fix headers for orthomatic
#adds sample name and leading zeros to contig fastas
 for fa in samples/*.fa ; do f=${fa%%.*}; f1=$(basename $f); echo $f1;
	awk '/^>/{printf ">" "%05d\n", ++i; next}{print}' < $fa > ${fa}2;
	awk -v new=">${f1}_" '{sub(/>/, new)}1'  ${fa}2 > tmp 2>/dev/null;
	mv tmp samples/${f1}.fa;
done



# run orthomatic
./orthomatic_nuc.sh -T 24 -i gene_sets -r Aoroides -t samples -e 1e-50 -g ref_genome

#post orthomatic
# if a fasta from reference species was not included in samples/ and is needed in final alignments:

seqfile="../Seed_OGs/*OGs_renamed.fa"
blastout="../ref_OGs/Aoroides_cds.fa_ref_blastout_2"

../pullOGsfromBlast.py $seqfile $blastout 2 1
parallel --jobs 24 'mafft --auto {} > {.}.aln' ::: OGfastasets/*.fa
mv OGfastasets/*.aln OGfastasALN/

#prep alignments for tree inference
mkdir aligned_ready/

for file in OGfastasALN/*aln
do
if [ -s $file ];
then
echo "$file" >> partitionfilelist
cat $file | sed -e 's/|/_/' | awk -F_ '{print $1}' > ${file%.aln}.fas
mv ${file%.aln}.fas aligned_ready/
 fi;
done

#filter alignments using miniumum number of sequences
MINSEQ=83; grep -c ">" aligned_ready/g*.fas | awk -v minseq="$MINSEQ" -F: '$2 > minseq' | cut -f1 -d':' > fastalist_${MINSEQ}

mkdir iq_aln_${MINSEQ}
while read line; do cp $line iq_aln_${MINSEQ}/ ; done < fastalist_${MINSEQ}
#create table for partion occupancy
while read line; do base=${line##*/}; grep '>' $line | cut -f3 -d"/" |sed -e "s|>|$base,|g"   >> rinput; done < fastalist_${MINSEQ}
R --slave -e 'dat<-read.csv("rinput",header=F);library(reshape2); dat$V3<-rep(1,nrow(dat));dat2<-recast(dat,V1~V2);write.csv(dat2,"partition_table.csv")'
####
# run IQ-tree on dir with alignements (iqtree-2.2.6)
module purge
module load anaconda/colsa
conda activate iqtree-2.2.6
iqtree -s iq_aln_${MINSEQ}/ -m MFP -B 1000  -T 24


