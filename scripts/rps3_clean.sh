#Install my working conda environments :) Hopefully these work for you too!
conda create -n hmmer --file hmmer.txt
conda create -n megahit --file megahit.txt
##!! OPTIONAL
#This is just where I keep bowtie2 and samtools. It's not project specific.
#If you already have an environment that has this, feel free to use your own.
#I don't think I'm doing anything version specific in this code for bowtie or samtools.
conda create -n map_software --file map_software.txt

conda activate megahit
#run megahit to get contigs from your files. This is the part that will probably require the most tweaking
#here, I had unmerged forward and reverse reads. Your command might look different
#../coral_forward_files is the path to a list of forward fastqs with one fastq per line
#../coral_reverse_files is the path to a list of reverse fastqs with one fastq per line, matching the order of the forward files.
{ echo "megahit --mem-flag 2 -1 "; cat ../coral_forward_files; echo " -2 "; cat ../coral_reverse_files; echo " -o coralMAG --out-prefix coral"; } > coralmag
#you might not need 12 hours... I think the job actually ran in about 2 on ~300 million reads. Can't remember though.
#And this is the most computationally limiting step in the whole pipeline.
python3 $HOME/bin/ls6_launcher_creator.py -j coralmag -n coralmag -a IBN21018 -e cbscott@utexas.edu -t 12:00:00 -N 1 -w 1 -q normal

#Note! I had to do this for all of my substrates separately, and then merged the contig output files into one file
#called "allcontigs.fa" - this is what I'm using in the rest of the program.


#cd into where your contigs are.
cd /scratch/06909/cbscott/STX_SpawnMeta/rps3
#I have a file here hat is allcontigs.fa - the result of
conda activate hmmer
#clean up your contig names, if you've combined across sources,
#there are probably duplicate names. Either way, it doesn't hurt to make sure they're:
#a- machine readable for later steps in the pipeline where I grep, and,
#b- unique
reformat.sh in=allcontigs.fa out=allcontigs.trimws.fa trd
reformat.sh in=allcontigs.trimws.fa out=allcontigs.clean.fa uniquenames


#to run the hmm model, we want to translate our DNA sequences to proten sequences.
#Use prodigal for this:
#need bam files aligned to our all.contigs.fa?
echo "prodigal -p meta -i allcontigs.clean.fa -a allcontigs.proteinonly" >testprod
python3 $HOME/bin/ls6_launcher_creator.py -j testprod -n testprod -a IBN21018 -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 1 -q development


#Use a hidden markov model developed by someone else (Diamond 2019) to quickly
#predict which contigs contian rps3 sequences

#1. Get the hmm from Diamond:
#this will clone the github repo into your current directory.
git clone https://github.com/SDmetagenomics/Angelo2019_Paper.git
cp Angelo2019_Paper/HMMs/rps3/rpS3_Diamond2019.hmm . #copy the hmm into your current directory for convenvience

#2. Run hmmer to search your protein sequences! This was super speedy for me.
echo "hmmsearch --cpu 24 --pfamtblout allcontigs.pfam.rps3.tbl --tblout allcontigs.rps3.tbl rpS3_Diamond2019.hmm allcontigs.proteinonly" > searchtest
python3 $HOME/bin/ls6_launcher_creator.py -j searchtest -n searchtest -a IBN21018 -e cbscott@utexas.edu -t 00:10:00 -N 1 -w 1 -q development

#Great! You now know which contig names contain rps3 sequences.
#we now need to extract these contig ids from our nucleotide fasta and protein fasta files
#########GET RPS3 PROT SEQUENCES
awk '{print $1}' allcontigs.rps3.tbl > rps3.seqids
#!!!!!!!IF YOUR CONTIGS DON'T START WITH K THIS WILL NOT WORK.
#This is the DEFAULT behavior of megahit, so it should be fine. double check your output to be sure.
grep "^k" rps3.seqids > rps3.seqids.clean
seqtk subseq allcontigs.proteinonly rps3.seqids.clean > contig.rps3.hits.prot.fa
seqtk subseq allcontigs.clean.fa rps3.seqids.clean  > contig.rps3.hits.dna.fa


#Ok. Now we have contigs that we think are rps3. Since we're going to use these
#to calculate species abundance, we're going to cluster them by sequence similarity.
#We only want to keep rps3 sequences that ate *dissimilar* to other sequences,
#because they probably correspond to one species group.
#To do this, use usearch. Usearch has to be downloaded from source.

#1. Go to this website: https://www.drive5.com/usearch/download.html
#2. Right click on the link usearch11.0.667_i86linux32.gz (or whatever the current linux version is called)
#and copy the link address
#3. cd into your software folder, or wherever you want to put this (I recommend NOT in scratch)
#4. wget LINK ADDRESS
#5. gunzip NAME OF DOWNLOADED THING.gz
#6. Rename the DOWNLOADED THING to something more user friendly:
#mv usearch11.0.667_i86linux32 usearch
#7. Make it executable:
#chmod +x usearch

#export the path to the software
export usearch=/work/06909/cbscott/ls6/software/usearch

#this will cluster things with 99% similarilty. At some point, I recommend you
#play with this threshold (-id) to check the sensitivity of this clustering.
#I reduced mine to .95, and didn't get different results, so I kept it at .99

#This is fast! I ran it on the login node on 450 sequenes. If you have more, I recommend making a short dev job.
#Do it for proteins and fasta sequences. Make sure the results seem to align pretty well.
$usearch -cluster_fast contig.rps3.hits.prot.fa -id 0.99 -sort length -maxrejects 0 -maxaccepts 0 -centroids all_rpS3_centroids.prot.fa
$usearch -cluster_fast contig.rps3.hits.dna.fa -id 0.99 -sort length -maxrejects 0 -maxaccepts 0 -centroids all_rpS3_centroids.dna.fa

conda activate map_software #this is the conda environment where I keep samtools, bowtie, etc.
#if you already have these somewhere else, you don't need to install this conda environment,
#there's nothing special or project-specific about it.

#create a bowtie2 index from your clustered rps3 sequenes
export CONTIGS=/scratch/06909/cbscott/STX_SpawnMeta/rps3/all_rpS3_centroids.dna.fa
echo "bowtie2-build --large-index --threads 6 $CONTIGS $CONTIGS" >btb
#this took <30 mins for 450 seqs
python3 $HOME/bin/ls6_launcher_creator.py -j btb -n btb -a IBN21018 -e cbscott@utexas.edu -t 02:00:00 -N 1 -w 1 -q development

#######NOW MAP YOUR FASTQS BACK TO THESE CONTIGS! THIS IS HOW WE WILL GET SPECIES ABUNDANCE!#######
export mapdir=/scratch/06909/cbscott/STX_SpawnMeta/rps3/map2rps3 #where I"mputting my mapped files
#create a list of the fastq locations, for me, this is the full path to my forward read files
#you'll have to edit this mapping code to correctly deal with the structure of your reads.
#note, here I've used --end-to-end mapping instead of --local because I want to be more stringent
#with what I accept as a hit.
>map_contigs
for file in `cat fastqlocs`; do
echo "bowtie2 -x $CONTIGS -1 $file -2 ${file/_R1/_R2} --end-to-end --no-unal -p 24 -S ${file/_R1.merged.trim.fq/}.rps3.sam && \
samtools sort -o ${file/_R1.merged.trim.fq/}.rps3.bam ${file/_R1.merged.trim.fq/}.rps3.sam && \
samtools index ${file/_R1.merged.trim.fq/}.rps3.bam" >> map_contigs; done
python3 $HOME/bin/ls6_launcher_creator.py -j map_contigs -n map_contigs -a IBN21018 -e cbscott@utexas.edu -t 00:15:00 -N 1 -w 12 -q development
#I messed up some of the paths above, so I had to move my bams to my target directory here:
cp ../clean_fastqs/*.rps3* /scratch/06909/cbscott/STX_SpawnMeta/rps3/map2rps3


###!!!!! IMPORTANT!!!!!! ######
#We only want the reads which mapped with high quality.
#Otherwise, you will have both duplicate mappings in your data and lots of spurious hits.
#It's much better to keep fewer high quality mappings than more low quality mappings with low accuracy.
cd $mapdir
#I got frustrated and did this on login - is pretty fast for 100 files. If you have more,
#make a job.
for file in *.bam; do
samtools view -bSq 30 $file > ${file/_S*.rps3.bam}.q30.rps3.bam && \
samtools index ${file/_S*.rps3.bam}.q30.rps3.bam; done

#need to standardize these guys. Get average fold coverage per contig.
#I just did this on an idev node, but it's pretty fast.
conda activate hmmer

for file in *.q30.rps3.bam; do
pileup.sh twocolumn=t in=$file out=${file/.q30.rps3.bam}.coverage; done
mkdir coverage
mv *.coverage coverage #ONLY HAVE THESE .coverage FILES IN COVERAGE folder - will make your life easier later when we pull things local.

#for each sample, also need to standardize by total sequencing effort... e.g. the number of reads in the raw fasta.
#i'll just do this for the forward reads. It shouldn't matter, because I'm just going to scale
#coverage by this value. We'll divide it by 4 in R to get the actual read count (these are fastqs with quality info)

#this takes forever :( there's probably a faster way?
>rc
>fname
for file in *R1.merged.trim.fq; do
  wc -l $file >> rc &&
  echo ${file/_S*_R1.merged.trim.fq} >>fname; done

paste fname rc > sequencing_effort.txt #I know this duplicates the names, but is going to make your life a lot easier in R


##### even though we won't really use it, let's also get raw read count per contig.
#we can use this to sanity check our standardization, etc. later
for file in *q30.rps3.bam; do
samtools idxstats $file > ${file/_S*.q30.rps3.bam}.rps3cov; done
mv *.rps3cov contig_counts
cd contig_counts

#move sequencing effort, coverage files, and raw read count per contig files locally.
#Pro Tip:
scp -r username@tacc.ls6.tacc.utexas.edu:/path/to/FOLDER destination/path
#will recursively copy an entire FOLDER to your local machine rather than having to do each file one by one.

#for example,
scp -r cbscott@ls6.tacc.utexas.edu:mytaccdir/contig_counts /Projects/STX_Meta/data
#made a subdirectory called contig_counts within my local project directory.
