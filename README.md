**Streaming workflow for phylogenomic analysis and comparative proteomics**

**Contributors**

Qingxiang Guo


**About**

The following contents contains source code for the analyses and plots in my Biology paper, 2022, “Quantitative insights into the contribution of nematocysts to the adaptive success of cnidarians based on proteomic analysis”

**Abstract**

Cnidarians (such as corals, anemones, jellyfish) are ancient and successful animals. Previous studies have offered a variety of explanations for the adaptive success of certain cnidarian taxa. However, common strategies for the long-term persistence of cnidarians have not been identified. One factor that may contribute to their evolutionary success of the lineage is the nematocyst, a sting organelle able to deliver venom into prey or enemies. Using bioinformatics analyses, we aimed to quantitatively investigate the role of nematocysts in cnidarian adaptation. We identified the extensive species specific adaptation in nematocyst proteins (NEMs) and demonstrate that both a unique evolutionary pattern of NEMs and the long evolutionary lag between nematocysts and cnidarians support their key adaptive role. Further, we find NEMs experience approximately 50% more adaptive changes on average compared to non-NEMs, and positively selected cnidarian-conserved proteins are enriched in NEMs. These results support a key role of nematocysts in successful cnidarian adaptation and provide a general quantitative framework for assessing the role of a phenotypic novelty in adaptation. Moreover, the findings will be critical for reassessing the evolutionary history of many established models, enhancing our understanding of both the mechanisms and evolutionary preference of adaptive evolution.

**Notes**

Written by Qingxiang Guo, qingxiang.guo@outlook.com, distributed without any guarantees or restrictions.

**Codes**

**1. Orthologous gene family analysis by OrthoMCL**

**1.1 Disposing the MySQL**

vim /etc/my.cnf

myisam\_sort\_buffer\_size=128000M

read\_buffer\_size=2000M

mysql> show variables LIKE 'myisam\_max\_sort\_file\_size';          

mysql> show variables LIKE 'myisam\_sort\_buffer\_size';

mysql> show variables LIKE 'read\_buffer\_size';

\# Restart

service mysqld restart

**1.2 Install MCL**

tar zxf /opt/software/mcl-14-137.tar.gz -C /opt/biosoft/

./configure --prefix=/opt/biosoft/mcl-14-137/ && make -j 4 && make install

echo 'PATH=$PATH:/opt/biosoft/mcl-14-137/bin/' >> ~/.bashrc 

source ~/.bashrc

**1.3 Configuration**

/home/gqx/proteome/tmp\_test/ortho\_test

cp /opt/biosoft/orthomclSoftware-v2.0.9/doc/OrthoMCLEngine/Main/orthomcl.config.template ./orthomcl.config

vim orthomcl.config

\# Editing the text

dbConnectString=dbi:mysql:orthomcl

dbLogin=gqx

dbPassword=123456

**1.4 Install the MySQL tables required for OrthoMCL**

mysql -u gqx -p123456

CREATE DATABASE IF NOT EXISTS orthomcl;

orthomclInstallSchema orthomcl.config

**1.5 Prepare the input fasta data**

\# Edit the file name to make it recognizable for OrthoMCL

perl -p -e "s/\_/\|/g" Aurelia\_aurita\_Rachamim.fasta > fix.fasta

tar zxf /home/train/00.incipient\_data/data\_for\_genome\_comparison/compliantFasta.tar.gz -C ./

orthomclAdjustFasta ncra ncra\_proteins.fasta 1

rm ncra\_proteins.fasta

**1.6 Remove low-quality sequence**

\# Longer than 30 aa, and the stop codon ratio threshold is 20%

orthomclFilterFasta compliantFasta/ 10 20

**1.7 All-vs-all BLAST on the filtered result goodProteins.fasta**

makeblastdb -in goodProteins.fasta -dbtype prot -title orthomcl -parse\_seqids -out orthomcl -logfile orthomcl.log

nohup blastp -db orthomcl -query goodProteins.fasta -seg yes -out blastresult.tab -evalue 1e-5 -num\_threads 64 -outfmt 6 &

**1.8 Process the results of BLAST to get the similarity information between the sequences, so as to facilitate importing into MySQL**

orthomclBlastParser blastresult.tab compliantFasta/ > similarSequences.txt

\# Then set myisam\_max\_sort\_file\_size (five times larger) in mysql according to the size of the similarSequences.txt file, 

myisam\_max\_sort\_file\_size=2955M

perl -p -i -e 's/0\t0/1\t-181/' similarSequences.txt

**1.9 Import similarSequences.txt to MySQL database**

orthomclLoadBlast orthomcl.config similarSequences.txt

**1.10** **Find similar sequence pairs: pairs that best match each other between species (orthologs pairs); pairs that best match each other and are better than interspecies matches within species (in-paralogs pairs); the first two are combined into co-ortholog pairs**

orthomclPairs orthomcl.config orthomcl\_pairs.log cleanup=no

**1.11 Export similar sequence pairs from MySQL**

orthomclDumpPairsFiles orthomcl.config

\# output: mclInput, pairs file

**1.12 Use mcl to cluster pairs (Ortholog Cluster Groups) and number the classes**

mcl mclInput --abc -I 1.5 -o mclOutput

orthomclMclToGroups OCG 1 < mclOutput > groups.txt

\# groups.txt contains the final cluster information

**1.13 Manual extraction of cluster results**

orthomcl\_two.pl groups.txt mhon mwul

orthomcl\_two.pl groups.txt mhon tkit

orthomcl\_two.pl groups.txt mwul tkit

orthomcl\_three.pl groups.txt

**1.14 Plot three-set Venn diagram of the cluster results**

\# First, match each sequence ID to OG in the result

orthomcl\_result\_2\_header.pl groups.txt

cd ..

mv 0\_orthomcl/header ./

\# Create header file one by one for each species

grep ">" 0\_orthomcl/compliantFasta/cfle.fasta > cfle

perl -p -i -e "s/>//g" cfle

\# Then replace the original sequence ID according to the corresponding relationship

orthomcl\_header\_2\_OG.pl mhon header

Rscript three\_set\_venn.R

**1.15 Statistics of the orthologs**

orthomcl\_genes\_number\_stats.pl groups.txt compliantFasta > genes\_number\_stats.txt

**1.16 Extract the single-copy orthologs from the OrthoMCL results**

orthomcl\_findSingleCopyOrtholog.pl groups.txt compliantFasta/

**1.17 Batch alignment of sequences**

fastas2aligned\_by\_mafft.pl SingleCopyOrthologGroups 8

\# You will get the final allSingleCopyOrthologsAlign.fasta for downstream analysis (e. g., phylogenomics)

**2. Orthologous gene family analysis by HaMStR and phylogenomic data matrices construction**

**2.1 Run HaMStR to get orthologs**

\# Use model organisms core set

nohup hamstr -protein -sequence\_file=../ALALA.fasta -taxon=ALALA -hmmset=modelorganisms\_hmmer3 -refspec=DROME -relaxed -central -force -cpu 24 -eval\_blast=0.01 -eval\_hmm=0.01 &

\# Use basal metazoan core set

nohup hamstr -protein -sequence\_file=rename\_all.fasta -taxon=CEBOR -hmmset=basalmets\_hmmer3 -refspec=nemve\_2309  -relaxed -central -force -cpu 48 -eval\_blast=0.01 -eval\_hmm=0.01  &

\# To run many species in batch

produce\_MODEL\_hamstr\_cmd\_list.pl Ctenactis\_echinata

for dir in \*/; do produce\_MODEL\_hamstr\_cmd\_list.pl $dir; done

for dir in \*/; do produce\_BASELMET\_hamstr\_cmd\_list.pl $dir; done

**2.2 Combine different OGs**

grep --p "hydma" basalmets\_hmmer3.fa > basl\_head

fast\_extract\_seq\_from\_fasta.pl basalmets\_hmmer3.fa basl\_head > basl.fasta

cat ./\* > model.fasta

mv /home/gqx/proteome/10\_phylogenomic\_tree/species/Hydrozoa/Hydra\_vulgaris/model/fa\_dir\_HYVUL\_modelorganisms\_hmmer3\_DROME\_relaxed/model.fasta ./

grep --p "HYVUL\|" model.fasta > model\_head

fast\_extract\_seq\_from\_fasta.pl model.fasta model\_head > model2.fasta

makeblastdb -in model2.fasta -out model -dbtype prot

nohup blastp -query basl.fasta -db model -max\_target\_seqs 1 -outfmt 6 -out basl\_model -evalue 1e-10 -num\_threads 32 &

perl -lane 'print $\_ if ($F[2] >= 70)' basl\_model > basl\_model\_70

\# Take the species with high coverage as a bridge

\# Extracted correspondence

extract\_OGs\_relationship.pl basl\_model\_70

cat ./connect\_by\_human/rela ./connect\_by\_hydra/rela ./connect\_by\_nemve/rela ./connect\_by\_Trichoplax/rela > relationship

remove\_duplicate.pl relationship

\# Get the same genes from all species into the same file

cat ./\*/met/hamstrsearch\_\* >> /home/gqx/proteome/10\_phylogenomic\_tree/my\_alignment/met\_out

cat ./\*/\*/model/hamstrsearch\_\* > /home/gqx/proteome/10\_phylogenomic\_tree/my\_alignment/model\_out

perl -p -i -e "s/^(\d+)\|\w+?\|/\1\|/g" met\_out

perl -p -i -e "s/^(\d+)\|\w+?\|/\1\|/g" model\_out

\# The missing reference sequence is added

grep --p "nemve\_2309|hydma\_2|triad\_390" basalmets\_hmmer3.fa > head

fast\_extract\_seq\_from\_fasta.pl basalmets\_hmmer3.fa head > haha.fa

perl -p -i -e "s/nemve\_2309/NEVEC/g" haha.fa

perl -p -i -e "s/hydma\_2/HYVUL/g" haha.fa

perl -p -i -e "s/triad\_390/TRADH/g" haha.fa

perl -p -i -e "s/^>(\S+)/\1\|1\|/g" haha.fa

perl -p -i -e "s/\|\n/\|/g" haha.fa

cat haha.fa >> /home/gqx/proteome/10\_phylogenomic\_tree/my\_alignment/met\_out

\# Check species

cat model\_out | cut -f 2 -d "|" > 22

remove\_duplicate.pl 22

wc -l duplicate\_remove

\# Then split out sequences under the same OG

phylogenomic\_seperate\_OGs\_1.pl met\_out

phylogenomic\_seperate\_OGs\_2.pl out

cd OGs

perl -p -i -e 'tr/,/\n/s' ./\*

for DIR in ./\*txt; do phylogenomic\_seperate\_OGs\_3.pl $DIR; done;

rm ./\*txt

\# Then find the intersection between model core set and metazoan core set

cd met\_model\_inter\_OGs/

cp ../met/OGs/\* ./

cp ../model/OGs/\* ./

inter\_OG\_met\_model\_remove.pl relationship

rm relationship

**2.3 Quality control before alignment**

mkdir trim\_OGs

cd trim\_OGs

cp ../met\_model\_inter\_OGs/\* ./

\# Change a more convenient sequence header

perl -p -i -e "s/>\d+\|(\w+)\|(\w+)\|\d\|/>\1\|\2/g" ./\*fasta

\# Delete sequence length shorter than 100 aa

phylogenomic\_delete\_short\_seq.pl

\# Delete OGs with taxon coverage less than specific species (n = 50 in present species)

./phylogenomic\_few\_tax.sh

\# Delete X at the first and last 20 aa of the sequence

./phylogenomic\_remove\_X.sh

**2.4 Alignment**

/opt/scripts/fastas2aligned\_by\_mafft\_hamstr\_version.pl trim\_OGs 16

**2.5 Trimming**

for dir in ./\*aln; do phylogenomic\_trimal\_BJZ2X.pl $dir; done

rm ./\*aln

rename .aln.BJZ2X .aln \*.aln.BJZ2X

for dir in ./\*aln; do trimal -in $dir -out ${dir}.trim -automated1; done

perl -p -i -e "s/^(>\S+)\s+$/\1\n/" ./\*aln.trim

./phylogenomic\_few\_tax\_2.sh

cd /home/gqx/proteome/10\_phylogenomic\_tree/single\_gene\_tree/alignment

cp ../../my\_alignment/trim\_OGs/\*trim ./

rename .aln.trim .fa \*.aln.trim

cp /opt/biosoft/PhyloTreePruner/remove\_short\_seqs.sh ./

\# Change the fasta file to flat format

for dir in ./\*fa; do flat\_the\_fasta\_seq.pl $dir; mv flated $dir.flat; done;

rm ./\*fa

rename .fa.flat .fa \*.fa.flat

\# Removing sequences that are more than 50 percent gaps

./remove\_short\_seqs.sh

**2.6 Build single-gene trees using FastTree**

phylogenomic\_FastTree\_OG.pl alignment 16

FastTreeMP -slow -gamma $ARGV[0]/$\_ > $ARGV[0]/$1.tre

**2.7 Use TreSpEx to remove Long-branch attraction (LBA) and saturation analyses**

cd TreSpEx/

mkdir single\_tree

cd single\_tree

\# TreSpEx needs the RAXML output format

cp /home/gqx/proteome/10\_phylogenomic\_tree/ML\_analysis/76cov\_105tx\_146og/out.fasta ./

cp /home/gqx/proteome/10\_phylogenomic\_tree/ML\_analysis/76cov\_105tx\_146og/new\_partition.txt ./

cp charset new\_partition.txt

perl -p -i -e "s/\s;\s//g" new\_partition.txt

perl -p -i -e "s/\s-\s/-/g" new\_partition.txt

perl -p -i -e "s/charset/AUTO,/g" new\_partition.txt

perl -p -i -e "s/^\S+\s//g" new\_partition.txt

./AMAS.py split -l new\_partition.txt -i out.fasta -f fasta -d aa

rm out.fasta new\_partition.txt

for DIR in ./\*fas; do  mkdir ${DIR}\_file; done;

for DIR in ./\*fas; do  mv $DIR  ${DIR}\_file; done;

for DIR in ./\*file; do  cd $DIR; remove\_wenhao\_seq.pl \*fas; awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' out > output; rm out; cd ..;  done;

for DIR in ./\*file; do  cd $DIR; mv output out\_\*; cd ..;  done;

\# Run iqtree

nohup iqtree-omp -s ./\*fas -wbt -bb 10000 -nt 1 -st AA &

\# Prepare the species name list, tree name list and alignment name list

cp /home/gqx/proteome/10\_phylogenomic\_tree/concatenated/146/otu/otu-name.otu ./

ls -1 ./\*contree > tree\_name

\# Make sure the format of alignments is relax phylip

cd /home/gqx/proteome/10\_phylogenomic\_tree/ML\_analysis/76cov\_105tx\_146og/TreSpEx/Longbranch\_and\_Saturation

mkdir align

cd align

cp ../../single\_tree/\*/\*fas ./

for dir in ./\*fas; do convertFasta2Phylip.sh $dir > ${dir}.phy ; done

ls -1 ./\*phy > align\_name

mv align\_name ../

mv ./\*phy ../

cd ..

rm align -rf

\# Perform LBA analysis

TreSpEx.v1.1.pl -fun e -ipt tree\_name -tf otu-name.otu

\# Saturation analysis

TreSpEx.v1.1.pl -fun g -ipt tree\_name -ipa align\_name

**2.7.1 Use density map to remove LBA (in R)**

data <- read.table("LB\_scores\_summary\_perPartition.txt", header=T,sep="\t",row.names = 1)

d <- density(data$LB\_score\_upper\_quartile)

pdf("density\_upper.pdf",width=6,height=6)

plot(d, main=NA, xlab="Mean value of upper quartile of LB scores", lwd=1.5)

\# plot(d, main=NA, xlab="Standard deviation of LB scores", lwd=1.5)

segments(x0=80.81,y0=0,x1=80.81,y1=0.0125,lwd=1.5, lty=2)

abline(h=0,col="grey")

\# Filling detected areas with red

x1 <- min(which(d$x >= 80.81))

x2 <- max(which(d$x < 140))

with(d, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="red"))

dev.off()

**2.7.2 Use density map to do saturation analysis (in R)**

Correlation\_Results/Correlation\_Slope\_Summary.txt

data <- read.table("Correlation\_Slope\_Summary.txt", header=T,sep="\t",row.names = 1)

d <- density(data$Slope)

pdf("density\_slope.pdf",width=6,height=6)

plot(d, main=NA, xlab="Slope of linear regression", lwd=1.5)

segments(x0= 0.224,y0=0,x1= 0.224,y1=1.5,lwd=1.5, lty=2)

abline(h=0,col="grey")

x1 <- min(which(d$x >= -1))

x2 <- max(which(d$x < 0.224))

with(d, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="red"))

dev.off()

**2.7.3 Use heatmap to analyze species-specific genes (in R)**

library(gplots)

data <- read.table("9999.txt", header=T,sep="\t",row.names = 1)

data <- as.matrix(data)

colors = c(seq(min(data),0,length=100),seq(1,200,length=100),seq(201,500,length=100),seq(501,700,length=100),seq(701,1000,length=100),seq(1001,1200,length=100),seq(1201,1600,length=100),seq(9999,9999,1))

my\_palette <- colorRampPalette(c("orange","yellow","green", "cyan","blue","purple","red","white"))(n = 700)

pdf("9999.pdf",width=18,height=18)

heatmap.2(data, density.info="density", col=my\_palette, symkey=F, symm=F, symbreaks=F,scale="none", trace="none", breaks=colors, keysize=1)

dev.off()

**2.8 Use PhyloTreePruner to remove paralogs**

\# Change the header to be recognizable for PhyloTreePruner

perl -p -i -e "s/(\d+)\|//g" 1.aln.trim

\# Feed the result of RAxML (RAxML\_bipartitions.out) to PhyloTreePruner

java -cp /usr/local/bin/ PhyloTreePruner RAxML\_bipartitions.out 3 2.aln.trim 0.5 r

java -cp /usr/local/bin/ PhyloTreePruner 113684.tre 71 113684.aln 0.5 r

cd single\_gene\_tree/

mkdir paralog\_screen

cd paralog\_screen/

cp ../alignment/\*long ./

cp ../alignment/\*tre ./

\# Change the fasta file format to flat fasta

for dir in ./\*fa.long; do fasta\_flat\_to\_multiline.py $dir 60; done;

rm ./\*fa.long

rename .fa.long\_multi-line.fasta  .fa.long \*.fa.long\_multi-line.fasta

cd ..

\# Run PhyloTreePruner

phylogenomic\_PhyloTreePruner.pl paralog\_screen 8

**2.9 Filtering for unique sequences**

SCaFoS

export SCAFOS=/opt/biosoft/SCaFoS

cpan install Tk

\# Install Tree-Puzzle

perl -p -i -e "s/\|/@/g" ./\*fa

scafos in=fasta out=otu

scafos in=fasta out=out otu=otu/otu-name.otu format=fp t=10

\# out.fasta is the result

\# Change the style to MARE style

phylogenomic\_scafos2MARE\_charset.pl out.len

**2.10 Visualize the phylogenomic matrices**

MARE charset out.fasta -m -t

**2.11 Exclude species with gene missing rate > 85%**

cd /home/gqx/proteome/10\_phylogenomic\_tree/concatenated/50\_taxon\_coverage\_remove\_poor\_species/out

mkdir count

cd count

split\_fasta\_file\_averagely.pl ../out.fasta 118

rm out.fasta

for DIR in ./\*fasta; do  mkdir ${DIR}\_file; done;

for DIR in ./\*fasta; do  mv $DIR  ${DIR}\_file; done;

for DIR in ./\*; do cd $DIR; AMAS.py summary -i ./\*fasta -f fasta -d aa; cd ..; done;

for DIR in ./\*; do cd $DIR; cat summary.txt | cut -f6 > missing | grep ">" \*fasta >> missing; cd ..; done;

for DIR in ./\*; do cd $DIR; perl -p -i -e "s/\n/\t/g" missing; cd ..; done;

cat ./\*/missing >> all\_missing

perl -p -i -e "s/\tMissing/\nMissing/g" all\_missing

rm ./\*file -rf

perl -lane 'print $\_ if ($F[1] gt 85)' all\_missing > all\_missing\_85

\# Then delete species in the out.fasta file

**3. Phylogenetic tree construction**

**3.1 Run RAXML**

\# LG+G4+F is a site-homogeneous model, while CAT+G4 is a site-heterogeneous model. It is possible to use PartitionFinder with different parameters on each module separately, but this hybrid homogeneous is not equal to the heterogeneous model. The evolutionary context within each module is also different.

\# CAT+G4 also considers that the evolutionary context may be different even within a module.

\# The optimal model is calculated automatically; and the alpha parameters of gamma distribution and stable frequency differ according to each gene.

phylogenomic\_scafos2RAXML\_partition.pl out.len

nohup raxmlHPC-PTHREADS-SSE3 -T 48 -s out.fasta -n 50cov\_105tx\_474og -f a -N 100 -x 12345 -p 12345 -m PROTGAMMAAUTO -q new\_partition.txt &

**3.2 Run IQ-TREE**

\# Change the charst into raxml new\_partition.txt

cp charset new\_partition.txt

perl -p -i -e "s/\s;\s//g" new\_partition.txt

perl -p -i -e "s/\s-\s/-/g" new\_partition.txt

perl -p -i -e "s/charset/protein,/g" new\_partition.txt

nohup iqtree-omp -s out.fasta -wbt -bb 10000 -nt AUTO -st AA -spp new\_partition.txt  -m MFP &

**3.3 Run Phylobayes**

/opt/biosoft/RAxML\_v8.2.9/usefulScripts/convertFasta2Phylip.sh ../allSingleCopyOrthologsAlign.fasta-gb > allSingleCopyOrthologsAlign.fasta-gb.phy

nohup mpirun -np 16 /opt/biosoft/phylobayes-mpi-1.7b/data/pb\_mpi -d allSingleCopyOrthologsAlign.fasta-gb.phy -cat -gtr -x 10 20000 gqx &

nohup mpirun -np 16 /opt/biosoft/phylobayes-mpi-1.7b/data/pb\_mpi -d allSingleCopyOrthologsAlign.fasta-gb.phy -cat -gtr -x 10 20000 gqx2&

\# Check convergence and output results

nohup /opt/biosoft/phylobayes-mpi-1.7b/data/bpcomp -x 1000 10 gqx gqx2 -o test &

nohup /opt/biosoft/phylobayes-mpi-1.7b/data/tracecomp -x 1000 gqx gqx2 -o test &

\# effsize > 300, max rel\_diff < 0.1, converged

**4. Remove bilaterial contamination from myxozoan sequences**

mkdir contamination\_blast

cd contamination\_blast

\# Make seed database: bilaterian.faa and cnidarian.faa

replace\_header\_for\_cdhit.pl -c bilaterian\_100 -t bilaterian

replace\_header\_for\_cdhit.pl -c cnidarian\_100 -t cnidarian

cat ./\* > all.fasta

\# Remove the non-myxozoan sequences in each matrices

cp ../fasta/ ./ -R

cd fasta

for dir in ./\*fa ; do perl -p -i -e "s/(\d)\n/\1\t/g" $dir; done

for dir in ./\*fa ; do phylogenomic\_filter\_no\_myxozoan\_seq.pl $dir; done

cd ..

mkdir myxo\_fasta

mv fasta/\*filter myxo\_fasta/

cd ..

makeblastdb -in all.fasta -out all -dbtype prot

./phylogenomic\_parafly\_contamination\_blast.pl myxo\_fasta/

ParaFly -c blast\_cmds -CPU 4 -shuffle -v

cd myxo\_fasta/

for dir in ./\*out ; do export LANG=C; export LC\_ALL=C; sort -k1,1 -k12,12gr -k11,11g -k3,3gr $dir | sort -u -k1,1 --merge > ${dir}\_best ; done

\# Get the bad list

for dir in ./\*out\_best ; do perl -ne '/bilaterian/ && print' $dir > ${dir}\_delete ; done

cd ..

mkdir delete\_list

mv myxo\_fasta/\*delete delete\_list/

cd delete\_list/

#find ./\* -size  0 -print0 |xargs -0 rm

for dir in ./\*delete ; do perl -p -i -e "s/^(\S+).\*$/\1/g" $dir ; done

\# If you only want to remove certain species

for dir in ./\*out\_best\_delete ; do perl -ne '/HESPP/ && print' $dir > ${dir}\_HESPP ; done

cd /home/gqx/proteome/11\_big\_tree/ML\_analysis/optimized\_105tx\_232og\_3/fasta

cp ../contamination\_blast/delete\_list/\* ./

phylogenomic\_remove\_contaminant\_of\_no\_cnidarian.sh

rm ./\*fa

rm ./\*best\_delete

rename .fa\_revised .fa \*.fa\_revised

**5. Approximately unbiased test (AU-test) of phylogenetic tree selection**

**5.1 Make constraint tree**

perl -p -i -e "s/:0\.\d+//g" bestree.nwk

phylogenomic\_constraint\_tree\_generator.pl tree\_topology mono\_list bestree.nwk

\# constraint.tre is the result file

**5.2 AU-test using RAXML and Consel**

\# Provide the main tree

nohup raxmlHPC-PTHREADS-SSE3 -T 48 -s out.fasta -n 50cov\_105tx\_474og -f a -N 100 -p 12345 -x 12345 -m PROTGAMMAAUTO -g constraint.tre -q new\_partition.txt &

\# Provide constraint tree to RAXML

\# Produce likehood value from the result tree

cat ../normal/RAxML\_bipartitions.50cov\_105tx\_474og ../constraint\_1/RAxML\_bipartitions.50cov\_105tx\_474og ../constraint\_2/RAxML\_bipartitions.50cov\_105tx\_474og ../constraint\_3/RAxML\_bipartitions.50cov\_105tx\_474og > tree

nohup raxmlHPC-PTHREADS-SSE3 -s out.fasta -z tree -m PROTGAMMAAUTO -T 6 -f G -n  50cov\_105tx\_474og  -q new\_partition.txt &

\# AU-test using Consel

mv RAxML\_perSiteLLs.50cov\_105tx\_474og RAxML\_perSiteLLs.sitelh

makermt --puzzle RAxML\_perSiteLLs.sitelh

consel RAxML\_perSiteLLs.rmt

catpv -s 9 RAxML\_perSiteLLs.pv

**5.3 AU-test using IQ-tree**

\# The main difference between IQ-tree and RAXML is that IQ-tree integrates the functionality of Consel. The second and third steps of RAXML are combined into one step

\# Get the main tree

nohup iqtree-1.6 -s new.fasta -wbt -bb 10000 -nt 6 -st AA -spp new\_partition.txt &

\# Provide the constraint tree to IQ-tree and get the result tree

nohup iqtree-1.6 -s new.fasta -wbt -bb 10000 -nt 6 -st AA -spp new\_partition.txt -g constraint.tre &

\# Use the constraint tree in former step as guide tree and conduct PMSF optimization 

mv new\_partition.txt.contree guide.tre

\# Here is a question: when doing PMSF, do you need to use both partitioned and C60,C20 (similar to CAT, site-heterogeneous substitution models), IQ-treee allows to do these at the same time, but each partition and C20 are not linked, so for most data will be over-parameterized so no partition is used here

\# Use LG+G4+F+C2 model

nohup iqtree-1.6 -s new.fasta -nt 16 -st AA -m "LG+G4+FMIX{C20pi1:1:0.0559910600,C20pi2:1:0.0514824870,C20pi3:1:0.0812922124,C20pi4:1:0.0721976867,C20pi5:1:0.0556718858,C20pi6:1:0.0331003080,C20pi7:1:0.0589501763,C20pi8:1:0.0263756889,C20pi9:1:0.0307584220,C20pi10:1:0.0376701125,C20pi11:1:0.0303058290,C20pi12:1:0.0808775576,C20pi13:1:0.0263349134,C20pi14:1:0.0579101455,C20pi15:1:0.0371248064,C20pi16:1:0.0586867766,C20pi17:1:0.0561479138,C20pi18:1:0.0349810886,C20pi19:1:0.0544937394,C20pi20:1:0.0596471901}+F" -ft guide.tre -g constrain.tre &

\# Running AU-test

cat normal/new\_partition.txt.contree constraint\_1/PMSF/new.fasta.treefile constraint\_2/PMSF/new.fasta.treefile constraint\_3/PMSF/new.fasta.treefile > AU/tree

nohup iqtree-1.6 -s new.fasta -nt 16 -m "LG+G4+FMIX{C20pi1:1:0.0559910600,C20pi2:1:0.0514824870,C20pi3:1:0.0812922124,C20pi4:1:0.0721976867,C20pi5:1:0.0556718858,C20pi6:1:0.0331003080,C20pi7:1:0.0589501763,C20pi8:1:0.0263756889,C20pi9:1:0.0307584220,C20pi10:1:0.0376701125,C20pi11:1:0.0303058290,C20pi12:1:0.0808775576,C20pi13:1:0.0263349134,C20pi14:1:0.0579101455,C20pi15:1:0.0371248064,C20pi16:1:0.0586867766,C20pi17:1:0.0561479138,C20pi18:1:0.0349810886,C20pi19:1:0.0544937394,C20pi20:1:0.0596471901}+F" -z tree -n 0 -zb 10000 -au &

\# Results are summarized in new.fasta.iqtree

**6. Fast-evolving site removal analysis**

**6.1 Install DistEst**

tar zxf dist\_est.tar.gz

\# distEST needs two function library, id05-1.0.0.tar.gz and ve11-1.0.0.tar.gz

\# unzip these files, copy the id05i.f, ve11d.f in src directory into the dist\_estv1.1 directory

make

\# If you see information below, you succeed

gcc  -O3  -c dcdflib.c -o dcdflib.o

gcc  -O3  -c dist\_est.c -o dist\_est.o

gcc  -O3  -c dist\_est\_fn.c -o dist\_est\_fn.o

gcc  -O3  -c ipmpar.c -o ipmpar.o

gfortran -O3  -c id05i.f -o id05i.o

gfortran -O3  -c ve11d.f -o ve11d.o

gfortran -O3 -o dist\_est dcdflib.o dist\_est.o dist\_est\_fn.o ipmpar.o id05i.o ve11d.o -lm

**6.2 Calculating evolutionary rates in DistEst**

\# Prepare hsp70.ctl， hsp70.dat， hsp70.fitchml.tre

dist\_est hsp70.ctl

\# hst70.ctl is as follows

treefile =  hsp70.fitchml.tre  \* treefile

seqfile = hsp70.dat            \* sequence data

nchar = 20             \* amino acid data

model =  3             \* empirical + F

aaRatefile = lg.dat \* LG substitution model

nrate = 101            \* number of rates

ub = 10.0              \* upper bound for rates

\# The resule file: DE.dat  rate\_est.dat

**6.3 Remove sites from fast to slow**

AMAS.py convert -i hsp70.dat -f phylip -d aa -u fasta

flat\_the\_fasta\_seq.pl metazoan.phy-out.fas

mv flated alignment.fasta

phylogenomic\_increment\_remove\_fast\_evolving\_sites.pl <fasta\_alignment>  <rate\_est.dat>  <num\_sites\_per>

\# For instance, phylogenomic\_increment\_remove\_fast\_evolving\_sites.pl alignment.fasta rate\_est.dat 100

**6.4 Build phylogenetic tree for each matrices**

nohup iqtree-omp -s out.fasta -wbt -bb 10000 -nt 24 -st AA -m LG+F+G &

\# You can also use PMSF model

nohup iqtree-omp -s out.fasta -wbt -bb 10000 -nt 24 -st AA -m LG+C60+F+G -ft guide.tre &

**6.5 Plot fast-evolving site removal analysis using R**

library(ggplot2)

data <- read.csv("table.csv", header=T)

data <- data.frame(A=c(data[1]),B=c(data[2]),C=c(data[3]))

ggplot(data=data, aes(x=seq, y=bp, group=supp))+geom\_line(aes(color=supp))+scale\_color\_brewer(palette="Dark2")+scale\_x\_continuous(breaks = c(0, 3500, 7000,10500,14000,17500,21000,24500,28000,31500,35000,38500,42000,45500,49000))+theme\_set(theme\_bw())+theme(legend.position="bottom")

ggsave("plot2.pdf", width=8, height=4)

**7. Molecular clock analysis using r8s**

\# Strick model: A rigorous model that assumes that all branches evolve at the same rate.

\# Fixed local clock: a type of relaxed model, assuming that the taxa of one or more branches are distinct from the whole, but of course, they must be monophyletic.

\# Uncorrelated relaxed clock: a type of relaxed model, where each branch has its own evolutionary rate and is uncorrelated.

\# Random local clock: a model between the discrete model and the relaxed model.

\# To use r8s,aAt least one internal node needs to be fixed. The whole process includes the first validation to determine the smooth value and the formal run using the optimal smooth value. cvstart=0 cvinc=0.5 cvnum=10 means to start from 0 and go 10 at 0.5 intervals, taking values to smooth. It is better to set cvstart as zero, do not choose a negative number. In addition, if it is a large phylogenetic tree, you must have collapse parameter to merge branches. Otherwise there will be some branches because of very short and be treated as zero length, which can not be calculated.

**7.1 Prepare the configuration files r8s\_in.txt**

cp /opt/biosoft/r8s1.81/r8s\_in\_template.txt ./

mv r8s\_in\_template.txt r8s\_in.txt

\# Then edit the configuration files based on conditions

**7.1.1 If you have 2 or more fossil fix points (i.e. clear fossil time points)**

\# First prepare r8s\_in\_validation.txt for first validation to set the bet smooth value

#NEXUS

begin trees;

tree tree\_1 = [&R] ((umay:0.377714,(pchr:0.126947,((cput:0.143183,slac:0.089884):0.040327,(scom:0.181745,(abis:0.150637,(lbic:0.088905,ccin:0.132642):0.022419):0.025092):0.02715):0.035777):0.207622):0.142535,((scer:0.463674,psti:0.376001):0.331119,((anig:0.213676,snod:0.244829):0.042769,(tree:0.147333,(ncra:0.134451,cpar:0.146872):0.034971):0.111483):0.206145):0.142535);

end;

begin r8s;

blformat lengths=persite nsites=320272 ulrametric=no;

collapse;

MRCA AGA ccin lbic scom;

MRCA ASC psti anig cpar snod tree;

MRCA BOL slac cput;

fixage taxon=AGA age=122.74;

fixage taxon=ASC age=517.55;

fixage taxon=BOL age=104.23;

divtime method=PL crossv=yes fossilfixed=yes cvstart=0 cvinc=0.5 cvnum=10;

showage;

describe plot=cladogram;

describe plot=chrono\_description;

end;

\# First validation

r8s -b -f r8s\_in\_validation.txt > cross\_validation\_1

\# Find the line with the smallest raw error, and then find the smooth value in the same line

\# Then revise the configuration file as a new file: r8s\_in\_true\_run.txt

#NEXUS

begin trees;

tree tree\_1 = [&R] ((umay:0.377714,(pchr:0.126947,((cput:0.143183,slac:0.089884):0.040327,(scom:0.181745,(abis:0.150637,(lbic:0.088905,ccin:0.132642):0.022419):0.025092):0.02715):0.035777):0.207622):0.142535,((scer:0.463674,psti:0.376001):0.331119,((anig:0.213676,snod:0.244829):0.042769,(tree:0.147333,(ncra:0.134451,cpar:0.146872):0.034971):0.111483):0.206145):0.142535);

end;

begin r8s;

blformat lengths=persite nsites=320272 ulrametric=no;

collapse;

MRCA AGA ccin lbic scom;

MRCA ASC psti anig cpar snod tree;

MRCA BOL slac cput;

fixage taxon=AGA age=122.74;

fixage taxon=ASC age=517.55;

fixage taxon=BOL age=104.23;

set smoothing=4;

set num\_time\_guesses=3;

set checkGradient=yes;

divtime method=PL algorithm=TN;

showage;

describe plot=cladogram;

describe plot=chrono\_description;

end;

\# Then run r8s

r8s -b -f r8s\_in\_true\_run.txt > true\_run\_1

\# Extract the time tree from the result file

perl -e 'while (<>) { if (m/tree tree\_1 = (.\*);/) { my $tree = $1; $tree =~ s/(BOL|AGA|ASC)//g; print "$tree\n"; } }' true\_run\_1 > tree.txt

**7.1.2 If you have 2 or more constraint points (most of the time you are in this condition)**

#NEXUS

begin trees;

`	`tree tree\_1 = [&R] (((Cerianthus\_borealis:0.290575,(((((Edwarsiella\_lineata:0.066246,Nematostella\_vectensis:0.081163):0.055217,((Exaiptasia\_pallida:0.067063,Hormathia\_digitata:0.134163):0.034484,((Bolocera\_tuediae:0.051153,Megalactis\_griffithsi:0.06898):0.043791,(Heteractis\_crispa:0.056233,(Anthopleura\_elegantissima:0.052055,(Anemonia\_viridis:0.00563,Anemonia\_sulcata:0.004281):0.02116):0.011187):0.01741):0.038479):0.062074):0.064797,((Antipathes\_caribbeana:0.053896,Plumapathes\_pennacea:0.045453):0.108419,((Corynactis\_australis:0.066395,(Ricordea\_yuma:0.065164,((Amplexidiscus\_fenestrafer:0.053832,Rhodactis\_indosinensis:0.014454):0.01069,Discosoma\_spp:0.03959):0.040398):0.011425):0.028652,((((Acropora\_digitifera:0.057801,Acropora\_millepora:0.008548):0.005001,Acropora\_tenuis:0.021788):0.090596,(Astreopora\_spp:0.048089,(Porites\_lobata:0.002024,Porites\_australiensis:0.006654):0.048218):0.012396):0.014759,(((Lobactis\_scutaria:0.008759,Ctenactis\_echinata:0.006969):0.027244,(Montastraea\_cavernosa:0.024341,(Pseudodiploria\_strigosa:0.027266,(((Favia\_lizardensis:0.002407,Favia\_spp:0.10512):0.006193,Platygyra\_carnosus:0.16142):0.00544,Orbicella\_faveolata:0.03004):0.003371):0.005356):0.017989):0.011364,(Madracis\_auretenra:0.027553,(Pocillopora\_damicornis:0.035619,(Stylophora\_pistillata:0.030801,(Seriatopora\_hystrix:0.015099,Seriatopora\_spp:0.035017):0.010937):0.007567):0.034336):0.022343):0.017003):0.027674):0.07751):0.026608):0.020789,Protopalythoa\_variabilis:0.277447):0.063276,((Pennatula\_rubra:0.055175,Renilla\_reniformis:0.107084):0.072312,(Corallium\_rubrum:0.073471,(Nephthyigorgia\_spp:0.062532,((Gorgonia\_ventalina:0.051904,Leptogorgia\_sarmentosa:0.023995):0.013823,(Eunicella\_cavolinii:0.011368,Eunicella\_verrucosa:0.002738):0.025411):0.011668):0.049632):0.032084):0.188635):0.018531):0.036866,((Polypodium\_hydriforme:0.566426,(((Enteromyxum\_leei:0.648484,Kudoa\_iwatai:0.561639):0.680711,(Sphaeromyxa\_zaharoni:0.953538,(Henneguya:0.780518,((((Thelohanellus\_kitauei:0.307058,Myxobolus\_turpisrotundus:0.304868):0.180837,((Myxobolus\_wulii:0.22627,(Myxobolus\_honghuensis:0.015607,Myxobolus\_ampullicapsulatus:0.070408):0.17434):0.057403,Myxobolus\_pendula:0.245172):0.256762):0.21192,Myxobolus\_cerebralis:0.432823):0.219181,Myxobolus\_physophilus:0.60296):0.170284):0.129202):0.147595):0.130403,Tetracapsuloides\_bryosalmonae:0.80208):0.20069):0.142613,(((((Atolla\_vanhoeffeni:0.040317,Periphylla\_periphylla:0.065065):0.150184,((Aurelia\_aurita:0.079173,(Cassiopea\_xamachana:0.132511,(Stomolophus\_meleagris:0.05262,Rhopilema\_esculentum:0.037626):0.059135):0.020893):0.040423,((Cyanea\_capillata:0.013458,Cyanea\_nozakii:9.14E-4):0.112443,(Chrysaora\_fuscescens:0.054366,Pelagia\_noctiluca:0.058332):0.054913):0.032858):0.097924):0.032958,(Chironex\_fleckeri:0.098164,(Alatina\_alata:0.081036,Tripedalia\_cystophora:0.180489):0.02203):0.132863):0.024249,(Craterolophus\_convolvulus:0.100106,(Calvadosia\_campanulata:0.081895,(Haliclystus\_sanjuanensis:0.150622,Calvadosia\_cruxmelitensis:0.047193):0.020468):0.023085):0.204894):0.024348,((Liriope\_tetraphylla:0.144402,(Craspedacusta\_sowerbyi:0.163299,Aegina\_citrea:0.22533):0.02255):0.055681,((Ectopleura\_larynx:0.169364,(Hydra\_viridissima:0.041963,(Hydra\_oligactis:0.030741,Hydra\_vulgaris:0.028101):0.038597):0.092278):0.078655,(Millepora\_alcicornis:0.226375,((Turritopsis\_spp:0.15472,(Podocoryna\_carnea:0.102354,(Hydractinia\_polyclina:0.008327,Hydractinia\_symbiolongicarpus:0.007533):0.059407):0.051791):0.030935,(Physalia\_physalis:0.111875,((Nanomia\_bijuga:0.07543,Agalma\_elegans:0.042964):0.044384,(Craseoa\_lathetica:0.059619,Abylopsis\_tetragona:0.121102):0.037884):0.042855):0.075418):0.015359):0.016351):0.058833):0.111132):0.063606):0.055702):0.028066,((((Helobdella\_robusta:0.428146,Crassostrea\_gigas:0.25635):0.054388,Drosophila\_melanogaster:0.536965):0.050571,((Ciona\_intestinalis:0.459065,Taeniopygia\_guttata:0.299825):0.052585,Strongylocentrotus\_purpuratus:0.349709):0.026884):0.077267,(Trichoplax\_adhaerens:0.483558,(Amphimedon\_queenslandica:0.550768,Mnemiopsis\_leidyi:0.725216):0.059918):0.045348):0.028066);

end;

begin r8s;

blformat lengths=persite nsites=51598 ulrametric=no;

collapse;

MRCA OUT Helobdella\_robusta Mnemiopsis\_leidyi Trichoplax\_adhaerens;

MRCA MED Chironex\_fleckeri Agalma\_elegans Thelohanellus\_kitauei;

MRCA HEX Protopalythoa\_variabilis Edwarsiella\_lineata Lobactis\_scutaria;

MRCA HYD Liriope\_tetraphylla Agalma\_elegans Hydra\_vulgaris;

MRCA CNI Protopalythoa\_variabilis Chironex\_fleckeri Thelohanellus\_kitauei;

prune taxon=OUT;

fixage taxon=HYD age=500;

constrain taxon=MED min\_age=570;

constrain taxon=HEX min\_age=540;

constrain taxon=CNI max\_age=741;

set minRateFactor=0.000001;

set minDurFactor=0.000001;

divtime method=PL crossv=yes fossilconstrained=yes cvstart=0 cvinc=1 cvnum=15;

showage;

describe plot=cladogram;

describe plot=chrono\_description;

end;

\# Then conduct first validation

r8s -b -f r8s\_in\_validation.txt > cross\_validation\_1

\# Find the line with the smallest raw error, and then find the log10smoon value in the same line

\# Then revise the configuration file as a new file: r8s\_in\_true\_run.txt 

#NEXUS

begin trees;

`	`tree tree\_1 = [&R] (((Cerianthus\_borealis:0.290575,(((((Edwarsiella\_lineata:0.066246,Nematostella\_vectensis:0.081163):0.055217,((Exaiptasia\_pallida:0.067063,Hormathia\_digitata:0.134163):0.034484,((Bolocera\_tuediae:0.051153,Megalactis\_griffithsi:0.06898):0.043791,(Heteractis\_crispa:0.056233,(Anthopleura\_elegantissima:0.052055,(Anemonia\_viridis:0.00563,Anemonia\_sulcata:0.004281):0.02116):0.011187):0.01741):0.038479):0.062074):0.064797,((Antipathes\_caribbeana:0.053896,Plumapathes\_pennacea:0.045453):0.108419,((Corynactis\_australis:0.066395,(Ricordea\_yuma:0.065164,((Amplexidiscus\_fenestrafer:0.053832,Rhodactis\_indosinensis:0.014454):0.01069,Discosoma\_spp:0.03959):0.040398):0.011425):0.028652,((((Acropora\_digitifera:0.057801,Acropora\_millepora:0.008548):0.005001,Acropora\_tenuis:0.021788):0.090596,(Astreopora\_spp:0.048089,(Porites\_lobata:0.002024,Porites\_australiensis:0.006654):0.048218):0.012396):0.014759,(((Lobactis\_scutaria:0.008759,Ctenactis\_echinata:0.006969):0.027244,(Montastraea\_cavernosa:0.024341,(Pseudodiploria\_strigosa:0.027266,(((Favia\_lizardensis:0.002407,Favia\_spp:0.10512):0.006193,Platygyra\_carnosus:0.16142):0.00544,Orbicella\_faveolata:0.03004):0.003371):0.005356):0.017989):0.011364,(Madracis\_auretenra:0.027553,(Pocillopora\_damicornis:0.035619,(Stylophora\_pistillata:0.030801,(Seriatopora\_hystrix:0.015099,Seriatopora\_spp:0.035017):0.010937):0.007567):0.034336):0.022343):0.017003):0.027674):0.07751):0.026608):0.020789,Protopalythoa\_variabilis:0.277447):0.063276,((Pennatula\_rubra:0.055175,Renilla\_reniformis:0.107084):0.072312,(Corallium\_rubrum:0.073471,(Nephthyigorgia\_spp:0.062532,((Gorgonia\_ventalina:0.051904,Leptogorgia\_sarmentosa:0.023995):0.013823,(Eunicella\_cavolinii:0.011368,Eunicella\_verrucosa:0.002738):0.025411):0.011668):0.049632):0.032084):0.188635):0.018531):0.036866,((Polypodium\_hydriforme:0.566426,(((Enteromyxum\_leei:0.648484,Kudoa\_iwatai:0.561639):0.680711,(Sphaeromyxa\_zaharoni:0.953538,(Henneguya:0.780518,((((Thelohanellus\_kitauei:0.307058,Myxobolus\_turpisrotundus:0.304868):0.180837,((Myxobolus\_wulii:0.22627,(Myxobolus\_honghuensis:0.015607,Myxobolus\_ampullicapsulatus:0.070408):0.17434):0.057403,Myxobolus\_pendula:0.245172):0.256762):0.21192,Myxobolus\_cerebralis:0.432823):0.219181,Myxobolus\_physophilus:0.60296):0.170284):0.129202):0.147595):0.130403,Tetracapsuloides\_bryosalmonae:0.80208):0.20069):0.142613,(((((Atolla\_vanhoeffeni:0.040317,Periphylla\_periphylla:0.065065):0.150184,((Aurelia\_aurita:0.079173,(Cassiopea\_xamachana:0.132511,(Stomolophus\_meleagris:0.05262,Rhopilema\_esculentum:0.037626):0.059135):0.020893):0.040423,((Cyanea\_capillata:0.013458,Cyanea\_nozakii:9.14E-4):0.112443,(Chrysaora\_fuscescens:0.054366,Pelagia\_noctiluca:0.058332):0.054913):0.032858):0.097924):0.032958,(Chironex\_fleckeri:0.098164,(Alatina\_alata:0.081036,Tripedalia\_cystophora:0.180489):0.02203):0.132863):0.024249,(Craterolophus\_convolvulus:0.100106,(Calvadosia\_campanulata:0.081895,(Haliclystus\_sanjuanensis:0.150622,Calvadosia\_cruxmelitensis:0.047193):0.020468):0.023085):0.204894):0.024348,((Liriope\_tetraphylla:0.144402,(Craspedacusta\_sowerbyi:0.163299,Aegina\_citrea:0.22533):0.02255):0.055681,((Ectopleura\_larynx:0.169364,(Hydra\_viridissima:0.041963,(Hydra\_oligactis:0.030741,Hydra\_vulgaris:0.028101):0.038597):0.092278):0.078655,(Millepora\_alcicornis:0.226375,((Turritopsis\_spp:0.15472,(Podocoryna\_carnea:0.102354,(Hydractinia\_polyclina:0.008327,Hydractinia\_symbiolongicarpus:0.007533):0.059407):0.051791):0.030935,(Physalia\_physalis:0.111875,((Nanomia\_bijuga:0.07543,Agalma\_elegans:0.042964):0.044384,(Craseoa\_lathetica:0.059619,Abylopsis\_tetragona:0.121102):0.037884):0.042855):0.075418):0.015359):0.016351):0.058833):0.111132):0.063606):0.055702):0.028066,((((Helobdella\_robusta:0.428146,Crassostrea\_gigas:0.25635):0.054388,Drosophila\_melanogaster:0.536965):0.050571,((Ciona\_intestinalis:0.459065,Taeniopygia\_guttata:0.299825):0.052585,Strongylocentrotus\_purpuratus:0.349709):0.026884):0.077267,(Trichoplax\_adhaerens:0.483558,(Amphimedon\_queenslandica:0.550768,Mnemiopsis\_leidyi:0.725216):0.059918):0.045348):0.028066);

end;

begin r8s;

blformat lengths=persite nsites=51598 ulrametric=no;

collapse;

MRCA OUT Helobdella\_robusta Mnemiopsis\_leidyi Trichoplax\_adhaerens;

MRCA MED Chironex\_fleckeri Agalma\_elegans Thelohanellus\_kitauei;

MRCA HEX Protopalythoa\_variabilis Edwarsiella\_lineata Lobactis\_scutaria;

MRCA HYD Liriope\_tetraphylla Agalma\_elegans Hydra\_vulgaris;

MRCA CNI Protopalythoa\_variabilis Chironex\_fleckeri Thelohanellus\_kitauei;

prune taxon=OUT;

fixage taxon=HYD age=500;

constrain taxon=MED min\_age=570;

constrain taxon=HEX min\_age=540;

constrain taxon=CNI max\_age=741;

set minRateFactor=0.000001;

set minDurFactor=0.000001;

set smoothing=1e-06;

set num\_time\_guesses=3;

set checkGradient=yes;

divtime method=PL algorithm=TN;

showage;

describe plot=cladogram;

describe plot=chrono\_description;

end;

\# Run r8s

r8s -b -f r8s\_in\_true\_run.txt > true\_run\_1

\# Extract the time tree from the result file

perl -e 'while (<>) { if (m/tree tree\_1 = (.\*);/) { my $tree = $1; $tree =~ s/(MED|HEX|HYD|CNI)//g; print "$tree\n"; } }' true\_run\_1 > tree.txt

**8. Selection pressure analysis using Hyphy**

\# Make configuration file

generate\_busted\_hyphy.sh alignment.phy tree.nwk > script.bf

HYPHYMPI script.bf

\# The content of generate\_busted\_hyphy.sh is as follows:

inputRedirect = {};

inputRedirect["01"]="Universal"; // genetic code

inputRedirect["02"]="/home/gqx/proteome/16\_selection\_analysis/disintergrin/Disintergin\_nucl.aln"; // codon data

inputRedirect["03"]="/home/gqx/proteome/16\_selection\_analysis/disintergrin/1.nwk"; // tree

inputRedirect["04"]="All"; // Test for selection on a branch

inputRedirect["05"]=""; // complete selection

ExecuteAFile ("/usr/local/lib/hyphy/TemplateBatchFiles/SelectionAnalyses/BUSTED.bf", inputRedirect);

\# If you want to use multithreading

mpirun -np 12 /usr/local/bin/HYPHYMPI ./script.bf

**9. Orthologous gene family analysis in comparative proteomics and visualization using UpSetR**

\# For the orthologous gene family analysis by OrthoMCL, please see section 1.

\# Provide matrix format that is suitable for UpSetR

cat aaur\_OG avir\_OG cfle\_OG hvul\_OG mhon\_OG mwul\_OG tkit\_OG > OG

remove\_duplicate.pl OG

mv duplicate\_remove gene

produce\_matrix\_for\_UpSet.pl gene aaur\_OG

\# Then combine a csv file that contains the gene and species information, eight\_matrix\_strict.csv

\# The content of Rscript upset\_eight.R is as follows

library(UpSetR)

data <- read.csv("eight\_matrix.csv",header=T)

pdf("out.pdf",  width=11, height=6)

upset(data, nsets = 8, order.by="freq", sets.bar.color = "#56B4E9", mainbar.y.label = "OG Intersections", sets.x.label = "ProteinsPer Species",nintersects = 60)

dev.off()

**10. Functional comparison in comparative proteomics and visualization using UpSetR**

**10.1 Get the Interpro number for each proteins**

\# Use Blast2GO to do Interpro annotation, collect the Interpro number (IPR) and save it as list file.

\# Then run the bash file: interpro\_ID\_converter.sh, as follows:

grep "IPR" list  > out

split\_lines\_by\_semicolon.pl out

grep "IPR" out2 > out3

remove\_duplicate.pl out3

perl -p -i -e 's/\s+(\S+)$//g' duplicate\_remove

mv duplicate\_remove duplicate\_remove1

remove\_duplicate.pl duplicate\_remove1

rm duplicate\_remove1

mv duplicate\_remove final\_IPR

perl -p -i -e 's/^\s//g' final\_IPR

remove\_duplicate.pl final\_IPR

rm final\_IPR

mv duplicate\_remove final\_IPR

convert\_interpro\_2\_name.pl final\_IPR entry.list

rm out out2 out3

cat convert\_name | cut -f 3 > domain\_name

\# The final\_IPR contains the information for visualization

**10.2 Visualize the results using UpSetR**

\# Change the name of final\_IPR for each species

\# Use *Chironex* as an example

cp ../5\_Chironex\_domain/final\_IPR ./Chironex

cat Anemonia  Aurelia  Chironex  Hydra  honghuensis  kitauei  shata wulii > OG

remove\_duplicate.pl OG

mv duplicate\_remove gene

\# For each species

produce\_matrix\_for\_UpSet.pl gene Anemonia

rename \_matrix .xls \*matrix

\# Combine gene and OG\_matrix into a single csv: eight\_matrix\_strict.csv

Rscript upset\_eight.R

**11. Circos analysis and visualization**

**11.1 Install circos**

http://circos.ca/distribution/circos-0.69-6.tgz

tar zxf /opt/software/circos-0.69-6.tgz

cd /opt/biosoft/circos-0.69-6/bin

./list.modules

./list.modules | perl -p -e 's/\n/ /' | perl -p -e 's/^/sudo cpan -i /; s/$/\n/' | sh

./circos -modules | perl -e 'while (<>) { print "$1 " if m/missing\s+(\S+)/ }' | perl -p -e 's/^/sudo cpan -i /; s/$/\n/' | sh     # Install missing modules

#Install GD

cd /opt/sysoft/

tar Jxf  /opt/software/libgd-2.1.1.tar.xz

./configure --prefix=/usr/

make -j 4

sudo make install

gdlib-config

ln -s /usr/lib64/libgd.so.2.0.0 /usr/lib64/libgd.so.3

wget  http://circos.ca/distribution/lib/GD-2.53.tar.gz

tar  zxf /opt/software/GD-2.53.tar.gz

perl Makefile.PL

make

make install

./circos -modules

./gddiag

./circos --help

echo 'PATH=$PATH:/opt/biosoft/circos-0.69-6/bin/' >> ~/.bashr

source ~/.bashrc

\# Dowmload circos\_tools

wget http://circos.ca/distribution/circos-tools-0.22.tgz

tar zxf /opt/software/circos-tools-0.22.tgz

cd /opt/biosoft/circos-tools-0.22/tools/binlinks

./run

**11.2 Install Circoletto**

wget https://github.com/infspiredBAT/Circoletto/archive/master.zip

unzip master.zip

cd Circoletto-master

\# If you are dealing with more than 200 sequences, change the max\_ideograms parameter

vim /opt/biosoft/circos-0.69-6/etc/housekeeping.conf

\# Also change circoletto.pl, feed it with the path of other softwares

$path2circos      = '/opt/biosoft/circos-0.69-6';     #  path to Circos

$path2circostools = '/opt/biosoft/circos-tools-0.22'; #  path to Circos utilities

**11.3 Install blast-2.2.25**

wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy/2.2.25/blast-2.2.25-x64-linux.tar.gz

echo "PATH=$PATH:/opt/biosoft/blast-2.2.25/bin" >> ~/.bashrc

source ~/.bashrc

**11.4 Run Circoletto**

\# Use analysis between the protein of *Myxobolus honghuensis* and other six species as an example

cat aaur.fasta avir.fasta cfle.fasta cfus.fasta cnoz.fasta hvul.fasta > rest.fa

\# Circoletto can not deal with asterisks.

perl -p -i -e "s/\\*//g" rest.fa

\# Change the parameter in Circoletto.pl

$max\_sequences = 10000

circoletto.pl --cpus 12 --query mhon.fasta --database rest.fa --out\_type svg --e\_value 1e-20 --score2colour eval --scoreratio2colour max  --ribocolours2allow '(blue|green|orange|red)' --z\_by score --no\_labels

\# Use Perl to process link file to adjust the link details

perl -p -i -e 's/stroke\_color=(\S+?),//g' ./\*.blasted.links

perl -p -i -e 's/stroke\_thickness=(\d+?),//g' ./\*.blasted.links

perl -p -i -e 's/\_a3/\_a2/g' ./\*.blasted.links

perl -p -i -e 's/thickness=(\d+)/thickness=4/g' ./\*.blasted.links

paste -d ',' ./\*.blasted.links ff > new

rm ./\*.links

mv new  cl0011919395.blasted.links

\# Use Perl to process karyotype file to adjust the color

perl -p -i -e 's/(\S+)$/dgreen/ if ($\_ =~ /chr/ and $\_ =~ /aaur/)' ./\*.blasted.karyotype

perl -p -i -e 's/(\S+)$/lblue/ if ($\_ =~ /chr/ and $\_ =~ /avir/)' ./\*.blasted.karyotype

perl -p -i -e 's/(\S+)$/dred/ if ($\_ =~ /chr/ and $\_ =~ /cfle/)' ./\*.blasted.karyotype

perl -p -i -e 's/(\S+)$/dyellow/ if ($\_ =~ /chr/ and $\_ =~ /cfus/)' ./\*.blasted.karyotype

perl -p -i -e 's/(\S+)$/black/ if ($\_ =~ /chr/ and $\_ =~ /hvul/)' ./\*.blasted.karyotype

perl -p -i -e 's/(\S+)$/dblue/ if ($\_ =~ /chr/ and $\_ =~ /mhon/)' ./\*.blasted.karyotype

perl -p -i -e 's/(\S+)$/lpurple/ if ($\_ =~ /chr/ and $\_ =~ /cnoz/)' ./\*.blasted.karyotype

perl -p -i -e 's/(\S+)$/dpurple/ if ($\_ =~ /chr/ and $\_ =~ /tkit/)' ./\*.blasted.karyotype

perl -p -i -e 's/(\S+)$/dorange/ if ($\_ =~ /chr/ and $\_ =~ /mwul/)' ./\*.blasted.karyotype

perl -p -i -e 's/(\S+)$/dyellow/ if ($\_ =~ /chr/ and $\_ =~ /csha/)' ./\*.blasted.karyotype

\# Edit the conf file as follows:

<ideogram>

<spacing>

default             = 0.0003r

break               = 10u

axis\_break\_at\_edge  = no

axis\_break          = no

axis\_break\_style    = 1

</spacing>

thickness           = 35p

fill                = yes

radius              = 0.66r

show\_label          = no

label\_with\_tag      = yes

label\_font          = condensed

label\_color         = grey

label\_radius        = dims(ideogram,radius) + 0.104r

label\_size          =

band\_stroke\_thickness = 0

show\_bands          = no

fill\_bands          = yes

</ideogram>

show\_ticks          = no

show\_tick\_labels    = no

<ticks>

radius              = dims(ideogram,radius\_outer)

multiplier          = 1e-6

<tick>

spacing             = 100u

rspacing            = 0.1

spacing\_type        = relative

size                = 5p

thickness           = 1p

color               = vvvvdgrey

show\_label          = no

label\_size          = 10p

label\_offset        = 5p

format              = %d

</tick>

</ticks>

karyotype           = cl0006216655.blasted.karyotype

<image>

dir  = .

file                = cl0006216655.blasted.svg

radius              = 1000p

background          = white

angle\_offset        = -90

24bit               = yes

auto\_alpha\_colors   = yes

auto\_alpha\_steps    = 5

</image>

chromosomes\_units   = 1

chromosomes\_display\_default = yes

<links>

ribbon              = no

<link>

show                = yes

file                = cl0006216655.blasted.links

record\_limit        = 2000

</link>

</links>

\# Run circus

circos -conf ./\*.conf -noparanoid

**12. Identify ortholog gene manually**

\# RBBH is mainly used to find orthologs, while the common blast includes all homologs, and RBBH can distinguish orthologs from paralogs. The NC method is actually an all-against-all blast, which includes paralogs, so it is generally used for homology identification. Note: homology does not necessarily mean ortholog. homology is more stringent.

**12.1 Pairwise Reciprocal Best BLAST Hits (RBBH)**

blast\_rbbh.py fish.fasta human.fasta -o fish\_human.tab -t blastp -a prot -i 1 -c 1 --threads 48

\# If you want to do RBBH for three species

cat Avir\_Aaur.tab | cut -f 1,2 > 22

cat Avir\_Hydra.tab | cut -f 1,2 > 33

cat Aaur\_Hydra.tab | cut -f 1,2 > 44

run\_RBBH.pl ../Anemonia\_viridis\_Rachamim\_2014.fasta ../Aurelia\_aurita\_Rachamim\_2014.fasta Anemonia Aurelia

proteome\_RBBH\_get\_three\_inter.pl Anemonia\_Aurelia Anemonia\_Hydra Aurelia\_Hydra

**12.2 Neighborhood correlation (NC) method**

\# Compare Anemonia\_viridis\_Rachamim\_2014.fasta and Hydra\_vulgaris\_Rachamim\_2014.fasta

run\_neighborhood\_correlation.pl Anemonia\_viridis\_Rachamim\_2014.fasta Hydra\_vulgaris\_Rachamim\_2014.fasta Anemonia Hydra 0.6 

proteome\_NC\_result\_process.pl out

NC\_thresh.pl NC\_result 0.8

cat out2 | cut -d ' ' -f 1 > 1

remove\_duplicate.pl 1

wc -l duplicate\_remove

mv out2 NC\_0.8

**Please cite this paper as**

Guo, Q., Atkinson, S. D., Xiao, B., Zhai, Y., Bartholomew, J. L., & Gu, Z. (2022). A myxozoan genome reveals mosaic evolution in a parasitic cnidarian. BMC biology, 20(1), 1-19. 

**License**

All source code, i.e. scripts/\*.pl, scripts/\*.sh or scripts/\*.py are under the MIT license.




