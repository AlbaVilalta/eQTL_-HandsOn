## Author: Alba Vilalta
## Date: 23/11/2018
## eQTL Hands-On

sudo setxkbmap -layout es
sudo docker run -v $PWD:$PWD -w $PWD -it dgarrimar/eqtlmapping
PATH=$PATH:$PWD/bin

####################
## Task 1
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf{.gz,.gz.tbi} --directory-prefix input/unprocessed/1000g #em descarrego el genotpi i l'índex .tbi corresponent

###################
## Task 2
# Get GEUVADIS samples from the metadata
cut -f1 input/unprocessed/geuvadis/geuvadis.metadata.txt | sed '1d' | sort | uniq > tmp/geuvadis.samples.txt
# Subset the VCF (common samples, biallelic SNPs and indels, MAF >= 0.05, no duplicates)
bcftools view -v snps,indels -m 2 -M 2 -q 0.05:minor -S tmp/geuvadis.samples.txt -Ob input/unprocessed/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | bcftools norm -d all -Oz -o tmp/genotypes.chr22.vcf.gz
# Subset the VCF so that there are at least 10 individuals per genotype group and compress it (for indexing we require 'bgzip' compression)
filter.genotype.py -t 10 -g <(zcat tmp/genotypes.chr22.vcf.gz) | bgzip > input/processed/genotypes.chr22.vcf.gz
# Index the VCF
tabix -p vcf input/processed/genotypes.chr22.vcf.gz

# Q1: What do the bcftools options employed mean?
# view -v snps -- to only view biallelic snps
# indels -m2 -M2 -9 0.05:minor -S -- to exclude samples with MAF min 5
# -Ob -- to extract genotypes
# nomr -d -- to nomalize indels (-d)
# filter.genotype.py -t 10 -g -- to select those variants with more than 10 individuals
# tabix -p -- to index the VVF file
# Q2: How many variants do you get in input/processed/genotypes.chr22.vcf.gz?
zcat input/processed/genotypes.chr22.vcf.gz | grep -v "#" | wc -l
# 74656
# Q3: How many samples do you have before and after subsetting?
# before = 2504
bcftools stats input/unprocessed/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > tmp/stats.before 
less -S tmp/stats.before  answer: 2504
# after = 445
bcftools stats input/processed/genotypes.chr22.vcf.gz > tmp/stats.after 
less -S tmp/stats.after  answer: 445

###################
## Task3
# Q1: Which version of GENCODE is GEUVADIS using?
# Utilitza la versió 12 de Gencode
# Q2: To which genome assembly does this annotation correspond?
# Correspon a GRCh37
# Q3: How many protein coding genes are annotated in the last version (v29)?
# 83129 coding proteins genes
# Q4:Which command do you use to do this?
PATH=$PATH:$PWD/bin

release= 12
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$release/gencode.v$release.annotation.gtf.gz
mv gencode.v$release.annotation.gtf.gz input/unprocessed/gencode/gencode.annotation.gtf.gz
zcat input/unprocessed/gencode/gencode.annotation.gtf.gz | grep "gene_type \"protein_coding\"\|gene_type \"lincRNA\"" | gtf2bed.sh > tmp/gencode.annotation.bed
head tmp/gencode.annotation.bed

# Q5: But how to get the TSS positions and the gene lengths from it?
# TSS position: 
# +upstream --> 2n column
# -downstream --> 3rd column
# Lengths: hem de restar la tercera i la segona columna
# Q6: to which BED coordinates would correspond the GTF coordinates chr1 10 20? Why?
# BED --> chr1 9 20, perquè l'index és: 0-index
# Q7: Why do we need to use tmpfile below?
# Perquè no podem llegir i escriure alhora al mateix fitxer

awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$3-$2,$6}' tmp/gencode.annotation.bed > tmpfile; mv tmpfile tmp/gencode.annotation.bed
awk 'BEGIN{OFS="\t"}{if($6=="+"){print $1,$2,$2+1,$4,$5,$6}else{print $1,$3-1,$3,$4,$5,$6}}' tmp/gencode.annotation.bed > tmpfile; mv tmpfile tmp/gencode.annotation.bed
sed -i "s/^chr//" tmp/gencode.annotation.bed

###################
## Task 4
join -1 4 -2 1 -t $'\t' <(sort -k4,4 tmp/gencode.annotation.bed) <(zcat input/unprocessed/geuvadis/geuvadis.gene_expr_rpkm.tsv.gz | sort -k1,1) > tmp/joint.tsv
awk '$2==22' tmp/joint.tsv > tmp/joint.chr22.tsv
paste <(awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$5,$6}' tmp/joint.chr22.tsv) <(cut -f1-6 --complement tmp/joint.chr22.tsv) | sort -k1,1V -k2,2n > tmp/joint.chr22.bed
cat <(zcat input/unprocessed/geuvadis/geuvadis.gene_expr_rpkm.tsv.gz | head -1 | sed "s/TargetID/#chr\tstart\tend\tgene\tlength\tstrand/") tmp/joint.chr22.bed > tmp/genes.chr22.rpkm.bed

# Q1: Of all genes considered, which have lower expression levels, protein-coding or lincRNA?
# lincRNA
# Q2: Why do we need gene expression to be normal?
# Perquè cada gen tingui una distribució propera a la normal / estàndard per tal de poder comaprar-los millor
# Q3: How would you check that quantile normalization worked? 
# A partir d'un gràfic
# Q4: and that gene expression of a gene follows a normal distribution?
# Cada gen té una distribució similar a l'estàndard normal

normalize.R -i tmp/genes.chr22.rpkm.bed -o tmp/genes.chr22.norm.bed
bgzip tmp/genes.chr22.norm.bed
tabix -p bed tmp/genes.chr22.norm.bed.gz
mv tmp/genes.chr22.norm.bed.gz* input/processed

###################
## Task 5
# Before:
check.norm.R -i input/processed/genes.chr22.norm.bed.gz -o result/plots/check.norm.pdf
# After:
check.norm.R -i input/processed/genes.chr22.norm.bed.gz -o result/plots/check.norm.pdf
# Q1: What can you see?
# Es generen dos gràfics per cada un dels arxius. En el gràfic generat a partir de les dades normalitzades s'oberva que el rang dels valors d'expressió va de 3 a -3 i en el gràfic Q-Q plot s'observa una recta perfecta amb tots els valors. En canvi, en els gràfics amb les dades no normalitzades no s'observa un rang gaire específic, ja que hi ha bastanta dispersió pel que fa als valors d'expressió i en el Q-Q plot, els quartils teòrics són bastant diferents les mostres.  

###################
## Task 6
# Q1: Which ones would you select?
# geuvadis.metadata.txt i 1000g.phase3_metadata.txt

head -1 input/unprocessed/geuvadis/geuvadis.metadata.txt | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++) {print i, $i}}'
head -1 input/unprocessed/1000g/1000g.phase3_metadata.txt | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++) {print i, $i}}'

###################
## Task 7
QTLtools pca --bed input/processed/genes.chr22.norm.bed.gz --scale --center --out result/expression
QTLtools pca --vcf input/processed/genotypes.chr22.vcf.gz --scale --center --maf 0.05 --distance 50000 --out result/genotypes

# Q1: What do the parameters employed mean?
# QTLtools pca -- to perform Principal Component Analysis (PCA) on genotype or molecular phenotype data.
# Q2: Which information do the output files contain?
# expression.pca i genotypes.pca -- contain the individual coordinates on the Principal Components (PCs)
# expression.pca_stats i genotypes.pca_stats -- contain the percentages of the variance explained by each PC

pcaPlot.R -i result/expression -o result/plots/expression.pca.pdf
pcaPlot.R -i result/genotypes -o result/plots/genotypes.pca.pdf

# Q3: What can you observe in the plots?
# expression plot -- veiem que el PCs dels diferents grups presenten una distribució molt semblant, ja que es troben tots focalitzats a la mateixa zona i barrejats. 
# genotypes plot -- veiem que els PCs els troben agrupats en dos grups clarament distingibles.

pcaPlot.R -i result/genotypes --metadata input/unprocessed/1000g/1000g.phase3_metadata.txt --color super_pop --out result/plots/genotypes.pca.super_pop.pdf
pcaPlot.R -i result/expression --metadata input/unprocessed/geuvadis/geuvadis.metadata.txt --color super_pop --out result/plots/expression.pca.super_pop.pdf

# Q4: With this information, which covariates seem more relevant to explain the variability in the data?
# Aportem color a les variables del gràfic i ens indica que un grup és l'europeu i l'altre l'africà

# Generate a common metadata with info about the population, gender and laboratory.
join -j 1 -t $'\t' <(sort -k1,1 input/unprocessed/1000g/1000g.phase3_metadata.txt) <(cut -f1,20 input/unprocessed/geuvadis/geuvadis.metadata.txt | sort -k1,1 | uniq) > tmp/metadata.txt
# Set names for the new metadata
sed -i '1s/^/sampleID\tpop\tsuper_pop\tgender\tlab\n/' tmp/metadata.txt
# Build a linear model and plot the contribution of each factor in the metadata to the total variance
var_partition.R -i input/processed/genes.chr22.norm.bed.gz -m tmp/metadata.txt --formula "~ (1|gender) + (1|pop) + (1|lab)" -o result/plots/vp.pdf

# Q5: Which are the factors that explain more variance?
# lab and pop

###################
## Task 8
# Compute 10 PEER factors
peer.R -i input/processed/genes.chr22.norm.bed.gz -p 10 -o tmp/peer.tsv
# Check how much variance do the first 5 PEER explain in comparison with the known factors
var_partition.R -i input/processed/genes.chr22.norm.bed.gz -m <(paste tmp/peer.tsv tmp/metadata.txt) -f "~ (1|pop) + (1|lab) + PEER1 + PEER2 + PEER3 + PEER4 + PEER5" -o result/plots/vp.peer.pdf

# Q1: How much variance do they explain? On average is it more or less than the explained by the known factors?
# Expliquen fins a un 80% (aproximat) de variança. És més que l'explicada pels factors coneguts, ja que no arribaba al 70%. 

# 'Rscript -e' is just a trick to run an R script without opening an interactive R session in the console. ;)
join -j 1 -t $'\t' tmp/metadata.txt tmp/peer.tsv  | Rscript -e 'write.table(t(read.table(file("stdin", open = "r", blocking = T), h = F)), file = "input/processed/covariates.tsv", quote = F, sep = "\t", col.names = F, row.names = F)'
# Compress it
gzip input/processed/covariates.tsv
# 1.4. cis eQTL mapping (nominal pass)
QTLtools cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --nominal 0.01 --out result/nominals.txt

###################
## Task 9
# Which information contains each field in the ouptut file?
# 1. The phenotype ID
# A cada columna s'explica:
# 2. The chromosome ID of the phenotype
# 3. The start position of the phenotype
# 4. The end position of the phenotype
# 5. The strand orientation of the phenotype
# 6. The total number of variants tested in cis
# 7. The distance between the phenotype and the tested variant (accounting for strand orientation)
# 8. The ID of the tested variant
# 9. The chromosome ID of the variant
# 10. The start position of the variant
# 11. The end position of the variant
# 12. The nominal P-value of association between the variant and the phenotype
# 13. The corresponding regression slope
# 14. A binary flag equal to 1 is the variant is the top variant in cis
# Q1: Are there pairs genotype-phenotype with exactly the same p-value and effect size (β)? How is this possible?
# Les parelles no tenen ni el mateix p-value ni la mateixa effect size. Això és possible PERQUEEEEEEEE¿¿¿
# Q2: What do you observe?
pvdist.R -i result/nominals.txt --col 12 -o result/plots/pvdist.pdf
# Q3: Which SNPs did you select? What do you observe?
# SNPs selected: rs5746938 rs36084991
plink --ld rs5746938 rs36084991 --vcf input/processed/genotypes.chr22.vcf.gz --out tmp/ld2
#  R-sq = 1              D' = 1

#   Haplotype     Frequency    Expectation under LE
#   ---------     ---------    --------------------
#         AGT      0.246067                0.060549
#         TGT     -0                       0.185518
#          AG     -0                       0.185518
#          TG      0.753933                0.568414
#   In phase alleles are AGT/TG

QTLtools cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --permute 1000 --out result/permutations.txt

###################
## Task 10
for j in $(seq 1 16); do
  echo "cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --permute 1000 --chunk $j 16 --out result/permutations_$j.txt"
done | xargs -P4 -n14 QTLtools
cat result/permutations_*.txt > result/permutations.txt; rm result/permutations_*.txt 
R
p <- read.table("result/permutations.txt")                                                      # Read input file
pdf("result/plots/pv-correlation.pdf",  paper = 'a4r', width = 9, height = 6)                   # Open PDF device
plot(p[, 18], p[, 19], xlab = "pv (perm)", ylab = "pv (beta)")                                  # Plot p-values
abline(0, 1, col = "red")                                                                       # Add red line 1=1
plot(-log10(p[, 18]), -log10(p[, 19]), xlab = "-log10 pv (perm)", ylab = "-log10 pv (beta)")    # Repeat in -log10 space to check the behaviour of the small p-values.
abline(0, 1, col = "red")
dev.off()                                                                                       # Close device
quit("no")

###################
## Task 11
mtc.R -n result/nominals.txt -p result/permutations.txt --method 'bonferroni' --alpha 0.05 --out tmp/bonferroni.txt
mtc.R -n result/nominals.txt -p result/permutations.txt --method 'fdr' --alpha 0.05 --out tmp/fdr.txt
mtc.R -n result/nominals.txt -p result/permutations.txt --method 'perm-fdr' --alpha 0.05 --out result/eqtls.tsv
# Global empirical P-value threshold = 1.20e-02

# Q1: How many significant eQTLs do we find in each case in comparison with the nominal pass?
awk '{if ($12 < 0.05) print $12}' result/nominals.txt | wc -l
# 95206
awk '{if ($12 < 0.05) print $12}' tmp/FDR.txt | wc -l
# 17267
awk '{if ($12 < 0.05) print $12}' result/eqtls.tsv | wc -l
# 12044
awk '{if ($12 < 0.05) print $12}' tmp/bonferroni.txt | wc -l
# 6007

###################
## Task 12
eQTLviewer.R -i <(head -n 10 result/eqtls.tsv) -g input/processed/genotypes.chr22.vcf.gz -e input/processed/genes.chr22.norm.bed.gz -o result/plots/eQTLs_head.pdf --verbose

###################
## Task 13
# Download from ftp server
rsync -av rsync://ftp.ensembl.org/ensembl/pub/grch37/release-86/regulation/homo_sapiens/AnnotatedFeatures.gff.gz input/unprocessed/ensembl
# Get chr, start, end and feature name in BED format
zcat input/unprocessed/ensembl/AnnotatedFeatures.gff.gz | awk 'BEGIN{FS=OFS="\t"}{print $1, $4-1, $5, $9}' | sed -r 's/Name=([^;]+);.*/\1/' | grep -v '^GL' | sort -V > tmp/ERB.bed
# Merge overlapping features of the same type 
# e.g. chr1 100 200 feat1            chr1 100 300 feat1
#      chr1 150 300 feat1     =>     chr1 100 250 feat2
#      chr1 100 250 feat2
for feat in $(cut -f4 tmp/ERB.bed | sort | uniq); do
 bedtools merge -i <(grep -Fw $feat tmp/ERB.bed) -c 4 -o distinct
done > input/processed/ERB.collapsed.bed
# Remove 'chr' from chromosome names (-i option to modify the file 'in place')
sed -i "s/^chr//" input/processed/ERB.collapsed.bed
for feat in $(cut -f4 input/processed/ERB.collapsed.bed | sort | uniq); do
 QTLtools fenrich --qtl <(sed '1d' result/eqtls.tsv | awk '{print $9, $10-1, $10, $8, $1, "."}') --tss tmp/gencode.annotation.bed  --bed <(grep -Fw $feat input/processed/ERB.collapsed.bed) --out tmp/enrich.txt > /dev/null; echo "$(cat tmp/enrich.txt) $feat"
done | grep -Fwv inf | grep -Fwv nan > result/enrichments.txt
plot.enrich.R -i result/enrichments.txt -o result/plots/enrich.pdf

# Q1: Which are the top enriched features? Which kind of factors are they? 
# H3K36me3 i PolII. El primer és un factor que es diposita a les histones a mesura que són desplaçades per les ARN polimerasa II durant la transcripció i el segon és l'ARN polimerassa II, un enzim que catalitza la transcripció de l'ADN.
# Q2: What does an odds ratio lower than one mean?
# Un ratio menor de odds significa una menor possiblitat que tingui lloc. 

###################
## Task 14
sed '1d' result/eqtls.tsv | cut -f8 | sort | uniq > tmp/eqtls_snps.tsv

# Respostes resoltes a partit del link de l'Ari: http://grch37.ensembl.org/Homo_sapiens/Tools/VEP/Results?db=core;tl=315nTgfC4grzij8q-4700617
# Q1: Which kind of consequences have they, according to the VEP? In which proportion? 
#	intron_variant: 47%
#	upstream_gene_variant: 15%
#	downstream_gene_variant: 14%
#	non_coding_transcript_variant: 11%
#	NMD_transcript_variant: 6%
#	regulatory_region_variant: 2%
#	intergenic_variant: 2%
#	non_coding_transcript_exon_variant: 1%
#	3_prime_UTR_variant: 1%
#	Others

# Q2: How many eQTLs are high impact variants? Which consequences are related to those high impact variants?
wc -l result/IMPACT_is_HIGH.txt 
# 24
uniq <(cut -f4 result/IMPACT_is_HIGH.txt)
#	Consequence
#	splice_acceptor_variant,non_coding_transcript_variant
#	stop_gained
#	splice_donor_variant,frameshift_variant
#	frameshift_variant
#	stop_gained
#	splice_acceptor_variant
#	splice_acceptor_variant,NMD_transcript_variant
#	splice_acceptor_variant
#	splice_acceptor_variant,non_coding_transcript_variant
#	stop_gained
# Q3: Out of all high impact variants, how many of them are falling in acceptor splice sites of protein coding genes? 
wc -l <(awk '$4=="splice_acceptor_variant"' result/IMPACT_is_HIGH.txt) 
# 3 /dev/fd/63

###################
## Task 15
# Generate a list of sGenes
cut -f1 result/eqtls.tsv | sed '1d' | sed 's/\..\+//' | sort | uniq > tmp/egenes.txt
# We will use as background all the genes (PC and lincRNA) in chr22
awk '{if($1==22) print $4}' tmp/gencode.annotation.bed | sed 's/\..\+//' | sort | uniq > tmp/bg.txt

# Q1: In which biological processes are your eGenes enriched? Which molecular functions and components correspond to those processes?
#	Pocesses: response to lipopolysaccharide i response to molecule of bacterial origin
#	Molecular functions: Ras guanyl-nucleotide exchange factor activity
#	Components: endoplasmic reticulum

###################
## Task 16
# Generate input files for QTLtools rtc
grep -Fwf <(cut -f1 result/eqtls.tsv ) result/permutations.txt > tmp/rtc_input
cut -f4,7 input/unprocessed/gwas/gwas.catalog.hg19.bed > tmp/gwas_trait
# Download the file 'hotspots_b37_hg19.bed' from QTLtools website
wget http://jungle.unige.ch/QTLtools_examples/hotspots_b37_hg19.bed --directory-prefix tmp
# Remove 'chr' from chromosome names (-i option to modify the file 'in place')
sed -i 's/^chr//' tmp/hotspots_b37_hg19.bed
# Run RTC
QTLtools rtc --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --hotspot tmp/hotspots_b37_hg19.bed --gwas-cis tmp/gwas_trait tmp/rtc_input --out result/rtc.txt

# Q1: How many pairs of variants have a RTC value above 0.9? 
awk '{if ($20 > 0.9) print $20}' result/rtc.txt | wc -l
#	39 pairs of variants
# Q2: For each pair, we have a GWAS hit and an eQTL. Find one example so that the gene to which the eQTL is associated is relevant for the trait/disease to which the GWAS variant is associated. Explore the literature and the biological databases that you know to gather more information. 
awk '{if ($20 > 0.99) print $1,$2,$3,$20}' result/rtc.txt
other_variant our_variant RTC
rs5750673 rs55850024 0.994118
rs909685 rs909685 1
rs2069235 rs909685 0.995885
rs4822024 rs4820438 0.994118
rs9611565 rs4820438 0.99902
rs2234052 rs9607799 0.99902
rs138177673 rs47341 0.994186

# Escullo el primer
grep rs5750673 input/unprocessed/gwas/gwas.catalog.hg19.bed
22	39110123	39110124	rs5750673	0	+	Hand grip strength
# Q3: Which consequences, according to the variant effect predictor, do these co-localized eQTL variants have?
# intron_variant: 75%
# NMD_transcript_variant: 13%
# non_coding_transcript_variant: 

###################
## Task 17
# Generate the ID/Z-scores input file. Select your favourite gene (e.g. gene=ENS00000000000.0).
# Set k (number of variants) to 50
gene=ENSG00000100211.6 #escullo el gen de l'exemple anterior
compZscore.R --gene $gene --nominal result/nominals.txt -k 50 --output tmp/$gene.rs_z
# Generate the LD matrix 
plink --r square --snps $(cut -f1 tmp/$gene.rs_z) --vcf input/processed/genotypes.chr22.vcf.gz --out tmp/$gene
CAVIAR -z tmp/$gene.rs_z -l tmp/$gene.ld -o result/$gene
head result/ENSG00000100211.6.log 

# Q1: How many variants are there in the credible (ρ=0.95) set? For each of these variants, which is the probability to be causal?
wc -l result/ENSG00000100211.6_set
# 5
cat result/ENSG00000100211.6_set
# rs55850024
# rs3827356
# rs5750673
# rs5757257
# rs4820346
# Q2: Which are the p-values and effect sizes of these variants? How are they in comparison to the p-values and effect sizes of other variants tested for the same gene? 
grep -f result/ENSG00000100211.6_set result/ENSG00000100211.6_post
# rs55850024	0.357026	0.714052
# rs3827356	0.0634794	0.126959
# rs5750673	0.0282665	0.0565329
# rs5757257	0.0282665	0.0565329
# rs4820346	0.5	1
# Q4: Which consequences, according to the variant effect predictor, do these variants have?
# --> rs55850024
# downstream_gene_variant: 56%
# intron_variant: 22%
# NMD_transcript_variant: 11%
# upstream_gene_variant: 11%
# --> rs3827356
# intron_variant: 52%
# non_coding_transcript_variant: 33%
# downstream_gene_variant: 10%
# upstream_gene_variant: 5%
# --> rs5750673
# intron_variant: 75%
# NMD_transcript_variant: 13%
# non_coding_transcript_variant: 12%
# --> rs5757257
# ntron_variant: 75%
# non_coding_transcript_variant: 13%
# NMD_transcript_variant: 13%
# --> rs4820346
# upstream_gene_variant: 38%
# downstream_gene_variant: 25%
# intron_variant: 25%
# non_coding_transcript_exon_variant: 12%

###################
## Task 18
# Define the gene corresponding to the co-localized or fine-mapped variants of interest
gene=ENSG00000100211.6
cat <(echo "MarkerName P.value") <(grep $gene result/nominals.txt | cut -d " " -f8,12) > tmp/metal.$gene




