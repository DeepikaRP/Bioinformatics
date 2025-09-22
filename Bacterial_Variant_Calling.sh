###### RAW DATA ######

#Retrive reference fastq file from ENA browser
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR396/SRR396537/SRR396537.fastq.gz

#assemble reference file using SPAdes
spades.py -o contigs.fasta --only-assembler -k 33 --cov-cutoff auto -s SRR396537.fastq

#Quality check
fastqc *.fastq 

#Merge fastQ files
cat *.fastq > merged.fastq



###### MAPPING ######

#Index fasta file
bwa index ref_cotigs.fasta

#Align reads
bwa mem ref_contigs.fasta merged.fastq > map.sam

#SAM to BAM
samtools view -S -b map.sam > map.bam

#Sort BAM file
samtools sort map.bam > sorted_map.bam

#Check BAM file reads
samtools flagstat sorted_map.bam

#Index BAM file
samtools index sorted_map.bam



###### VARIANT CALLING ######

#Variant calling using BCFTools
bcftools mpileup -f ref_contigs.fasta sorted_map.bam | bcftools call -mv -Ov -o calls.vcf

#Check VCF file
cat calls.vcf | head -n 10

#Filter SNPs
bcftools view -v snps map/calls.vcf -Ov -o snps.vcf

#Check VCF file 
cat snps.vcf | head -n 10



###### SnpEff ANNOTATION ######

#Replace ID in VCF file
sed '/^[^#]/s/NODE_[0-9]\+_length_[0-9]\+_cov_[0-9.]\+/Chromosome/g' snps.vcf > modified_snps.vcf

#Annotate VCF file
snpEff Escherichia_coli_str_k_12_substr_mg1655 -s stats.html modified_snps.vcf  -o vcf > snp_ann.vcf
