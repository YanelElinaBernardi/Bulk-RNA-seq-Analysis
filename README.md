# Bulk RNA-seq Analysis

Bulk RNA sequencing (bulk RNA-seq) allows the analysis of the average gene expression in a population of cells. The transcriptomic profile obtained from bulk RNA-seq represents an average across all cells, unlike single-cell RNA sequencing (scRNA-seq), which enables the study of gene expression at the level of individual cells.

This pipeline includes the following key steps:

1. RNA Library Download
2. Quality Control (QC)
3. Read Trimming and Cleaning
4. Concatenation of Reads
5. Removal of Ribosomal RNA (rRNA)
6. Alignment to the Reference Genome
7. Gene Quantification

<br>

**Step 1: Download Bulk RNA Library**

Use the fastq-dump from the [*fastx_toolkit*](https://github.com/agordon/fastx_toolkit) package to download the degradome library from the NCBI Sequence Read Archive [SRA](https://www.ncbi.nlm.nih.gov/sra).

```bash
# For SINGLE-END sequencing
fastq-dump SRR000000                    
# For PAIRED-END sequencing
fastq-dump --split-files SRR000000      

# Replace SRR000000 with the correct identification number for each library.
```
<br>

If the files are compressed, decompress them as needed:

```bash
# Decompress GZ files:
gunzip SRR000000.fastq.gz
# Decompress ZIP files:
unzip SRR000000.fastq.zip
```
<br>

**Step 2: Quality control (QC)**

Evaluate read quality using [*FastQC*](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

```bash
fastqc SRR000000.fastq  # (version v0.10.1)
```
Review the generated reports to identify potential issues such as adapter contamination or low-quality regions.

<br>

**Step 3: Library cleaning**

We used the FastQC report to identify adapters and contaminating sequences. There are many trimming programs but in this analysis we chose
[Trimommatic](https://github.com/usadellab/Trimmomatic) because it allows us to trim the reads to the desired size, remove adapters and filter out low quality reads and contaminants.
For the code we prepared a FASTA file called remove_sequences.fasta with contaminating sequences and adapters, which looked like this example:

```
>contaminant_1
CTTATATTAGGCTTCTCCTCAGCGAAAATCACTGGCCGTCGTTTTACATG
>contaminant_2
AGAGTATTAGGCTTCTCCTCAGCGGAATTCACTGGCCGTCGTTTTACATG
>TruSeq Adapter
GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTATCTCGTAT
```

```bash
# For SINGLE-END reads:
java -jar trimmomatic-0.39.jar SE -phred33 \
    SRR000000.fastq.gz SRR000000_trimmed.fastq.gz \
    ILLUMINACLIP:remove_sequences.fasta:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36

# For PAIRED-END reads:
java -jar trimmomatic-0.39.jar PE -phred33 \
    SRR000000_1.fastq.gz SRR000000_2.fastq.gz \
    SRR000000_forward_paired.fastq.gz SRR000000_forward_unpaired.fastq.gz \
    SRR000000_reverse_paired.fastq.gz SRR000000_reverse_unpaired.fastq.gz \
    ILLUMINACLIP:remove_sequences.fasta:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36
```
<br>

**Step 3.3: Quality control**

Run FastQC again to re-evaluate quality with [*FastQC*](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

```bash
fastqc SRR000000_output.fastq.gz  
fastqc SRR000000_forward_paired.fastq.gz
fastqc SRR000000_reverse_paired.fastq.gz
```

<br>

#### Step 4: Concatenation of Reads

If the reads are split across multiple lanes or files, concatenate them:

```bash
# For SINGLE-END reads:
cat SRR000000_L1_R1_1.fastq.gz SRR000000_L1_R1_2.fastq.gz > SRR000000_concat.fastq.gz

# For PAIRED-END reads:
cat SRR000000_L1_R1.fastq.gz SRR000000_L2_R1.fastq.gz > SRR000000_forward.fastq.gz
cat SRR000000_L1_R2.fastq.gz SRR000000_L2_R2.fastq.gz > SRR000000_reverse.fastq.gz
```

<br>

#### Step 5: Removal of Ribosomal RNA (rRNA)

Filter rRNA sequences with [Bowtie](https://bowtie-bio.sourceforge.net/index.shtml). First, you need to create a FASTA file with all rRNA sequences or use an rRNA database. In this case, Sec_rRNA-all.fasta contained all the rRNA sequences from the organism from which the cells in the RNA isolation originated. -->

<br>

5.1: Index rRNA FASTA

```bash
bowtie-build Sec_rRNA-all.fasta INDEXrRNA
```

5.2: Alignment to Remove rRNA

```bash
# For SINGLE-END reads:
bowtie -v 1 --un SRR000000_Clean.fastq INDEXrRNA SRR000000_concat.fastq.gz --al SRR000000_Ribosomal.fastq.gz
# For PAIRED-END reads:
bowtie -v 1 --un SRR000000_Clean_Paired.fastq INDEXrRNA -1 SRR000000_forward.fastq.gz SRR000000_reverse.fastq.gz --al SRR000000_Ribosomales_Paired.fastq
```

**Key Parameters:**
- `-v 1`: Allows one mismatch in the alignment.
- `--un`: Outputs reads that do not align (keep these for further analysis).
- `--al`: Outputs reads that align to rRNA (discard these)


> **Expected report:**

```
# reads processed: 66368980
# reads with at least one reported alignment: 814686 (1.23%)
# reads that failed to align: 65554294 (98.77%)
Reported 814686 alignments
```
<br>

Please note!
After filtering ribosomal reads with Bowtie, the file you get for unaligned reads (SRR000000_Clean_Paired.fastq) is not separated into two files (_1 and _2) for forward and reverse reads. 
This is a problem when you try to use HISAT2, which requires paired-end reads to be provided as two separate files. To split the single file generated by Bowtie you can use:

```
seqtk seq -1 SRR000000_Clean_Paired.fastq > SRR000000_Clean_Paired_1.fastq
seqtk seq -2 SRR000000_Clean_Paired.fastq > SRR000000_Clean_Paired_2.fastq
```
<br>

#### Step 6: Alignment to the Reference Genome

Align reads to a reference genome using *[Hisat2](https://daehwankimlab.github.io/hisat2)*. Download the genome from [Ensembl](https://www.ensembl.org) or [UCSC](https://genome.ucsc.edu/).

6.1: Index the Genome

```bash
hisat2-build hg38.fa IndexGenome # For example: hg38.fa.gz (human genome)
```

6.2: Alignment

```bash
# For SINGLE-END reads:
hisat2 --no-unal --no-softclip -x IndexGenome -U SRR000000_Clean.fastq -S SRR000000.sam -summary-file SRR000000_summary.txt

# For PAIRED-END reads:
hisat2 --no-unal --no-softclip --dta -x IndexGenome -1 SRR000000_Clean_Paired_1.fastq -2 SRR000000_Clean_Paired_2.fastq -S SRR000000_Paired.sam --summary-file SRR000000_Paired_summary.txt
```
> **Expected report:**

```
# Single End:
67606074 reads; of these:
  67606074 (100.00%) were unpaired; of these:
    6461149 (9.56%) aligned 0 times
    60061930 (88.84%) aligned exactly 1 time
    1082995 (1.60%) aligned >1 times
90.44% overall alignment rate

# Paired End:
38153328 reads; of these:
  38153328 (100.00%) were paired; of these:
    2000000 (5.24%) aligned concordantly 0 times
    35100000 (92.00%) aligned concordantly exactly 1 time
    533328 (1.40%) aligned concordantly >1 times
    ----
    2000000 pairs aligned concordantly 0 times; of these:
      50000 (0.25%) aligned discordantly 1 time
    ----
    1950000 pairs aligned 0 times concordantly or discordantly; of these:
      3900000 mates make up the pairs; of these:
        1000000 (25.64%) aligned 0 times
        2500000 (64.10%) aligned exactly 1 time
        400000 (10.26%) aligned >1 times
95.76% overall alignment rate
```

#### Step 7: Gene Quantification

7.1 Sorting SAM and convert to a BAM format

```bash
samtools sort SRR000000.sam > SRR000000.bam
samtools sort SRR000000_Paired.sam > SRR000000_Paired.bam
```
7.2 Indexing the BAM file for counting

```bash
samtools index SRR000000.bam
samtools index SRR000000_Paired.bam
```

7.3 Counting

Use *[featureCounts](https://subread.sourceforge.net/featureCounts.html)* for gene-level quantification with gene annotation file (GTF format). Download the filefrom [Ensembl](https://www.ensembl.org) or[UCSC](https://genome.ucsc.edu/).

```bash
# We aligned with human genome, so now we use human gene annotation.
# Make sure all BAM files stay in the same directory because we will generate a TAB file with the counts of all imputs (*).

# For SINGLE-END reads:
featureCounts -t exon -g gene_id -a hg38.ensGene.gtf -o Matrix_single.tab *.bam 

# For PAIRED-END reads:
featureCounts -p -t exon -g gene_id -a hg38.ensGene.gtf -o Matrix_paired.tab *.bam
```
