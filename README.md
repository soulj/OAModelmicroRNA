# ACL model of osteoarthritis mRNA and miRNA analysis

The following R notebooks can be used to generate the bioinformatics figures and tables shown in the paper:

* 01_ACLmRNA.Rmd - DESeq2 analysis of the ACL rupture model mRNA-seq
* 02_ACLmiRNA.Rmd - DESeq2 analysis of the ACL rupture model smallRNA-seq
* 03_mir199KD.Rmd - RNA-seq Differential expression, gene ontology and target analysis of mir199 inhibited HACs
* 04_DMMDiffExp.Rmd - DESeq2 analysis of the DMM OA model mRNA-seq

For the smallRNA-seq data the nextflow core smrnaseq v1.1.0
was run using:

```
nextflow run nf-core/smrnaseq -r 1.1.0 --input "fastqFiles/*.fastq.gz" --genome GRCm38 --protocol 'custom' --three_prime_adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --mirtrace_protocol illumina --max_cpus 6 -profile singularity 
```

Skeletalvis pipeline was used to process the RNA-seq data (github.com/soulj/SkeletalVis-Pipeline)


