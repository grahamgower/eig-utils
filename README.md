# eig-utils
 
Utilities for manipulation of EIGENSTRAT format files, as used by
[EIGENSOFT](https://github.com/DReichLab/EIG) and
[ADMIXTOOLS](https://github.com/DReichLab/AdmixTools).

vcf2eig - call pseudohaploid genotypes from vcfs and output EIGENSTRAT format  
eigreduce - filter existing EIGENSTRAT files (e.g. remove individuals/sites)  
eig2phylip - convert EIGENSTRAT files to phylip format (e.g. for RAxML)  
eig2multifasta - convert EIGENSTRAT files to multifasta (e.g. for Saguaro)  
distance - print pairwise distance between all individuals  
abba-baba - D stats reimplementation  

# Prerequisites
`vcf2eig` uses **htslib** to parse vcf/bcf files.

# Installation
Clone the git repository, then build with `make`.

# vcf2eig
`vcf2eig` calls pseudohaploid genotypes from the read pileup.  This program
accepts vcf/bcf files as input, not sam/bam, in order that indels may be
excluded more easily.  The genotype call in the vcf (if any) is ignored.
A haploid REF or ALT genotype is called, based on the read pileup, either
randomly sampling an allele in proportion to the allele depths (default),
or taking the majority allele (-j flag).
An example pipeline might be as follows.

1. Call variants (making sure to output required field FORMAT/AD)
```
for b in s1.bam s2.bam s3.bam; do
  out=${b%.bam}.bcf.gz
  samtools mpileup -f ref.fa -t AD -g -u $b \
    | bcftools call -c -O b -o $out
  bcftools index $out
done
```

2. Mark snps within 3 bp of indels with the `SnpGap` filter.
```
for b in s1.bcf.gz s2.bcf.gz s3.bcf.gz; do
  out=SnpGap.$b
  bcftools filter --SnpGap 3 -O b -o $out $b
  bcftools index $out
```

3. Convert to eigenstrat format, ignoring sites with the SnpGap filter
   (-F SnpGap), doing majority allele pseudohaploidisation (-j),
   outputting monomorphic/invariant (-m) and singleton (-s) sites.
   This can be slow, so splitting by chromosome (-r) and working in parallel is
   recommended for large vcfs.
```
nchrom=34 # number of autosomes for your organism
for c in $(seq $nchrom); do
  chr="chr$c" # chr prefix used in most references
  echo vcf2eig -r $chr -F SnpGap -m -s -j -o eigdata.$chr
done | parallel -j 8
```

4. Merge output files from individual chromosomes.
   The *.snp and *.geno files contain a single line for each locus,
   and can simply be catenated.
```
for c in $(seq $nchrom); do cat eigdata.chr$c.snp; done > all_chr.snp
for c in $(seq $nchrom); do cat eigdata.chr$c.geno; done > all_chr.geno
```

# Further info
More information about usage and options can be obtained from the utilities
by running them with no parameters, or inspection of the code.
```
$ ./vcf2eig
Convert vcfs to Eigenstrat/AdmixTools format.

 Only biallelic SNPs are considered, and genotypes are haploidised.
 REF or ALT is called by randomly sampling in proportion to depth
 using FORMAT/DPR or FORMAT/AD fields.

usage: ./vcf2eig [...] file1.vcf [... fileN.vcf]
   -r STR           Comma separated list of regions to include []
   -R FILE          Bed file listing regions to include []
   -t               Exclude transitions [no]
   -m               Output monomorphic sites [no]
   -s               Output singleton sites [no]
   -u               Output uninformative sites (all alleles missing) [no]
   -F STR           Ignore sites with STR in any file's FILTER column []
   -a INT           Number of autosomes [29]
   -o STR           Output file prefix [out]
   -j               Use majority allele for genotype call [no]
   -d FILE          Max depth for samples (file format: Sample  MaxDepth) []
```

```
$ ./eigreduce
usage: ./eigreduce [...] file.geno file.snp [file.ind]
   -m               Output monomorphic sites [no]
   -s               Output singleton sites [no]
   -i NEW.IND       Output only individuals (and reorder) from NEW.IND []
   -R BED           Output only the regions specified in the bed file []
                        (intervals must be non-overlapping)
   -t INT           Thin the output to no more than one site per INT bp. [1]
   -o STR           Output file prefix [eigreduce.out]
```

```
$ ./eig2phylip
usage: ./eig2phylip -o OUTPREFIX f.ind f.geno f.snp
   -n               Output numbers (0,1,2,?) instead of IUPAC codes [no]
   -m               Output monomorphic sites [no]
   -o STR           Output file prefix [phylip.out]
```
