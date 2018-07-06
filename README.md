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
Clone the git repository, then build with `make`.  Some minor changes to the
Makefile may be requred, e.g. to specify the path to htslib.

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

3. Convert to EIGENSTRAT format, ignoring sites with the SnpGap filter
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
ln eigdata.chr1.ind all_chr.ind
```

# eig2phylip

`eig2phylip` was written to provide input for
[RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html).
The EIGENSTRAT snp/geno files have one locus per line with columns
corresponding to individuals; in contrast, the phylip format is
transposed, with one inidividual per line and each column corresponding
to one locus.
`eig2phylip` takes the most memory efficient approach to conversion, and
outputs one file per individual, writing genotypes as the input files
are parsed.  The result is that the output files need to be catenated
once `eig2phylip` completes.

```
tmppfx=temp.phylip
eig2phylip -o $tmppfx all_chr.ind all_chr.geno all_chr.snp

# catenate output into single file
(   cat ${tmppfx}.HEADER.txt
    for f in ${tmppfx}.*; do
       case "$f" in
         "${tmppfx}.HEADER.txt") continue ;;
         "${tmppfx}.invariant.txt") continue ;;
       esac
      cat $f
   done
) > all_chr.phy

mv ${tmppfx}.invariant.txt all_chr.invariant.txt
rm ${tmppfx}.*
```

Addionally, `eig2phylip` outputs a file `*.invariant.txt`, which
contains counts of monomorphic/invariant sites for each of the four
nucleotides.  This can be used with the RAxML option `--asc-corr=stamatakis`
in conjuction with RAxML's `ASC_GTRGAMMA` model.  Note that the input
data to `eig2phylip` must contain monomorphic sites (use `vcf2eig`'s -m flag)
for the counts to be obtained.  It is also important to maintain singleton
sites (`vcf2eig`'s -s flag) for terminal branch lengths in RAxML's output to
have meaning.  This is in contrast to use with EIGENSOFT and ADMIXTOOLS,
where monomorphic and singleton sites contribute nothing, but dramatically
increases file sizes.  `eigreduce` can be used to remove monomorphic and
singleton sites from existing EIGENSTRAT files, which will reduce memory
consumption and run times when using EIGENSOFT and ADMIXTOOLS.

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
