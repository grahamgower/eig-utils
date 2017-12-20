/*
 * Copyright (c) 2016 Graham Gower <graham.gower@gmail.com>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include <unistd.h>
#include <getopt.h>
#include <errno.h>
#include <ctype.h>

#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>


typedef struct {
	int ignore_transitions;
	int ignore_monomorphic;
	int ignore_singleton;
	int ignore_uninformative;
	int haploidise_majority_allele;
	int num_autosomes;
	int filter_pass;
	char *depth_fn;
	char *filter_str;
	char *regions_fn;
	char *regions;
	char *oprefix;
} opt_t;

/*
 * Map chromosome name to an integer.
 */
int
chrmap(opt_t *opt, const char *chr)
{
	const char *c = chr;
again:
	switch (tolower(*c)) {
		case 'c':
			if (tolower(c[1]) == 'h' && tolower(c[2]) == 'r') {
				c += 3;
				goto again;
			}
			break;
		case 'x':
			if (tolower(c[1]) == 'y')
				// 'XY': pseudoautosomal region
				return 91;
			else
				return opt->num_autosomes+1;
		case 'y':
			return opt->num_autosomes+2;
		case 'm':
			if (c[1] == 0 || (tolower(c[1]) == 't' && c[2] == 0))
				return 90;
	}

	int d = atoi(c);
	if (d < 1 || d > 87) {
		fprintf(stderr, "Unsupported chromosome name ``%s''.\n"
				"Eigensoft/AdmixTools have a hardcoded limit of 89 chromosome numbers,\n"
				"of which two are reserved for sex chromosomes.  If you're using contigs\n"
				"or scaffolds, try renaming the 87 biggest to numbers.\n", chr);
		return -1;
	}

	return d;
}

int
vcf2eig(opt_t *opt, char **vcflist, int n)
{
	bcf_srs_t *sr;
	bcf_hdr_t *hdr;
	bcf1_t *rec;
	int ret;
	int i, j;
	int nr;
	int nsamples = 0;
       
	int *gt; // genotypes (number of REF alleles for each sample)
	int ndpr_arr = 0;
	int32_t *dpr_arr = NULL; // DPR or AD array
	int *max_dp;

	char ref, alt;

	char buf[strlen(opt->oprefix) + 10];
	FILE *snp_fp, *geno_fp, *ind_fp;

	sr = bcf_sr_init();
	if (sr == NULL) {
		ret = -1;
		goto err0;
	}

	sr->require_index = 1;
	sr->collapse = COLLAPSE_ANY;

	if (opt->regions_fn && bcf_sr_set_regions(sr, opt->regions_fn, 1) < 0) {
		ret = -2;
		goto err0;
	} else if (opt->regions && bcf_sr_set_regions(sr, opt->regions, 0) < 0) {
		ret = -2;
		goto err0;
	}

	for (i=0; i<n; i++) {
		if (bcf_sr_add_reader(sr, vcflist[i]) == 0) {
			fprintf(stderr, "%s: %s\n", vcflist[i], bcf_sr_strerror(sr->errnum));
			ret = -3;
			goto err0;
		}
		nsamples += bcf_hdr_nsamples(bcf_sr_get_header(sr, i));
	}

	gt = calloc(nsamples, sizeof(int));
	if (gt == NULL) {
		perror("calloc");
		ret = -4;
		goto err1;
	}

	sprintf(buf, "%s.snp", opt->oprefix);
	snp_fp = fopen(buf, "w");
	if (snp_fp == NULL) {
		fprintf(stderr, "%s: %s\n", buf, strerror(errno));
		ret = -5;
		goto err2;
	}

	sprintf(buf, "%s.geno", opt->oprefix);
	geno_fp = fopen(buf, "w");
	if (geno_fp == NULL) {
		fprintf(stderr, "%s: %s\n", buf, strerror(errno));
		ret = -6;
		goto err3;
	}

	sprintf(buf, "%s.ind", opt->oprefix);
	ind_fp = fopen(buf, "w");
	if (ind_fp == NULL) {
		fprintf(stderr, "%s: %s\n", buf, strerror(errno));
		ret = -7;
		goto err4;
	}

	for (i=0; i<sr->nreaders; i++) {
		hdr = bcf_sr_get_header(sr, i);
		for (j=0; j<bcf_hdr_nsamples(hdr); j++) {
			const char *sample = bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, j);
			fprintf(ind_fp, "%.39s\tU\t%.39s\n", sample, sample);
		}
	}
	fclose(ind_fp);

	max_dp = calloc(nsamples, sizeof(int));
	if (max_dp == NULL) {
		perror("calloc");
		ret = -8;
		goto err4;
	}

	if (opt->depth_fn) {
		FILE *fp;
		char *linebuf = NULL;
		size_t linebuflen = 0;
		ssize_t nbytes;

		fp = fopen(opt->depth_fn, "r");
		if (fp == NULL) {
			fprintf(stderr, "%s: %s\n", opt->depth_fn, strerror(errno));
			ret = -9;
			goto err5;
		}

		// parse max depths into max_dp; this is O(nsamples^2)
		for (;;) {
			int dp;
			int x; // sample index
			char *p;
			char *sample;

			errno = 0;
			if ((nbytes = getline(&linebuf, &linebuflen, fp)) == -1) {
				if (errno) {
					fprintf(stderr, "getline: %s: %s\n",
							opt->depth_fn, strerror(errno));
					fclose(fp);
					ret = -10;
					goto err5;
				}
				break;
			}

			if (linebuf[nbytes-1] == '\n')
				linebuf[nbytes-1] = 0;

			p = linebuf;

			// columns are: sample maxDepth
			sample = p;

			// skip to end of string and null terminate
			while (*p != ' ' && *p != '\t')
				p++;
			*(p++) = '\0';

			// next column
			while (*p == ' ' || *p == '\t')
				p++;

			dp = atoi(p);

			for (i=0, x=0; i<sr->nreaders; i++) {
				hdr = bcf_sr_get_header(sr, i);
				for (j=0; j<bcf_hdr_nsamples(hdr); j++, x++) {
					if (!strcmp(sample, bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, j)))
						break;
				}
				if (j != bcf_hdr_nsamples(hdr))
					break;
			}
			if (x < nsamples) {
				// got one
				max_dp[x] = dp;
			}

		}

		/* manually verify samples names were as expected
		int x;
		for (i=0, x=0; i<sr->nreaders; i++) {
			hdr = bcf_sr_get_header(sr, i);
			for (j=0; j<bcf_hdr_nsamples(hdr); j++, x++) {
				fprintf(stderr, "%s\t%d\n", bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, j), max_dp[x]);
			}
		}*/

		fclose(fp);
	}

	unsigned short xsubi[3] = {31,41,59}; // random state
	char *chr = NULL;

	while ((nr = bcf_sr_next_line(sr))) {

		ref = '.';
		alt = '.';

		int x = 0; // current sample index
		int chrid = -1;
		int pos = -1;

		for (i=0; i<n; i++) {
			hdr = bcf_sr_get_header(sr, i);
			rec = bcf_sr_get_line(sr, i);

			if (rec == NULL || rec->n_allele > 2
					|| (opt->filter_pass && !bcf_has_filter(hdr, rec, "."))) {
				// set missing
				for (j=0; j<bcf_hdr_nsamples(hdr); j++)
					gt[x++] = 9;
				continue;
			}

			if (!bcf_is_snp(rec) || bcf_get_info(hdr, rec, "INDEL"))
				break;

			if (opt->filter_str && bcf_has_filter(hdr, rec, opt->filter_str) > 0)
				break;

			//printf("[%s:%d, %d] nr=%d, qual=%lf, n_allele=%d\n", bcf_seqname(hdr, rec), rec->pos+1, i, nr, rec->qual, rec->n_allele);

			if (bcf_get_format_int32(hdr, rec, "DPR", &dpr_arr, &ndpr_arr) < 1) {
				// DPR is deprecated, try AD
				if (bcf_get_format_int32(hdr, rec, "AD", &dpr_arr, &ndpr_arr) < 1) {
					fprintf(stderr, "Error: %s: no FORMAT/DPR nor FORMAT/AD field at %s:%d.\n",
						vcflist[i], bcf_seqname(hdr, rec), rec->pos+1);
					ret = -11;
					goto err6;
				}
			}

			if (ref == '.')
				ref = rec->d.allele[0][0];

			for (j=0; j<bcf_hdr_nsamples(hdr); j++) {
				if (rec->n_allele == 1) {
					// we have variant type VCF_REF
					if (dpr_arr[j] == 0)
						// no reads covering locus
						gt[x++] = 9;
					else
						// only REF alleles observed
						gt[x++] = 2;
					continue;
				}

				// rec->n_allele == 2

				// only take SNPs
				if (bcf_get_variant_type(rec, 1) != VCF_SNP)
					break;

				char alt_i = rec->d.allele[1][0];
				if (alt == '.')
					alt = alt_i;
				else if (alt != alt_i) {
					// different ALT genotypes called in different samples
					break;
				}

				int32_t dp_ref = dpr_arr[j*2];
				int32_t dp_alt = dpr_arr[j*2+1];
				if (dp_ref == 0) {
					if (dp_alt == 0)
						// no reads covering locus
						gt[x++] = 9;
					else
						// only ALT alleles observed
						gt[x++] = 0;
					continue;
				} else {
					if (max_dp[x] && dp_ref+dp_alt > max_dp[x]) {
						// depth outlier, ignore
						gt[x++] = 9;
						continue;
					}

					if (opt->haploidise_majority_allele) {
						if (dp_ref > dp_alt) {
							gt[x++] = 2;
							continue;
						} else if (dp_ref < dp_alt) {
							gt[x++] = 0;
							continue;
						}
					}
					// randomly select a read
					gt[x++] = (nrand48(xsubi) % (dp_ref + dp_alt) < dp_ref) ? 2 : 0;
				}
			}

			if (j != bcf_hdr_nsamples(hdr))
				// filtered
				break;

			if (chrid == -1) {
				if (chr) {
					free(chr);
					chr = NULL;
				}
				chr = strdup(bcf_seqname(hdr, rec));
				chrid = chrmap(opt, bcf_seqname(hdr, rec));
				if (chrid == -1) {
					ret = -12;
					goto err6;
				}
				pos = rec->pos+1;
			}
		}

		if (i != n)
			// filtered
			continue;

		if (opt->ignore_monomorphic && alt == '.')
			// skip monomorphic (REF) sites
			continue;

		{
			int gt_sum = 0;
			int mono_ref = 0;
			for (x=0; x<nsamples; x++) {
				if (gt[x] != 9) {
					gt_sum += gt[x];
					mono_ref += 2;
				}
			}

			// skip sites with no information
			if (opt->ignore_uninformative && mono_ref == 0)
				continue;

			// skip monomorphic sites (all samples share an allele)
			if (opt->ignore_monomorphic && (gt_sum == mono_ref || gt_sum == 0))
				continue;

			// skip singleton sites (all samples but one share an allele)
			if (opt->ignore_singleton && (gt_sum == mono_ref-2 || gt_sum == 2))
				continue;
		}

		// filter transitions
		if (opt->ignore_transitions &&
		   ((ref == 'C' && alt == 'T') || (alt == 'C' && ref == 'T') ||
		   (ref == 'G' && alt == 'A') || (alt == 'G' && ref == 'A')))
			continue;

		// write out genotypes
		for (x=0; x<nsamples; x++)
			putc('0'+gt[x], geno_fp);
		putc('\n', geno_fp);

		// write out SNP
		fprintf(snp_fp, "%s_%d\t%d\t0.0\t%d\t%c\t%c\n",
				chr, pos, chrid, pos, ref, alt=='.'?'X':alt);
	}

	ret = 0;
err6:
	if (dpr_arr != NULL)
		free(dpr_arr);
	if (chr)
		free(chr);
err5:
	free(max_dp);
err4:
	fclose(geno_fp);
err3:
	fclose(snp_fp);
err2:
	free(gt);
err1:
        bcf_sr_destroy(sr);
err0:
	return ret;
}

void
usage(char *argv0)
{
	fprintf(stderr, "Convert vcfs to Eigenstrat/AdmixTools format.\n\n"
			" Only biallelic SNPs are considered, and genotypes are haploidised.\n"
			" REF or ALT is called by randomly sampling in proportion to depth\n"
			" using FORMAT/DPR or FORMAT/AD fields.\n\n");

	fprintf(stderr, "usage: %s [...] file1.vcf [... fileN.vcf]\n", argv0);
	fprintf(stderr, "   -r STR           Comma separated list of regions to include []\n");
	fprintf(stderr, "   -R FILE          Bed file listing regions to include []\n");
	fprintf(stderr, "   -t               Exclude transitions [no]\n");
	fprintf(stderr, "   -m               Output monomorphic sites [no]\n");
	fprintf(stderr, "   -s               Output singleton sites [no]\n");
	fprintf(stderr, "   -u               Output uninformative sites (all alleles missing) [no]\n");
	//fprintf(stderr, "   -f               Treat calls not PASSing the FILTER as missing [no]\n"); // Don't use for QUAL filtering, it biases low coverage samples towards the reference
	fprintf(stderr, "   -F STR           Ignore sites with STR in any file's FILTER column []\n");
	fprintf(stderr, "   -a INT           Number of autosomes [29]\n");
	fprintf(stderr, "   -o STR           Output file prefix [out]\n");
	fprintf(stderr, "   -j               Use majority allele for genotype call [no]\n");
	fprintf(stderr, "   -d FILE          Max depth for samples (file format: Sample\tMaxDepth) []\n");
	exit(1);
}

int
main(int argc, char **argv)
{
	opt_t opt;
	int c;

	opt.depth_fn = NULL;
	opt.filter_str = NULL;
	opt.filter_pass = 0;
	opt.regions_fn = NULL;
	opt.regions = NULL;
	opt.ignore_transitions = 0;
	opt.ignore_monomorphic = 1;
	opt.ignore_singleton = 1;
	opt.ignore_uninformative = 1;
	opt.haploidise_majority_allele = 0;
	opt.num_autosomes = 29;
	opt.oprefix = "out";

	while ((c = getopt(argc, argv, "a:d:r:R:o:F:tmsufj")) != -1) {
		switch (c) {
			case 't':
				opt.ignore_transitions = 1;
				break;
			case 'm':
				opt.ignore_monomorphic = 0;
				break;
			case 's':
				opt.ignore_singleton = 0;
				break;
			case 'u':
				opt.ignore_uninformative = 0;
				break;
			case 'R':
				opt.regions_fn = optarg;
				break;
			case 'r':
				opt.regions = optarg;
				break;
			case 'o':
				opt.oprefix = optarg;
				break;
			case 'a':
				opt.num_autosomes = strtoul(optarg, NULL, 0);
				if (opt.num_autosomes < 1 || opt.num_autosomes > 87) {
					fprintf(stderr, "num_autosomes=`%s' out of range\n", optarg);
					return -1;
				}
				break;
			case 'f':
				opt.filter_pass = 1;
				break;
			case 'F':
				opt.filter_str = optarg;
				break;
			case 'j':
				opt.haploidise_majority_allele = 1;
				break;
			case 'd':
				opt.depth_fn = optarg;
				break;
			default:
				usage(argv[0]);
		}
	}

	if (argc-optind < 1) {
		usage(argv[0]);
	}

	if (opt.regions && opt.regions_fn) {
		fprintf(stderr, "-R and -r are mutually exclusive.");
		usage(argv[0]);
	}

	if (vcf2eig(&opt, argv+optind, argc-optind) < 0)
		return -1;

	return 0;
}

