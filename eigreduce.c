/*
 * Copyright (c) 2016,2017 Graham Gower <graham.gower@gmail.com>
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
#include <getopt.h>
#include <errno.h>
#include <stdint.h>

#include "kbtree.h"

typedef struct {
	uint64_t left; // (chrom_id << 32) | left_pos
	uint32_t right;
} interval_t;

#define left_cmp(a, b) (((b).left < (a).left) - ((a).left < (b).left))
KBTREE_INIT(bed, interval_t, left_cmp);

// skip to next column
#define next(x) \
		while (*x != ' ' && *x != '\t') x++; \
		while (*x == ' ' || *x == '\t') x++;

typedef struct {
	int ignore_monomorphic;
	int ignore_singleton;
	char *geno_fn;
	char *snp_fn;
	char *regions_fn;
	char *oprefix;
} opt_t;

/*
 * Parse bed intervals and place into a b-tree.
 * This assumes non-overlapping intervals.
 */
kbtree_t(bed) *
parse_bed(char *fn)
{
	kbtree_t(bed) *bt;
	FILE *fp;
	char *p;
	char *buf = NULL;
	size_t buflen = 0;
	ssize_t nbytes;

	uint32_t chrom, from, to;

	fp = fopen(fn, "r");
	if (fp == NULL) {
		fprintf(stderr, "fopen: %s: %s\n", fn, strerror(errno));
		bt = NULL;
		goto err0;
	}

	bt = kb_init(bed, KB_DEFAULT_SIZE);

	for (;;) {
		if ((nbytes = getline(&buf, &buflen, fp)) == -1) {
			if (errno) {
				fprintf(stderr, "getline: %s: %s\n",
						fn, strerror(errno));
				kb_destroy(bed, bt);
				bt = NULL;
				goto err1;
			}
			break;
		}

		if (buf[nbytes-1] == '\n')
			buf[nbytes-1] = 0;

		p = buf;

		// columns are: chrom from to ...
		chrom = atoi(p);
		next(p);
		from = atoi(p);
		next(p);
		to = atoi(p);

		interval_t in = {(uint64_t)chrom<<32 | from, to};
		kb_putp(bed, bt, &in);
	}

err1:
	fclose(fp);
err0:
	return bt;
}

int
parse_eig(opt_t *opt)
{
	int i;
	int ret;
	FILE *geno_fp, *snp_fp;
	FILE *out_geno_fp, *out_snp_fp;
	char *gbuf, *sbuf;
	size_t gbuflen, sbuflen;
	ssize_t gnbytes, snbytes;
	int64_t n_sites;
	char tmpfn[4096];
	kbtree_t(bed) *bt = NULL;

	geno_fp = fopen(opt->geno_fn, "r");
	if (geno_fp == NULL) {
		fprintf(stderr, "fopen: %s: %s\n", opt->geno_fn, strerror(errno));
		ret = -1;
		goto err0;
	}

	snprintf(tmpfn, 4096, "%s.geno", opt->oprefix);
	out_geno_fp = fopen(tmpfn, "w");
	if (out_geno_fp == NULL) {
		fprintf(stderr, "fopen: %s: %s\n", tmpfn, strerror(errno));
		ret = -2;
		goto err1;
	}

	snp_fp = fopen(opt->snp_fn, "r");
	if (snp_fp == NULL) {
		fprintf(stderr, "fopen: %s: %s\n", opt->snp_fn, strerror(errno));
		ret = -3;
		goto err2;
	}

	snprintf(tmpfn, 4096, "%s.snp", opt->oprefix);
	out_snp_fp = fopen(tmpfn, "w");
	if (out_snp_fp == NULL) {
		fprintf(stderr, "fopen: %s: %s\n", tmpfn, strerror(errno));
		ret = -4;
		goto err3;
	}

	if (opt->regions_fn) {
		bt = parse_bed(opt->regions_fn);
		if (bt == NULL) {
			ret = -5;
			goto err4;
		}
	}

	n_sites = 0;
	gbuf = sbuf = NULL;
	gbuflen = sbuflen = 0;

	for (;;) {
		gnbytes = getline(&gbuf, &gbuflen, geno_fp);
		snbytes = getline(&sbuf, &sbuflen, snp_fp);
		if (gnbytes == -1 || snbytes == -1)
			break;

		char ref;
		char alt;
		char *c = sbuf;
		int n_gts = gnbytes;
		int32_t chrom, pos;
		while (gbuf[n_gts-1] == '\n' || gbuf[n_gts-1] == '\r')
			n_gts--;

		// columns are: snpid chrom gpos pos ref alt
		next(c);
		chrom = atoi(c);
		next(c);
		next(c);
		pos = atoi(c);
		next(c);
		ref = *c;
		next(c);
		alt = *c;

		if (ref == '\n' || ref == '\r' || alt == '\n' || alt == '\r') {
			fprintf(stderr, "%s: line %jd: missing ref/alt field(s)\n",
					opt->snp_fn, (intmax_t)n_sites+1);
			ret = -6;
			goto err5;
		}

		// check for invariant or singletons sites
		int ref_i = 0, alt_i = 0;
		for (i=0; i<n_gts; i++) {
			switch (gbuf[i]) {
				case '2':
					ref_i++;
					break;
				case '0':
					alt_i++;
					break;
			}
		}
		if (opt->ignore_monomorphic && (alt_i == 0 || ref_i == 0))
			continue;
		if (opt->ignore_singleton && (alt_i == 1 || ref_i == 1))
			continue;

		if (opt->regions_fn) {
			interval_t *l = NULL, *u = NULL;
			interval_t x = {(uint64_t)chrom<<32 | pos, 0};

			// get nearest lower and upper entries to pos in the b-tree
			kb_interval(bed, bt, x, &l, &u);

			/*
			printf("[%d,%d], l=[%d,%d], u=[%d,%d], skip=%d\n",
					chrom, pos,
					l ? (uint32_t)l->left : -1, l ? l->right: -1,
					u ? (uint32_t)u->left : -1, u ? u->right: -1,
					l==NULL || (l && pos > l->right));
			*/

			if (l == NULL || (l && pos > l->right))
				// not in the intervals
				continue;
		}

		fputs(gbuf, out_geno_fp);
		fputs(sbuf, out_snp_fp);

		n_sites++;
	}

	if (errno) {
		fprintf(stderr, "getline: %s: %s\n",
				gnbytes==-1 ? opt->geno_fn : opt->snp_fn,
				strerror(errno));
		ret = -7;
		goto err5;
	} else if (gnbytes != -1 || snbytes != -1) {
		fprintf(stderr, "%s has more entries than %s -- truncated file?\n",
				gnbytes!=-1 ? opt->geno_fn : opt->snp_fn,
				gnbytes!=-1 ? opt->snp_fn : opt->geno_fn);
		ret = -8;
		goto err5;
	}

	ret = 0;
err5:
	if (opt->regions_fn)
		kb_destroy(bed, bt);
err4:
	if (gbuflen)
		free(gbuf);
	if (sbuflen)
		free(sbuf);

	fclose(out_snp_fp);
err3:
	fclose(snp_fp);
err2:
	fclose(out_geno_fp);
err1:
	fclose(geno_fp);
err0:
	return ret;
}

void
usage(char *argv0)
{
	fprintf(stderr, "usage: %s [...] file.geno file.snp\n", argv0);
	fprintf(stderr, "   -m               Output monomorphic sites [no]\n");
	fprintf(stderr, "   -s               Output singleton sites [no]\n");
	fprintf(stderr, "   -R BED           Output only the regions specified in the bed file []\n"
			"                        (intervals must be non-overlapping)\n");
	fprintf(stderr, "   -o STR           Output file prefix [eigreduce.out]\n");
	exit(1);
}

int
main(int argc, char **argv)
{
	opt_t opt = {0,};
	int c;

	opt.oprefix = "eigreduce.out";
	opt.ignore_monomorphic = 1;
	opt.ignore_singleton = 1;

	while ((c = getopt(argc, argv, "mso:R:")) != -1) {
		switch (c) {
			case 'R':
				opt.regions_fn = optarg;
				break;
			case 'o':
				opt.oprefix = optarg;
				break;
			case 'm':
				opt.ignore_monomorphic = 0;
				break;
			case 's':
				opt.ignore_singleton = 0;
				break;
			default:
				usage(argv[0]);
		}
	}

	if (argc-optind != 2)
		usage(argv[0]);

	opt.geno_fn = argv[optind];
	opt.snp_fn = argv[optind+1];

	return parse_eig(&opt);
}

