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
#include <getopt.h>
#include <errno.h>
#include <stdint.h>

typedef struct {
	int ignore_monomorphic;
	int ignore_singleton;
	char *geno_fn;
	char *snp_fn;
	char *oprefix;
} opt_t;

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
		while (gbuf[n_gts-1] == '\n' || gbuf[n_gts-1] == '\r')
			n_gts--;

// skip to next column
#define next(x) \
		while (*x != ' ' && *x != '\t') x++; \
		while (*x == ' ' || *x == '\t') x++;

		// columns are: snpid chrom gpos pos ref alt
		next(c);
		next(c);
		next(c);
		next(c);

		ref = *c;
		next(c);
		alt = *c;

		if (ref == '\n' || ref == '\r' || alt == '\n' || alt == '\r') {
			fprintf(stderr, "%s: line %jd: missing ref/alt field(s)\n",
					opt->snp_fn, (intmax_t)n_sites+1);
			ret = -5;
			goto err4;
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

		fputs(gbuf, out_geno_fp);
		fputs(sbuf, out_snp_fp);

		n_sites++;
	}

	if (errno) {
		fprintf(stderr, "getline: %s: %s\n",
				gnbytes==-1 ? opt->geno_fn : opt->snp_fn,
				strerror(errno));
		ret = -6;
		goto err4;
	} else if (gnbytes != -1 || snbytes != -1) {
		fprintf(stderr, "%s has more entries than %s -- truncated file?\n",
				gnbytes!=-1 ? opt->geno_fn : opt->snp_fn,
				gnbytes!=-1 ? opt->snp_fn : opt->geno_fn);
		ret = -7;
		goto err4;
	}

	ret = 0;
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
	fprintf(stderr, "   -o STR           Output file prefix [eigreduce.out]\n");
	exit(1);
}

int
main(int argc, char **argv)
{
	opt_t opt;
	int c;

	opt.oprefix = "eigreduce.out";
	opt.ignore_monomorphic = 1;
	opt.ignore_singleton = 1;

	while ((c = getopt(argc, argv, "mso:")) != -1) {
		switch (c) {
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

