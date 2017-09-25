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
#include <ctype.h>

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
	int thinning_interval;
	char *geno_fn;
	char *snp_fn;
	char *ind_fn;
	char *new_ind_fn;
	char *regions_fn;
	char *oprefix;
} opt_t;

typedef struct ind {
	char *s; // name
	struct ind *next;
} indlist_t;

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

void
free_indlist(indlist_t *head)
{
	indlist_t *ind;
	for (ind=head; ind!=NULL; ) {
		indlist_t *tmp = ind;
		ind = ind->next;
		free(tmp->s);
		free(tmp);
	}
}

int
parse_ind(char *fn, indlist_t **indlist, int *n)
{
	int ret;
	FILE *fp;
	char *buf = NULL;
	size_t buflen = 0;
	int n_indivs = 0;
	indlist_t *head = NULL, *cur = NULL;

	*indlist = NULL;
	*n = 0;

	fp = fopen(fn, "r");
	if (fp == NULL) {
		fprintf(stderr, "fopen: %s: %s\n", fn, strerror(errno));
		ret = -1;
		goto err0;
	}

	while (getline(&buf, &buflen, fp) != -1) {
		char *c, *c0 = buf;

		/* skip leading spaces created by EIGENSOFT */
		while (*c0 == ' ' || *c0 == '\t')
			c0++;

		c = c0;
		while (*c != ' ' && *c != '\t' && *c != '\n' && *c != '\r')
			c++;
		*c = '\0';

		indlist_t *ind = calloc(1, sizeof(indlist_t));
		if (ind == NULL) {
			fprintf(stderr, "calloc: %s\n", strerror(errno));
			ret = -2;
			goto err1;
		}
		ind->s = strdup(c0);
		ind->next = NULL;
		if (head == NULL)
			head = cur = ind;
		else {
			cur->next = ind;
			cur = ind;
		}
		n_indivs++;
	}

	if (errno) {
		fprintf(stderr, "getline: %s: %s", fn, strerror(errno));
		ret = -3;
		goto err1;
	}

	*indlist = head;
	*n = n_indivs;
	ret = 0;

err1:
	if (ret)
		free_indlist(head);
	if (buflen)
		free(buf);
	fclose(fp);
err0:
	return ret;
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
	int64_t lineno;
	char tmpfn[4096];
	kbtree_t(bed) *bt = NULL;
	int *ind_map = NULL;
	int n_old, n_new;

	int32_t last_chrom, last_pos;

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

	if (opt->new_ind_fn) {
		indlist_t *x, *y, *old_ind, *new_ind;
		int oi, ni;

		if (parse_ind(opt->ind_fn, &old_ind, &n_old)) {
			ret = -5;
			goto err4;
		}
		if (parse_ind(opt->new_ind_fn, &new_ind, &n_new)) {
			free_indlist(old_ind);
			ret = -6;
			goto err4;
		}

		ind_map = malloc(n_new*sizeof(int));
		if (ind_map == NULL) {
			perror("malloc");
			free_indlist(old_ind);
			free_indlist(new_ind);
			ret = -7;
			goto err4;
		}

		for (x=new_ind, ni=0; x!=NULL; x=x->next, ni++) {
			for (y=old_ind, oi=0; y!=NULL; y=y->next, oi++) {
				if (!strcmp(x->s, y->s))
					break;
			}
			if (oi == n_old) {
				fprintf(stderr, "%s: %s not found in %s\n",
						opt->new_ind_fn, x->s, opt->ind_fn);
				free_indlist(old_ind);
				free_indlist(new_ind);
				ret = -8;
				goto err5;
			}
			ind_map[ni] = oi;
		}
	}

	if (opt->regions_fn) {
		bt = parse_bed(opt->regions_fn);
		if (bt == NULL) {
			ret = -9;
			goto err5;
		}
	}

	lineno = 0;
	gbuf = sbuf = NULL;
	gbuflen = sbuflen = 0;
	last_chrom = last_pos = -1;

	for (;;) {
		lineno++;
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

		// skip leading spaces
		while(*c == ' ' || *c == '\t')
			c++;

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
					opt->snp_fn, (intmax_t)lineno);
			ret = -10;
			goto err6;
		}

		if (chrom == last_chrom && pos < last_pos+opt->thinning_interval)
			continue;

		if (opt->new_ind_fn)
			n_gts = n_new;

		/*
		 * check for missing rows, invariant, or singletons sites
		 */
		int ref_i = 0, alt_i = 0;
		for (i=0; i<n_gts; i++) {
			int j = opt->new_ind_fn ? ind_map[i] : i;
			switch (gbuf[j]) {
				case '2':
					ref_i++;
					break;
				case '0':
					alt_i++;
					break;
			}
		}
		if (alt_i+ref_i == 0)
			continue;
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

		for (i=0; i<n_gts; i++) {
			int j = opt->new_ind_fn ? ind_map[i] : i;
			fputc(gbuf[j], out_geno_fp);
		}
		fputc('\n', out_geno_fp);
		fputs(sbuf, out_snp_fp);

		last_chrom = chrom;
		last_pos = pos;
	}

	if (errno) {
		fprintf(stderr, "getline: %s: %s\n",
				gnbytes==-1 ? opt->geno_fn : opt->snp_fn,
				strerror(errno));
		ret = -11;
		goto err6;
	} else if (gnbytes != -1 || snbytes != -1) {
		fprintf(stderr, "%s has more entries than %s -- truncated file?\n",
				gnbytes!=-1 ? opt->geno_fn : opt->snp_fn,
				gnbytes!=-1 ? opt->snp_fn : opt->geno_fn);
		ret = -12;
		goto err6;
	}

	ret = 0;
err6:
	if (gbuflen)
		free(gbuf);
	if (sbuflen)
		free(sbuf);

	if (opt->regions_fn)
		kb_destroy(bed, bt);
err5:
	if (ind_map)
		free(ind_map);
err4:
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

int
parse_bp(char *s)
{
	int x;

	switch (tolower(s[strlen(s)-1])) {
		case 'k':
			x = 1000;
			break;
		case 'm':
			x = 1000000;
			break;
		default:
			x = 1;
			break;
	}

	return x*atoi(s);
}

void
usage(char *argv0)
{
	fprintf(stderr, "usage: %s [...] file.geno file.snp [file.ind]\n", argv0);
	fprintf(stderr, "   -m               Output monomorphic sites [no]\n");
	fprintf(stderr, "   -s               Output singleton sites [no]\n");
	fprintf(stderr, "   -i NEW.IND       Output only individuals (and reorder) from NEW.IND []\n");
	fprintf(stderr, "   -R BED           Output only the regions specified in the bed file []\n"
			"                        (intervals must be non-overlapping)\n");
	fprintf(stderr, "   -t INT           Thin the output to no more than one site per INT bp. [1]\n");
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
	opt.thinning_interval = 1;

	while ((c = getopt(argc, argv, "i:mo:R:st:")) != -1) {
		switch (c) {
			case 'i':
				opt.new_ind_fn = optarg;
				break;
			case 'm':
				opt.ignore_monomorphic = 0;
				break;
			case 'o':
				opt.oprefix = optarg;
				break;
			case 'R':
				opt.regions_fn = optarg;
				break;
			case 's':
				opt.ignore_singleton = 0;
				break;
			case 't':
				opt.thinning_interval = parse_bp(optarg);
				if (opt.thinning_interval < 0 || opt.thinning_interval > 100*1000*1000) {
					fprintf(stderr, "Thinning interval `%s' is out of range\n", optarg);
					usage(argv[0]);
				}
				break;
			default:
				usage(argv[0]);
		}
	}

	if (argc-optind != 2 && argc-optind != 3)
		usage(argv[0]);

	opt.geno_fn = argv[optind];
	opt.snp_fn = argv[optind+1];
	if (argc-optind == 3) {
		opt.ind_fn = argv[optind+2];
	}

	if (opt.new_ind_fn && opt.ind_fn == NULL) {
		fprintf(stderr, "The .ind file matching the input data has "
				"not been provided, but is required to use "
				"the -i option.\n");
		usage(argv[0]);
	}

	return parse_eig(&opt);
}

