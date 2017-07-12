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
#include <ctype.h>
#include <math.h>

#include "kmath.h"

// skip to next column
#define next(x) \
		while (*x != ' ' && *x != '\t') x++; \
		while (*x == ' ' || *x == '\t') x++;

typedef struct {
	size_t binsize;
	char *quadlist_fn;
	krand_t *seed;
} opt_t;

typedef struct ind {
	char *s; // name
	char *g; // group
	struct ind *next;
} indlist_t;

typedef struct group {
	char *g;
	int *ii; // indexes of individuals
	int n_inds;
	struct group *next;
} grouplist_t;

typedef struct quad {
	char *s[4];
	grouplist_t *gp[4];
	struct quad *next;
} quadpop_t;

void
free_quadlist(quadpop_t *head)
{
	quadpop_t *quad;
	int i;
	for (quad=head; quad!=NULL; ) {
		quadpop_t *tmp = quad;
		quad = quad->next;
		for (i=0; i<4; i++) {
			if (tmp->s[i])
				free(tmp->s[i]);
		}
		free(tmp);
	}
}

int
parse_quadlist(char *fn, quadpop_t **quadlist, int *n)
{
	FILE *fp;
	int ret;
	char *buf = NULL;
	size_t buflen = 0;
	int n_quads = 0;
	int lineno = 0;
	int i;
	quadpop_t *head = NULL, *cur = NULL;

	*quadlist = NULL;
	*n = 0;

	fp = fopen(fn, "r");
	if (fp == NULL) {
		fprintf(stderr, "fopen: %s: %s\n", fn, strerror(errno));
		ret = -1;
		goto err0;
	}

	while (getline(&buf, &buflen, fp) != -1) {
		char *c0, *c;
		lineno++;

		quadpop_t *quad = calloc(1, sizeof(quadpop_t));
		if (quad == NULL) {
			fprintf(stderr, "calloc: %s\n", strerror(errno));
			ret = -2;
			goto err1;
		}
		quad->next = NULL;
		if (head == NULL)
			head = cur = quad;
		else {
			cur->next = quad;
			cur = quad;
		}
		n_quads++;

		c = buf;
		for (i=0; i<4; i++) {
			c0 = c;
			while (*c != ' ' && *c != '\t' && *c != '\n' && *c != '\r')
				c++;
			if (i!=3 && (*c == '\n' || *c == '\r')) {
				fprintf(stderr, "%s:%d: invalid, expected 4 columns, found %d\n",
						fn, lineno, i);
				ret = -3;
				goto err1;
			}
			*(c++) = '\0';

			quad->s[i] = strdup(c0);
			while (*c == ' ' || *c == '\t')
				c++;
		}

		//printf("[%d] `%s'\t`%s'\t`%s'\t`%s'\n", n_quads-1, quad->s[0], quad->s[1], quad->s[2], quad->s[3]);
	}

	*quadlist = head;
	*n = n_quads;
	ret = 0;
err1:
	if (ret)
		free_quadlist(head);
	if (buflen)
		free(buf);
	fclose(fp);
err0:
	return ret;
}

void
free_indlist(indlist_t *head)
{
	indlist_t *ind;
	for (ind=head; ind!=NULL; ) {
		indlist_t *tmp = ind;
		ind = ind->next;
		free(tmp->s);
		free(tmp->g);
		free(tmp);
	}
}

void
free_grouplist(grouplist_t *head)
{
	grouplist_t *group;
	for (group=head; group!=NULL;) {
		grouplist_t *tmp = group;
		group = group->next;
		free(tmp->g);
		if (tmp->ii)
			free(tmp->ii);
		free(tmp);
	}
}

int
parse_ind(char *fn, indlist_t **indlist, int *_n_inds, grouplist_t **grouplist, int *_n_groups)
{
	int ret;
	FILE *fp;
	char *buf = NULL;
	size_t buflen = 0;
	int n_inds = 0, n_groups = 0;
	indlist_t *head = NULL, *cur = NULL;
	grouplist_t *gr, *ghead = NULL, *gcur = NULL;
	void *tmp;

	*indlist = NULL;
	*_n_inds = 0;
	*grouplist = NULL;
	*_n_groups = 0;

	fp = fopen(fn, "r");
	if (fp == NULL) {
		fprintf(stderr, "fopen: %s: %s\n", fn, strerror(errno));
		ret = -1;
		goto err0;
	}

	while (getline(&buf, &buflen, fp) != -1) {
		char *c0, *c = buf;

		indlist_t *ind = calloc(1, sizeof(indlist_t));
		if (ind == NULL) {
			fprintf(stderr, "calloc: %s\n", strerror(errno));
			ret = -2;
			goto err1;
		}
		ind->next = NULL;
		if (head == NULL)
			head = cur = ind;
		else {
			cur->next = ind;
			cur = ind;
		}
		n_inds++;

		c0 = c;
		while (*c != ' ' && *c != '\t' && *c != '\n' && *c != '\r')
			c++;
		*(c++) = '\0';
		while (*c == ' ' || *c == '\t')
			c++;

		ind->s = strdup(c0);
		next(c);

		c0 = c;
		while (*c != ' ' && *c != '\t' && *c != '\n' && *c != '\r')
			c++;
		*(c++) = '\0';
		ind->g = strdup(c0);

		// find group
		for (gr=ghead; gr!=NULL; gr=gr->next) {
			if (!strcmp(ind->g, gr->g))
				break;
		}

		if (gr == NULL) {
			gr = calloc(1, sizeof(grouplist_t));
			if (gr == NULL) {
				fprintf(stderr, "calloc: %s\n", strerror(errno));
				ret = -3;
				goto err1;
			}
			gr->next = NULL;
			if (ghead == NULL)
				ghead = gcur = gr;
			else {
				gcur->next = gr;
				gcur = gr;
			}
			gr->g = strdup(ind->g);
			n_groups++;
			//printf("adding group[%d]: `%s'\n", n_groups-1, gr->g);
		}

		// append individual to group's index list
		gr->n_inds++;
		tmp = realloc(gr->ii, gr->n_inds*sizeof(*gr->ii));
		if (tmp == NULL) {
			perror("realloc");
			ret = -4;
			goto err1;
		}
		gr->ii = tmp;
		gr->ii[gr->n_inds-1] = n_inds-1;
	}
	
	if (errno) {
		fprintf(stderr, "getline: %s: %s", fn, strerror(errno));
		ret = -3;
		goto err1;
	}

	*indlist = head;
	*_n_inds = n_inds;
	*grouplist = ghead;
	*_n_groups = n_groups;
	ret = 0;

err1:
	if (ret) {
		free_indlist(head);
		free_grouplist(ghead);
	}
	if (buflen)
		free(buf);
	fclose(fp);
err0:
	return ret;
}

int
parse_eig(opt_t *opt, char *ind_fn, char *geno_fn, char *snp_fn)
{
	int i, j;
	int ret;
	FILE *geno_fp, *snp_fp;
	char *gbuf, *sbuf;
	size_t gbuflen, sbuflen;
	ssize_t gnbytes, snbytes;
	uint64_t lineno;
	indlist_t *indlist;
	grouplist_t *grouplist, *group;
	quadpop_t *quadlist, *quad;
	int n_indivs, n_groups, n_quads;
	int max_group_size;

	int chrom = -1, pos = -1;
	int last_chrom = -1, last_pos = -1;

	struct count_t {
		uint64_t aaaa;
		uint64_t aaab;
		uint64_t aaba;
		uint64_t abaa;
		uint64_t baaa;
		uint64_t bbaa;
		uint64_t abba;
		uint64_t baba;
		uint64_t nsites;
	} *counts;


	if (parse_ind(ind_fn, &indlist, &n_indivs, &grouplist, &n_groups) < 0) {
		ret = -1;
		goto err0;
	}
	max_group_size = 0;
	for (group=grouplist; group!=NULL; group=group->next) {
		if (group->n_inds > max_group_size)
			max_group_size = group->n_inds;
	}

	if (parse_quadlist(opt->quadlist_fn, &quadlist, &n_quads) < 0) {
		ret = -2;
		goto err1;
	}

	for (quad=quadlist; quad!=NULL; quad=quad->next) {
		for (i=0; i<4; i++) {
			char *s = quad->s[i];
			for (group=grouplist; group!=NULL; group=group->next) {
				if (!strcmp(s, group->g)) {
					quad->gp[i] = group;
					break;
				}
			}
			if (group == NULL) {
				fprintf(stderr, "Error: %s: group `%s' not found in %s\n",
						opt->quadlist_fn, s, ind_fn);
				ret = -3;
				goto err2;
			}
		}
	}

	counts = calloc(n_quads, sizeof(*counts));
	if (counts == NULL) {
		perror("calloc");
		ret = -4;
		goto err2;
	}

	geno_fp = fopen(geno_fn, "r");
	if (geno_fp == NULL) {
		fprintf(stderr, "fopen: %s: %s\n", geno_fn, strerror(errno));
		ret = -5;
		goto err3;
	}

	snp_fp = fopen(snp_fn, "r");
	if (snp_fp == NULL) {
		fprintf(stderr, "fopen: %s: %s\n", snp_fn, strerror(errno));
		ret = -6;
		goto err4;
	}

	gbuf = sbuf = NULL;
	gbuflen = sbuflen = 0;

	int gtmap[256];
	for (i=0; i<256; i++)
		gtmap[i] = 0xf;
	gtmap['0'] = 0;
	gtmap['2'] = 1;

	// header
	printf("chr\tblockstart\tP1\tP2\tP3\tP4\tAAAA\tAAAB\tAABA\tABAA\tBAAA\tBBAA\tABBA\tBABA\tnsites\n");

	lineno = 0;

	for (;;) {
		gnbytes = getline(&gbuf, &gbuflen, geno_fp);
		snbytes = getline(&sbuf, &sbuflen, snp_fp);
		if (gnbytes == -1 || snbytes == -1)
			break;

		lineno++;

		char ref;
		char alt;
		char *c = sbuf;
		int n_gts = gnbytes;
		while (gbuf[n_gts-1] == '\n' || gbuf[n_gts-1] == '\r')
			n_gts--;

		if (n_gts != n_indivs) {
			fprintf(stderr, "%s has %d individuals, but %s: line %jd supplies %d genotypes\n",
					ind_fn, n_indivs, geno_fn, lineno, n_gts);
			ret = -8;
			goto err5;
		}

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
					snp_fn, lineno);
			ret = -9;
			goto err5;
		}

		//

		if (chrom != last_chrom || last_pos+opt->binsize < pos) {
			if (last_chrom != -1) {
				for (i=0, quad=quadlist; quad!=NULL; quad=quad->next, i++) {
					printf("%d\t%d\t%s\t%s\t%s\t%s\t%jd\t%jd\t%jd\t%jd\t%jd\t%jd\t%jd\t%jd\t%jd\n",
							last_chrom,
							last_pos,
							quad->s[0],
							quad->s[1],
							quad->s[2],
							quad->s[3],
							counts[i].aaaa,
							counts[i].aaab,
							counts[i].aaba,
							counts[i].abaa,
							counts[i].baaa,
							counts[i].bbaa,
							counts[i].abba,
							counts[i].baba,
							counts[i].nsites
							);
				}
				memset(counts, 0, n_quads*sizeof(*counts));
			}
			last_chrom = chrom;
			last_pos = pos;
		}

		for (i=0, quad=quadlist; quad!=NULL; quad=quad->next, i++) {
			int gt[4];
			int gtlist[max_group_size];
			int n_gts;
			int k;

			for (j=0; j<4; j++) {
				grouplist_t *gr = quad->gp[j];
				if (gr->n_inds == 1) {
					int x = gtmap[(int)gbuf[gr->ii[0]]];
					if (x == 0xf)
						// no data
						break;
					gt[j] = x;
					continue;
				}

				n_gts = 0;
				for (n_gts=0, k=0; k<gr->n_inds; k++) {
					int x = gtmap[(int)gbuf[gr->ii[k]]];
					if (x == 0xf)
						continue;
					gtlist[n_gts++] = x;
				}

				if (n_gts == 0)
					break;

				if (n_gts == 1)
					gt[j] = gtlist[0];
				else
					gt[j] = gtlist[kr_rand(opt->seed) % n_gts];
			}

			if (j != 4)
				// missing genotype(s)
				continue;
			
			int pattern = gt[0]<<3 | gt[1]<<2 | gt[2]<<1 | gt[3];
			//printf("[%x] %d%d%d%d\n", pattern, gt[0], gt[1], gt[2], gt[3]);

			switch (pattern) {
				case 0xf:
				case 0x0:
					counts[i].aaaa++;
					break;
				case 0xe:
				case 0x1:
					counts[i].aaab++;
					break;
				case 0xd:
				case 0x2:
					counts[i].aaba++;
					break;
				case 0xb:
				case 0x4:
					counts[i].abaa++;
					break;
				case 0x8:
				case 0x7:
					counts[i].baaa++;
					break;
				case 0xc:
				case 0x3:
					counts[i].bbaa++;
					break;
				case 0x6:
				case 0x9:
					counts[i].abba++;
					break;
				case 0xa:
				case 0x5:
					counts[i].baba++;
					break;
			}
			counts[i].nsites++;

		}
	}

	if (errno) {
		fprintf(stderr, "getline: %s: %s\n",
				gnbytes==-1 ? geno_fn : snp_fn,
				strerror(errno));
		ret = -10;
		goto err5;
	} else if (gnbytes != -1 || snbytes != -1) {
		fprintf(stderr, "%s has more entries than %s -- truncated file?\n",
				gnbytes!=-1 ? geno_fn : snp_fn,
				gnbytes!=-1 ? snp_fn : geno_fn);
		ret = -11;
		goto err5;
	}

	if (chrom != last_chrom || last_pos != pos) {
		for (i=0, quad=quadlist; quad!=NULL; quad=quad->next, i++) {
			printf("%d\t%d\t%s\t%s\t%s\t%s\t%jd\t%jd\t%jd\t%jd\t%jd\t%jd\t%jd\t%jd\t%jd\n",
					last_chrom,
					last_pos,
					quad->s[0],
					quad->s[1],
					quad->s[2],
					quad->s[3],
					counts[i].aaaa,
					counts[i].aaab,
					counts[i].aaba,
					counts[i].abaa,
					counts[i].baaa,
					counts[i].bbaa,
					counts[i].abba,
					counts[i].baba,
					counts[i].nsites
					);
		}
	}

	ret = 0;
err5:
	if (gbuflen)
		free(gbuf);
	if (sbuflen)
		free(sbuf);
	fclose(snp_fp);
err4:
	fclose(geno_fp);
err3:
	free(counts);
err2:
	free_quadlist(quadlist);
err1:
	free_indlist(indlist);
	free_grouplist(grouplist);
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
	fprintf(stderr, "usage: %s [-b BLOCKSIZE] -p 4poplist.txt f.ind f.geno f.snp\n", argv0);
	exit(1);
}

int
main(int argc, char **argv)
{
	opt_t opt;
	int c;
	int ret;

	memset(&opt, 0, sizeof(opt));
	opt.binsize = 5*1000*1000; // 5mb

	while ((c = getopt(argc, argv, "b:p:")) != -1) {
		switch (c) {
			case 'b':
				opt.binsize = parse_bp(optarg);
				break;
			case 'p':
				opt.quadlist_fn = optarg;
				break;
			default:
				usage(argv[0]);
		}
	}

	if (argc-optind != 3 || opt.quadlist_fn == NULL)
		usage(argv[0]);

	opt.seed = kr_srand(31415);

	ret = parse_eig(&opt, argv[optind], argv[optind+1], argv[optind+2]);

	free(opt.seed);

	return ret;
}

