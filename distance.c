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
#include <math.h>

typedef struct ind {
	char *s; // name
	struct ind *next;
} indlist_t;

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
		char *c = buf;
		while (*c != ' ' && *c != '\t' && *c != '\n' && *c != '\r')
			c++;
		*c = '\0';

		indlist_t *ind = calloc(1, sizeof(indlist_t));
		if (ind == NULL) {
			fprintf(stderr, "calloc: %s\n", strerror(errno));
			ret = -2;
			goto err1;
		}
		ind->s = strdup(buf);
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
parse_eig(char *ind_fn, char *geno_fn, char *snp_fn)
{
	int i, j, x;
	int ret;
	FILE *geno_fp, *snp_fp;
	char *gbuf, *sbuf;
	size_t gbuflen, sbuflen;
	ssize_t gnbytes, snbytes;
	indlist_t *indlist, *ind1,  *ind2;
	int n_indivs;
	int64_t n_sites;
	int64_t *sumsq, *nsites;

	if (parse_ind(ind_fn, &indlist, &n_indivs) < 0) {
		ret = -1;
		goto err0;
	}

	geno_fp = fopen(geno_fn, "r");
	if (geno_fp == NULL) {
		fprintf(stderr, "fopen: %s: %s\n", geno_fn, strerror(errno));
		ret = -2;
		goto err1;
	}

	snp_fp = fopen(snp_fn, "r");
	if (snp_fp == NULL) {
		fprintf(stderr, "fopen: %s: %s\n", snp_fn, strerror(errno));
		ret = -3;
		goto err2;
	}

	nsites = calloc(sizeof(int64_t), n_indivs*(n_indivs-1)/2);
	sumsq = calloc(sizeof(int64_t), n_indivs*(n_indivs-1)/2);
	if (nsites == NULL || sumsq == NULL) {
		perror("calloc");
		ret = -4;
		goto err3;
	}

	n_sites = 0;
	gbuf = sbuf = NULL;
	gbuflen = sbuflen = 0;

	int gtmap[] = {['0']=0, ['1']=-0x1000, ['2']=1, ['9']=-0x1000};

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

		if (n_gts != n_indivs) {
			fprintf(stderr, "%s has %d individuals, but %s: line %jd supplies %d genotypes\n",
					ind_fn, n_indivs, geno_fn, (intmax_t)n_sites+1, n_gts);
			ret = -5;
			goto err4;
		}

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
					snp_fn, (intmax_t)n_sites+1);
			ret = -6;
			goto err4;
		}

		//

		for (x=0, i=0; i<n_indivs; i++) {
			for (j=i+1; j<n_indivs; x++, j++) {
				if (gbuf[i] == '9' || gbuf[j] == '9')
					continue;
				int ci = gtmap[(int)gbuf[i]];
				int cj = gtmap[(int)gbuf[j]];
				sumsq[x] += (ci-cj)*(ci-cj);
				nsites[x]++;
			}
		}

		n_sites++;
	}

	if (errno) {
		fprintf(stderr, "getline: %s: %s\n",
				gnbytes==-1 ? geno_fn : snp_fn,
				strerror(errno));
		ret = -7;
		goto err4;
	} else if (gnbytes != -1 || snbytes != -1) {
		fprintf(stderr, "%s has more entries than %s -- truncated file?\n",
				gnbytes!=-1 ? geno_fn : snp_fn,
				gnbytes!=-1 ? snp_fn : geno_fn);
		ret = -8;
		goto err4;
	}

	for (x=0, ind1=indlist; ind1!=NULL; ind1=ind1->next) {
		for (ind2=ind1->next; ind2!=NULL; x++, ind2=ind2->next) {
			printf("%s\t%s\t%g\n",
					ind1->s, ind2->s,
					(double)sumsq[x]/nsites[x]);
		}
	}

	ret = 0;
err4:
	if (gbuflen)
		free(gbuf);
	if (sbuflen)
		free(sbuf);
err3:
	if (nsites)
		free(nsites);
	if (sumsq)
		free(sumsq);
	fclose(snp_fp);
err2:
	fclose(geno_fp);
err1:
	free_indlist(indlist);
err0:
	return ret;
}

void
usage(char *argv0)
{
	fprintf(stderr, "usage: %s f.ind f.geno f.snp\n", argv0);
	exit(1);
}

int
main(int argc, char **argv)
{
	/*
	int c;
	while ((c = getopt(argc, argv, "")) != -1) {
		switch (c) {
			default:
				usage(argv[0]);
		}
	}*/

	if (argc-optind != 3)
		usage(argv[0]);

	return parse_eig(argv[optind], argv[optind+1], argv[optind+2]);
}

