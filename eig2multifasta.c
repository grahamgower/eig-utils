/*
 * Convert eig format to multifasta, written for input to Saguaro.
 *
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

// ambiguity codes
static char iupac[] = {
	['A'<<8|'A']='A',
	['A'<<8|'C']='M',
	['A'<<8|'G']='R',
	['A'<<8|'T']='W',
	['C'<<8|'A']='M',
	['C'<<8|'C']='C',
	['C'<<8|'G']='S',
	['C'<<8|'T']='Y',
	['G'<<8|'A']='R',
	['G'<<8|'C']='S',
	['G'<<8|'G']='G',
	['G'<<8|'T']='K',
	['T'<<8|'A']='W',
	['T'<<8|'C']='Y',
	['T'<<8|'G']='K',
	['T'<<8|'T']='T',

	['N'<<8|'A']='N',
	['N'<<8|'C']='N',
	['N'<<8|'G']='N',
	['N'<<8|'T']='N',
	['A'<<8|'N']='N',
	['C'<<8|'N']='N',
	['G'<<8|'N']='N',
	['T'<<8|'N']='N',
	['N'<<8|'N']='N',
};

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
parse_eig(char *ind_fn, char *geno_fn, char *snp_fn, char *oprefix)
{
	int i;
	int ret;
	FILE *geno_fp, *snp_fp,
	     **indfplist;
	char *gbuf, *sbuf;
	size_t gbuflen, sbuflen;
	ssize_t gnbytes, snbytes;
	indlist_t *indlist, *curind;
	int n_indivs;
	int64_t n_sites;

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

	indfplist = calloc(n_indivs, sizeof(FILE *));
	if (indfplist == NULL) {
		fprintf(stderr, "calloc: %s\n", strerror(errno));
		ret = -4;
		goto err3;
	}

	for (i=0, curind=indlist; curind!=NULL; i++, curind=curind->next) {
		char tmpfn[4096];
		snprintf(tmpfn, 4096, "%s.%s.singlefasta", oprefix, curind->s);
		indfplist[i] = fopen(tmpfn, "w");
		if (indfplist[i] == NULL) {
			fprintf(stderr, "fopen: %s: %s\n", snp_fn, strerror(errno));
			ret = -5;
			goto err4;
		}
		fprintf(indfplist[i], ">%s\n", curind->s);
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

		if (n_gts != n_indivs) {
			fprintf(stderr, "%s has %d individuals, but %s: line %jd supplies %d genotypes\n",
					ind_fn, n_indivs, geno_fn, (intmax_t)n_sites+1, n_gts);
			ret = -6;
			goto err4;
		}

		// convertf prefixes lines with spaces, skip them
		while (*c == ' ' || *c == '\t') c++;

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
			ret = -7;
			goto err4;
		}

		if (ref != 'A' && ref != 'C' && ref != 'G' && ref != 'T' && ref != 'N') {
			fprintf(stderr, "%s: line %jd: do not understand ref allele '%c'\n", snp_fn, (intmax_t)n_sites+1, ref);
			ret = -8;
			goto err4;
		}
		if (alt == 'X')
			alt = 'N';
		else if (alt != 'A' && alt != 'C' && alt != 'G' && alt != 'T' && alt != 'N') {
			fprintf(stderr, "%s: line %jd: do not understand alt allele '%c'\n", snp_fn, (intmax_t)n_sites+1, ref);
			ret = -9;
			goto err4;
		}

		if ((n_sites%60) == 0 && n_sites > 0) {
			for (i=0; i<n_indivs; i++)
				fputc('\n', indfplist[i]);
		}


		for (i=0; i<n_indivs; i++) {
			switch (gbuf[i]) {
				case '2':
					fputc(ref, indfplist[i]);
					break;
				case '1':
					fputc(iupac[ref<<8|alt], indfplist[i]);
					break;
				case '0':
					fputc(alt, indfplist[i]);
					break;
				default:
					fputc('N', indfplist[i]);
			}
		}

		n_sites++;
	}

	for (i=0; i<n_indivs; i++)
		fputc('\n', indfplist[i]);

	if (errno) {
		fprintf(stderr, "getline: %s: %s\n",
				gnbytes==-1 ? geno_fn : snp_fn,
				strerror(errno));
		ret = -10;
		goto err4;
	} else if (gnbytes != -1 || snbytes != -1) {
		fprintf(stderr, "%s has more entries than %s -- truncated file?\n",
				gnbytes!=-1 ? geno_fn : snp_fn,
				gnbytes!=-1 ? snp_fn : geno_fn);
		ret = -11;
		goto err4;
	}


	ret = 0;
err4:
	if (gbuflen)
		free(gbuf);
	if (sbuflen)
		free(sbuf);
	for (i=0; i<n_indivs; i++) {
		if (indfplist[i])
			fclose(indfplist[i]);
	}
	free(indfplist);
err3:
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
	fprintf(stderr, "usage: %s -o OUTPREFIX f.ind f.geno f.snp\n", argv0);
	exit(1);
}

int
main(int argc, char **argv)
{
	int c;
	char *oprefix = "phylip.out";

	while ((c = getopt(argc, argv, "o:")) != -1) {
		switch (c) {
			case 'o':
				oprefix = optarg;
				break;
			default:
				usage(argv[0]);
		}
	}

	if (argc-optind != 3)
		usage(argv[0]);

	return parse_eig(argv[optind], argv[optind+1], argv[optind+2], oprefix);
}

