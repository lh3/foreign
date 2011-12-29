#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include "awk.h"
#include "addon.h"

const char *lh3_col_defn = NULL;
const char *valid_coldefs[] = {"header", "bed", "bedgraph", "sam", "vcf", "gff", "gtf", NULL};

static void set_colnm_aux(const char *p, int col)
{
	const char *q;
	Cell *x;
	for (q = p; *q; ++q)
		if (!isdigit(*q)) break;
	if (*q == 0) return; /* do not set if string p is an integer */
	if ((x = lookup(p, symtab)) != NULL)
		x->tval = NUM, x->fval = col;
}

int isvalid_coldef(const char *request)
{
    int i;
    for (i = 0; valid_coldefs[i] != NULL; ++i)
        if (strcmp(request, valid_coldefs[i]) == 0)
            return 1;
    return 0;
}

void print_valid_coldefs() 
{
    printf("valid -c options include:\n");
    int i;
    for (i = 0; valid_coldefs[i] != NULL; ++i) 
        printf("  %s\n", valid_coldefs[i]);
}

void lh3_set_colnm()
{
    if (lh3_col_defn == NULL) return;
    
    if (strcmp(lh3_col_defn, "header") == 0) {
    	char *p, *q, c;
    	int i;
    	for (p = record; *p && isspace(*p); ++p); /* skip leading spaces */
    	for (i = 1, q = p; *q; ++q) {
    		if (!isspace(*q)) continue;
    		c = *q; /* backup the space */
    		*q = 0; /* terminate the field */
    		set_colnm_aux(p, i);
    		*q = c; /* change back */
    		++i;
    		for (p = q + 1; *p && isspace(*p); ++p); /* skip contiguous spaces */
    		q = p;
    	}
    	set_colnm_aux(p, i); /* the last column */
    } 
    else if (strcmp(lh3_col_defn, "bed") == 0) {
        set_colnm_aux("chrom", 1);
        set_colnm_aux("start", 2);
        set_colnm_aux("end", 3);
        set_colnm_aux("name", 4);
        set_colnm_aux("score", 5);
        set_colnm_aux("strand", 6);
        // force tab delimited input and output
        *FS = "\t";
        *OFS = "\t";
    }
    else if (strcmp(lh3_col_defn, "bedgraph") == 0) {
        set_colnm_aux("chrom", 1);
        set_colnm_aux("start", 2);
        set_colnm_aux("end", 3);
        set_colnm_aux("score", 4);
        // force tab delimited input and output
        *FS = "\t";
        *OFS = "\t";
    }
    else if (strcmp(lh3_col_defn, "sam") == 0) {
        set_colnm_aux("qname", 1);
        set_colnm_aux("flag", 2);
        set_colnm_aux("rname", 3);
        set_colnm_aux("pos", 4);
        set_colnm_aux("mapq", 5);
        set_colnm_aux("cigar", 6);
        set_colnm_aux("rnext", 7);
        set_colnm_aux("pnext", 8);
        set_colnm_aux("tlen", 9);
        set_colnm_aux("seq", 10);
        set_colnm_aux("qual", 11);
        // todo: any intelligent way to handle tags?
        // force tab delimited input and output
        *FS = "\t";
        *OFS = "\t";
        // auto-report any header lines
        while (getrec(&record, &recsize, 1) > 0 && record[0] == '@') {
            printf("%s\n", record);
        }
    }
    else if (strcmp(lh3_col_defn, "vcf") == 0) {
        set_colnm_aux("chrom", 1);
        set_colnm_aux("pos", 2);
        set_colnm_aux("id", 3);
        set_colnm_aux("ref", 4);
        set_colnm_aux("alt", 5);
        set_colnm_aux("qual", 6);
        set_colnm_aux("filter", 7);
        set_colnm_aux("info", 8);
        // todo: any intelligent way to handle genotypes?
        // force tab delimited input and output
        *FS = "\t";
        *OFS = "\t";
        // auto-report any header lines
        while (getrec(&record, &recsize, 1) > 0 && record[0] == '#') {
            printf("%s\n", record);
        }
    }
    else if (strcmp(lh3_col_defn, "gff") == 0 || strcmp(lh3_col_defn, "gtf") == 0) {
        set_colnm_aux("seqname", 1);
        set_colnm_aux("source", 2);
        set_colnm_aux("feature", 3);
        set_colnm_aux("start", 4);
        set_colnm_aux("end", 5);
        set_colnm_aux("score", 6);
        set_colnm_aux("filter", 7);
        set_colnm_aux("strand", 8);
        set_colnm_aux("group", 9);     // allow group or attribute, as
        set_colnm_aux("attribute", 9); // GFF v1 used group, v2 uses attribute
        // force tab delimited input and output
        *FS = "\t";
        *OFS = "\t";
        // auto-report any header lines
        while (getrec(&record, &recsize, 1) > 0 && record[0] == '#') {
            printf("%s\n", record);
        }
    }
}
