#include <ctype.h>
#include <stdio.h>
#include "awk.h"
#include "addon.h"

int lh3_has_colnm = 0;

void lh3_set_colnm()
{
	char *p, *q, c;
	Cell *x;
	int i;
	if (lh3_has_colnm == 0) return;
	for (p = record; *p && isspace(*p); ++p); /* skip leading spaces */
	for (i = 1, q = p; *q; ++q) {
		if (!isspace(*q)) continue;
		c = *q; /* backup the space */
		*q = 0; /* terminate the field */
		x = lookup(p, symtab);
		if (x != NULL) x->fval = i, x->tval = NUM;
		else setsymtab(p, "", i, NUM, symtab);
		*q = c; /* change back */
		++i;
		for (p = q + 1; *p && isspace(*p); ++p); /* skip contiguous spaces */
		q = p;
	}
	if (*p) { /* the last column name */
		x = lookup(p, symtab);
		if (x != NULL) x->fval = i, x->tval = NUM;
		else setsymtab(p, "", i, NUM, symtab);
	}
}
