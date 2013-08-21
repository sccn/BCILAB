/*
 *  mex_interface.c
 *  mex_svmperf
 *
 *  Created by Oscar Luaces on 22/02/08.
 *  Copyright 2008 AIC-University of Oviedo at Gijón. All rights reserved.
 *
 */

#include "mex_interface.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#ifdef COMPILE_STRSEP

/*-
 * Copyright (c) 1990, 1993
 *  The Regents of the University of California.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *  This product includes software developed by the University of
 *  California, Berkeley and its contributors.
 * 4. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */


/*
 * Get next token from string *stringp, where tokens are possibly-empty
 * strings separated by characters from delim.
 *
 * Writes NULs into the string at *stringp to end tokens.
 * delim need not remain constant from call to call.
 * On return, *stringp points past the last NUL written (if there might
 * be further tokens), or is NULL (if there are definitely no more tokens).
 *
 * If *stringp is NULL, strsep returns NULL.
 */
char *
strsep(stringp, delim)
register char **stringp;
register const char *delim;
{
	register char *s;
	register const char *spanp;
	register int c, sc;
	char *tok;
	
	if ((s = *stringp) == NULL)
		return (NULL);
	for (tok = s;;) {
		c = *s++;
		spanp = delim;
		do {
			if ((sc = *spanp++) == c) {
				if (c == 0)
					s = NULL;
				else
					s[-1] = 0;
				*stringp = s;
				return (tok);
			}
		} while (sc != 0);
	}
	/* NOTREACHED */
}
#endif


void create_argc_argv(const mxArray *options,int *argc,char **argv[]) {

	
	mwSize buflen = mxGetN(options)*sizeof(mxChar)+1;
	char *buf = (char*)mxMalloc(buflen);
	char **ap, **argv_ptr=(char **)my_malloc(MAX_ARGVS*sizeof(char *));
		mxGetString(options, buf, buflen);

		*argc=1;
	argv_ptr[0]="OLR";
	for (ap = (argv_ptr+1); (*ap = strsep(&buf, " \t")) != NULL;)
		if (**ap != '\0') {
			(*argc)++;
			if (++ap >= &argv_ptr[MAX_ARGVS])
				break;
		}
	
	*argv=argv_ptr;
}

void my_fprintf(FILE *dummy,const char *fmt, ...) {
	if (dummy==stdout) {
		char p[512];
		va_list ap;
		va_start(ap, fmt);
		(void) vsnprintf(p, 512, fmt, ap);
		va_end(ap);
		mexPrintf(p);
	}
}

void my_free(void *ptr) {
	if (ptr!=NULL)
		mxFree(ptr);
}

void *my_realloc(void *ptr, size_t size) {
	if (size==0)
		return(ptr);
	if (ptr==NULL)
		return (my_malloc(size));
	return mxRealloc(ptr,size);
}

void nop() {
}