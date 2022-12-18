/* parser.h */

/*
	$Id: parser.h 9932 2006-11-12 14:26:53Z katayama $
	$Log$
	Revision 1.1  1998/06/17 08:43:12  higuchit
	*** empty log message ***

*/

#ifndef PARSER_H
#define PARSER_H

/* This is for C only, no namespace needed */

typedef struct{
	int pline, cline;
	char *master, *buf, *ptr;
} PRS;

extern PRS *prsopen();
extern void prsrewind();
extern int prstell();
extern int prsnext();
extern void prsclose();

#endif /* PARSER_H */

