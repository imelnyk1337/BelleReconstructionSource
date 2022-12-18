/* qq2panther.h */

/*
	$Id: qq2panther.h 9932 2006-11-12 14:26:53Z katayama $
	$Log$
	Revision 1.3  2003/05/16 23:01:05  katayama
	For gcc3.3

	Revision 1.2  1998/10/14 12:50:24  jtanaka
	uses qq2panther to get Ptype information.
	
	Revision 1.1  1998/09/25 13:56:41  katayama
	New function to write panther table from qq decay table

	Revision 1.2  1998/09/09 08:07:01  jtanaka
	added Ptype(idhep).

	Revision 1.1  1998/07/16 07:03:07  higuchit
	QQ2PDT -> QQREAD

*/

/* This is for C only, no namespace needed */

#ifndef PARTICLE_QQ2PANTHER_H
#define PARTICLE_QQ2PANTHER_H

#include <stdarg.h>

#ifndef __cplusplus

int qq2panther_(const char *, ...);
struct pdtent *qqgetpdt();

#else

extern "C" int qq2panther_(const char *, ...);
extern "C" struct pdtent *qqgetpdt(const char *);

#endif /* !__cplusplus */

#endif /* PARTICLE_QQ2PANTHER_H */

