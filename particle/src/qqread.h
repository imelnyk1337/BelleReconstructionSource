/* qqread.h */

/*
	$Id: qqread.h 9932 2006-11-12 14:26:53Z katayama $
	$Log$
	Revision 1.3  1998/10/14 12:50:25  jtanaka
	uses qq2panther to get Ptype information.

	Revision 1.2  1998/09/09 08:07:01  jtanaka
	added Ptype(idhep).

	Revision 1.1  1998/07/16 07:03:07  higuchit
	QQ2PDT -> QQREAD

*/

/* killed on 10/8/98 */
#if 0

#ifndef PARTICLE_QQREAD_H
#define PARTICLE_QQREAD_H
#include "belle.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


typedef char pname[16];

#define MAXNPDT 1024

struct pdtent{
	pname name;
	int pdgid;
	int stable;
	float mass /* in GeV */, charge, spin, ctau;
	float width /* FWHM in GeV */, mass_min, mass_max;
};

#define MAXNDECAY 4096

struct decay{
	pname mother;
	int matrix;
	float br;
	size_t ndaughter;
	pname daughter[8];
};

#ifndef __cplusplus

int qqread();
struct pdtent  *qqgetpdt_with_name();
struct pdtent **qqgetpdt_with_pdgid();

#else

extern "C" int qqread(const char *file_qq[], const size_t nfile_qq);
extern "C" struct pdtent  *qqgetpdt_with_name(const char *name);
extern "C" struct pdtent **qqgetpdt_with_pdgid(const int pdgid);

#endif /* !__cplusplus */

#endif /* PARTICLE_QQREAD_H */

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
