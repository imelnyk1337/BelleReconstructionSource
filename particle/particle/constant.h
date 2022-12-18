//
// $Id: constant.h 9932 2006-11-12 14:26:53Z katayama $
//
// $Log$
// Revision 1.7  2000/05/16 14:46:29  jtanaka
// added constructor of Mdst_ecl.
//
// Revision 1.6  2000/04/14 12:40:57  jtanaka
// updated for new table "mdst_vee2".
//
// Revision 1.5  2000/01/05 06:10:39  jtanaka
// major updates:please see BELLE whiteboard.
//
// Revision 1.4  1999/04/09 15:01:10  jtanaka
// Added new member functions, "const", "&"...etc. Removed some member functions.
//
// Revision 1.3  1998/07/23 12:51:45  katayama
// conform ANSI
//
// Revision 1.2  1998/07/02 09:29:10  higuchit
// null flag -> `usable' flag
//
// Revision 1.1  1998/07/01 11:56:53  jtanaka
// new header file
//
//
//
#ifndef PARTICLE_CLASS_CONSTANT_H
#define PARTICLE_CLASS_CONSTANT_H
#include "belle.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


//concerning "usable" of particle calss
//#define UNUSABLE   0
//#define USABLE     1
const unsigned UNUSABLE = 0;
const unsigned USABLE   = 1;

//MDST_TABLES
const unsigned PC_CHARGED =  1;
const unsigned PC_GAMMA   =  2;
const unsigned PC_VEE     =  4;
const unsigned PC_KLONG   =  8;
const unsigned PC_PI0     =  16;
const unsigned PC_TRK     =  32;
const unsigned PC_ELID    =  64;
const unsigned PC_HEPEVT  =  128;
const unsigned PC_VEE2    =  256;
const unsigned PC_ECL     =  512;
const unsigned PC_ALL     =  PC_CHARGED+PC_GAMMA+PC_VEE+PC_KLONG+PC_PI0+PC_TRK+PC_ELID+PC_HEPEVT+PC_VEE2+PC_ECL;

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* PARTICLE_CLASS_CONSTANT_H */
