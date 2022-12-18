//
// $Id: kfittererror.h 9932 2006-11-12 14:26:53Z katayama $
//
// $Log$
// Revision 1.5  2002/03/27 23:42:30  jtanaka
// Add new fitter to the mass-constraint fit.
//
// Revision 1.4  1999/05/26 12:02:29  jtanaka
// speed up etc.
//
// Revision 1.3  1999/03/29 05:38:16  jtanaka
// new class structure, new+old interfaces
//
// Revision 1.2  1998/07/23 11:58:01  jtanaka
// endif HEADER --> endif /* HEADER */
//
// Revision 1.1  1998/01/22 03:21:39  jtanaka
// Updated from Tanaka san. New Interface etc.
//
//
// kfittererror of KFitter
//
// ver1.0  : 1998/01
// ver3.0  : 1999/03
//
// author  : jtanaka
// e-mail  : jtanaka@hep.phys.s.u-tokyo.ac.jp
//
#ifndef KFITTERERROR_H
#define KFITTERERROR_H
#include "belle.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//No error
const unsigned int KF_NO_ERROR   = 0;
//
//Main
//
const unsigned int KF_DIV_ZERO   = 1;
const unsigned int KF_TRACK_SIZE = 2;
const unsigned int KF_INVERSE    = 3;
const unsigned int KF_INIT_CHISQ = 4;
const unsigned int KF_ARCSIN     = 5;
//
//Input Interface
//
const unsigned int KF_INPUT_DIV_ZERO    = 11;
const unsigned int KF_INPUT_TRACK_SIZE  = 12;
const unsigned int KF_INPUT_MATRIX_SIZE = 13;
const unsigned int KF_INPUT_CORRE_SIZE  = 14;
//
//Output Interface
//
const unsigned int KF_OUTPUT_DIV_ZERO   = 21;
const unsigned int KF_OUTPUT_TRACK_NUM  = 22;
const unsigned int KF_OUTPUT_OUT_RANGE  = 23;
const unsigned int KF_OUTPUT_INVERSE    = 24;
//
//
//
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif  /* KFITTERERROR_H */
