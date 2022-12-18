//
// $Id: kfitter_ini.h 9932 2006-11-12 14:26:53Z katayama $
//
// $Log$
// Revision 1.7  2002/03/27 23:42:29  jtanaka
// Add new fitter to the mass-constraint fit.
//
// Revision 1.6  1999/12/10 15:50:10  jtanaka
// add refitter for slowpion to kvertexfitter and BF is changed sligthly
//
// Revision 1.5  1999/07/19 15:58:35  jtanaka
// added Mdst_charged
//
// Revision 1.4  1999/05/26 12:02:29  jtanaka
// speed up etc.
//
// Revision 1.3  1999/03/29 05:38:14  jtanaka
// new class structure, new+old interfaces
//
// Revision 1.2  1998/07/23 11:58:00  jtanaka
// endif HEADER --> endif /* HEADER */
//
// Revision 1.1  1998/01/22 03:21:38  jtanaka
// Updated from Tanaka san. New Interface etc.
//
//
// kfitter_ini of KFitter
//
// ver1.0  : 1998/01
// ver3.0  : 1999/03
//
// author  : jtanaka
// e-mail  : jtanaka@hep.phys.s.u-tokyo.ac.jp
//
#ifndef KFITTER_INI_H
#define KFITTER_INI_H
#include "belle.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//
//
//
const unsigned KF_PARAMETER_NUMBER    =  6;
const unsigned KF_SUBPARAMETER_NUMBER =  7;
const unsigned KF_NUM6                =  6;
const unsigned KF_NUM7                =  7;
const unsigned KF_BEFORE_FIT          =  0;
const unsigned KF_AFTER_FIT           =  1;
const unsigned KF_AT_DECAY_POINT      =  0;
const unsigned KF_NOT_DECAY_POINT     =  1;
const unsigned KF_NO_OVER_ITERATION   =  0;
const unsigned KF_OVER_ITERATION      =  1;
const double   KF_PHOTON_VELOCITY     =  0.00299792458;

const int      KF_FIX_MASS            =  1;
const int      KF_UNFIX_MASS          =  0;

const unsigned KF_MODE_WO_CORRELATION   = 0;
const unsigned KF_MODE_WITH_CORRELATION = 1;

const unsigned KF_MAX_TRACK_NUMBER     = 10;
const unsigned KF_MAX_TRACK_NUMBER2    = 30;
const unsigned KF_MAX_ITERATION_NUMBER = 15;
const double   KF_INIT_CHI2            = 1.0e+30;

//const double   KF_MAGNETIC_FIELD       = 1.5;
// Bfield 21 : ave of (0,0,10cm), (0,0,0), (0,0,-10cm) = (0,0,14.5361)
//const double   KF_MAGNETIC_FIELD       = 1.45361;
// Bfield 21 : ave of (0,0,1cm), (0,0,0), (0,0,-1cm) = (0,0,14.5638)
const double   KF_MAGNETIC_FIELD       = 1.45638;

#define KF_WITH_OLD_INTERFACE 1
#define KF_WITH_BELLE_INTERFACE 1
const unsigned KF_DUMP_MEASUREMENT = 0;
const unsigned KF_DUMP_CORE_MATRIX = 1;
const unsigned KF_DUMP_FITTED      = 2;
//
//
//
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif  /* KFITTER_INI_H */
