//
// 1999/05/28 16:00:00  skkim
// Added good_gamma 
/*
 * $Id: mdst.h 10613 2008-09-03 11:36:24Z katayama $
 *
 * $Log$
 * Revision 1.23  2003/04/15 02:48:15  katayama
 * pi0resol added
 *
 * Revision 1.22  2002/03/14 01:19:31  katayama
 * correction functions are mvoed to fix_mdst
 *
 * Revision 1.21  2001/12/25 18:55:13  hitoshi
 * added mdst_good_event function.
 *
 * Revision 1.20  2001/12/15 06:39:03  hitoshi
 * updated tracing back to hepevt (by Kakuno).
 *
 * Revision 1.19  2001/12/12 02:03:57  hitoshi
 * updated for mdst_ecl/trk (by Kakuno).
 *
 * Revision 1.18  2001/12/06 07:48:45  hitoshi
 * added void scale_error(void);
 *
 * Revision 1.17  2001/12/05 07:39:42  hitoshi
 * added make_pi0 and correct_ecl declarations.
 *
 * Revision 1.16  2001/12/04 12:09:33  hitoshi
 * added scale_error, scale gamma energy, benergy.
 *
 * Revision 1.15  2001/11/05 12:39:27  hitoshi
 * updated so that it can handle factor per 100 runs.
 *
 * Revision 1.14  2001/01/31 11:48:30  hitoshi
 * added a function mdst2mdst (by Hamasaki).
 *
 * Revision 1.13  2000/12/10 03:18:14  hitoshi
 * added a function remove_duplicates.
 *
 * Revision 1.12  2000/05/19 09:27:43  hitoshi
 * added a function copy_vee_to_vee2
 *
 * Revision 1.11  2000/05/02 10:09:20  hitoshi
 * added momentum scaling function.
 *
 * Revision 1.10  1999/07/13 03:16:27  katayama
 * added rectrk->mdst_charged (temporary)
 *
 * Revision 1.9  1999/07/10 12:17:39  katayama
 * added trk->charged
 *
 * Revision 1.8  1999/06/29 13:25:17  katayama
 * from pcs
 *
 * Revision 1.7  1999/05/29 05:05:40  katayama
 * from skim san
 *
 *
 * Revision 1.6  1999/04/27 23:45:14  katayama
 * Use const keyword
 *
 * Revision 1.5  1999/01/16 10:31:10  katayama
 * clean up includes
 *
 * Revision 1.4  1999/01/13 15:24:25  hitoshi
 * added default cut values.
 *
 * Revision 1.3  1999/01/13 00:15:51  katayama
 * added good_charged
 *
 * Revision 1.2  1997/09/19 01:07:45  katayama
 * Added Id and Log
 *
 */
#if !defined(ANAL_MDST_H)
#define ANAL_MDST_H
#include "belle.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


class Mdst_charged;
class Mdst_ecl;
class Mdst_gamma;
class Gen_hepevt;
class Mdst_trk;
class Rectrk;
class Mdst_vee2;
class Mdst_pi0;
  

/// returns correponding hepevt
const Gen_hepevt &get_hepevt(const Mdst_trk&, int ith=0);

/// returns correponding hepevt
const Gen_hepevt &get_hepevt(const Mdst_charged&, int ith=0);

/// returns correponding hepevt(This hepevt might be made at inside of ECL).
/// the following example will return what most of you want!
///   Mdst_ecl & e = ...
///   const Gen_hepevt & h(gen_level(get_hepevt(e)));
const Gen_hepevt &get_hepevt(const Mdst_ecl&, int ith=0);

/// returns correponding hepevt(This hepevt might be made at inside of ECL).
/// the following example will return what most of you want!
///   Mdst_gamma & g = ...
///   const Gen_hepevt & h(gen_level(get_hepevt(g)));
const Gen_hepevt &get_hepevt(const Mdst_gamma&, int ith=0);

/// returns correponding hepevt.
/// the following example will return what most of you want!
///   Mdst_vee2 & g = ...
///   const Gen_hepevt & h(gen_level(get_hepevt(g)));
const Gen_hepevt &get_hepevt(const Mdst_vee2&, int ith=0);

/// returns correponding hepevt.
/// the following example will return what most of you want!
///   Mdst_vee2 & g = ...
///   const Gen_hepevt & h(gen_level(get_hepevt(g)));
const Gen_hepevt &get_hepevt(const Mdst_pi0&, int ith=0);


/// trace back to mother hepevt which has isthep>0 or 
/// isthep =-10(in case of daughters of Ks/Lambda)
const Gen_hepevt &gen_level(const Gen_hepevt &);

Mdst_charged &mdst_charged(const Mdst_trk&, int ith=0);

Mdst_charged &mdst_charged(const Rectrk&, int ith=0);

bool good_charged(const Mdst_charged&, float cl_cut=1.e-25, 
float dz_cut=40., float dr_cut=30.);

bool good_gamma(const Mdst_gamma&, float ecut=0.02, float e925cut=0.75, 
float widcut=6.0, int nhcut=0, float nncut=0.0);

bool good_gamma(const Mdst_ecl&, float ecut=0.02, float e925cut=0.75, 
float widcut=6.0, int nhcut=0, float nncut=0.0);

void copy_vee_to_vee2 (void);

void remove_duplicates (void);

void mdst2mdst ( bool fillMdstCharged=true, bool fillMdstKlong=true,
		 bool fillMdstGamma=true,   bool fillMdstPi0=true,
		 bool fillMdstEclCr=true );

void mdst2xref ( bool fillTrk=true, bool fillEcl=true, bool fillKlong=true );

double Benergy();

// pi0 mass resolution function.
double Pi0resol(double p, double theta, char* side, bool mcdata, 
		int exp, int option );

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif





