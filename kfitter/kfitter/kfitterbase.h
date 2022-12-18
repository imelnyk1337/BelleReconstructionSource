//
// $Id: kfitterbase.h 10002 2007-02-26 06:56:17Z katayama $
//
// $Log$
// Revision 1.12  2004/04/25 16:39:00  hitoshi
// make vertices from rectrk (by Karim)
//
// Revision 1.11  2002/04/22 11:07:44  jtanaka
// bug fix for virtual functions (Kataoka-san's info.).
//
// Revision 1.10  2002/03/27 23:42:30  jtanaka
// Add new fitter to the mass-constraint fit.
//
// Revision 1.9  2002/01/03 11:04:36  katayama
// Point3D and other header files are cleaned
//
// Revision 1.8  2001/12/23 09:58:22  katayama
// removed Strings.h
//
// Revision 1.7  2000/06/07 12:31:45  jtanaka
// add m_errorMsgFlag to control "dump error message".
// Default is ON. errorMsg(0) --> OFF.
//
// Revision 1.6  2000/03/07 17:09:12  jtanaka
// bug fixed, use "similarity" etc.
//
// Revision 1.5  2000/03/07 10:52:37  katayama
// compatibility with CC5.0
//
// Revision 1.4  1999/07/19 18:06:56  jtanaka
// bug fix : chisq
//
// Revision 1.3  1999/07/19 15:58:35  jtanaka
// added Mdst_charged
//
// Revision 1.2  1999/05/26 12:02:29  jtanaka
// speed up etc.
//
// Revision 1.1  1999/03/29 05:38:15  jtanaka
// new class structure, new+old interfaces
//
//
//
// kfitterbase of KFitter
//
// ver3.0  : 1999/03
//
// author  : jtanaka
// e-mail  : jtanaka@hep.phys.s.u-tokyo.ac.jp
//
#ifndef KFITTERBASE_H
#define KFITTERBASE_H
#include "belle.h"
#include "kfitter/kfitterparticle.h"
#include "kfitter/kfittererror.h"
#include <iostream>
#include <vector>
#include "belleCLHEP/Matrix/Matrix.h"
#include "belleCLHEP/Matrix/SymMatrix.h"
#include "belleCLHEP/Vector/ThreeVector.h"
#include "belleCLHEP/Vector/LorentzVector.h"
#ifndef CLHEP_POINT3D_H
#include "belleCLHEP/Geometry/Point3D.h"
#endif
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
#if KF_WITH_BELLE_INTERFACE
class Mdst_charged;
class Rectrk;
#endif

class kfitterbase {
public:
  kfitterbase(void);
  virtual ~kfitterbase(void);

  //input interface
  void addTrack(const HepLorentzVector &p, 
                const HepPoint3D &x,
                const HepSymMatrix &e, 
                const double q);
#if KF_WITH_OLD_INTERFACE
  void addTrack(const HepLorentzVector&,
		const HepPoint3D&,
		const HepSymMatrix&, 
		const double,
		const double); // obsolete
#endif
  void addTrack(const kfitterparticle&);
#if KF_WITH_BELLE_INTERFACE
  void addTrack(const Mdst_charged&, 
		const double mass = 0.139570,
		const unsigned id = 2);
  void addTrack(const Rectrk&, 
		const double mass = 0.139570,
		const unsigned id = 2);
#endif
  void magneticField(const double);
  virtual void correlation(const HepMatrix&);
  virtual void correlation(void);
  void errorMsg(unsigned int);
  
  //fit
  virtual unsigned fit(void);

  //output interface
  unsigned          error(void);
  unsigned          tracks(void);
  unsigned          iteration(void);
  virtual unsigned  dgf(void);
  virtual double    chisq(void);
  virtual double    cl(void);
  virtual double    chisq(const unsigned);
  double            magneticField(void);
  HepLorentzVector  momentum(const unsigned);
  HepPoint3D        position(const unsigned);
  HepSymMatrix      error(const unsigned);
  kfitterparticle   track(const unsigned);
  virtual HepMatrix         correlation(const unsigned, 
					const unsigned, 
					const unsigned = KF_AFTER_FIT);

  virtual void dump(const unsigned = KF_DUMP_MEASUREMENT);

#if KF_WITH_OLD_INTERFACE
  void add_track(const HepLorentzVector &p, const Hep3Vector &x,
                 const HepSymMatrix &e, double q, double m) { addTrack(p,static_cast<const HepPoint3D&>(x),e,q,m); }
  void add_track(const kfitterparticle &kp) { addTrack(kp); }
  void set_magneticfield(double m){ magneticField(m); }
  void set_correlation(const HepMatrix &e) { correlation(e); }
  void set_correlation() { correlation(); }

  int start() { return static_cast<int>(fit()); }

  int               get_dgf() { return static_cast<int>(dgf()); }
  double            get_chisq() { return chisq(); }
  double            get_CL() { return cl(); }
  HepLorentzVector  get_4momentum(int i) { return momentum(i); }
  Hep3Vector        get_position(int i) { return static_cast<Hep3Vector>(position(i)); }
  HepSymMatrix      get_err_track(int i) { return error(i); }
  kfitterparticle   get_track(int i) { return track(i); }
  HepMatrix         get_err_track(int i, int j) { return correlation(i,j); }
#endif //KF_WITH_OLD_INTERFACE

protected:
  std::vector<kfitterparticle> m_plist;
  double                       m_magneticField;
  std::vector<HepMatrix>       m_correlation_b;
  
  //statistic parameter
  unsigned m_dgf;
  double   m_chisq;
  double   m_cl;
  
  //other parameters
  unsigned m_trackNum;
  unsigned m_necessaryTrackNum;
  unsigned m_errorFlag;
  unsigned m_correlationFlag;
  unsigned m_overIterationFlag;
  unsigned m_errorMsgFlag;

  //input functions
  virtual unsigned m_setInputMatrix(void) = 0;
  virtual unsigned m_setCorrelation(void);

  //output functions
  virtual unsigned m_setOutputMatrix(void) = 0;
  HepSymMatrix m_makeError(const HepLorentzVector&,
			   const HepMatrix&);
  HepMatrix    m_makeError(const HepLorentzVector&,
			   const HepLorentzVector&,
			   const HepMatrix&);
  HepMatrix    m_makeError2(const HepLorentzVector&,
			    const HepMatrix&);

  HepSymMatrix m_makeError3(const HepLorentzVector &,
                            const HepMatrix &,
                            const int isFixMass);
  HepMatrix m_makeError3(const HepLorentzVector &p1,
                         const HepLorentzVector &p2,
                         const HepMatrix &e,
                         const int isFixMass1,
                         const int isFixMass2);
  HepMatrix    m_makeError4(const HepLorentzVector &p,
                            const HepMatrix &e);

  //fit functions
  virtual unsigned m_makeCoreMatrix(void) = 0;
  virtual unsigned m_calDgf(void) = 0;

  //matrix
  HepSymMatrix m_V_al_0;
  HepMatrix    m_al_0;
  HepMatrix    m_al_1;
  HepMatrix    m_al_a;
  HepMatrix    m_property;

  HepMatrix    m_D;
  HepMatrix    m_d;

  HepMatrix    m_V_D;

  HepMatrix    m_V_al_1;
  HepMatrix    m_lam;
};

double chisq2Confi(const int, const double);

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* KFITTERBASE_H */
