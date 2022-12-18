//
// $Id: kfitterparticle.h 10002 2007-02-26 06:56:17Z katayama $
//
// $Log$
// Revision 1.7  2002/03/27 23:42:30  jtanaka
// Add new fitter to the mass-constraint fit.
//
// Revision 1.6  2002/01/03 11:04:36  katayama
// Point3D and other header files are cleaned
//
// Revision 1.5  2000/03/07 17:09:13  jtanaka
// bug fixed, use "similarity" etc.
//
// Revision 1.4  1999/03/29 05:38:16  jtanaka
// new class structure, new+old interfaces
//
// Revision 1.3  1998/09/08 11:52:07  jtanaka
// non-const --> const
//
// Revision 1.2  1998/07/23 11:58:01  jtanaka
// endif HEADER --> endif /* HEADER */
//
// Revision 1.1  1998/01/22 03:21:39  jtanaka
// Updated from Tanaka san. New Interface etc.
//
//
// kfitterparticle of KFitter
//
// ver1.0  : ?
// ver3.0  : 1999/03
//
// author  : jtanaka
// e-mail  : jtanaka@hep.phys.s.u-tokyo.ac.jp
//
#ifndef KFITTERPARTICLE_H
#define KFITTERPARTICLE_H
#include "belle.h"
#include "kfitter/kfitter_ini.h"
#include <iostream>
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

class kfitterparticle{
public:
  //constructor
  kfitterparticle(void);
  kfitterparticle(const HepLorentzVector &p, 
                  const HepPoint3D &x, 
                  const HepSymMatrix &e,
                  const double q,
                  const unsigned = KF_BEFORE_FIT);
#if KF_WITH_OLD_INTERFACE
  kfitterparticle(const HepLorentzVector&, 
		  const HepPoint3D&, 
		  const HepSymMatrix&,
		  const double,
		  const double,
		  const unsigned = KF_BEFORE_FIT); // This is obsolete
  // since mass is obtained from HepLorentzVector's information.
#endif
  kfitterparticle(const kfitterparticle&);
  //destructor
  virtual ~kfitterparticle(void){};
  
  //input
  void momentum(const HepLorentzVector&,
		const unsigned = KF_BEFORE_FIT);
  void position(const HepPoint3D&,
		const unsigned = KF_BEFORE_FIT);
  void error(const HepSymMatrix&,
	     const unsigned = KF_BEFORE_FIT);
  void charge(const double);
  //void mass(const double);
  void vertex(const HepPoint3D&);
  void errVertex(const HepSymMatrix&);

  //output
  HepLorentzVector momentum(const unsigned = KF_AFTER_FIT);
  HepLorentzVector momentum(const unsigned = KF_AFTER_FIT) const;
  HepPoint3D       position(const unsigned = KF_AFTER_FIT);
  HepPoint3D       position(const unsigned = KF_AFTER_FIT) const;
  HepSymMatrix     error(const unsigned = KF_AFTER_FIT);
  HepSymMatrix     error(const unsigned = KF_AFTER_FIT) const;
  double           charge(void);
  double           charge(void) const;
  double           mass(void);
  double           mass(void) const;
  HepPoint3D       vertex(void);
  HepPoint3D       vertex(void) const;
  HepSymMatrix     errVertex(void);
  HepSymMatrix     errVertex(void) const;
  HepMatrix        mompos(const unsigned); // momentum+position
  HepMatrix        mompos(const unsigned) const;

  //fitter
  double           getFitParameter(const unsigned,
				   const unsigned);
  double           getFitParameter(const unsigned,
				   const unsigned) const;
  HepMatrix        getFitParameter(const unsigned);
  HepMatrix        getFitParameter(const unsigned) const;
  HepSymMatrix     getFitError(const unsigned);
  HepSymMatrix     getFitError(const unsigned) const;

  //operator
  kfitterparticle & operator = (const kfitterparticle &);

#if KF_WITH_OLD_INTERFACE
  void set_4momentum(const HepLorentzVector &p) { momentum(p); }
  void set_position(const Hep3Vector &x) { position(static_cast<const HepPoint3D&>(x)); }
  void set_error(const HepSymMatrix &e) { error(e); }
  void set_4momentum(int f, const HepLorentzVector &p) { f==0 ? momentum(p) : momentum(p,KF_AFTER_FIT); }
  void set_position(int f, const Hep3Vector &x) { f==0 ? position(static_cast<const HepPoint3D&>(x)) : position(static_cast<const HepPoint3D&>(x),KF_AFTER_FIT); }
  void set_error(int f, const HepSymMatrix &e) { f==0 ? error(e) : error(e,KF_AFTER_FIT); }
  //void set_charge_mass(double c,double m) { charge(c); mass(m); }
  void set_vertex(const Hep3Vector &v) { vertex(static_cast<const HepPoint3D&>(v)); }
  void set_err_vertex(const HepSymMatrix &e) { errVertex(e); }
   
  HepLorentzVector get_4momentum() const { return momentum(); }
  Hep3Vector       get_position() const { return static_cast<Hep3Vector>(position()); }
  HepSymMatrix     get_error() const { return error(); }
  HepLorentzVector get_4momentum(int f) const { return f==0 ? momentum(KF_BEFORE_FIT) : momentum(); }
  Hep3Vector       get_position(int f) const { return f==0 ? static_cast<Hep3Vector>(position(KF_BEFORE_FIT)) : static_cast<Hep3Vector>(position()); }
  HepSymMatrix     get_error(int f) const { return f==0 ? error(KF_BEFORE_FIT) : error(); }
  double           get_charge() const { return charge(); }
  double           get_mass() const { return mass(); }
  double           get_fit_parameter(int i) const { return getFitParameter(i,KF_BEFORE_FIT); }
  HepSymMatrix     get_fit_error() const { return getFitError(KF_BEFORE_FIT); }
  Hep3Vector       get_vertex() const { return static_cast<Hep3Vector>(vertex()); }
  HepSymMatrix     get_err_vertex() const { return errVertex(); }
#endif // KF_WITH_OLD_INTERFACE

private:
  //before fit
  HepLorentzVector   m_momentum_b;
  HepPoint3D         m_position_b;
  HepSymMatrix       m_error_b;
  //after fit
  HepLorentzVector   m_momentum_a;
  HepPoint3D         m_position_a;
  HepSymMatrix       m_error_a;
  //properties
  double             m_charge;
  double             m_mass;
  //vertex
  HepPoint3D         m_vertex;
  HepSymMatrix       m_errVertex;
};
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* KFITTERPARTILCE_H */
