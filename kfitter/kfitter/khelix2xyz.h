//
// $Id: khelix2xyz.h 10002 2007-02-26 06:56:17Z katayama $
//
// $Log$
// Revision 1.9  2002/03/27 23:42:30  jtanaka
// Add new fitter to the mass-constraint fit.
//
// Revision 1.8  2002/01/03 11:04:37  katayama
// Point3D and other header files are cleaned
//
// Revision 1.7  1999/03/29 05:38:17  jtanaka
// new class structure, new+old interfaces
//
// Revision 1.6  1998/07/23 11:58:02  jtanaka
// endif HEADER --> endif /* HEADER */
//
// Revision 1.5  1998/06/30 21:06:05  jtanaka
// modification : this class uses helix class.
//
// Revision 1.4  1998/05/07 04:11:16  jtanaka
// version 2.1
//
// Revision 1.3  1998/01/22 03:21:40  jtanaka
// Updated from Tanaka san. New Interface etc.
//
// Revision 1.2  1997/10/24 07:05:08  katayama
// Updated from Tanaka san. Bug fixes etc.
//
//
// khelix2xyz of KFitter
//
// ver1.5  : 1997/11/24
// ver1.6  : 1998/06/17 We do not need this class 
//                      because we can use Helix class in belle soft.
//                      In this version1.6, I use Helix class.
//                      If we use this class, we waste CPU Time.
//                      I recommend to use Helix class.
// ver3.0  : 1999/03
//
// author  : jtanaka
// e-mail  : jtanaka@hep.phys.s.u-tokyo.ac.jp
//
#ifndef KHELIX2XYZ_H
#define KHELIX2XYZ_H 
#include "belle.h"
#include "helix/Helix.h"
#include "belleCLHEP/Matrix/Vector.h"
#include "belleCLHEP/Matrix/Matrix.h"
#include "belleCLHEP/Matrix/SymMatrix.h"
#include "belleCLHEP/Vector/LorentzVector.h"
#ifndef CLHEP_POINT3D_H
#include "belleCLHEP/Geometry/Point3D.h"
#endif
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class helix2Xyz{
public:
  helix2Xyz(void);
  ~helix2Xyz(void){};
  unsigned update(void);

  //interface
  void helix(const double&, 
	     const double&, 
	     const double&, 
	     const double&, 
	     const double&);
  void helix(const unsigned&, 
	     const double&);
  void helix(const double*);
  void errHelix(const unsigned&,
		const double&);
  void errHelix(const double*);
  void mass(const double&);
  void pivot(const HepPoint3D&);

  HepLorentzVector momentum(void);
  HepPoint3D       position(void);
  HepSymMatrix     errMomentumPosition(void);
  double           xyz(const unsigned&);
  double           errXyz(const unsigned&, const unsigned&);

private:
  unsigned         m_errorFlag;
  double           m_mass;
  HepVector        m_helix;
  HepSymMatrix     m_errHelix;
  HepPoint3D       m_pivot;
  HepLorentzVector m_momentum;
  HepPoint3D       m_position;
  HepSymMatrix     m_errMomentumPosition;
};
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* KHELIX2XYZ_H */
