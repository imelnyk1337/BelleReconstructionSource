//
// $Id: kmassfitter.h 9932 2006-11-12 14:26:53Z katayama $
//
// $Log$
// Revision 1.12  2002/04/22 11:07:44  jtanaka
// bug fix for virtual functions (Kataoka-san's info.).
//
// Revision 1.11  2002/03/27 23:42:30  jtanaka
// Add new fitter to the mass-constraint fit.
//
// Revision 1.10  2000/03/08 07:44:00  jtanaka
// add new member function.
//
// Revision 1.9  2000/03/07 17:09:13  jtanaka
// bug fixed, use "similarity" etc.
//
// Revision 1.8  2000/03/07 10:52:38  katayama
// compatibility with CC5.0
//
// Revision 1.7  1999/03/29 05:38:18  jtanaka
// new class structure, new+old interfaces
//
// Revision 1.6  1998/09/28 13:26:08  jtanaka
// bug fix in kmakemother.h and <vector.h> -> <vector>
//
// Revision 1.5  1998/09/08 11:52:07  jtanaka
// non-const --> const
//
// Revision 1.4  1998/07/23 11:58:03  jtanaka
// endif HEADER --> endif /* HEADER */
//
// Revision 1.3  1998/01/22 03:21:40  jtanaka
// Updated from Tanaka san. New Interface etc.
//
// Revision 1.2  1997/10/24 07:05:09  katayama
// Updated from Tanaka san. Bug fixes etc.
//
//
// kmassfitter of KFitter
//
// ver2.0  : 1998/01
// ver3.0  : 1999/03
//
// author  : jtanaka
// e-mail  : jtanaka@hep.phys.s.u-tokyo.ac.jp
//
#ifndef KMASSFITTER_H
#define KMASSFITTER_H
#include "belle.h"
#include "kfitter/kfitterbase.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class kmassfitter : public kfitterbase 
{
public:
  kmassfitter(void);
  ~kmassfitter(void);
  
  //input interface
  void vertex(const HepPoint3D&);
  void errVertex(const HepSymMatrix&);
  void errVertexTrack(const HepMatrix&);
  void errVertexTrack(void);
  void invariantMass(const double);
  void atDecayPoint(void);
  void notDecayPoint(void);
  void fixMass(void);
  void unfixMass(void);
  void correlation(const HepMatrix&);
  void correlation(void);

  //output interface
  HepPoint3D   vertex(const unsigned = KF_AFTER_FIT);
  HepSymMatrix errVertex(const unsigned = KF_AFTER_FIT);
  HepMatrix    errVertexTrack(const unsigned, 
			      const unsigned = KF_AFTER_FIT);
  double       invariantMass(void);
  unsigned     decayPoint(void);
  unsigned     fitWithVertex(void);
  double       chisq(void);
  double       chisq(const unsigned n);
  HepMatrix    correlation(const unsigned n, const unsigned m, const unsigned flag = KF_AFTER_FIT);

#if KF_WITH_OLD_INTERFACE
  void set_vertex(const Hep3Vector &v){ vertex(static_cast<const HepPoint3D&>(v)); }
  void set_vertex(double x,double y,double z)
    { HepPoint3D t(x,y,z); vertex(t); }
  void set_invariant_mass(double mass){ invariantMass(mass); }
  void set_at_decaypoint() { atDecayPoint(); }
  void set_not_decaypoint(){ notDecayPoint(); }
#endif //KF_WITH_OLD_INTERFACE

private:
  //parameter before fitting
  HepPoint3D   m_vertex_b;
  HepSymMatrix m_errVertex_b;
  std::vector<HepMatrix> m_errVertexTrack_b;
  unsigned m_errVertexTrackFlag;

  //parameter after fitting
  HepPoint3D   m_vertex_a;
  HepSymMatrix m_errVertex_a;
  std::vector<HepMatrix> m_errVertexTrack_a;

  //other parameters
  double   m_invariantMass;
  unsigned m_fitIncludingVertex;
  unsigned m_atDecayPoint;
  std::vector<int> m_isFixMass;

  //input functions
  unsigned m_setInputMatrix(void);
  unsigned m_setCorrelation(void);

  //output functions
  unsigned m_setOutputMatrix(void);

  //fit functions
  unsigned m_makeCoreMatrix(void);
  unsigned m_calDgf(void);
};
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* KMASSFITTER_H */
