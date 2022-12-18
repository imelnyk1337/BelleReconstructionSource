//
// $Id: kmassvertexfitter.h 9932 2006-11-12 14:26:53Z katayama $
//
// $Log$
// Revision 1.13  2002/04/22 11:07:44  jtanaka
// bug fix for virtual functions (Kataoka-san's info.).
//
// Revision 1.12  2002/03/27 23:42:30  jtanaka
// Add new fitter to the mass-constraint fit.
//
// Revision 1.11  2001/12/24 12:03:27  katayama
// gcc 3.0 and headers are cleaned up
//
// Revision 1.10  2000/03/07 17:09:14  jtanaka
// bug fixed, use "similarity" etc.
//
// Revision 1.9  2000/03/07 10:52:38  katayama
// compatibility with CC5.0
//
// Revision 1.8  1999/03/29 05:38:18  jtanaka
// new class structure, new+old interfaces
//
// Revision 1.7  1998/09/28 13:26:08  jtanaka
// bug fix in kmakemother.h and <vector.h> -> <vector>
//
// Revision 1.6  1998/09/08 11:52:07  jtanaka
// non-const --> const
//
// Revision 1.5  1998/07/23 11:58:03  jtanaka
// endif HEADER --> endif /* HEADER */
//
// Revision 1.4  1998/05/07 04:11:16  jtanaka
// version 2.1
//
// Revision 1.3  1998/01/22 03:21:41  jtanaka
// Updated from Tanaka san. New Interface etc.
//
// Revision 1.2  1997/10/24 07:05:10  katayama
// Updated from Tanaka san. Bug fixes etc.
//
//
// kmassvertexfitter of KFitter
//
// ver2.0  : 1998/01
// ver3.0  : 1999/03
//
// author  : jtanaka
// e-mail  : jtanaka@hep.phys.s.u-tokyo.ac.jp
//
#ifndef KMASSVERTEXFITTER_H
#define KMASSVERTEXFITTER_H
#include "belle.h"
#include "kfitter/kfitterbase2.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class kmassvertexfitter : public kfitterbase2 
{
public:
  kmassvertexfitter(void);
  ~kmassvertexfitter(void);

  //interface
  void initialVertex(const HepPoint3D&);
  void invariantMass(const double);
  void fixMass(void);
  void unfixMass(void);
  void correlation(const HepMatrix&);
  void correlation(void);

  //output interface
  HepPoint3D   vertex(const unsigned = KF_AFTER_FIT);
  HepSymMatrix errVertex(void);
  HepMatrix    errVertexTrack(const unsigned);
  double       invariantMass(void);
  double       chisq(void);
  double       chisq(const unsigned n);
  HepMatrix    correlation(const unsigned n, const unsigned m, const unsigned flag = KF_AFTER_FIT);

#if KF_WITH_OLD_INTERFACE
  void set_initial_vertex(const Hep3Vector &v){ initialVertex(static_cast<const HepPoint3D&>(v)); }
  void set_initial_vertex(double x,double y,double z)
    { HepPoint3D t(x,y,z); initialVertex(t); }
  void set_invariant_mass(double mass){ invariantMass(mass); }

  Hep3Vector        get_vertex(){ return static_cast<Hep3Vector>(vertex()); }
  double            get_vx()    { return vertex().x(); }
  double            get_vy()    { return vertex().y(); }
  double            get_vz()    { return vertex().z(); }
  HepSymMatrix      get_err_vertex(){ return errVertex(); }
  double            get_err_vertex(int i,int j){ return errVertex()[i][j]; }
  HepMatrix         get_err_vertex_track(int i){ return errVertexTrack(i); }
#endif //KF_WITH_OLD_INTERFACE

private:
  //before fitting
  HepPoint3D m_vertex_b;

  //after fitting
  HepPoint3D   m_vertex_a;
  HepSymMatrix m_errVertex_a;
  std::vector<HepMatrix> m_errVertexTrack_a;

  //other parameters
  double     m_invariantMass;

  //input functions
  unsigned m_setInputMatrix(void);
  unsigned m_setInputSubMatrix(void);
  unsigned m_setCorrelation(void);

  //output functions
  unsigned m_setOutputMatrix(void);
  std::vector<int> m_isFixMass;

  //fit functions
  unsigned m_makeCoreMatrix(void);
  unsigned m_calDgf(void);
};
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* KMASSVERTEXFITTER_H */
