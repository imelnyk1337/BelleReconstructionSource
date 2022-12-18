//
// $Id: kmakemother.h 9932 2006-11-12 14:26:53Z katayama $
//
// $Log$
// Revision 1.11  2002/03/27 23:42:30  jtanaka
// Add new fitter to the mass-constraint fit.
//
// Revision 1.10  2001/12/24 12:03:27  katayama
// gcc 3.0 and headers are cleaned up
//
// Revision 1.9  2000/03/07 17:09:13  jtanaka
// bug fixed, use "similarity" etc.
//
// Revision 1.8  2000/03/07 10:52:37  katayama
// compatibility with CC5.0
//
// Revision 1.7  1999/03/29 05:38:17  jtanaka
// new class structure, new+old interfaces
//
// Revision 1.6  1998/09/28 13:26:07  jtanaka
// bug fix in kmakemother.h and <vector.h> -> <vector>
//
// Revision 1.5  1998/09/08 11:52:07  jtanaka
// non-const --> const
//
// Revision 1.4  1998/07/23 11:58:02  jtanaka
// endif HEADER --> endif /* HEADER */
//
// Revision 1.3  1998/01/22 03:21:40  jtanaka
// Updated from Tanaka san. New Interface etc.
//
// Revision 1.2  1997/10/24 07:05:09  katayama
// Updated from Tanaka san. Bug fixes etc.
//
//
// kmakemother of KFitter
//
// ver2.0  : 1998/01
// ver3.0  : 1999/03
//
// author  : jtanaka
// e-mail  : jtanaka@hep.phys.s.u-tokyo.ac.jp
//
#ifndef KMAKEMOTHER_H
#define KMAKEMOTHER_H
#include "belle.h"
#include "kfitter/kfitterparticle.h"
#include <vector>
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class kmakemother
{
public:
  kmakemother(void);
  ~kmakemother(void);

  //input interface
  void addTrack(const HepLorentzVector&, 
		const HepPoint3D&,
		const HepSymMatrix&,
		const double,
		const unsigned = KF_AFTER_FIT);
  void addTrack(const kfitterparticle&);
  void magneticField(const double);
  void vertex(const HepPoint3D&);
  void errVertex(const HepSymMatrix&);
  void correlation(const HepMatrix&);
  void correlation(void);
  void errVertexTrack(const HepMatrix&);
  void errVertexTrack(void);
  void atDecayPoint(void);
  void notDecayPoint(void);
  void beforeAfter(const unsigned);

  //output interface
  kfitterparticle   track(void);
  HepLorentzVector  momentum(void);
  HepPoint3D        position(void);
  HepSymMatrix      error(void);

  //main
  unsigned make(void);

#if KF_WITH_OLD_INTERFACE
  void add_track(const HepLorentzVector &p, const Hep3Vector &x,
                 const HepSymMatrix &e, double c) { addTrack(p, static_cast<const HepPoint3D&>(x),e,c);}
  void add_track(const kfitterparticle &kp) { addTrack(kp); }
  void set_magneticfield(double m) { magneticField(m); }
  void set_vertex(const Hep3Vector &v) { vertex(static_cast<const HepPoint3D&>(v)); }
  void set_err_vertex(const HepSymMatrix& e) { errVertex(e); }
  void set_correlation(const HepMatrix &e) { correlation(e); }
  void set_correlation() { correlation(); }
  void set_vertex_track_correlation(const HepMatrix &e) { errVertexTrack(e); }
  void set_vertex_track_correlation() { errVertexTrack(); }
  void set_at_decaypoint() { atDecayPoint(); }
  void set_not_decaypoint(){ notDecayPoint(); }

  kfitterparticle   get_mother_track(){ return track(); }
  HepLorentzVector  get_4momentum()   { return momentum(); }
  Hep3Vector        get_position()    { return static_cast<Hep3Vector>(position()); }
  HepSymMatrix      get_err_track()   { return error(); }

  int start() { return static_cast<int>(make()); }
#endif //KF_WITH_OLD_INTERFACE

private:
  unsigned m_errorFlag;
  unsigned m_atDecayPoint;
  unsigned m_trackNum;
  unsigned m_errVertexFlag;
  unsigned m_correlationFlag;
  unsigned m_errVertexTrackFlag;
  unsigned m_beforeAfterFlag;
  double   m_magneticField;

  //...children
  std::vector<kfitterparticle> m_plist;
  std::vector<HepMatrix>       m_correlation;
  std::vector<HepMatrix>       m_errVertexTrack;
  
  //...vertex
  Hep3Vector              m_vertex;
  HepSymMatrix            m_errVertex;
  
  //...utilities
  void m_delMdelC(HepMatrix&);
  void m_error(HepSymMatrix&);

  //...mother
  double           m_charge;
  kfitterparticle  m_mother;
};
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* KMAKEMOTHER_H */
