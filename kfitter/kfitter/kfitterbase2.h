//
// $Id: kfitterbase2.h 9932 2006-11-12 14:26:53Z katayama $
//
// $Log$
// Revision 1.2  2002/03/27 23:42:30  jtanaka
// Add new fitter to the mass-constraint fit.
//
// Revision 1.1  1999/03/29 05:38:15  jtanaka
// new class structure, new+old interfaces
//
//
//
// kfitterbase2 of KFitter
//
// ver3.0  : 1999/03
//
// author  : jtanaka
// e-mail  : jtanaka@hep.phys.s.u-tokyo.ac.jp
//
#ifndef KFITTERBASE2_H
#define KFITTERBASE2_H
#include "belle.h"
#include "kfitter/kfitterbase.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class kfitterbase2 : public kfitterbase 
{
public:
  kfitterbase2(void){};
  virtual ~kfitterbase2(void){};

  virtual unsigned fit(void);

  virtual void dump(const unsigned = KF_DUMP_MEASUREMENT);

protected:
  virtual unsigned m_setInputSubMatrix(void) = 0;

  //matrix
  HepMatrix    m_E;
  HepMatrix    m_V_E;
  HepMatrix    m_lam0;
  HepMatrix    m_v;
  HepMatrix    m_v_a;
  HepMatrix    m_V_Dt;
  HepMatrix    m_Cov_v_al_1;
};
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* KFITTERBASE2_H */
