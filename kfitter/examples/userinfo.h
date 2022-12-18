#include "belle.h"
#include "particle/Particle.h"
#include "particle/ParticleUserInfo.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


// For Interface to Set UserInfo Class

void setUserInfo(Particle &p);
void setUserInfo(vector<Particle> &p);

// UserInfo Class

class UserInfo : public ParticleUserInfo
{
public:
  /// Default constructor
  UserInfo();

  /// Copy constructor
  UserInfo(const UserInfo &);

  /// Destructor
  virtual ~UserInfo();

  /// constructs self object.
  UserInfo * clone(void) const;

  /// Copy operator
  UserInfo & operator = (const UserInfo &);

public:
  void chisq(const double &v) { m_chisq = v; }
  const double & chisq(void) const { return m_chisq; }

  void cl(const double &v) { m_cl = v; }
  const double & cl(void) const { return m_cl; }

  void ndf(const unsigned &v) { m_ndf = v; }
  const unsigned & ndf(void) const { return m_ndf; }

  void mass(const double &v, const int i) { m_mass[i] = v; }
  const double & mass(const int i) const { return m_mass[i]; }

  void px(const double &v, const int i) { m_px[i] = v; }
  const double & px(const int i) const { return m_px[i]; }

  void py(const double &v, const int i) { m_py[i] = v; }
  const double & py(const int i) const { return m_py[i]; }

  void pz(const double &v, const int i) { m_pz[i] = v; }
  const double & pz(const int i) const { return m_pz[i]; }

  void epx(const double &v, const int i) { m_epx[i] = v; }
  const double & epx(const int i) const { return m_epx[i]; }

  void epy(const double &v, const int i) { m_epy[i] = v; }
  const double & epy(const int i) const { return m_epy[i]; }

  void epz(const double &v, const int i) { m_epz[i] = v; }
  const double & epz(const int i) const { return m_epz[i]; }

  void type(const unsigned &v) { m_type = v; }
  const unsigned & type(void) const { return m_type; }

  void vx(const double &v, const int i) { m_vx[i] = v; }
  const double & vx(const int i) const { return m_vx[i]; }

  void vy(const double &v, const int i) { m_vy[i] = v; }
  const double & vy(const int i) const { return m_vy[i]; }

  void vz(const double &v, const int i) { m_vz[i] = v; }
  const double & vz(const int i) const { return m_vz[i]; }

  void evx(const double &v, const int i) { m_evx[i] = v; }
  const double & evx(const int i) const { return m_evx[i]; }

  void evy(const double &v, const int i) { m_evy[i] = v; }
  const double & evy(const int i) const { return m_evy[i]; }

  void evz(const double &v, const int i) { m_evz[i] = v; }
  const double & evz(const int i) const { return m_evz[i]; }

  void CL(const double &v, const int i) { m_CL[i] = v; }
  const double & CL(const int i) const { return m_CL[i]; }

  void dM(const double &v, const int i) { m_dM[i] = v; }
  const double & dM(const int i) const { return m_dM[i]; }

private:
  double m_chisq;
  double m_cl;
  unsigned m_ndf;

  double m_mass[3];
  double m_px[3];
  double m_py[3];
  double m_pz[3];
  double m_epx[3];
  double m_epy[3];
  double m_epz[3];
  double m_vx[3];
  double m_vy[3];
  double m_vz[3];
  double m_evx[3];
  double m_evy[3];
  double m_evz[3];
  double m_CL[3];
  double m_dM[3];

  unsigned m_type;
};
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
