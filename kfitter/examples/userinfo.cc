#include "belle.h"
#include "userinfo.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


// For Interface to Set UserInfo Class

void
setUserInfo(Particle &p)
{
  p.userInfo(UserInfo());
}

void 
setUserInfo(vector<Particle> &p)
{
  for(int i=0;i<p.size();++i)setUserInfo(p[i]);
}

// UserInfo Class

UserInfo::UserInfo()
  : m_chisq(-1.), m_cl(-1.), m_ndf(0), 
    m_type(0)
{
  for(int i=0;i<3;++i){
    m_mass[i] = 0.;
    m_px[i]   = 0.;
    m_py[i]   = 0.;
    m_pz[i]   = 0.;
    m_epx[i]   = 0.;
    m_epy[i]   = 0.;
    m_epz[i]   = 0.;
    m_vx[i]   = 0.;
    m_vy[i]   = 0.;
    m_vz[i]   = 0.;
    m_evx[i]  = 0.;
    m_evy[i]  = 0.;
    m_evz[i]  = 0.;
    m_CL[i]   = 0.;
    m_dM[i]   = 0.;
  }
}

UserInfo::~UserInfo()
{
}

UserInfo::UserInfo(const UserInfo &x)
  : m_chisq(x.m_chisq), m_cl(x.m_cl), m_ndf(x.m_ndf),
    m_type(x.m_type)
{
  for(int i=0;i<3;++i){
    m_mass[i] = x.m_mass[i];
    m_px[i]   = x.m_px[i];
    m_py[i]   = x.m_py[i];
    m_pz[i]   = x.m_pz[i];
    m_epx[i]  = x.m_epx[i];
    m_epy[i]  = x.m_epy[i];
    m_epz[i]  = x.m_epz[i];
    m_vx[i]   = x.m_vx[i];
    m_vy[i]   = x.m_vy[i];
    m_vz[i]   = x.m_vz[i];
    m_evx[i]  = x.m_evx[i];
    m_evy[i]  = x.m_evy[i];
    m_evz[i]  = x.m_evz[i];
    m_CL[i]   = x.m_CL[i];
    m_dM[i]   = x.m_dM[i];
  }
}

UserInfo*
UserInfo::clone(void) const
{
  UserInfo *x = new UserInfo(*this);
  return x;
}

UserInfo &
UserInfo::operator = (const UserInfo &x)
{
  m_chisq = x.m_chisq;
  m_cl    = x.m_cl;
  m_ndf   = x.m_ndf;
  for(int i=0;i<3;++i){
    m_mass[i] = x.m_mass[i];
    m_px[i]   = x.m_px[i];
    m_py[i]   = x.m_py[i];
    m_pz[i]   = x.m_pz[i];
    m_epx[i]  = x.m_epx[i];
    m_epy[i]  = x.m_epy[i];
    m_epz[i]  = x.m_epz[i];
    m_vx[i]   = x.m_vx[i];
    m_vy[i]   = x.m_vy[i];
    m_vz[i]   = x.m_vz[i];
    m_evx[i]  = x.m_evx[i];
    m_evy[i]  = x.m_evy[i];
    m_evz[i]  = x.m_evz[i];
    m_CL[i]   = x.m_CL[i];
    m_dM[i]   = x.m_dM[i];
  }
  m_type = x.m_type;
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
