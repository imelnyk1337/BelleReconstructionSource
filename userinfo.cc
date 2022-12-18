#include "userinfo.h"

// #if defined(BELLE_NAMESPACE)
namespace Belle {
// #endif

// UserInfo Class

void UserInfo::init() {
    m_msComb = -1.;
    m_msKvf = -1.;
    m_msKmv = -1;
    m_chisq = -1.;
    m_chisqKvf = -1.;
    m_cl = -1.;
    m_clKvf = -1.;
    m_dist2IP = -1.;
    m_dist2IPmvf = -1.;
    m_dist2Mother = -1.;
    m_useTube = false;
    m_useKmvf = false;
    m_isAdoptCut = true;
    double m_wMass = 1000.;
    double m_maxChi2 = 1.e6;
    double m_helicity = -1.;
    double m_pid = -1.;
}
    
    
UserInfo::UserInfo() {
    init();
}

UserInfo::~UserInfo() {}

UserInfo::UserInfo(Particle& a_particle) {
    init();
    particle = &a_particle;
}

UserInfo::UserInfo(const UserInfo& x) : ParticleUserInfo(x) {
    m_msComb = x.m_msComb;
    m_msKvf = x.m_msKvf;
    m_msKmv = x.m_msKmv; 
    m_chisq = x.m_chisq; 
    m_chisqKvf = x.m_chisqKvf; 
    m_cl = x.m_cl;
    m_clKvf = x.m_clKvf;
    m_dist2IP = x.m_dist2IP;
    m_dist2IPmvf = x.m_dist2IPmvf;
    m_dist2Mother = x.m_dist2Mother;
    m_useTube = x.m_useTube;
    m_useKmvf = x.m_useKmvf;
    m_isAdoptCut = x.m_isAdoptCut;
    m_wMass = x.m_wMass;
    m_maxChi2 = x.m_maxChi2;
    m_helicity = x.m_helicity;
    m_pid = x.m_pid;
    particle = x.particle;
}

UserInfo* UserInfo::clone(void) const {
    UserInfo* x = new UserInfo(*this);
    return x;
}

UserInfo& UserInfo::operator = (const UserInfo &x) {
    return *(new UserInfo(x));
}

void UserInfo::msComb(double v) {m_msComb = v;}
double UserInfo::msComb() const {return m_msComb;}

void UserInfo::msKvf(double v) {m_msKvf = v;}
double UserInfo::msKvf() const {return m_msKvf;}

void UserInfo::msKmv(double v) {m_msKmv = v;};
double UserInfo::msKmv() const {return m_msKmv;};

void UserInfo::chisq(double v) {m_chisq = v;}
double UserInfo::chisq() const {return m_chisq;}

void UserInfo::chisqKvf(double v) {m_chisqKvf = v;}
double UserInfo::chisqKvf() const {return m_chisqKvf;}

void UserInfo::cl(double v) {m_cl = v;}
double UserInfo::cl() const {return m_cl;}

void UserInfo::clKvf(double v) {m_clKvf = v;}
double UserInfo::clKvf() const {return m_clKvf;}

void UserInfo::dist2IP(double v) {m_dist2IP = v;}
double UserInfo::dist2IP() const {return m_dist2IP;}

void UserInfo::dist2IPmvf(double v) {m_dist2IPmvf = v;}
double UserInfo::dist2IPmvf() const {return m_dist2IPmvf;}

void UserInfo::dist2Mother(double v) {m_dist2Mother = v;}
double UserInfo::dist2Mother() const {return m_dist2Mother;}

void UserInfo::useTube(bool v) {m_useTube = v;}
bool UserInfo::useTube() const {return m_useTube;}

void UserInfo::useKmvf(bool v) {m_useKmvf = v;}
bool UserInfo::useKmvf() const {return m_useKmvf;}

void UserInfo::isAdoptCut(bool v) {m_isAdoptCut = v;}
bool UserInfo::isAdoptCut() const {return m_isAdoptCut;}

void UserInfo::wMass(double v) {m_wMass = v;}
double UserInfo::wMass() const {return m_wMass;}

void UserInfo::maxChi2(double v) {m_maxChi2 = v;}
double UserInfo::maxChi2() const {return m_maxChi2;}

void UserInfo::helicity(double v) {m_helicity = v;}
double UserInfo::helicity() const {return m_helicity;}

void UserInfo::probpid(double v) {m_pid = v;}
double UserInfo:: probpid() const {return m_pid;}


// #if defined(BELLE_NAMESPACE)
} // namespace Belle
// #endif
