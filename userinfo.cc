#include "userinfo.h"

namespace Belle {

    // Object initialization
    void UserInfo::init() {
        m_msComb          = -1.;
        m_msKvf           = -1.;
        m_msKmvf          = -1.;
        m_chisqKvf        = -1.;
        m_chisqKvfNdf     = -1.;
        m_chisqKmvf       = -1.;
        m_chisqKmvfNdf    = -1.;
        m_clKvf           = -1.;
        m_clKmvf          = -1.;
        m_dist2IP         = -1.;
        m_dist2IPKmvf     = -1.;
        m_dist2Mother     = -1.;
        m_useTube         = false;
        m_useKmvf         = false;
        m_isAdoptCut      = true;
        m_wMass    = 1000.;
        m_helicity = -1.;
        m_pid      = -1.;
        m_vertexMode = "vertex-fit";
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
        m_msComb        = x.m_msComb;
        m_msKvf         = x.m_msKvf;
        m_msKmvf        = x.m_msKmvf;
        m_chisqKvf      = x.m_chisqKvf;
        m_chisqKvfNdf   = x.m_chisqKvfNdf;
        m_chisqKmvfNdf = x.m_chisqKmvfNdf;
        m_chisqKmvf     = x.m_chisqKmvf;
        m_clKvf         = x.m_clKvf;
        m_clKmvf        = x.m_clKmvf;
        m_dist2IP       = x.m_dist2IP;
        m_dist2IPKmvf   = x.m_dist2IPKmvf;
        m_dist2Mother   = x.m_dist2Mother;
        m_useTube       = x.m_useTube;
        m_useKmvf       = x.m_useKmvf;
        m_isAdoptCut    = x.m_isAdoptCut;
        m_wMass         = x.m_wMass;
        m_helicity      = x.m_helicity;
        m_pid           = x.m_pid;
        m_vertexMode    = x.m_vertexMode;
        particle        = x.particle;
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

    void UserInfo::msKmvf(double v) {m_msKmvf = v;}
    double UserInfo::msKmvf() const {return m_msKmvf;}

    void UserInfo::chisqKvf(double v) {m_chisqKvf = v;}
    double UserInfo::chisqKvf() const {return m_chisqKvf;}

    void UserInfo::chisqKvfNdf(double v) {m_chisqKvfNdf = v;}
    double UserInfo::chisqKvfNdf() const {return m_chisqKvfNdf;}

    void UserInfo::chisqKmvf(double v) {m_chisqKmvf = v;}
    double UserInfo::chisqKmvf() const {return m_chisqKmvf;}

    void UserInfo::chisqKmvfNdf(double v) {m_chisqKmvfNdf = v;}
    double UserInfo::chisqKmvfNdf() const {return m_chisqKmvfNdf;}

    void UserInfo::clKvf(double v) {m_clKvf = v;}
    double UserInfo::clKvf() const {return m_clKvf;}

    void UserInfo::clKmvf(double v) {m_clKmvf = v;}
    double UserInfo::clKmvf() const {return m_clKmvf;};

    void UserInfo::dist2IP(double v) {m_dist2IP = v;}
    double UserInfo::dist2IP() const {return m_dist2IP;}

    void UserInfo::dist2IPKmvf(double v) {m_dist2IPKmvf = v;}
    double UserInfo::dist2IPKmvf() const {return m_dist2IPKmvf;}

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

    void UserInfo::helicity(double v) {m_helicity = v;}
    double UserInfo::helicity() const {return m_helicity;}

    void UserInfo::probPid(double v) { m_pid = v;}
    double UserInfo:: probPid() const {return m_pid;}

    void UserInfo::vertexMode(std::string s)  {m_vertexMode = s;}
    std::string UserInfo::vertexMode() const {return m_vertexMode;}



    } // namespace Belle

