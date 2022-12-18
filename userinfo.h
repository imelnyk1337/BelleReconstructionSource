#include "belle.h"
#include "particle/Particle.h"
#include "particle/ParticleUserInfo.h"
#include <iostream>
#include <string>

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

using namespace std;

class UserInfo : public ParticleUserInfo {
private:
    double m_msComb; // mass after combination
    double m_msKvf;  // mass after vertex fitter
    double m_msKmv;
    double m_chisq;     // chi2 after last vertexing
    double m_chisqKvf;  // chi2 after Kvf vertexing
    double m_cl;        // cl after last vertexing
    double m_clKvf;     // cl after Kvf vertexing
    double m_dist2IP;       // distance to IP after vertex
    double m_dist2IPmvf;    // distance to IP after mass-vertex fitter
    double m_dist2Mother;   // distance to Mother
    bool m_useTube;  // use tube in vertexing
    bool m_useKmvf;  // use mass-vertex fitter to correct momentum
    bool m_isAdoptCut; // pass through cuts
    double m_wMass;
    double m_maxChi2;
    double m_helicity;
    double m_pid;   // probability of particle identification for Dss children, def value = -1
    Particle* particle;

    void init();

public:
    // Default constructor
    UserInfo();

    // Particle constructor
    UserInfo(Particle&);

    // Copy constructor
    UserInfo(const UserInfo&);

    // Destructor
    virtual ~UserInfo();

    // constructs self object.
    UserInfo* clone(void) const;

    // Copy operator
    UserInfo& operator = (const UserInfo &);

//    void SetParticlePointer(Particle *p) { particle = p; }

    void msComb(double v);
    double msComb() const;

    void msKvf(double v);
    double msKvf() const;

    void msKmv(double v);
    double msKmv() const;

    void chisq(double v);
    double chisq() const;

    void chisqKvf(double v);
    double chisqKvf() const;

    void cl(double v);
    double cl() const;

    void clKvf(double v);
    double clKvf() const;

    void dist2IP(double v);
    double dist2IP() const;

    void dist2IPmvf(double v);
    double dist2IPmvf() const;

    void dist2Mother(double v);
    double dist2Mother() const;
  
    void useTube(bool v);
    bool useTube() const;
  
    void useKmvf(bool v);
    bool useKmvf() const;
  
    void isAdoptCut(bool v);
    bool isAdoptCut() const;

    void wMass(double v);
    double wMass() const;

    void maxChi2(double v);
    double maxChi2() const;

    void helicity(double v);
    double helicity() const;

    void probpid(double v);
    double probpid() const;
  
};
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
