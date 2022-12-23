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
    double m_msComb;             // invariant mass after combination
    double m_msKvf;              // invariant mass after Kvf vertex fitting
    double m_msKmvf;             // invariant mass after Kmvf vertex fitting (mass-constraint fit)
    double m_chisq;              // chi2 after the last vertexing (vertex fitting)
    double m_chisqKvf;           // chi2 after Kvf vertex fitting
    double m_chisqProbKvf;       // Prob_chi2 after Kvf vertex fitting
    double m_chisqKmvf;          // chi2 after Kmvf vertex fitting (mass-constraint fit)
    double m_chisqProbKmvf;      // Prob_chi2 after Kmvf vertex fitting (mass-constraint fit)
    double m_cl;                 // cl after the last vertexing
    double m_clKvf;              // cl after Kvf vertexing
    double m_clKmvf;             // cl after Kmvf vertexing
    double m_dist2IP;            // distance to IP after vertex
    double m_dist2IPKmvf;        // distance to IP after Kmvf vertex fitting (mass-constraint fit)
    double m_dist2Mother;        // distance to Mother
    bool   m_useTube;            // use tube in vertexing
    bool   m_useKmvf;            // use mass-constraint fit (Kmvf) to correct momentum
    bool   m_isAdoptCut;         // pass through cuts
    double m_wMass;              // mass window for particle pre-selection (soft cuts)
    double m_maxChi2;            // maximum chi2 of vertexing
    double m_helicity;           // cos of helicity angle Theta for Dss children (which one exactly .. ??)
    double m_pid;                // probability of particle identification for Dss children, def value = -1.
    Particle* particle;          // pointer to the source instance (Particle class object)

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

    // void SetParticlePointer(Particle* p) {particle = p;}

    void msComb(double v);
    double msComb() const;

    void msKvf(double v);
    double msKvf() const;

    void msKmvf(double v);
    double msKmvf() const;

    void chisq(double v);
    double chisq() const;

    void chisqKvf(double v);
    double chisqKvf() const;

    void probChi2Kvf(double v);
    double probChi2Kvf() const;

    void chisqKmvf(double v);
    double chisqKmvf() const;

    void probChi2Kmvf(double v);
    double probChi2Kmvf() const;

    void cl(double v);
    double cl() const;

    void clKvf(double v);
    double clKvf() const;

    void clKmvf(double v);
    double clKmvf() const;

    void dist2IP(double v);
    double dist2IP() const;

    void dist2IPKmvf(double v);
    double dist2IPKmvf() const;

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
