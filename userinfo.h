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
        double m_chisqKvf;           // chi2 after Kvf vertex fitting
        double m_chisqKvfNdf;        // chi2 ndf after Kvf vertex fitting
        double m_chisqKmvf;          // chi2 after Kmvf vertex fitting (mass-constraint fit)
        double m_chisqKmvfNdf;       // chi2 ndf after Kmvf vertex fitting (mass-constraint fit)
        double m_clKvf;              // confidence level after Kvf vertexing
        double m_clKmvf;             // confidence level after Kmvf vertexing
        double m_dist2IP;            // distance to IP after vertex
        double m_dist2IPKmvf;        // distance to IP after Kmvf vertex fitting (mass-constraint fit)
        double m_dist2Mother;        // distance to Mother
        bool   m_useTube;            // use tube in vertexing
        bool   m_useKmvf;            // use mass-constraint fit (Kmvf) to correct momentum
        bool   m_isAdoptCut;         // pass through cuts
        double m_wMass;              // mass window for particle pre-selection (soft cuts)
        double m_helicity;           // cos of helicity angle Theta for Dss children (which one exactly .. ??)
        double m_pid;                // probability of particle identification for Dss children, def value = -1.
        std::string m_vertexMode;    // vertexing mode
        Particle* particle;          // pointer to the source instance (Particle class object)

        void init();

    public:
        // Default constructor
        UserInfo();

        // Particle constructor
        explicit UserInfo(Particle&);

        // Copy constructor
        UserInfo(const UserInfo&);

        // Destructor
        virtual ~UserInfo();

        // constructs self object.
        UserInfo* clone() const;

        // Copy operator
        UserInfo& operator = (const UserInfo &);

        void msComb(double v);
        double msComb() const;

        void msKvf(double v);
        double msKvf() const;

        void msKmvf(double v);
        double msKmvf() const;

        void chisqKvf(double v);
        double chisqKvf() const;

        void chisqKvfNdf(double v);
        double chisqKvfNdf() const;

        void chisqKmvf(double v);
        double chisqKmvf() const;

        void chisqKmvfNdf(double v);
        double chisqKmvfNdf() const;

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

        void helicity(double v);
        double helicity() const;

        void probPid(double v);
        double probPid() const;

        void vertexMode(std::string s);
        std::string vertexMode() const;

        };

#if defined(BELLE_NAMESPACE)
    } // namespace Belle
#endif
