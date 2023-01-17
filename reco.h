#include <iostream>
#include <cstdio>
#include <ctime>
#include <vector>
#include <string>
#include <list>

#include "belle.h"

#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"

#ifdef HAVE_BOOST_FOREACH
#undef HAVE_BOOST_FOREACH
#endif

#include "panther/panther.h"
#include MDST_H
#include BELLETDF_H
#include HEPEVT_H
#include EVTCLS_H

#include "particle/utility.h"
#include "particle/combination.h"
#include "particle/ParticleUserInfo.h"
#include "kid/atc_pid.h"
#include "mdst/mdst.h"
#include "mdst/Muid_mdst.h"
#include "eid/eid.h"
#include "ip/IpProfile.h"
#include "benergy/BeamEnergy.h"
#include "kfitter/kvertexfitter.h"

#include "tagv/TagV.h"
#include "hamlet/Hamlet.h"


#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

    using namespace std;

    std::clock_t start;
    double duration;
    // Beam parameters
    static const double E_HER               = BeamEnergy::E_HER();
    static const double E_LER               = BeamEnergy::E_LER();
    static const double CROSS_ANGLE         = BeamEnergy::Cross_angle();

    // Parameters for track selection
    static const double
                             minProbPID_Kn  = 0.2,
                             minProbPID_Kp  = 0.6,
                             minProbProtPID = 0.0,
                             maxProbEl      = 1.0,
                             maxProbPion    = 0.9,
                             dRcut          = 0.5,
                             dZcut          = 3.0;
    // Parameters for candidate selection
    static const double
                             B_MASS         = 5.279,
                             KP_MASS        = 0.4937,
                             K0_MASS        = 0.4977,
                             DP_MASS        = 1.869,
                             D0_MASS        = 1.864,
                             JPSI_MASS      = 3.097,
                             PI_MASS        = 0.1396,
                             dM_phi         = 0.020,
                             dM_2317        = 0.270,
                             M_2317_min     = 2.25,
                             M_2317_max     = 2.40,
                             dM_Bs0         = 0.100;
    // Photon and pi0 selection
    static const double
                            eGammaMin       = 0.030, // 0.100;
                            wMassPi0GG      = 0.030, // 0.020;
                            minPi0GammP     = 0.030,
                            minPi0pStar     = 0.040;

    static const double
                           dM_Dgr           = 0.150,  // GeV, mass window for D ground tag;
                           dM_V0            = 0.150,  // GeV, peak mass window for V0 tag;
                           dM_Ksr0          = 0.150,  // GeV, peak mass window for K*0 tag;
                           dM_Rho           = 0.175,  // 0.275; // 0.130; // GeV, peak mass window for Rho tag;
                           dM_Dss           = 0.150,  // GeV, mass window for DS+;
                           dM_Dsst          = 0.150,  // GeV, mass window for D*S+;
                           wDst             = 0.050,  // 0.008; // wWideMassDs/2.;
                           wB               = 0.150,  // 0.120; // wWideMassB/2.;
                           maxR2            = 0.5,
                           maxChisq         = 200.;


    static bool            McFlag           = false;  // determined w/ getEventInfo
    
    static const bool      stDumpBs0        = true;
    static const bool      stDump2317       = false;
    static const bool      stDumpDss        = false;
    // Magnetic field
    static const double    BF               = 1.5;
    bool                   useBF            = true;
    bool                   useVTX           = true;
    bool                   debugHel         = false;




    static const int nTrk = 4;
    std::vector< vector<Particle> > trkV(nTrk);
    std::string trkTit[nTrk] = {"pi+", "pi-", "K+", "K-"};
    
    std::vector<Particle> gammaV, pi0, phi0, Ksr0, Ksr0bar, Dss_p, Dss_m, 
            Dss_p_2317, Dss_m_2317, Bs0, Bs0bar; // BsStar0, BsStar0bar, Upsilon_5S;

    class Reco : public Module {
    public:
        Reco(void);

        ~Reco(void) {};

        void init(int *) {
            std::cout << "---- Reco initialization -----" << std::endl;
            start = std::clock();
            /* initialization of Ptype is done inside Hamlet */
            // Hamlet::init();
            
        };

        void term(void) {};
        void disp_stat(const char*) {};
        void hist_def(void);
        void event(BelleEvent*, int*);
        void begin_run(BelleEvent*, int*);
        void end_run(BelleEvent*, int*) {};
        void other(int*, BelleEvent*, int*) {};

        BelleTuple* TP_Bs0;
        // BelleTuple* TP_Dss, *TP_Dss_2317, *TP_Bs0;
        // BelleTuple *TP_phi0, *TP_Ksr0, *TP_Dss, *TP_Dss_2317, *TP_Bs0;
    };

    extern "C" Module_descr* mdcl_Reco() { /* main */
        Reco* module = new Reco;
        Module_descr* dscr = new Module_descr("Reco", module);
        IpProfile::define_global(dscr);
        BeamEnergy::define_global(dscr);
        return dscr;
    }

    Reco::Reco(void) {
        std::cout << "---- Reco constructor -----" << std::endl;
    };


    void Reco::begin_run(BelleEvent*, int*) {

        std::cout << std::endl << "---- Reco's begin_run function -----" << std::endl;
        eid::init_data();
        // Get IP profile data from $BELLE_POSTGRES_SERVER
        duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

        // Hamlet::begin_run();
        IpProfile::begin_run();
        // Dump IP profile data to stdout (optional)
        IpProfile::dump();
        
        BeamEnergy::begin_run();
        
        std::cout << "duration: "     << duration << " sec"        << std::endl;
        std::cout << " Benergy(): "   << Benergy()                 << std::endl;
        std::cout << " E_beam_corr: " << BeamEnergy::E_beam_corr() << std::endl;
        std::cout << " E_beam_err: "  << BeamEnergy::E_beam_err()  << std::endl;
        std::cout << " E_LER: "       << BeamEnergy::E_LER()       << std::endl;
        std::cout << " E_HER: "       << BeamEnergy::E_HER()       << std::endl;
        std::cout << " E_beam_orig: " << BeamEnergy::E_beam_orig() << std::endl;
        std::cout << " E_LER_orig: "  << BeamEnergy::E_LER_orig()  << std::endl;
        std::cout << " E_HER_orig: "  << BeamEnergy::E_HER_orig()  << std::endl;
        std::cout << " E_beam2: "     << BeamEnergy::E_beam2()     << std::endl;
        std::cout << " Cross_angle: " << BeamEnergy::Cross_angle() << std::endl;
    }

    void Reco::hist_def(void) {

        extern BelleTupleManager* BASF_Histogram;

        std::cout << "----- Reco's hist_def function -----" << std::endl;

        std::string s_info        = " expn runn evtn r2 ipx ipy ipz evtgen ";
        std::string s_B0s         = " gen_bs chg_bs mbc_bs de_bs de2_bs msV_bs msM_bs msC_bs chiV_bs chiM_bs prbV_bs"
                                    " prbM_bs cl_bs clV_bs clM_bs px_bs py_bs pz_bs p_bs pVx_bs pVy_bs pVz_bs xx_bs xy_bs xz_bs"
                                    " dVx_bs dVy_bs dVz_bs pt_bs ph_bs th_bs eta_bs cth_bs ";
        std::string s_Dss         = " gen_ds chg_ds ms_ds chi_ds pt_ds ph_ds th_ds ";
        std::string s_Dss_1       = " gen_ds1 chg_ds1 msV_ds1 msM_ds1 msC_ds1 chiV_ds1 chiM_ds1 prbV_ds1 prbM_ds1 cl_ds1"
                                    " clV_ds1 clM_ds1 pVx_ds1 pVy_ds1 pVz_ds1 xx_ds1 xy_ds1 xz_ds1 dVx_ds1 dVy_ds1 dVz_ds1"
                                    " psr_ds1 px_ds1 py_ds1 pz_ds1 p_ds1 pt_ds1 ph_ds1 th_ds1 id0_ds1 id1_ds1 ";
        std::string s_Dss_2       = " gen_ds2 chg_ds2 msV_ds2 msM_ds2 msC_ds2 chiV_ds2 chiM_ds2 prbV_ds2 prbM_ds2 cl_ds2 clV_ds2"
                                    " clM_ds2 pVx_ds2 pVy_ds2 pVz_ds2 xx_ds2 xy_ds2 xz_ds2 dVx_ds2 dVy_ds2 dVz_ds2 psr_ds2"
                                    " px_ds2 py_ds2 pz_ds2 p_ds2 pt_ds2 ph_ds2 th_ds2 id0_ds2 id1_ds2 ";
        std::string s_Dss_gen     = " ms_ds_t px_ds_t py_ds_t pz_ds_t e_ds_t th_ds_t ct_ds_t ph_ds_t rh_ds_t ";
        std::string s_Dss_gen_1   = " ms_ds_t1 px_ds_t1 py_ds_t1 pz_ds_t1 e_ds_t1 th_ds_t1 ct_ds_t1 ph_ds_t1 rh_ds_t1 ";
        std::string s_Dss_gen_2   = " ms_ds_t2 px_ds_t2 py_ds_t2 pz_ds_t2 e_ds_t2 th_ds_t2 ct_ds_t2 ph_ds_t2 rh_ds_t2 ";

        std::string s_child       = " gen_ch chg_ch ms_ch chi_ch pt_ch hel_ch ind_ch "; // ind_ch=1 - phipi, 2 - KsrK
        // Add ind_ch
        std::string s_child_1     = " gen_ch1 chg_ch1 msV_ch1 msM_ch1 msC_ch1 chiV_ch1 chiM_ch1 prbV_ch1 prbM_ch1 cl_ch1"
                                    " clV_ch1 clM_ch1 pVx_ch1 pVy_ch1 pVz_ch1 xx_ch1 xy_ch1 xz_ch1 dVx_ch1 dVy_ch1 dVz_ch1"
                                    " psr_ch1 px_ch1 py_ch1 pz_ch1 p_ch1 pt_ch1 ph_ch1 th_ch1 ind_ch1 hel_ch1 ";
        std::string s_child_2     = " gen_ch2 chg_ch2 msV_ch2 msM_ch2 msC_ch2 chiV_ch2 chiM_ch2 prbV_ch2 prbM_ch2 cl_ch2"
                                    " clV_ch2 clM_ch2 pVx_ch2 pVy_ch2 pVz_ch2 xx_ch2 xy_ch2 xz_ch2 dVx_ch2 dVy_ch2 dVz_ch2"
                                    " psr_ch2 px_ch2 py_ch2 pz_ch2 p_ch2 pt_ch2 ph_ch2 th_ch2 ind_ch2 hel_ch2 ";

        std::string s_2317        = " gen_d17 chg_d17 msV_d17 msM_d17 msC_d17 chiV_d17 chiM_d17 prbV_d17 prbM_d17"
                                    " cl_d17 clV_d17 clM_d17 pVx_d17 pVy_d17 pVz_d17 xx_d17 xy_d17 xz_d17 dVx_d17 dVy_d17"
                                    " dVz_d17 psr_d17 px_d17 py_d17 pz_d17 p_d17 pt_d17 ph_d17 th_d17 hel_d17 ";
        std::string s_2317_gen    = " ms_d17_t px_d17_t py_d17_t pz_d17_t e_d17_t th_d17_t ct_d17_t ph_d17_t rh_d17_t ";

        std::string s_pi_d17      = " gen_p0_d eg1_p0_d eg2_p0_d psr_p0_d mgg_p0_d gg1_p0_d gg2_p0_d ";
        std::string s_pi_B0s      = " gen_p0_b eg1_p0_b eg2_p0_b psr_p0_b mgg_p0_b gg1_p0_b gg2_p0_b ";

        std::string s_2317_wPi    = s_2317  + s_2317_gen  + s_pi_d17;
        std::string s_Dss1_wCh    = s_Dss_1 + s_Dss_gen_1 + s_child_1;
        std::string s_Dss2_wCh    = s_Dss_2 + s_Dss_gen_2 + s_child_2;

        std::string s_Dss_sum     = s_info  + s_Dss + s_Dss_gen  + s_child;
        std::string s_2317_sum    = s_Dss_sum       + s_2317_wPi;
        std::string s_B0s_sum     = s_info  + s_B0s + s_Dss1_wCh + s_2317_wPi + s_Dss2_wCh + s_pi_B0s;

        // TP_Dss               = BASF_Histogram->ntuple("Dss",    s_Dss_sum);
        // TP_Dss_2317          = BASF_Histogram->ntuple("Ds2317", s_2317_sum);
        TP_Bs0               = BASF_Histogram->ntuple("Bs0", s_B0s_sum);
    }

#if defined(BELLE_NAMESPACE)
}//namespace Belle
#endif

