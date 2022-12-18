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
        

    const float B_MASS    = 5.279;
    const float KP_MASS   = 0.4937;
    const float K0_MASS   = 0.4977;
    const float DP_MASS   = 1.869;
    const float D0_MASS   = 1.864;
    const float JPSI_MASS = 3.097;
    const float PI_MASS   = 0.1396;

//     double dM_Ksr = 0.100;
    double dM_phi = 0.020;
//     double dM_Dss  = 0.050;
    double dM_2317 = 0.270;
    double M_2317_min = 2.25, M_2317_max = 2.40;
    double dM_Bs0  = 0.080;
    
    double eGammaMin   = 0.030;
    double wMassPi0GG  = 0.020; // 0.020;
    double minPi0GammP = 0.030;
    double minPi0pStar = 0.040;

    double dM_Dgr  = 0.100;  // GeV, mass window for D ground tag
    double dM_V0   = 0.060;  // GeV, peak mass window for V0 tag
    double dM_Ksr0 = 0.120;  // GeV, peak mass window for K*0 tag
    double dM_Rho  = 0.175; // 0.275; // 0.130;  // GeV, peak mass window for Rho tag
    double dM_Dss  = 0.100;  // GeV, mass window for DS+ 
    double dM_Dsst = 0.100;  // GeV, mass window for D*S+
    double wDst  = 0.020; //0.008; //wWideMassDs/2.;
    double wB    = 0.110;//0.120;//wWideMassB/2.;

    double maxR2 = 0.5;
    double maxChisq = 200.;


    bool McFlag = false;  // determined w/ getEventInfo
    
    bool stDumpBs0  = true;
    bool stDump2317 = false;
    bool stDumpDss  = false;
  
    const double BF = 1.5;
    bool useBF = true;
    bool useVTX = true;
    bool debugHel = false;

    const int nTrk = 4;
    std::vector< vector<Particle> > trkV(nTrk);
    string trkTit[nTrk] = {"pi+", "pi-", "K+", "K-"};
    
    std::vector<Particle> gammaV, pi0, phi0, Ksr0, Ksr0bar, Dss_p, Dss_m, 
            Dss_p_2317, Dss_m_2317, Bs0, Bs0bar, BsStar0, BsStar0bar, Upsilon_5S;

    class Reco : public Module {
    public:
        Reco(void);

        ~Reco(void) {};

        void init(int *) {
            std::cout << "---- Reco initialization -----" << std::endl << std::endl;
            start = std::clock();
            /* initialization of Ptype is done inside Hamlet */
    //         Hamlet::init();
            
        };

        void term(void) {};
        void disp_stat(const char*) {};
        void hist_def(void);
        void event(BelleEvent*, int*);
        void begin_run(BelleEvent*, int*);
        void end_run(BelleEvent*, int*) {};
        void other(int*, BelleEvent*, int*) {};

        BelleTuple *TP_phi0, *TP_Ksr0, *TP_Dss, *TP_Dss_2317, *TP_Bs0;
    };

    extern "C" Module_descr* mdcl_Reco() { /* main */
        Reco* module = new Reco;
        Module_descr* dscr = new Module_descr("Reco", module);
        IpProfile::define_global(dscr);
        BeamEnergy::define_global(dscr);
        return dscr;
    }

    Reco::Reco ( void ) {	
        std::cout << "---- Reco constructor -----" << std::endl << std::endl; 
    };


    void Reco::begin_run(BelleEvent*, int *) {
        std::cout << std::endl << "---- Reco's begin_run function -----" << std::endl;
        eid::init_data();
        //   Get IP profile data from $BELLE_POSTGRES_SERVER
        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    //     Hamlet::begin_run();
        IpProfile::begin_run();
        //  Dump IP profile data to STDOUT (optional) 
        IpProfile::dump();
        
        BeamEnergy::begin_run();
        
        std::cout << "duration: "<< duration << " sec" << std::endl << std::endl;    
        std::cout << " Benergy(): " << Benergy() << std::endl;
        
        std::cout << " E_beam_corr: " << BeamEnergy::E_beam_corr() << std::endl;
        std::cout << " E_beam_err: " << BeamEnergy::E_beam_err() << std::endl;
        std::cout << " E_LER: " << BeamEnergy::E_LER() << std::endl;
        std::cout << " E_HER: " << BeamEnergy::E_HER() << std::endl;
        std::cout << " E_beam_orig: " << BeamEnergy::E_beam_orig() << std::endl;
        std::cout << " E_LER_orig: " << BeamEnergy::E_LER_orig() << std::endl;
        std::cout << " E_HER_orig: " << BeamEnergy::E_HER_orig() << std::endl;
        std::cout << " E_beam2: " << BeamEnergy::E_beam2() << std::endl;
        std::cout << " Cross_angle: " << BeamEnergy::Cross_angle() << std::endl;
    }

    void Reco::hist_def(void) {

        extern BelleTupleManager* BASF_Histogram;

        std::cout << "---- Reco's hist_def function -----" << std::endl; 

        string s_info        = " expn runn evtn r2 ipx ipy ipz evtgen "; 

        string s_B0s         = " gen_bs chg_bs mbc_bs de_bs  ms_bs chi_bs vx_bs vy_bs vz_bs pt_bs ph_bs th_bs ";
        string s_Dss         = " gen_ds chg_ds  ms_ds chi_ds pt_ds ph_ds th_ds "; 
        string s_Dss_1       = " gen_ds1 chg_ds1  ms_ds1 chi_ds1 pt_ds1 ph_ds1 th_ds1 "; 
        string s_Dss_2       = " gen_ds2 chg_ds2  ms_ds2 chi_ds2 pt_ds2 ph_ds2 th_ds2 ";
        string s_Dss_gen     = " ms_ds_t  px_ds_t  py_ds_t  pz_ds_t  e_ds_t ";
        string s_Dss_gen_1   = " ms_ds_t1 px_ds_t1 py_ds_t1 pz_ds_t1 e_ds_t1 ";
        string s_Dss_gen_2   = " ms_ds_t2 px_ds_t2 py_ds_t2 pz_ds_t2 e_ds_t2 ";

        string s_phi0        = " gen_ph0 chg_ph0 ms_ph0 chi_ph0 pt_ph0 hel_ph0 ihp_ph0 "; 
        string s_Ksr0        = " gen_ks0 chg_ks0 ms_ks0 chi_ks0 pt_ks0 hel_ks0 ihp_ks0 "; 
        string s_child       = " gen_ch  chg_ch  ms_ch  chi_ch  pt_ch  hel_ch  ind_ch  "; // ind_ch=1 - phipi, 2 - KsrK
        string s_child_1     = " gen_ch1 chg_ch1 ms_ch1 chi_ch1 pt_ch1 hel_ch1 ind_ch1 "; 
        string s_child_2     = " gen_ch2 chg_ch2 ms_ch2 chi_ch2 pt_ch2 hel_ch2 ind_ch2 "; 

        string s_2317        = " gen_d17 chg_d17  ms_d17 chi_d17 pt_d17 psr_d17 ph_d17 th_d17 hel_d17 px_d17 py_d17 pz_d17 vx_d17 vy_d17 vz_d17 ";
        string s_2317_gen    = " ms_d17_t px_d17_t py_d17_t pz_d17_t e_d17_t ";

        string s_pi_d17      = " gen_p0_d eg1_p0_d  eg2_p0_d  psr_p0_d  mgg_p0_d gg1_p0_d gg2_p0_d ";
        string s_pi_B0s      = " gen_p0_b eg1_p0_b  eg2_p0_b  psr_p0_b  mgg_p0_b gg1_p0_b gg2_p0_b ";
        
        string s_2317_wPi    = s_2317  + s_2317_gen  + s_pi_d17;
        string s_Dss1_wCh    = s_Dss_1 + s_Dss_gen_1 + s_child_1;
        string s_Dss2_wCh    = s_Dss_2 + s_Dss_gen_2 + s_child_2;

        string s_Dss_sum     = s_info  + s_Dss + s_Dss_gen  + s_child;
        string s_2317_sum    = s_Dss_sum       + s_2317_wPi;
        string s_B0s_sum     = s_info  + s_B0s + s_Dss1_wCh + s_2317_wPi + s_Dss2_wCh + s_pi_B0s;

        TP_Dss               = BASF_Histogram->ntuple("dss",    s_Dss_sum);
        TP_Dss_2317          = BASF_Histogram->ntuple("ds2317", s_2317_sum);
        TP_Bs0               = BASF_Histogram->ntuple("Bs0",    s_B0s_sum);
    }

#if defined(BELLE_NAMESPACE)
}//namespace Belle
#endif

