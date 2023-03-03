#include <iostream>
#include <cstdio>
#include <ctime>
#include <vector>
#include <string>


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
// necessary for dealing with systematic errors on data,
// idle and shoud be commented out for MC simulations
// modity depanding on you analysis requirements
// look here http://belle.kek.jp/group/pid_joint/
//# include "kid_eff_o6.h"


#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

    using namespace std;

    std::clock_t start;
    double duration;

    // Parameters for track selection
    static const double
                             minProbPID_Kn  = 0.2,
                             minProbPID_Kp  = 0.6,
                             maxProbPion    = 0.9,
                             dRcut          = 0.5,
                             dZcut          = 3.0;
    // Parameters for candidate selection
    static const double

                             dM_2317        = 0.270,
                             dM_Bs0         = 0.250;
    // Photon and pi0 selection
    static const double
                            eGammaMin       = 0.030, // 0.100;
                            wMassPi0GG      = 0.030, // 0.020;
                            minPi0GammP     = 0.030;

    static const double
                           dM_Dgr           = 0.150,  // GeV, mass window for D ground tag;
                           dM_V0            = 0.150,  // GeV, peak mass window for V0 tag;
                           dM_Ksr0          = 0.150,  // GeV, peak mass window for K*0 tag;
                           dM_Dss           = 0.150;  // GeV, mass window for DS+;



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


    // :::::::::::::::::: Main class Reco :::::::::::::::
    // defines the behavior of the entire Reco module
    class Reco : public Module {
    public:
        Reco();

        ~Reco() = default;

        void init(int*);
        void term();
        void disp_stat(const char*);
        void hist_def();
        void event(BelleEvent*, int*);
        void begin_run(BelleEvent*, int*);
        void end_run(BelleEvent*, int*);
        void other(int*, BelleEvent*, int*);

        BelleTuple* TP_Bs0;
        // BelleTuple* TP_Dss, *TP_Dss_2317, *TP_Bs0;
        // BelleTuple *TP_phi0, *TP_Ksr0, *TP_Dss, *TP_Dss_2317, *TP_Bs0;

        /* dealing with systematic errors in the dataset,
        * not applicable to a MC simulation;
        * uncomment the chunk of code below and
        * modify considering your analysis requirements
        *
        private:
            KID_eff_06 kideff_k_1, kideff_k_2, kideff_k_3;
            KID_eff_06 kideff_pi_1, kideff_pi_2, kideff_pi_3;
        */
    };

extern "C" Module_descr* mdcl_Reco() {  /* main */
        Reco* module = new Reco;
        Module_descr* dscr = new Module_descr("Reco", module);
        IpProfile::define_global(dscr);
        BeamEnergy::define_global(dscr);
        return dscr;
    }

    Reco::Reco() {
        std::cout << "---- Reco constructor -----" << std::endl;
    };

    void Reco::init(int *) {
        std::cout << "---- Reco initialization -----" << std::endl;
        start = std::clock();
        /* initialization of Ptype is done inside Hamlet */
        // Hamlet::init();
        /* dealing with systematic errors in the dataset,
         * not applicable to a MC simulation;
         * uncomment the chunk of code below and
         * modify considering your analysis requirements
         *
         * kideff_k_1.init (0.9, 1, "kaon_1", "data/kideff-2006-svd1-all.dat");
         * kideff_k_2.init (0.9, 1, "kaon_2", "data/kideff-2006-svd2-all.dat");
         * kideff_k_3.init (0.9, 1, "kaon_3", "data/kideff-2009-newtrk-all.dat");
         * kideff_pi_1.init(0.4, 3, "pion_1", "data/kideff-2006-svd1-all.dat");
         * kideff_pi_2.init(0.4, 3, "pion_2", "data/kideff-2006-svd2-all.dat");
         * kideff_pi_3.init(0.4, 3, "pion_3", "data/kideff-2009-newtrk-all.dat");
         */
    };

    void Reco::term() {
        /*
         * dealing with systematic erors,
         * uncomment the chunk of code
         * below and modify considering
         * your analysis requirements
         *
         * kideff_k_1.calculate();
         * kideff_k_2.calculate();
         * kideff_k_3.calculate();
         * kideff_pi_1.calculate();
         * kideff_pi_2.calculate();
         * kideff_pi_3.calculate();
         * kideff_k_1.dump();
         * kideff_k_2.dump();
         * kideff_k_3.dump();
         * kideff_pi_1.dump();
         * kideff_pi_2.dump();
         * kideff_pi_3.dump();
         *
        */
    }

    void Reco::disp_stat(const char*) {};

    void Reco::hist_def() {

        extern BelleTupleManager* BASF_Histogram;

        std::cout << "----- Reco's hist_def function -----" << std::endl;

        std::string s_info        = " expn runn evtn r2 ipx ipy ipz evtgen ";

        std::string s_B0s         = " chg_bs gen_bs ms_bs msm_bs chi_bs chim_bs chn_bs chmn_bs"
                                    " cl_bs clm_bs dvx_bs dvy_bs dvz_bs"
                                    " px_bs py_bs pz_bs p_bs pt_bs"
                                    " ph_bs th_bs mbc_bs de_bs de2_bs";

        std::string s_Bs0_gen     = " ms_bs_t px_bs_t py_bs_t pz_bs_t e_bs_t th_bs_t ct_bs_t ph_bs_t rh_bs_t ";

        std::string s_Dss         = " gen_ds chg_ds ms_ds chi_ds pt_ds ph_ds th_ds ";

        std::string s_Dss_1       = " chg_ds1 gen_ds1 c00_ds1 c01_ds1 c10_ds1"
                                    " ms_ds1 msm_ds1 chi_ds1 chim_ds1 chn_ds1 chmn_ds1"
                                    " cl_ds1 clm_ds1 dvx_ds1 dvy_ds1 dvz_ds1"
                                    " psr_ds1 px_ds1 py_ds1 pz_ds1 p_ds1 pt_ds1 ph_ds1 th_ds1 p00_ds1 p01_ds1 p10_ds1 ";

        std::string s_Dss_2       = " chg_ds2 gen_ds2 c00_ds2 c01_ds2 c10_ds2"
                                    " ms_ds2 msm_ds2 chi_ds2 chim_ds2 chn_ds2 chmn_ds2"
                                    " cl_ds2 clm_ds2 dvx_ds2 dvy_ds2 dvz_ds2"
                                    " psr_ds2 px_ds2 py_ds2 pz_ds2 p_ds2 pt_ds2 ph_ds2 th_ds2 p00_ds2 p01_ds2 p10_ds2 ";

        std::string s_Dss_gen     = " ms_ds_t px_ds_t py_ds_t pz_ds_t e_ds_t th_ds_t ct_ds_t ph_ds_t rh_ds_t ";

        std::string s_Dss_gen_1   = " ms_ds_t1 px_ds_t1 py_ds_t1 pz_ds_t1 e_ds_t1 th_ds_t1 ct_ds_t1 ph_ds_t1 rh_ds_t1 ";

        std::string s_Dss_gen_2   = " ms_ds_t2 px_ds_t2 py_ds_t2 pz_ds_t2 e_ds_t2 th_ds_t2 ct_ds_t2 ph_ds_t2 rh_ds_t2 ";

        std::string s_child       = " gen_ch chg_ch ms_ch chi_ch pt_ch hel_ch ind_ch "; // ind_ch=1 - phipi, 2 - KsrK

        std::string s_child_1     = " chg_ch1 gen_ch1 ind_ch1 ms_ch1 msm_ch1 chi_ch1 chim_ch1 chn_ch1 chmn_ch1"
                                    " cl_ch1 clm_ch1 dvx_ch1 dvy_ch1 dvz_ch1"
                                    " pt_ch1 ph_ch1 th_ch1 psr_ch1 px_ch1 py_ch1 pz_ch1 ";

        std::string s_child_2     = " chg_ch2 gen_ch2 ind_ch2 ms_ch2 msm_ch2 chi_ch2 chim_ch2 chn_ch2 chmn_ch2"
                                    " cl_ch2 clm_ch2 dvx_ch2 dvy_ch2 dvz_ch2"
                                    " pt_ch2 ph_ch2 th_ch2 psr_ch2 px_ch2 py_ch2 pz_ch2 ";

        std::string s_2317        = " chg_d17 gen_d17 ms_d17 msm_d17 chi_d17 chim_d17 chn_d17 chmn_d17"
                                    " cl_d17 clm_d17 dVx_d17 dVy_d17 dVz_d17 psr_d17"
                                    " px_d17 py_d17 pz_d17 p_d17 pt_d17 ph_d17 th_d17 hel_d17 ";

        std::string s_2317_gen    = " ms_d17_t px_d17_t py_d17_t pz_d17_t e_d17_t th_d17_t ct_d17_t ph_d17_t rh_d17_t ";

        std::string s_pi_d17      = " gen_p0_d gg1_p0_d gg2_p0_d eg1_p0_d eg2_p0_d psr_p0_d p_p0_d mgg_p0_d ms_p0_d msm_p0_d"
                                    " chi_p0_d chm_p0_d chn_p0_d cmn_p0_d cl_p0_d clm_p0_d ";
        std::string s_pi_B0s      = " gen_p0_b gg1_p0_b gg2_p0_b eg1_p0_b eg2_p0_b psr_p0_b p_p0_b mgg_p0_b ms_p0_b msm_p0_b"
                                    " chi_p0_b chm_p0_b chn_p0_b cmn_p0_b cl_p0_b clm_p0_b ";


        std::string s_2317_wPi    = s_2317  + s_2317_gen  + s_pi_d17;
        std::string s_Dss1_wCh    = s_Dss_1 + s_Dss_gen_1 + s_child_1;
        std::string s_Dss2_wCh    = s_Dss_2 + s_Dss_gen_2 + s_child_2;

        std::string s_Dss_sum     = s_info  + s_Dss + s_Dss_gen  + s_child;
        std::string s_2317_sum    = s_Dss_sum       + s_2317_wPi;
        std::string s_B0s_sum     = s_info  + s_B0s + s_Bs0_gen  + s_Dss1_wCh + s_2317_wPi + s_Dss2_wCh + s_pi_B0s;

        // TP_Dss               = BASF_Histogram->ntuple("Dss",    s_Dss_sum);
        // TP_Dss_2317          = BASF_Histogram->ntuple("Ds2317", s_2317_sum);
        TP_Bs0               = BASF_Histogram->ntuple("Bs0", s_B0s_sum);
    }

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
        
        std::cout << "duration: "        << duration << " sec"        << std::endl;
        std::cout << " Run_flag(): "     << BeamEnergy::Run_flag()    << std::endl;
        std::cout << " Benergy(): "      << Benergy()                 << std::endl;
        std::cout << " E_beam_corr: "    << BeamEnergy::E_beam_corr() << std::endl;
        std::cout << " E_beam_err: "     << BeamEnergy::E_beam_err()  << std::endl;
        std::cout << " E_LER: "          << BeamEnergy::E_LER()       << std::endl;
        std::cout << " E_HER: "          << BeamEnergy::E_HER()       << std::endl;
        std::cout << " E_beam_orig: "    << BeamEnergy::E_beam_orig() << std::endl;
        std::cout << " E_LER_orig: "     << BeamEnergy::E_LER_orig()  << std::endl;
        std::cout << " E_HER_orig: "     << BeamEnergy::E_HER_orig()  << std::endl;
        std::cout << " E_beam2: "        << BeamEnergy::E_beam2()     << std::endl;
        std::cout << " Cross_angle: "    << BeamEnergy::Cross_angle() << std::endl;
        std::cout << " Ecm(), sqrt(s): " << BeamEnergy::Ecm()         << std::endl;

        std::cout << "BeamEnergy::dump(): " << std::endl;
        BeamEnergy::dump();
    }




    void Reco::end_run(BelleEvent*, int*) {};
    void Reco::other(int*, BelleEvent*, int*) {}





#if defined(BELLE_NAMESPACE)
}//namespace Belle
#endif

