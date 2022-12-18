//
// $Id: frec_ana.cc 9944 2006-11-29 07:36:07Z katayama $
//
// $Log$
// Revision 1.10  2004/11/22 05:24:35  matumot
// updates for new kid_statistics
// module parameter, "LogFileName" is introduced
//
// Revision 1.9  2004/07/29 06:53:10  matumot
// minor modification on counting method
//
// Revision 1.8  2004/05/28 10:05:54  matumot
// adjusted pi0 mass cut to 0.1178 < Mpi0 < 0.1502 GeV/c^2
// ( consistent with pi0 mass cuts in current fix_mdst module )
//
// Revision 1.7  2004/05/25 13:58:15  matumot
// add offset in ID of kid_statistics to avoid
// duplication for using frec_skim0 and frec_ana
// modules at the same time.
//
// Revision 1.6  2004/05/04 21:22:16  matumot
// fixed for CC complier
// - namespace was changed : fullrecon -> fullreconNS
//   to avoid overlap for fullrecon struct used by fullrecon.tdf
// - use new operator for dynamic allocation of array
// - add std:: for sort function
//
// Revision 1.4  2004/04/19 09:27:31  matumot
// revival version
//
//
// ==========================================================
//  File Name : frec_ana.cc
// ----------------------------------------------------------
//  Creation    ; 2004.04.07
//  Description ; Implimentation for the frec_ana class
//  Author      ; T.Matsumoto (matumot@bmail.kek.jp)
// ----------------------------------------------------------

#include <fstream>

#include "belle.h"
#include "fullrecon/frec_ana.h"
#include "fullrecon/frec_util.h"

#include "basf/module_descr.h"
#include "particle/Particle.h"
#include "kid/atc_pid.h"
#include "kid/kid_statistics.h"
#include "ip/IpProfile.h"
#include "benergy/BeamEnergy.h"
#include "mdst/mdst.h"

#include BELLETDF_H
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


using namespace fullreconNS;

// ----------------
// BASF interface
// ----------------
extern "C" Module_descr  *mdcl_frec_ana(void){

  frec_ana  *m = new frec_ana;
  Module_descr *d = new Module_descr("frec_ana", m);

  BeamEnergy::define_global(d);

  d->define_param( "summaryFile", "summary file", 256, m->m_LogFileName );

  d->define_param( "LogFileName", "log file name", 256, m->m_LogFileName );

  d->define_param( "Monitor", "Monitor (1=on, 0=off)" , &m->m_Monitor );

  d->define_param( "flagHist", "Hist output (1=on, 0=off)" , &m->m_flagHist );

  d->define_param( "MC_flag", "MC flag (0=off, 1=on(for old), 2=on(for new))", &m->mm_mc );

  d->define_param( "SetSave_Brecon", "Save Brecon table (1=on 0=off)", &m->m_SetSave_Brecon);

  d->define_param( "HadronB_cut", "HadronB cut", &m->m_HadronB_cut );

  d->define_param( "useIpProfile", "IpProfile (1=on, 0=off)", &m->m_useIpProfile);

  d->define_param( "useBenergy", "Benergy (2=on(database), 1=on(Benergy func.), 0=off)", 
		   &m->m_useBenergy);

  d->define_param( "dEdxSel",  "dEdx sel. (1=on, 0=off)", &m->m_dEdxSel );

  d->define_param( "CorrEpi0Drho",  "Epi0 corr for Drho (1=on, 0=off)", &m->m_corrEpi0Drho );

  d->define_param( "VetoDa1",  "Da1 veto (1=on, 0=off)",       &m->m_vetoDa1 );

  d->define_param( "reconDDs",  "DDs recon (1=on, 0=off)",     &m->m_reconDDs );

  d->define_param( "selDDs",   "DDs selection (1=on, 0=off)",  &m->m_selDDs );

  d->define_param( "CutCont",  "Continuum cut (1=on, 0=off)",  &m->m_cutCont );

  d->define_param( "BestB", "Best B selection (1=on, 0=off)",  &m->m_bestB );

  d->define_param( "Dr_cut", "dr cut value", &m->m_Dr_cut);
  d->define_param( "Dz_cut", "dz cut value", &m->m_Dz_cut);

  d->define_param( "Pt_cut", "Pt cut value", &m->m_Pt_cut);

  d->define_param( "Ene_cut",         "Energy cut value",             &m->m_Ene_cut);
  d->define_param( "Ene_cut_slowPI0", "Energy cut value for slowPI0", &m->m_Ene_cut_slowPI0);

  d->define_param( "cutPi0", "pi0 mass cut (1=on, 0=off)", &m->m_cutPi0 );
  d->define_param( "Mdst_pi0_lower_limit", "Set Mdst_pi0 lower mass limit", &m->m_Mdst_pi0_lower_limit );
  d->define_param( "Mdst_pi0_upper_limit", "Set Mdst_pi0 upper mass limit", &m->m_Mdst_pi0_upper_limit );

  d->define_param( "dEdxSigma_cut", "dEdx sigma cut value", &m->m_dEdxSigma_cut);

  d->define_param( "D_Pstar_cut_low", "D(*) P* cut value", &m->m_D_Pstar_cut_low);
  d->define_param( "H_Pstar_cut_low", "H P* cut value",    &m->m_H_Pstar_cut_low);

  d->define_param( "ACCq", "PID value ACC", &m->m_ACCq);
  d->define_param( "TOFq", "PID value TOF", &m->m_TOFq);
  d->define_param( "CDCq", "PID value CDC", &m->m_CDCq);

  d->define_param( "PID_cut_K",        "PID cut kaon prob(K:PI) value",        &m->m_PID_cut_K);
  d->define_param( "PID_cut_K_loose",  "PID cut kaon prob(K:PI) value(loose)", &m->m_PID_cut_K_loose);
  d->define_param( "PID_cut_PI",       "PID cut pion prob(K:PI) value",        &m->m_PID_cut_PI);
  d->define_param( "PID_cut_PI_loose", "PID cut pion prob(K:PI) value(loose)", &m->m_PID_cut_PI);
  d->define_param( "PID_cut_PI_D",   "PID cut pion from D   prob(K:PI) value", &m->m_PID_cut_PI_D);
  d->define_param( "PID_cut_PI_rho", "PID cut pion from rho prob(K:PI) value", &m->m_PID_cut_PI_rho);
  d->define_param( "PID_cut_PI_a1",  "PID cut pion from a1  prob(K:PI) value", &m->m_PID_cut_PI_a1);

  d->define_param( "Ks_fanfan_cut", "turn on goodKs cut (1=on, 0=off)", &m->m_Ks_fanfan_cut);
  d->define_param( "Ks_cut", "Ks mass cut",                             &m->m_Ks_cut);

  d->define_param( "RHO_cut", "RHO mass cut",                 &m->m_RHO_cut);
  d->define_param( "RHO_Pstar_cut_low", "RHO P* low value",   &m->m_RHO_Pstar_cut_low);
  d->define_param( "RHO_Pstar_cut_high", "RHO P* high value", &m->m_RHO_Pstar_cut_high);

  d->define_param( "RHO0_cut", "RHO0 mass cut",                 &m->m_RHO0_cut);
  d->define_param( "RHO0_Pstar_cut_low",  "RHO0 P* low value",  &m->m_RHO0_Pstar_cut_low);
  d->define_param( "RHO0_Pstar_cut_high", "RHO0 P* high value", &m->m_RHO0_Pstar_cut_high);

  d->define_param( "A1_cut_low",  "A1 mass cut low",        &m->m_A1_cut_low);
  d->define_param( "A1_cut_high", "A1 mass cut high",       &m->m_A1_cut_high);
  d->define_param( "A1_Pstar_cut_low", "A1 P* low value",   &m->m_A1_Pstar_cut_low);
  d->define_param( "A1_Pstar_cut_high", "A1 P* high value", &m->m_A1_Pstar_cut_high);
  d->define_param( "A1_Chisq_cut", "A1 vertex chisq",       &m->m_A1_Chisq_cut);

  d->define_param( "Phi_cut", "Phi mass cut",               &m->m_Phi_cut);
  d->define_param( "Kst_cut", "Kst mass cut",               &m->m_Kst_cut);

  d->define_param( "D_fitmass1",   "D fit mass1",   &m->m_D_fitmass[1] );
  d->define_param( "D_fitmass2",   "D fit mass2",   &m->m_D_fitmass[2] );
  d->define_param( "D_fitmass3",   "D fit mass3",   &m->m_D_fitmass[3] );
  d->define_param( "D_fitmass4",   "D fit mass4",   &m->m_D_fitmass[4] );
  d->define_param( "D_fitmass5",   "D fit mass5",   &m->m_D_fitmass[5] );
  d->define_param( "D_fitmass6",   "D fit mass6",   &m->m_D_fitmass[6] );
  d->define_param( "D_fitmass7",   "D fit mass7",   &m->m_D_fitmass[7] );
  d->define_param( "D_fitmass101", "D fit mass101", &m->m_D_fitmass[101] );
  d->define_param( "D_fitmass102", "D fit mass102", &m->m_D_fitmass[102] );
  d->define_param( "D_fitmass103", "D fit mass103", &m->m_D_fitmass[103] );
  d->define_param( "D_fitmass104", "D fit mass104", &m->m_D_fitmass[104] );
  d->define_param( "D_fitmass105", "D fit mass105", &m->m_D_fitmass[105] );
  d->define_param( "D_fitmass106", "D fit mass106", &m->m_D_fitmass[106] );
  d->define_param( "D_fitmass301", "D fit mass301", &m->m_D_fitmass[301] );
  d->define_param( "D_fitmass302", "D fit mass302", &m->m_D_fitmass[302] );
  d->define_param( "D_fitmass303", "D fit mass303", &m->m_D_fitmass[303] );

  d->define_param( "D_cut1",   "D mass cut1",   &m->m_D_cut[1]);
  d->define_param( "D_cut2",   "D mass cut2",   &m->m_D_cut[2]);
  d->define_param( "D_cut3",   "D mass cut3",   &m->m_D_cut[3]);
  d->define_param( "D_cut4",   "D mass cut4",   &m->m_D_cut[4]);
  d->define_param( "D_cut5",   "D mass cut5",   &m->m_D_cut[5]);
  d->define_param( "D_cut6",   "D mass cut6",   &m->m_D_cut[6]);
  d->define_param( "D_cut7",   "D mass cut7",   &m->m_D_cut[7]);
  d->define_param( "D_cut101", "D mass cut101", &m->m_D_cut[101]);
  d->define_param( "D_cut102", "D mass cut102", &m->m_D_cut[102]);
  d->define_param( "D_cut103", "D mass cut103", &m->m_D_cut[103]);
  d->define_param( "D_cut104", "D mass cut104", &m->m_D_cut[104]);
  d->define_param( "D_cut105", "D mass cut105", &m->m_D_cut[105]);
  d->define_param( "D_cut106", "D mass cut106", &m->m_D_cut[106]);
  d->define_param( "D_cut301", "D mass cut301", &m->m_D_cut[301]);
  d->define_param( "D_cut302", "D mass cut302", &m->m_D_cut[302]);
  d->define_param( "D_cut303", "D mass cut303", &m->m_D_cut[303]);

  d->define_param( "Dalitz_cut", "Dalitz plot weight cut",             &m->m_Dalitz_cut);
  d->define_param( "D_vertexfit", "turn on D vertexfit (1=on, 0=off)", &m->m_D_vertexfit);
  d->define_param( "D_Chisq_cut", "D vertex chisq",                    &m->m_D_Chisq_cut);
  d->define_param( "D_massfit", "turn on D(*) massfit (1=on, 0=off)",  &m->m_D_massfit);

  d->define_param( "Dm_cut1",  "D* - D mass diff. cut1", &m->m_Dm_cut[1] );
  d->define_param( "Dm_cut2",  "D* - D mass diff. cut2", &m->m_Dm_cut[2] );
  d->define_param( "Dm_cut3",  "D* - D mass diff. cut3", &m->m_Dm_cut[3] );
  d->define_param( "Dm_cut4",  "D* - D mass diff. cut4", &m->m_Dm_cut[4] );
  d->define_param( "Dm_cut5",  "D* - D mass diff. cut5", &m->m_Dm_cut[5] );

  d->define_param( "MB_cut_low",  "B mass cut  (low)",   &m->m_MB_cut_low  );
  d->define_param( "MB_cut_high", "B mass cut  (high)",  &m->m_MB_cut_high );
  d->define_param( "DE_cut_low",  "Delta E cut (low)",   &m->m_DE_cut_low  );
  d->define_param( "DE_cut_low",  "Delta E cut (high)",  &m->m_DE_cut_high );

  d->define_param( "B_Thrust_cut",        "cos theta_T",       &m->m_B_Thrust_cut );
  d->define_param( "B_Thrust_cut_loose",  "cos theta_T_loose", &m->m_B_Thrust_cut_loose );
  d->define_param( "B_Thrust_cut_tight",  "cos theta_T_tight", &m->m_B_Thrust_cut_tight );

  d->define_param( "B_Direction_cut",     "cos theta_B",       &m->m_B_Direction_cut );
  d->define_param( "R2_cut", "R2",                             &m->m_R2_cut);

  d->define_param( "signal_MB_low",  "signal box Mb (low)",  &m->m_signal_MB_low );
  d->define_param( "signal_MB_high", "signal box Mb (high)", &m->m_signal_MB_high );
  d->define_param( "signal_DE_low",  "signal box DE (low)",  &m->m_signal_DE_low );
  d->define_param( "signal_DE_high", "signal box DE (high)", &m->m_signal_DE_high );

  return d;

}

// -----------------------
// Default constructor
// -----------------------
frec_ana::frec_ana() {

  m_name  = "frec_ana";

  m_Monitor        = 0;
  m_flagHist       = 1;

  mm_mc = 2;
  m_mc = 0;

  m_SetSave_Brecon = 0;
  m_HadronB_cut    = 1;
  m_useIpProfile   = 1;
  m_useBenergy     = 2;
  m_dEdxSel        = 0;
  m_corrEpi0Drho   = 0;
  m_vetoDa1        = 1;
  m_reconDDs       = 1;
  m_selDDs         = 1;
  m_cutCont        = 1;    
  m_bestB          = 1;

  // charged K & PI cut
  m_Dr_cut = 2.0;           // cm
  m_Dz_cut = 5.0;           // cm
  m_Pt_cut = 0.0;           // GeV/c

  // PI0 & gamma cut
  m_Ene_cut         = 0.05; // GeV
  m_Ene_cut_slowPI0 = 0.03; // (for slow pi0)

  m_cutPi0 = 1;
  m_Mdst_pi0_lower_limit = 0.1178; // GeV/c^2
  m_Mdst_pi0_upper_limit = 0.1502;

  // dEdx sigma cut
  m_dEdxSigma_cut   = 5.0;

  // P(*) cut
  m_D_Pstar_cut_low =  1.0; // for D(*) from B-->D(*)X
  m_H_Pstar_cut_low  = 1.5; // for prompt hadron from B-->D(*)pi

  // Ks cut
  m_Ks_fanfan_cut = 1;      // 0:OFF 1:ON 
  m_Ks_cut  = 0.03;         // GeV/c^2

  // PID cut
  m_ACCq = 3;
  m_TOFq = 1;
  m_CDCq = 5;

  m_PID_cut_PI       = 1.1;
  m_PID_cut_PI_loose = 1.1;
  m_PID_cut_K        = 0.6;
  m_PID_cut_K_loose  = 0.3;

  m_PID_cut_PI_D     = 1.1;
  m_PID_cut_PI_rho   = 0.99;
  m_PID_cut_PI_a1    = 0.99;

  // RHO cut
  m_RHO_cut              = 0.225;
  m_RHO_Pstar_cut_low    = 1.8;
  m_RHO_Pstar_cut_high   = 2.6;

  m_RHO0_cut             = 0.225;
  m_RHO0_Pstar_cut_low   = 0.0;
  m_RHO0_Pstar_cut_high  = 2.6;

  // a1 cut
  m_A1_cut_low         =  0.7;
  m_A1_cut_high        =  1.6;
  m_A1_Pstar_cut_low   =  1.5;
  m_A1_Pstar_cut_high  =  2.6;
  m_A1_Chisq_cut       = -1.;

  // phi, K* cut from Ds
  m_Phi_cut  = 0.020;
  m_Kst_cut  = 0.075;

  // D cut
  m_D_fitmass[1]   = 1.865;
  m_D_fitmass[2]   = 1.863;
  m_D_fitmass[3]   = 1.865;
  m_D_fitmass[4]   = 1.863;
  m_D_fitmass[5]   = 1.865;
  m_D_fitmass[6]   = 1.863;
  m_D_fitmass[7]   = 1.865;
  m_D_fitmass[101] = 1.870;
  m_D_fitmass[102] = 1.868;
  m_D_fitmass[103] = 1.870;
  m_D_fitmass[104] = 1.868;
  m_D_fitmass[105] = 1.870;
  m_D_fitmass[106] = 1.870;
  m_D_fitmass[301] = 1.969;
  m_D_fitmass[302] = 1.969;
  m_D_fitmass[303] = 1.969;

  m_D_cut[1]   = 0.006*5;
  m_D_cut[2]   = 0.015*3;
  m_D_cut[3]   = 0.006*5;
  m_D_cut[4]   = 0.030*2;
  m_D_cut[5]   = 0.006*5;
  m_D_cut[6]   = 0.006*5;
  m_D_cut[7]   = 0.006*5;

  m_D_cut[101] = 0.006*5;
  m_D_cut[102] = 0.015*3;
  m_D_cut[103] = 0.006*5;
  m_D_cut[104] = 0.015*3;
  m_D_cut[105] = 0.006*5;
  m_D_cut[106] = 0.006*5;

  m_D_cut[301] = 0.006*5;
  m_D_cut[302] = 0.006*5;
  m_D_cut[303] = 0.006*5;

  m_Dalitz_cut  = -1;
  m_D_vertexfit =  1;
  m_D_Chisq_cut = -1.;
  m_D_massfit   =  0;

  m_Dm[1] = Ptype_Dst0.mass()     - Ptype_D0.mass();
  m_Dm[2] = m_Dm[1];
  m_Dm[3] = Ptype_Dstplus.mass()  - Ptype_D0.mass();
  m_Dm[4] = Ptype_Dstplus.mass()  - Ptype_Dplus.mass();
  m_Dm[5] = Ptype_DSstplus.mass() - Ptype_DSplus.mass();

  m_Dm_cut[1] = 0.001*5;
  m_Dm_cut[2] = 0.010*2; 
  m_Dm_cut[3] = 0.0005*10;
  m_Dm_cut[4] = 0.001*5;
  m_Dm_cut[5] = 0.010*2; 

  // B cut
  m_MB_cut_low  =  5.2;
  m_MB_cut_high =  5.3;
  m_DE_cut_low  = -0.3;
  m_DE_cut_high =  0.3;

  m_B_Thrust_cut       = 0.8; 
  m_B_Thrust_cut_loose = 1.1;
  m_B_Thrust_cut_tight = 0.8;
  m_B_Direction_cut    = 1.1;
  m_R2_cut             = 0.5;

  m_signal_MB_low  =  5.27;
  m_signal_MB_high =  5.29;
  m_signal_DE_low  = -0.08;
  m_signal_DE_high =  0.06;

  // counter
  int offset = 1000;
  m_cid_evt       = offset +  0;
  m_cid_hadb      = offset +  1;
  m_cid_ip        = offset +  2;
  m_cid_benergy   = offset +  3;
  m_cid_b         = offset + 10;
  m_cid_b2dpi     = offset + 11;
  m_cid_b2drho    = offset + 12;
  m_cid_b2da1     = offset + 13;
  m_cid_b2dds     = offset + 14;
  m_cid_d         = offset + 20;
  m_cid_dst       = offset + 21;
  m_cid_a1        = offset + 22;
  m_cid_rho       = offset + 23;
  m_cid_ks        = offset + 24;
  m_cid_pi0       = offset + 25;
  m_cid_b_mc      = offset + 30;
  m_cid_b2dpi_mc  = offset + 31;
  m_cid_b2drho_mc = offset + 32;
  m_cid_b2da1_mc  = offset + 33;
  m_cid_b2dds_mc  = offset + 34;
  m_cid_d_mc      = offset + 40;
  m_cid_dst_mc    = offset + 41;
  m_cid_a1_mc     = offset + 42;
  m_cid_rho_mc    = offset + 43;
  m_cid_ks_mc     = offset + 44;
  m_cid_pi0_mc    = offset + 45;

  for ( unsigned i = 0; i<256; i++ ) m_LogFileName[i] = 0;

}

// ------------------
// Init function
// ------------------
void  frec_ana::init(int  *status){

  dout(Debugout::INFO,"frec_ana") << "[0m[01;35m"                                         << "\n"
	    << "Hello BASF"                                             << "\n"
            << "##### frec_ana, basf_if parameters #####" << "\n"
            << "Monitor          = " << m_Monitor 
	    << " : (0=off, 1=on)" << "\n"
            << "flagHist         = " << m_flagHist
	    << " : (0=off, 1=on)" << "\n"
            << "MC_flag          = " << mm_mc 
	    << " : (0=off, 1=on(for old), 2=on(for new))" << "\n"
            << "SetSave_Brecon   = " << m_SetSave_Brecon 
	    << " : (0=off, 1=on)" << "\n"
            << "HadronB_cut      = " << m_HadronB_cut
	    << " : (0=off, 1=on)" << "\n"
            << "useIpProfile     = " << m_useIpProfile 
	    << " : (0=off, 1=on)" << "\n"
            << "useBenergy       = " << m_useBenergy 
	    << " : (0=off, 1=on(Benergy func.), 2=on(database))" << "\n"
            << "dEdxSel          = " << m_dEdxSel
	    << " : (0=off, 1=on)" << "\n"
            << "pi0 mass cut     = " << m_cutPi0
            << " : (0=off, 1=on)" ;
  if ( m_cutPi0 ){
    dout(Debugout::INFO,"frec_ana") << " : [ "       << m_Mdst_pi0_lower_limit
              << " < Mpi0 < "  << m_Mdst_pi0_upper_limit      
	      << " GeV/c^2 ]";
  }
  dout(Debugout::INFO,"frec_ana") << std::endl;
  dout(Debugout::INFO,"frec_ana") << "Ks_fanfan_cut    = " << m_Ks_fanfan_cut
	    << " : (0=off, 1=on)" << "\n"
            << "D_vertexfit      = " << m_D_vertexfit
	    << " : (0=off, 1=on)" << "\n"
            << "D_massfit        = " << m_D_massfit
	    << " : (0=off, 1=on)" << "\n"
            << "corrEpi0Drho     = " << m_corrEpi0Drho 
	    << " : (0=off, 1=on)" << "\n"
            << "vetoDa1          = " << m_vetoDa1 
	    << " : (0=off, 1=on)" << "\n"
            << "reconDDs         = " << m_reconDDs
	    << " : (0=off, 1=on)" << "\n"
            << "selDDs           = " << m_selDDs
	    << " : (0=off, 1=on)" << "\n"
            << "cutCont          = " << m_cutCont
	    << " : (0=off, 1=on)" << "\n"
            << "bestB            = " << m_bestB
            << " : (0=off, 1=on) " << "\n"
	    << "[0m" << std::endl;

  if ( m_Monitor ){

  m_PID_cut_PI      = 1.1;
  m_PID_cut_PI_D    = 1.1;
  m_PID_cut_PI_rho  = 1.1;
  m_PID_cut_PI_a1   = 1.1;

  // RHO cut
  m_RHO_cut              = 0.3;
  m_RHO_Pstar_cut_low    = 1.8;
  m_RHO_Pstar_cut_high   = 2.6;

  m_RHO0_cut             = 0.25;
  m_RHO0_Pstar_cut_low   = 0.0;
  m_RHO0_Pstar_cut_high  = 2.6;

  // a1 cut
  m_A1_cut_low         = 0.7;
  m_A1_cut_high        = 1.8;
  m_A1_Pstar_cut_low   = 1.5;
  m_A1_Pstar_cut_high  = 2.6;
  m_A1_Chisq_cut       = -1.;

  m_D_fitmass[1]   = Ptype_D0.mass();
  m_D_fitmass[2]   = m_D_fitmass[1];
  m_D_fitmass[3]   = m_D_fitmass[1];
  m_D_fitmass[4]   = m_D_fitmass[1];
  m_D_fitmass[5]   = m_D_fitmass[1];
  m_D_fitmass[6]   = m_D_fitmass[1];
  m_D_fitmass[7]   = m_D_fitmass[1];
  m_D_fitmass[101] = Ptype_Dplus.mass();
  m_D_fitmass[102] = m_D_fitmass[101];
  m_D_fitmass[103] = m_D_fitmass[101];
  m_D_fitmass[104] = m_D_fitmass[101];
  m_D_fitmass[105] = m_D_fitmass[101];
  m_D_fitmass[106] = m_D_fitmass[101];
  m_D_fitmass[301] = Ptype_DSplus.mass();
  m_D_fitmass[302] = m_D_fitmass[301];
  m_D_fitmass[303] = m_D_fitmass[301];

  m_D_cut[1]   = 0.05;
  m_D_cut[2]   = m_D_cut[1];
  m_D_cut[3]   = m_D_cut[1];
  m_D_cut[4]   = 0.07;
  m_D_cut[5]   = m_D_cut[1];
  m_D_cut[6]   = m_D_cut[1];
  m_D_cut[7]   = m_D_cut[1];

  m_D_cut[101] = m_D_cut[1];
  m_D_cut[102] = m_D_cut[1];
  m_D_cut[103] = m_D_cut[1];
  m_D_cut[104] = m_D_cut[1];
  m_D_cut[105] = m_D_cut[1];
  m_D_cut[106] = m_D_cut[1];

  m_D_cut[301] = m_D_cut[1];
  m_D_cut[302] = m_D_cut[1];
  m_D_cut[303] = m_D_cut[1];

  m_Dalitz_cut  = -1.;
  m_D_Chisq_cut = -1.;

  m_Dm_cut[1] = 0.02;
  m_Dm_cut[2] = 0.03; 
  m_Dm_cut[3] = 0.02;
  m_Dm_cut[4] = 0.02;
  m_Dm_cut[5] = 0.03;

  }

  kid_statistics::init_sta( m_name );

  return;

}

// --------------------
// Term function
// -------------------
void  frec_ana::term(){

  dout(Debugout::INFO,"frec_ana") << "[0m[01;35m" << "\n"
	    << "Good bye BASF"  << "\n"
	    << "[0m"          << std::endl;
  return;

}

// --------------------
// Begin run
// -------------------
void  frec_ana::begin_run( BelleEvent  *evptr, int  *status ){

  // Judge data is Real Data or MC
  Belle_runhead_Manager&  rhd_mgr = Belle_runhead_Manager::get_manager();
  std::vector<Belle_runhead>::const_iterator  rhd = rhd_mgr.begin();

  if (rhd == rhd_mgr.end()) {
    dout(Debugout::ERR,"frec_ana") << "!!!!!Can not Access to Belle_RunHead!!!!!" << std::endl;
    dout(Debugout::ERR,"frec_ana") << "!!!!!Regard the data as MC sample!!!!!" << std::endl;
  }
  else if (rhd->ExpMC() == 1) {
    m_mc = 0;

    dout(Debugout::INFO,"frec_ana") << "[0m[01;34m"           << "\n"
	      << "Exp Analysis (frec_ana)"  << "\n"
	      << "ExpNo>> " << rhd->ExpNo() << "\n"
    	      << "RunNo>> " << rhd->RunNo() << "\n" 
	      << "Event>> " << rhd->NEvt()  << "\n"
	      << "[0m"    << std::endl; 
  }
  else {

    // set MC flag
    m_mc = mm_mc;             

    dout(Debugout::INFO,"frec_ana") << "[0m[01;34m"           << "\n"
	      << "MC Analysis (frec_ana, MC flag = "  << m_mc  << ")\n"
	      << "ExpNo>> " << rhd->ExpNo() << "\n"
    	      << "RunNo>> " << rhd->RunNo() << "\n" 
	      << "Event>> " << rhd->NEvt()  << "\n"
	      << "[0m"    << std::endl; 
  }

  // get IP profile data from database
  if( m_useIpProfile ){
    IpProfile::begin_run();

  dout(Debugout::INFO,"frec_ana") << "[0m[01;34m"         << "\n"
	    << "IP correct>> "  << IpProfile::position() << "\n"
	    << "[0m"          << std::endl;

  }

  // get Benergy data from database
  if( m_useBenergy == 2) {
    BeamEnergy::begin_run();
    dout(Debugout::INFO,"frec_ana") << "[0m[01;34m"      << "\n"
	      << "Benergy correct>> "  << BeamEnergy::E_beam2() << "\n"
	      << "[0m"               << std::endl;
  }else if ( m_useBenergy == 1 ){
    dout(Debugout::INFO,"frec_ana") << "[0m[01;34m"      << "\n"
	      << "Benergy correct>> "  << Benergy() << "\n"
	      << "[0m"               << std::endl;
  }

  return;

}

// -----------------
// display status
// -----------------
void  frec_ana::disp_stat( const char *disp ) {

  print_summary( dout(Debugout::RESULT,"frec_ana") );

  std::string sfile(m_LogFileName);
  if( sfile.length() == 0 ) return;
  dout(Debugout::INFO,"frec_ana") << "frec_ana> log file = [" << sfile << "]" << std::endl;
  std::ofstream sout(sfile.c_str());
  
  if( !sout ){
    dout(Debugout::ERR,"frec_ana") << "!!! can't open the file" << std::endl;
    return;
  }

    print_summary( sout );

}

void  frec_ana::print_summary( std::ostream & out ){

    int evt          = (int) kid_statistics::out_n( m_name, m_cid_evt );
    int hadb_evt     = (int) kid_statistics::out_n( m_name, m_cid_hadb );
    int ip_evt       = (int) kid_statistics::out_n( m_name, m_cid_ip );
    int benergy_evt  = (int) kid_statistics::out_n( m_name, m_cid_benergy );
    int b_evt        = (int) kid_statistics::out_n( m_name, m_cid_b );
    int b2dpi_evt    = (int) kid_statistics::out_n( m_name, m_cid_b2dpi );
    int b2drho_evt   = (int) kid_statistics::out_n( m_name, m_cid_b2drho );
    int b2da1_evt    = (int) kid_statistics::out_n( m_name, m_cid_b2da1 );
    int b2dds_evt    = (int) kid_statistics::out_n( m_name, m_cid_b2dds );
    int d_evt        = (int) kid_statistics::out_n( m_name, m_cid_d );
    int dst_evt      = (int) kid_statistics::out_n( m_name, m_cid_dst );
    int a1_evt       = (int) kid_statistics::out_n( m_name, m_cid_a1 );
    int rho_evt      = (int) kid_statistics::out_n( m_name, m_cid_rho );
    int ks_evt       = (int) kid_statistics::out_n( m_name, m_cid_ks );
    int pi0_evt      = (int) kid_statistics::out_n( m_name, m_cid_pi0 );
    int b_mcevt      = (int) kid_statistics::out_n( m_name, m_cid_b_mc );
    int b2dpi_mcevt  = (int) kid_statistics::out_n( m_name, m_cid_b2dpi_mc );
    int b2drho_mcevt = (int) kid_statistics::out_n( m_name, m_cid_b2drho_mc );
    int b2da1_mcevt  = (int) kid_statistics::out_n( m_name, m_cid_b2da1_mc );
    int b2dds_mcevt  = (int) kid_statistics::out_n( m_name, m_cid_b2dds_mc );
    int d_mcevt      = (int) kid_statistics::out_n( m_name, m_cid_d_mc );
    int dst_mcevt    = (int) kid_statistics::out_n( m_name, m_cid_dst_mc );
    int a1_mcevt     = (int) kid_statistics::out_n( m_name, m_cid_a1_mc );
    int rho_mcevt    = (int) kid_statistics::out_n( m_name, m_cid_rho_mc );
    int ks_mcevt     = (int) kid_statistics::out_n( m_name, m_cid_ks_mc );
    int pi0_mcevt    = (int) kid_statistics::out_n( m_name, m_cid_pi0_mc );

    out << std::endl;
    out << "frec_ana> ********** frec_ana diagnostics **********"
        << std::endl;

    out << std::endl;
    out << "---------- Event Counter ----------" << std::endl;
    out << "number of processed events : " << evt << std::endl;
    out << "number of HadronB   events : " << hadb_evt << std::endl;
    out << "number of IP used          : " << ip_evt << std::endl;
    out << "number of BeamEnergy used  : " << benergy_evt << std::endl;
    out << "number of B events         : " << b_evt << std::endl;
    out << "number of B(D pi )  events : " << b2dpi_evt << std::endl;
    out << "number of B(D rho)  events : " << b2drho_evt << std::endl;
    out << "number of B(D a1 )  events : " << b2da1_evt << std::endl;
    out << "number of B(D Ds )  events : " << b2dds_evt << std::endl;
    out << "number of good D    events : " << d_evt << std::endl;
    out << "number of good D*   events : " << dst_evt << std::endl;
    out << "number of good a1   events : " << a1_evt << std::endl;
    out << "number of good rho  events : " << rho_evt << std::endl;
    out << "number of good Ks   events : " << ks_evt << std::endl;
    out << "number of good PI0  events : " << pi0_evt << std::endl;
    out << std::endl;
  if( m_mc ) {
    out << "---------- Event Counter (MC true) ----------" << std::endl;
    out << "number of processed events : " << evt << std::endl;
    out << "number of B events         : " << b_mcevt << std::endl;
    out << "number of B(D pi )  events : " << b2dpi_mcevt << std::endl;
    out << "number of B(D rho)  events : " << b2drho_mcevt << std::endl;
    out << "number of B(D a1 )  events : " << b2da1_mcevt << std::endl;
    out << "number of B(D Ds )  events : " << b2dds_mcevt << std::endl;
    out << "number of good D    events : " << d_mcevt << std::endl;
    out << "number of good D*   events : " << dst_mcevt << std::endl;
    out << "number of good a1   events : " << a1_mcevt << std::endl;
    out << "number of good rho  events : " << rho_mcevt << std::endl;
    out << "number of good Ks   events : " << ks_mcevt << std::endl;
    out << "number of good PI0  events : " << pi0_mcevt << std::endl;
 }

}

// ---------------------
//   Event function
// ---------------------
void  frec_ana::event( BelleEvent  *evptr, int  *status ) {

  *status = 0;

  // --- count processed events
  kid_statistics::fill_sta( m_name, m_cid_evt, (float)1. ); 
  if ( m_flagHist ) H[1]->accumulate( (float)1.0, 1.0 );

  // -- count HadronB events
  bool  flag_HadronB = frec_util::HadronB();
  if ( flag_HadronB ) kid_statistics::fill_sta( m_name, m_cid_hadb, (float)1. ); 

  if ( m_HadronB_cut && flag_HadronB == false ) return;

  /// --- count IP used events
  if ( IpProfile::usable() ) kid_statistics::fill_sta( m_name, m_cid_ip, (float)1. ); 

  /// --- count BeamEnergy used events
  if ( BeamEnergy::usable() ) kid_statistics::fill_sta( m_name, m_cid_benergy, (float)1. ); 

  // create Mdst_sim_xref table
  if ( m_mc == 1 ) mdst2xref();

  // vector to store Candidate
  std::vector<Particle> All_ptcle;
  std::vector<Particle> Charged_all;
  std::vector<Particle> Pi_all;
  std::vector<Particle> PiLoose_all;
  std::vector<Particle> K_all;
  std::vector<Particle> Gamma;
  std::vector<Particle> Pi_zero, slowPi_zero;
  std::vector<Particle> K_short;
  std::vector<Particle> Rho_all;
  std::vector<Particle> Rho_zero;
  std::vector<Particle> A1_all;

  std::vector<Particle> D_zero;
  std::vector<Particle> D_pm;
  std::vector<Particle> Dst_zero;
  std::vector<Particle> Dst_pm;
  std::vector<Particle> DS_pm;
  std::vector<Particle> DSst_pm;

  std::vector<Particle> B_zero;
  std::vector<Particle> B_pm;
  std::vector<Particle> select_B_zero;
  std::vector<Particle> select_B_pm;

  // fundamental particles
  // ~~~~~~~~~~~~~~~~~~~~~~
  recon_Charged( Charged_all, Pi_all, PiLoose_all, K_all );
  recon_Gam( Gamma );
  recon_Pi0( Pi_zero, slowPi_zero );
  recon_Ks( K_short );
  recon_ALL( All_ptcle, Charged_all, Gamma );

  m_R2 = frec_util::R2( All_ptcle );
  if ( m_flagHist ) H[11]->accumulate( m_R2,1. );
  if ( m_cutCont && m_R2 > m_R2_cut ) return;

  recon_Rho( Rho_all, Pi_all, Pi_zero );
  recon_Rho0( Rho_zero, Pi_all );
  recon_A1( A1_all, Rho_zero, Pi_all );

  // D(*) reconstruction
  //~~~~~~~~~~~~~~~~~~~~
  recon_D( D_zero, D_pm, DS_pm, 
	   K_all, Pi_all, Pi_zero, K_short, Gamma );

  recon_Dst( Dst_zero, Dst_pm, DSst_pm,
	     D_zero, D_pm, DS_pm, Charged_all, slowPi_zero, Gamma ); 

  // B reconstruction
  // ~~~~~~~~~~~~~~~~~
  recon_B( B_zero, D_pm,     PiLoose_all, Rho_all, A1_all );
  recon_B( B_zero, Dst_pm,   PiLoose_all, Rho_all, A1_all );
  recon_B( B_pm,   D_zero,   PiLoose_all, Rho_all, A1_all );
  recon_B( B_pm,   Dst_zero, PiLoose_all, Rho_all, A1_all );

  if ( m_reconDDs ){
    recon_B_dds( B_zero, D_pm,     DS_pm   );
    recon_B_dds( B_zero, D_pm,     DSst_pm );
    recon_B_dds( B_zero, Dst_pm,   DS_pm   );
    recon_B_dds( B_zero, Dst_pm,   DSst_pm );
    recon_B_dds( B_pm,   D_zero,   DS_pm   );
    recon_B_dds( B_pm,   D_zero,   DSst_pm );
    recon_B_dds( B_pm,   Dst_zero, DS_pm   );
    recon_B_dds( B_pm,   Dst_zero, DSst_pm );
  }

  if ( B_zero.size() ) sel_B(B_zero, All_ptcle, select_B_zero);
  if ( B_pm.size() )   sel_B(B_pm,   All_ptcle, select_B_pm);

  if ( m_flagHist && (select_B_zero.size() || select_B_pm.size()) ) H[2]->accumulate((float)1.0, 1.0);

  if ( select_B_zero.size() ) set_bestB( select_B_zero );
  if ( select_B_pm.size() )   set_bestB( select_B_pm );

   if ( m_Monitor ){

     if ( select_B_zero.size() ) FillHist_Monitor( select_B_zero );
     if ( select_B_pm.size() )   FillHist_Monitor( select_B_pm );

   }else{

     if ( select_B_zero.size() ) FillHist_B( select_B_zero );
     if ( select_B_pm.size() )   FillHist_B( select_B_pm );

   }

   count_B( select_B_pm, select_B_zero );

    // save Brecon table, set *status
    if ( select_B_zero.size() + select_B_pm.size() ){
       *status = 1;
       if ( m_SetSave_Brecon ) save_Brecon( select_B_zero, select_B_pm );
    }

  return;

}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
