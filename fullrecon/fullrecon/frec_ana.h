//
// $Id: frec_ana.h 9932 2006-11-12 14:26:53Z katayama $
//
// $Log$
// Revision 1.10  2004/11/22 05:24:36  matumot
// updates for new kid_statistics
// module parameter, "LogFileName" is introduced
//
// Revision 1.9  2004/07/29 06:53:11  matumot
// minor modification on counting method
//
// Revision 1.8  2004/05/28 10:05:54  matumot
// adjusted pi0 mass cut to 0.1178 < Mpi0 < 0.1502 GeV/c^2
// ( consistent with pi0 mass cuts in current fix_mdst module )
//
// Revision 1.7  2004/05/04 21:22:37  matumot
// fixed for CC complier
// - namespace was changed : fullrecon -> fullreconNS
//   to avoid overlap for fullrecon struct used by fullrecon.tdf
// - use new operator for dynamic allocation of array
// - add std:: for sort function
//
// Revision 1.5  2004/04/19 09:35:21  matumot
// revival version
//
//
// ===================================================
//  File Name : frec_ana.h
// ---------------------------------------------------
//  Creation    ; 2004.04.07
//  Description ; Definition for the frec_ana class
//  Author      ; T.Matsumoto (matumot@bmail.kek.jp)
// ---------------------------------------------------

#ifndef  __FREC_ANA__
#define  __FREC_ANA__

#include <vector>
#include <iosfwd>

#include "belle.h"

#include  "basf/module.h"
#include  "particle/Particle.h"
#include  "kid/atc_pid.h"
#include  "benergy/BeamEnergy.h"

#include  "tuple/BelleTupleManager.h"

#include  HEPEVT_H               
#include  MDST_H                 
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// --------------------------
//  frec_ana class
// --------------------------
class frec_ana : public Module {

public :

  frec_ana ();
 ~frec_ana (){}

public :

  void  init( int* );
  void  term();

  void  begin_run( BelleEvent*, int* );
  void  end_run( BelleEvent*, int* ){}

  void  event( BelleEvent*, int* );

  void  hist_def( void );

  void  disp_stat( const  char* );
  void  other( int*, BelleEvent*, int* ) {}
  void  print_summary( std::ostream& );

  char m_LogFileName[256];

protected:

  // member functions
  // ~~~~~~~~~~~~~~~~~
  void  recon_Charged( std::vector<Particle>&, std::vector<Particle>&, 
 		       std::vector<Particle>&, std::vector<Particle>& );

  void  recon_Gam( std::vector<Particle>& );

  void  recon_Pi0( std::vector<Particle>&, std::vector<Particle>& );

  void  recon_Ks( std::vector<Particle>& );

  void  recon_Rho( std::vector<Particle>&, std::vector<Particle>&,
		   std::vector<Particle>& );

  void  recon_Rho0( std::vector<Particle>&, std::vector<Particle>& );

  void  recon_A1( std::vector<Particle>&, std::vector<Particle>&, 
		  std::vector<Particle>& );

  void  recon_D( std::vector<Particle>&, std::vector<Particle>&,
		 std::vector<Particle>&, std::vector<Particle>&,
		 std::vector<Particle>&, std::vector<Particle>&,
		 std::vector<Particle>&, std::vector<Particle>& );

  void  recon_Dst( std::vector<Particle>&, std::vector<Particle>&,
		   std::vector<Particle>&, std::vector<Particle>&,
		   std::vector<Particle>&, std::vector<Particle>&,
		   std::vector<Particle>&, std::vector<Particle>&,
		   std::vector<Particle>& );

  void  recon_B( std::vector<Particle>&, std::vector<Particle>&,
		 std::vector<Particle>&, std::vector<Particle>&,
		 std::vector<Particle>& );

  void  recon_B_dds( std::vector<Particle>&, std::vector<Particle>&,
		     std::vector<Particle>& );

  void  recon_ALL( std::vector<Particle>&, std::vector<Particle>&,
		   std::vector<Particle>& );

  bool sel_D( Particle&, const int );

  bool sel_Dst( Particle&, const int );

  bool presel_B( Particle&, const int );

  void sel_B( std::vector<Particle>&, std::vector<Particle>&,
	      std::vector<Particle>& );

  void set_bestB( std::vector<Particle>& );

  void FillHist_B( std::vector<Particle>& );
  void FillHist_Monitor( std::vector<Particle>& );
  
  void count_B( std::vector<Particle>&, std::vector<Particle>& );

  int  save_Brecon( std::vector<Particle>&, std::vector<Particle>& );

  int get_b_mode( Particle& );

private:

  int m_mc;               

public :                

  int m_Monitor;
  int m_flagHist;
 
  int mm_mc;

  int m_SetSave_Brecon;

  int m_HadronB_cut;
  int m_useIpProfile;
  int m_useBenergy;
  int m_dEdxSel;
  int m_corrEpi0Drho;
  int m_vetoDa1;
  int m_reconDDs;
  int m_selDDs;
  int m_cutCont;
  int m_bestB;

  double m_Dr_cut;
  double m_Dz_cut;

  double m_Pt_cut;
  double m_Ene_cut;
  double m_Ene_cut_slowPI0;

  int    m_cutPi0;
  double m_Mdst_pi0_lower_limit;
  double m_Mdst_pi0_upper_limit;

  double m_dEdxSigma_cut;
  double m_D_Pstar_cut_low;
  double m_H_Pstar_cut_low;
  
  int  m_ACCq;
  int  m_TOFq;
  int  m_CDCq;

  double m_PID_cut_K;
  double m_PID_cut_K_loose;
  double m_PID_cut_PI;
  double m_PID_cut_PI_loose;
  double m_PID_cut_PI_D;
  double m_PID_cut_PI_rho;
  double m_PID_cut_PI_a1;

  double m_Ks_fanfan_cut;
  double m_Ks_cut;

  double m_RHO_cut;
  double m_RHO_Pstar_cut_low;
  double m_RHO_Pstar_cut_high;

  double m_RHO0_cut;
  double m_RHO0_Pstar_cut_low;
  double m_RHO0_Pstar_cut_high;

  double m_A1_cut_low;
  double m_A1_cut_high;
  double m_A1_Pstar_cut_low;
  double m_A1_Pstar_cut_high;
  double m_A1_child_Pstar_cut;
  double m_A1_Chisq_cut;

  double m_Kst_cut;
  double m_Phi_cut;

  double m_D_fitmass[1000];
  double m_D_cut[1000];
  double m_Dalitz_cut;
  int    m_D_vertexfit;
  double m_D_Chisq_cut;
  int    m_D_massfit;

  double m_Dm[10];
  double m_Dm_cut[10];

  double m_MB_cut_low;
  double m_MB_cut_high;
  double m_DE_cut_low;
  double m_DE_cut_high;

  double m_B_Thrust_cut;
  double m_B_Thrust_cut_loose;
  double m_B_Thrust_cut_tight;
  double m_B_Direction_cut;
  double m_R2_cut;

  double m_signal_MB_low;
  double m_signal_MB_high;
  double m_signal_DE_low;
  double m_signal_DE_high;

private :               

  double m_R2;

  BelleHistogram  *H[990];

  BelleTuple  *nt_b;
  BelleTuple  *nt_dpi, *nt_drho, *nt_da1;

  // counter
  int m_cid_evt;
  int m_cid_hadb;
  int m_cid_ip;
  int m_cid_benergy;
  int m_cid_b;
  int m_cid_b2dpi;
  int m_cid_b2drho;
  int m_cid_b2da1;
  int m_cid_b2dds;
  int m_cid_d;
  int m_cid_dst;
  int m_cid_a1;
  int m_cid_rho;
  int m_cid_ks;
  int m_cid_pi0;
  int m_cid_b_mc;
  int m_cid_b2dpi_mc;
  int m_cid_b2drho_mc;
  int m_cid_b2da1_mc;
  int m_cid_b2dds_mc;
  int m_cid_d_mc;
  int m_cid_dst_mc;
  int m_cid_a1_mc;
  int m_cid_rho_mc;
  int m_cid_ks_mc;
  int m_cid_pi0_mc;

  const char* m_name;

};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
