//
// $Id: hist_def.cc 9932 2006-11-12 14:26:53Z katayama $
//
// $Log$
// Revision 1.5  2004/05/04 21:23:20  matumot
// fixed for CC complier
// - namespace was changed : fullrecon -> fullreconNS
//   to avoid overlap for fullrecon struct used by fullrecon.tdf
// - use new operator for dynamic allocation of array
// - add std:: for sort function
//
// Revision 1.3  2004/04/19 09:27:34  matumot
// revival version
//
//
// =====================================================
//  File Name : hist_def.cc
// -----------------------------------------------------
//  Creation    ; 2004.04.07
//  Description ; Implimentation for the frec_ana class
//                ( definition of histograms )
//  Author      ; T.Matsumoto (matumot@bmail.kek.jp)
// -----------------------------------------------------

#include "belle.h"
#include  "fullrecon/frec_ana.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


// ------------------------------------
//    histgram & ntuple define
// ------------------------------------
void  frec_ana::hist_def(void){

  extern BelleTupleManager  *BASF_Histogram;
  BelleTupleManager  *tm = BASF_Histogram;

  // Define 1D & 2D histogram
  // ~~~~~~~~~~~~~~~~~~~~~~~~
  // Event Counter
  H[1] = tm->histogram("Event Counter (All)", 6, -0.5, 5.5, 1);
  H[2] = tm->histogram("Event Counter (B)",   6, -0.5, 5.5);

  // Several basic distributions
  H[11] = tm->histogram("R2", 100, 0., 1., 11);
  H[12] = tm->histogram("dr (charged)", 300,  -3.,  3.);
  H[13] = tm->histogram("dz (charged)", 500, -10., 10.);
  H[14] = tm->histogram("Pt (charged)", 100,   0.,  5.);
  H[15] = tm->histogram("charge (charged)", 5, -2.5, 2.5);
  H[16] = tm->histogram("PID prob(K:PI)", 200,0.,1.);

  // Track numbers
  H[21] = tm->histogram("track number (charged)", 51, -0.5, 50.5, 21);
  H[22] = tm->histogram("track number (pi)",      51, -0.5, 50.5);
  H[23] = tm->histogram("track number (K)",       51, -0.5, 50.5);
  H[24] = tm->histogram("track number (gamma)",   51, -0.5, 50.5);
  H[25] = tm->histogram("track number (pi0)",     51, -0.5, 50.5);
  H[26] = tm->histogram("track number (Ks)",      51, -0.5, 50.5);

  // Momentum distributions
  H[31] = tm->histogram( "P(pi)", 200,0.,5.,31);
  H[32] = tm->histogram( "P(K)" , 200,0.,5.);
  H[33] = tm->histogram( "E(gamma)", 1000,0.,5.);
  H[34] = tm->histogram( "E(gamma) from pi0", 1000,0.,5.);
  H[35] = tm->histogram( "P(pi0)", 200,0.,5.);
  H[36] = tm->histogram( "P(Ks)", 200,0.,5.);

  H[41] = tm->histogram( "P(pi), MC", 200,0.,5.,41);
  H[42] = tm->histogram( "P(K),  MC" , 200,0.,5.);
  H[43] = tm->histogram( "E(gamma), MC", 1000,0.,5.);
  H[44] = tm->histogram( "E(gamma) from pi0, MC", 1000,0.,5.);
  H[45] = tm->histogram( "P(pi0), MC", 200,0.,5.);
  H[46] = tm->histogram( "P(Ks), MC", 200,0.,5.);

  // Mass distributions
  H[51] = tm->histogram( "M(pi0)", 200,0.11,0.16,51);
  H[52] = tm->histogram( "M(Ks)", 160,0.47,0.53 );

  H[61] = tm->histogram( "M(pi0), MC", 200,0.11,0.16,61);
  H[62] = tm->histogram( "M(Ks), MC", 160,0.47,0.53 );

  // Define Ntuple
  // ~~~~~~~~~~~~~
  nt_b = tm->ntuple("B_cand",    "Mbc deltaE best best2 best3 b_lund b_mode dst_mode d_mode dst2mode ds_mode pi_p pi_pid rhomass rho_p a1mass a1_p a1chisq dm_dst dmass dm_dst2 dsmass dvchisq dsvchisq d_p ds_p r2 cosT cosB bnum mc chisq purity", 1000);

  nt_dpi  = tm->ntuple("B_Dpi",  "Mbc deltaE b_lund b_mode dst_mode d_mode pi_p pi_pid dm_dst ms_dst ps_dst dmass dvchisq d_p kd_pid pid_pid mpi0_d ppi0_d mks_d pks_d mkpi_d mkpi0_d dalitz_d r2 cosT cosB mc chisq", 1001);
  nt_drho = tm->ntuple("B_Drho", "Mbc deltaE b_lund b_mode dst_mode d_mode rhomass rho_p rpi_p rpi0_p rpi0_m pi_pid dm_dst dmass dvchisq d_p kd_pid pid_pid r2 cosT cosB mc chisq", 1002);
  nt_da1  = tm->ntuple("B_Da1",   "Mbc deltaE b_lund b_mode dst_mode d_mode a1mass a1chisq a1_p rhomass rho_p rpi1_p rpi2_p api_p pi_pid dm_dst dmass dvchisq d_p kd_pid pid_pid r2 cosT cosB mc chisq", 1003);

  return;

}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
