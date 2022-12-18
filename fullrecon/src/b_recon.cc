//
// $Id: b_recon.cc 10548 2008-07-24 13:13:11Z hitoshi $
//
// $Log$
// Revision 1.10  2004/11/22 05:24:36  matumot
// updates for new kid_statistics
// module parameter, "LogFileName" is introduced
//
// Revision 1.9  2004/07/29 06:53:11  matumot
// minor modification on counting method
//
// Revision 1.8  2004/05/04 21:23:19  matumot
// fixed for CC complier
// - namespace was changed : fullrecon -> fullreconNS
//   to avoid overlap for fullrecon struct used by fullrecon.tdf
// - use new operator for dynamic allocation of array
// - add std:: for sort function
//
// Revision 1.7  2004/04/29 21:19:20  katayama
// doesn't compile with gcc-2.95.3
//
// Revision 1.6  2004/04/29 12:29:42  katayama
// for gcc33
//
// Revision 1.4  2004/04/19 09:27:33  matumot
// revival version
//
// ===================================================
//  File Name :  b_recon.cc
// ---------------------------------------------------
// Creation    ; 2004.04.07
// Description ; Implimentaion for the frec_ana class
//               ( B reconstruction )
// Author      ; T.Matsumoto (matumot@bmail.kek.jp)
// ----------------------------------------------------
#include "belle.h"
#include <algorithm>

#include  "fullrecon/frec_ana.h"
#include  "fullrecon/frec_util.h"

#include "particle/Particle.h"
#include "particle/ParticleUserInfo.h"

#include "belleCLHEP/Matrix/SymMatrix.h"
#include "kid/kid_statistics.h"
#include "benergy/BeamEnergy.h"

#include "hamlet/AnaBrecon.h"
#include  BRECON_H         
#include  FULLRECON_H

#include "belleutil/debugout.h"

#include "mdst/mdst.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


using namespace fullreconNS;

inline bool b_lowerChisq( Particle b1, Particle b2 ){
 return  dynamic_cast<UserInfo_B&>(b1.userInfo()).masschisq() < 
         dynamic_cast<UserInfo_B&>(b2.userInfo()).masschisq();
}

inline bool b_higherPurity( Particle b1, Particle b2 ){
 return  dynamic_cast<UserInfo_B&>(b1.userInfo()).purity() >
         dynamic_cast<UserInfo_B&>(b2.userInfo()).purity();
}

inline bool b_withBest( Particle b1, Particle b2 ){
 return  dynamic_cast<UserInfo_B&>(b1.userInfo()).flag_best() >
         dynamic_cast<UserInfo_B&>(b2.userInfo()).flag_best();
}

// Dalitz weight function
// ~~~~~~~~~~~~~~~~~~~~~~~
extern "C" {
  extern float weight_(float*, float*, float*);
};


// ---------------------------
//   B reconstruction
// ---------------------------
void  frec_ana::recon_B( std::vector<Particle>& B,  std::vector<Particle>& D,
			 std::vector<Particle>& PI, std::vector<Particle>& RHO,
			 std::vector<Particle>& A1 ){

  Ptype ptype;
  int b_mode(0);
  //++++++ b__mode +++++++++++
  //   1: B0B -->D+pi-
  //   2: B0B -->D*+pi- 
  //   3: B0B -->D+rho-
  //   4: B0B -->D*+rho-
  //   5: B0B -->D+a1-
  //   6: B0B -->D*+a1- 
  // 101: B-  -->D0pi-
  // 102: B-  -->D*0pi- 
  // 103: B-  -->D0rho-
  // 104: B-  -->D*0rho-
  // 105: B-  -->D0a1-
  // 106: B-  -->D*0a1-
  //+++++++++++++++++++++++++

  for( std::vector<Particle>::iterator d_i = D.begin();
       d_i != D.end(); d_i++ ) {    

    bool flag_monitor = false;

    int dst_mode = dynamic_cast<UserInfo_B&>(d_i->userInfo()).sub0_mode();
    int d_mode = dynamic_cast<UserInfo_B&>(d_i->userInfo()).sub1_mode();

    HepLorentzVector p_d = frec_util::p_cm(d_i->p());
    if ( p_d.vect().mag() < m_D_Pstar_cut_low ) continue;

    // monitor when D0-->K-pi+  from D0X-, or D*+(D0pi+)X-
    if ( d_mode == 1 && (dst_mode == 0 || dst_mode == 3) ) flag_monitor = true;

    // B -> D(*)pi-
    // ~~~~~~~~~~~~~~ 
    for( std::vector<Particle>::iterator pi_j = PI.begin();
	 pi_j != PI.end(); pi_j++ ) {

      if ( d_i->lund() * pi_j->lund() > 0 ) continue;

      HepLorentzVector p_pi = frec_util::p_cm(pi_j->p());
      if ( p_pi.vect().mag() < m_H_Pstar_cut_low ) continue;

      // remove duplication
      bool chk_flag = frec_util::prohib_dupli( *d_i, *pi_j );
      if( chk_flag == true ) continue;

      int charge = (int)d_i->charge() + (int)pi_j->charge();
      if ( charge == 0 ){
	if ( pi_j->charge() < 0 ) ptype = Ptype_B0B;
	else			  ptype = Ptype_B0;
      }else if ( charge > 0 )	  ptype = Ptype_Bplus;
      else			  ptype = Ptype_Bminus;

      Particle b = frec_util::make_2body( ptype, *d_i, *pi_j, m_mc );
      b_mode = get_b_mode( b );
      if ( presel_B( b , b_mode ) == true ) B.push_back(b);

    }// end of loop pi_j

    if ( m_Monitor && flag_monitor == false ) continue;

    // B -> D(*)rho- 
    // ~~~~~~~~~~~~~~~
    for( std::vector<Particle>::iterator rho_j = RHO.begin();
	 rho_j != RHO.end(); rho_j++ ) {

      if ( d_i->lund() * rho_j->lund() > 0 ) continue;

      // remove duplication
      bool chk_flag = frec_util::prohib_dupli( *d_i, *rho_j );
      if( chk_flag == true ) continue;
       
      int charge = (int)d_i->charge() + (int)rho_j->charge();
      if ( charge == 0 ){
	if ( rho_j->charge() < 0 ) ptype = Ptype_B0B;
	else			   ptype = Ptype_B0;
      }else if ( charge > 0 )	   ptype = Ptype_Bplus;
      else			   ptype = Ptype_Bminus;

      Particle b = frec_util::make_2body( ptype, *d_i, *rho_j, m_mc );
      b_mode = get_b_mode( b );
      if ( presel_B( b , b_mode ) == true ) B.push_back(b);

    } // end of loop rho_j

    // B -> D(*)a1-
    // ~~~~~~~~~~~~~~
    bool flag_da1 = true;

    if ( m_vetoDa1 ){

      // Modes for Veto
      // B-  --> D0a1-             : except for D0 -->K-pi+
      // B0B --> D+a1-             : except for D+ -->K-pi+pi+, Kspi+
      // B-  --> D*0(D0pi0)a1-     : D0 --> Kspi+pi-pi0
      // B-  --> D*0(D0gam)a1-     : except for D0 -->K-pi+
      // B0B --> D*+(D+pi0)a1-     : except for D+ -->K-pi+pi+, Kspi+

      if ( dst_mode == 0 ){
	if ( !(d_mode == 1 || d_mode == 101 || d_mode == 103) ) flag_da1 = false;
      }else if ( dst_mode == 1 ){
	if ( d_mode == 6 ) flag_da1 = false;
      }else if ( dst_mode == 2 ){
	if ( d_mode != 1 ) flag_da1 = false;
      }else if ( dst_mode == 4 ){
	if ( !(d_mode == 101 || d_mode == 103) ) flag_da1 = false;
      }

    }

    if ( flag_da1 == true ){

    for( std::vector<Particle>::iterator a1_j = A1.begin();
	 a1_j != A1.end(); a1_j++ ) {

      if ( d_i->lund() * a1_j->lund() > 0 ) continue;
      
      // remove duplication
      bool chk_flag = frec_util::prohib_dupli( *d_i, *a1_j );
      if( chk_flag == true ) continue;
       
      int charge = (int)d_i->charge() + (int)a1_j->charge();
      if ( charge == 0 ){
	if ( a1_j->charge() < 0 ) ptype = Ptype_B0B;
	else			  ptype = Ptype_B0;
      }else if ( charge > 0 )	  ptype = Ptype_Bplus;
      else			  ptype = Ptype_Bminus;

      Particle b = frec_util::make_2body( ptype, *d_i, *a1_j, m_mc );
      b_mode = get_b_mode( b );
      if ( presel_B( b , b_mode ) == true ) B.push_back(b);

     } // end of loop a_j

    }

  } // end of loop d_i
 
  return;
}

void  frec_ana::recon_B_dds( std::vector<Particle>& B,  
			     std::vector<Particle>& D, std::vector<Particle>& DS ){

  bool flag_monitor = false;
  if ( m_Monitor && flag_monitor == false ) return;

  Ptype ptype;
  int b_mode(0);
  //++++++ b__mode +++++++++++
  //  11: B0B -->D+Ds-
  //  12: B0B -->D+Ds*- 
  //  13: B0B -->D*+Ds-
  //  14: B0B -->D*+Ds*-
  // 111: B-  -->D0Ds-
  // 112: B-  -->D0Ds*- 
  // 113: B-  -->D*0Ds-
  // 114: B-  -->D*0Ds*-
  //+++++++++++++++++++++++++

  for( std::vector<Particle>::iterator d_i = D.begin();
       d_i != D.end(); d_i++ ) {    

    HepLorentzVector p_d = frec_util::p_cm(d_i->p());
    if ( p_d.vect().mag() < m_D_Pstar_cut_low ) continue;

    const int dst_mode = dynamic_cast<UserInfo_B&>(d_i->userInfo()).sub0_mode();
    const int d_mode = dynamic_cast<UserInfo_B&>(d_i->userInfo()).sub1_mode();

    for( std::vector<Particle>::iterator ds_j = DS.begin();
	 ds_j != DS.end(); ds_j++ ) {

      if ( d_i->lund() * ds_j->lund() > 0 ) continue;

      const int ds_mode = dynamic_cast<UserInfo_B&>(ds_j->userInfo()).sub1_mode();

      HepLorentzVector p_ds = frec_util::p_cm(ds_j->p());
      if ( p_ds.vect().mag() < m_D_Pstar_cut_low ) continue;

      // remove duplication
      bool chk_flag = frec_util::prohib_dupli( *d_i, *ds_j );
      if( chk_flag == true ) continue;

     bool flag_dds = true;

     if ( m_selDDs ){

      // *** DDs selection ***
      // - Ds ->KsK   : select all
      // - Ds ->KKpi  : apply phi mass cut except for
      //       D*+ --> D0pi+, D*0 --> D0pi0 (except for D0 --> Kspipipi0)
      //       D0  --> Kpi, D+ ->Kpipi, Kspi
      // - Ds ->Kpipi : don't select 

      if ( ds_mode == 302 ){

	if ( !ds_j ->usable() ) flag_dds = false;
	if ( dst_mode == 3 ) flag_dds = true;
	else if ( dst_mode == 1 && d_mode != 6 ) flag_dds = true;
	else if ( d_mode == 1 || d_mode == 101 || d_mode == 103 ) flag_dds = true;

      }else if ( ds_mode == 303 ) flag_dds = false;

    }

    if ( flag_dds == true ){

      int charge = (int)d_i->charge() + (int)ds_j->charge();
      if ( charge == 0 ){
	if ( ds_j->charge() < 0 ) ptype = Ptype_B0B;
	else			  ptype = Ptype_B0;
      }else if ( charge > 0 )	  ptype = Ptype_Bplus;
      else			  ptype = Ptype_Bminus;

      Particle b = frec_util::make_2body( ptype, *d_i, *ds_j, m_mc );
      b_mode = get_b_mode( b );
      if ( presel_B( b , b_mode ) == true ) B.push_back(b);

    }

    }// end of loop ds_j
  } // end of loop d_i

}

// -------------------------------
//   B selection
// -------------------------------
bool frec_ana::presel_B( Particle &b, const int b_mode ){

 bool preselect = false;

  const HepLorentzVector p_b = frec_util::p_cm( b.p() );

  double Ebeam;
  if ( m_useBenergy == 2 )      Ebeam = BeamEnergy::E_beam2();
  else if ( m_useBenergy == 1 ) Ebeam = Benergy();
  else                          Ebeam = frec_util::Ebeam();

  const double Delta_E = p_b.t() - Ebeam;

  double Mass_B2;
  // pi0 energy correction for B-->D(*)rho
  if ( m_corrEpi0Drho && abs(b.child(1).pType().lund()) == Ptype_RHOplus.lund() ){

    HepLorentzVector p_d = frec_util::p_cm(b.child(0).p());
    HepLorentzVector p_pirho = frec_util::p_cm(b.child(1).child(0).p());
    HepLorentzVector p_pi0rho = frec_util::p_cm(b.child(1).child(1).p());

    double epi0_new = Ebeam - p_d.e() - p_pirho.e();
    double pi0mass = Ptype_PI0.mass();

    if ( epi0_new > pi0mass ){

      double ppi0     = p_pi0rho.vect().mag();
      double ppi0_new = sqrt(epi0_new*epi0_new - pi0mass*pi0mass);

      HepLorentzVector pcorr_pi0( (ppi0_new/ppi0)*p_pi0rho.vect(), epi0_new );
      HepLorentzVector pcorr_b = p_d + p_pirho + pcorr_pi0;
      Mass_B2 =  Ebeam* Ebeam - pcorr_b.vect().mag2();
    }else {
      Mass_B2 = -1.;
    }

  }else {
    Mass_B2 =  Ebeam* Ebeam - p_b.vect().mag2();
  }

  if( Mass_B2 < 0.) return preselect;
  const double Mass_B  = sqrt(Mass_B2);
  if(!( Mass_B >  m_MB_cut_low && Mass_B < m_MB_cut_high && 
        Delta_E > m_DE_cut_low && Delta_E < m_DE_cut_high) ) return preselect;

  preselect = true;
  
  UserInfo_B x;  x.set(b);
  const int dst_mode  = dynamic_cast<UserInfo_B&>(b.child(0).userInfo()).sub0_mode();
  const int d_mode    = dynamic_cast<UserInfo_B&>(b.child(0).userInfo()).sub1_mode();
  int dsst_mode(0), ds_mode(0);  

  if ( abs(b.child(1).lund()) == Ptype_DSplus.lund() ||
       abs(b.child(1).lund()) == Ptype_DSstplus.lund() ){
    dsst_mode = dynamic_cast<UserInfo_B&>(b.child(1).userInfo()).sub0_mode();
    ds_mode   = dynamic_cast<UserInfo_B&>(b.child(1).userInfo()).sub1_mode();
  }

  dynamic_cast<UserInfo_B&>(b.userInfo()).Mass_B(Mass_B);
  dynamic_cast<UserInfo_B&>(b.userInfo()).Delta_E(Delta_E);
  dynamic_cast<UserInfo_B&>(b.userInfo()).b_mode(b_mode);
  dynamic_cast<UserInfo_B&>(b.userInfo()).sub0_mode(dst_mode);
  dynamic_cast<UserInfo_B&>(b.userInfo()).sub1_mode(d_mode);
  dynamic_cast<UserInfo_B&>(b.userInfo()).sub2_mode(dsst_mode);
  dynamic_cast<UserInfo_B&>(b.userInfo()).sub3_mode(ds_mode);
  if( m_mc && b.genHepevt()) dynamic_cast<UserInfo_B&>(b.userInfo()).flag_mc(1);

  return preselect;

}

void frec_ana::sel_B( std::vector<Particle>& B, std::vector<Particle>& ALL,
                      std::vector<Particle>& selB )
{

  if ( B.size() == 0 ) return;

  for(std::vector<Particle>::iterator b_i = B.begin();
      b_i != B.end(); b_i++) {

    bool select = false;

    const int b_mode    = dynamic_cast<UserInfo_B&>(b_i->userInfo()).b_mode();
    const int dst_mode  = dynamic_cast<UserInfo_B&>(b_i->userInfo()).sub0_mode();
    const int d_mode    = dynamic_cast<UserInfo_B&>(b_i->userInfo()).sub1_mode();
    const int dsst_mode = dynamic_cast<UserInfo_B&>(b_i->userInfo()).sub2_mode();
    const int ds_mode   = dynamic_cast<UserInfo_B&>(b_i->userInfo()).sub3_mode();

    double purity = frec_util::Bpurity_hadronic_DX( b_mode, dst_mode, d_mode, dsst_mode, ds_mode );
    if ( m_Monitor == 0 && m_bestB && purity <= 0. ) continue;

    // Thrust angle
    // ~~~~~~~~~~~~
//    bool good_thrust = true;
    std::vector<HepLorentzVector> plist_other;
    std::vector<HepLorentzVector> plist_bcand;

    frec_util::division( *b_i, ALL, plist_bcand, plist_other );
    thrust bcand( plist_bcand );  // Thrust for B candidate tracks
    thrust other( plist_other );  // Thrust for remaining tracks
    const float cos_th   = other.thr_axis().dot(bcand.thr_axis());
    if(plist_other.size()<2){
      dout(Debugout::DDEBUG,"frec_ana")
	<<" plist_other() size is "<<plist_other.size()<<" !!!!!!!" <<" thr: "<< other.thr()
	  <<" cos_th: "<<cos_th<<std::endl;
    }
    
//    if( bcand.thr() < 0.5 || other.thr() < 0.5 ) good_thrust = false;
    
//   if ( good_thrust == false ) {
//     std::cout<<"WARNING!!! rejected by good_thrust cut: "<<good_thrust
//       <<" bcand.thr(): "<<bcand.thr()<<" other.thr(): "<<other.thr()<<std::endl;
//     continue;
//   }
    
    // B flight derection
    //~~~~~~~~~~~~~~~~~~~~
    const HepLorentzVector p_b = frec_util::p_cm(b_i->p());
    const double cosB    = p_b.vect().cosTheta();

    // Select B candidates
    //~~~~~~~~~~~~~~~~~~~~~
    double thrust_cut = m_B_Thrust_cut;

    // for B- -->D0(Kpi)pi-, B- -->D*0(D0pi0)pi- and B0B -->D*+(D0pi+)pi-, give loose thrust cut
    if ( ( b_mode == 101 && dst_mode == 0 && d_mode == 1 ) ||
	 ( b_mode == 102 && dst_mode == 1 ) ||
	 ( b_mode ==   2 && dst_mode == 3 ) ) thrust_cut = m_B_Thrust_cut_loose;

    // for B --> D(*)a1, give tighter thrust cut
    if ( abs(b_i->child(1).lund()) == Ptype_A1plus.lund() ) thrust_cut = m_B_Thrust_cut_tight;

    if( m_cutCont == 0 ) select = true;
    else if( abs(cos_th)< thrust_cut && abs(cosB) < m_B_Direction_cut ) select = true;

    if( select == true ) {
      dynamic_cast<UserInfo_B&>(b_i->userInfo()).purity(purity);
      dynamic_cast<UserInfo_B&>(b_i->userInfo()).cos_thetaT(cos_th);
      dynamic_cast<UserInfo_B&>(b_i->userInfo()).cos_thetaB(cosB);
      selB.push_back(*b_i);
    }
  
   } // end of loop b_i

  return;

}

void frec_ana::set_bestB( std::vector<Particle>& B ){

  if( B.size() == 0 ) return;

  // calculate reduced chisq for each event
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for( std::vector<Particle>::iterator b_i = B.begin();
       b_i != B.end(); b_i++ ) {
 
    // selection with mass chisq
    const double Delta_E = dynamic_cast<UserInfo_B&>(b_i->userInfo()).Delta_E();

    // assume DeltaE_err is equal to it's err of invariant mass
    // ( effect of beam energy spread is neglected )
    double DeltaE_err(sqrt(b_i->momentum().dMass()));
    if ( !(DeltaE_err <100 ) || DeltaE_err == 0 ) DeltaE_err = 0.02;
    const double de_chisq = (Delta_E/DeltaE_err)*(Delta_E/DeltaE_err);

    double dof(1.), chisq(de_chisq);

     if( ( abs(b_i->child(0).lund()) == Ptype_D0.lund() ) || 
	  ( abs(b_i->child(0).lund()) == Ptype_Dplus.lund() ) ) {
	const double d_chisq = dynamic_cast<UserInfo_B&>(b_i->child(0).userInfo()).masschisq();
	chisq += d_chisq; dof   += 1.;
      }else if( ( abs(b_i->child(0).lund()) == Ptype_Dst0.lund() ) ||
	        ( abs(b_i->child(0).lund()) == Ptype_Dstplus.lund() ) ) {
	double dst_chisq   = 
	  dynamic_cast<UserInfo_B&>(b_i->child(0).userInfo()).masschisq();
	double d_chisq    = 
	  dynamic_cast<UserInfo_B&>(b_i->child(0).child(0).userInfo()).masschisq();
	chisq += (dst_chisq + d_chisq); dof   += 2.;
      }

     if( abs(b_i->child(1).lund()) == Ptype_DSplus.lund() ) {
	const double ds_chisq = dynamic_cast<UserInfo_B&>(b_i->child(1).userInfo()).masschisq();
	chisq += ds_chisq; dof   += 1.;
      }else if( abs(b_i->child(1).lund()) == Ptype_DSstplus.lund() ){
	double dsst_chisq   = 
	  dynamic_cast<UserInfo_B&>(b_i->child(1).userInfo()).masschisq();
	double ds_chisq    = 
	  dynamic_cast<UserInfo_B&>(b_i->child(1).child(0).userInfo()).masschisq();
	chisq += (dsst_chisq + ds_chisq); dof   += 2.;
      }

     dynamic_cast<UserInfo_B&>(b_i->userInfo()).masschisq(chisq/dof);

  } // end of loop b_i

  // select the best reduced chisq for each mode
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // int flag_veto[B.size()];
  int* flag_veto = new int[B.size()];

  for ( int i=0; i< B.size(); i++) flag_veto[i] = 0;

  int count = - 1;
  for( std::vector<Particle>::iterator b_i = B.begin();
       b_i != B.end(); b_i++ ) {

    count += 1;
    if ( flag_veto[count] ) continue;

    const int    b_mode_i    = dynamic_cast<UserInfo_B&>(b_i->userInfo()).b_mode();
    const int    dst_mode_i  = dynamic_cast<UserInfo_B&>(b_i->userInfo()).sub0_mode();
    const int    d_mode_i    = dynamic_cast<UserInfo_B&>(b_i->userInfo()).sub1_mode();
    const int    dsst_mode_i = dynamic_cast<UserInfo_B&>(b_i->userInfo()).sub2_mode();
    const int    ds_mode_i   = dynamic_cast<UserInfo_B&>(b_i->userInfo()).sub3_mode();
    const double chisq_i     = dynamic_cast<UserInfo_B&>(b_i->userInfo()).masschisq();

    int count2 = count;
    for( std::vector<Particle>::iterator b_j = b_i + 1;
	 b_j != B.end(); b_j++ ) {

    count2 += 1;

    const int    b_mode_j    = dynamic_cast<UserInfo_B&>(b_j->userInfo()).b_mode();
    const int    dst_mode_j  = dynamic_cast<UserInfo_B&>(b_j->userInfo()).sub0_mode();
    const int    d_mode_j    = dynamic_cast<UserInfo_B&>(b_j->userInfo()).sub1_mode();
    const int    dsst_mode_j = dynamic_cast<UserInfo_B&>(b_j->userInfo()).sub2_mode();
    const int    ds_mode_j   = dynamic_cast<UserInfo_B&>(b_j->userInfo()).sub3_mode();
    const double chisq_j     = dynamic_cast<UserInfo_B&>(b_j->userInfo()).masschisq();

      if ( b_mode_j    == b_mode_i &&
	   dst_mode_j  == dst_mode_i &&
	   d_mode_j    == d_mode_i && 
	   dsst_mode_j == dsst_mode_i &&
	   ds_mode_j   == ds_mode_i ){

	  if ( chisq_j  <  chisq_i ) flag_veto[count]  = 1;
          else                       flag_veto[count2] = 1;

      }

    } // end of loop b_j
  } // end of loop b_i

  count = -1;
  for( std::vector<Particle>::iterator b_i = B.begin();
       b_i != B.end(); b_i++ ) {
      count += 1;
      if ( flag_veto[count] == 0 ) dynamic_cast<UserInfo_B&>(b_i->userInfo()).flag_bestEach(1);
  }

  delete [] flag_veto;

  if ( m_bestB == 0 ) return;

    // select the best candidate for each event with best purity
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    int best_id(-1), best_id2(-1);
    double max_purity(-1.), max_purity2(-1.);

    count  = -1;
    for( std::vector<Particle>::iterator b_i = B.begin();
       b_i != B.end(); b_i++ ) {

    count +=1;

      if ( dynamic_cast<UserInfo_B&>(b_i->userInfo()).flag_bestEach() == 0 ) continue;
      if ( dynamic_cast<UserInfo_B&>(b_i->userInfo()).purity() <= 0. )       continue;

      const double purity =  dynamic_cast<UserInfo_B&>(b_i->userInfo()).purity();

      if( purity > max_purity2 || best_id2 < 0 ){
          max_purity2 = purity; best_id2   = count;
      }

      const double deltaE = dynamic_cast<UserInfo_B&>(b_i->userInfo()).Delta_E();
      if ( deltaE < m_signal_DE_low || deltaE > m_signal_DE_high ) continue;

      if( purity > max_purity || best_id < 0 ){
	   max_purity = purity; best_id   = count;
        }

    } // end of loop i

   if ( best_id  >= 0 ) dynamic_cast<UserInfo_B&>(B[best_id].userInfo()).flag_best(1);
   if ( best_id2 >= 0 ) dynamic_cast<UserInfo_B&>(B[best_id2].userInfo()).flag_best2(1);

    // select the best candidate for each event with best reduced chisq
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    int best_id3(-1);

    count  = -1;
    double min_chisq(-1.);
    for( std::vector<Particle>::iterator b_i = B.begin();
       b_i != B.end(); b_i++ ) {

      count +=1;

      if ( dynamic_cast<UserInfo_B&>(b_i->userInfo()).flag_bestEach() == 0 ) continue;
      if ( dynamic_cast<UserInfo_B&>(b_i->userInfo()).purity() <= 0. )       continue;

      const double chisq =  dynamic_cast<UserInfo_B&>(b_i->userInfo()).masschisq();

      if( chisq < min_chisq || best_id3 < 0 ){
          min_chisq = chisq; best_id3  = count;
      }

    } // end of loop i

   if ( best_id3 >= 0 ) dynamic_cast<UserInfo_B&>(B[best_id3].userInfo()).flag_best3(1);

  return;

}

// -------------------------------
//   Fill B info.
// -------------------------------
void frec_ana::FillHist_B( std::vector<Particle>& B ){


  const atc_pid selkpi( m_ACCq, m_TOFq, m_CDCq, KAON_CODE, PION_CODE );
  const atc_pid selkpi_noTOF( m_ACCq, -1, m_CDCq, KAON_CODE, PION_CODE );

  int bnum(0);
  for( std::vector<Particle>::iterator b_i = B.begin();
       b_i != B.end(); b_i++ ) {
       const int flag_bestEach = dynamic_cast<UserInfo_B&>(b_i->userInfo()).flag_bestEach();
       if ( flag_bestEach == 0 ) continue;
       const double purity  = dynamic_cast<UserInfo_B&>(b_i->userInfo()).purity();
       if ( m_bestB && purity <=0. ) continue;
       bnum += 1;
  }

  for( std::vector<Particle>::iterator b_i = B.begin();
       b_i != B.end(); b_i++ ) {

    const int flag_bestEach = dynamic_cast<UserInfo_B&>(b_i->userInfo()).flag_bestEach();
    if ( flag_bestEach == 0 ) continue;
    const int flag_best     = dynamic_cast<UserInfo_B&>(b_i->userInfo()).flag_best();
    const int flag_best2    = dynamic_cast<UserInfo_B&>(b_i->userInfo()).flag_best2();
    const int flag_best3    = dynamic_cast<UserInfo_B&>(b_i->userInfo()).flag_best3();

    const int b_mode        = dynamic_cast<UserInfo_B&>(b_i->userInfo()).b_mode();
    const int dst_mode      = dynamic_cast<UserInfo_B&>(b_i->userInfo()).sub0_mode();
    const int d_mode        = dynamic_cast<UserInfo_B&>(b_i->userInfo()).sub1_mode();
    const int dsst_mode     = dynamic_cast<UserInfo_B&>(b_i->userInfo()).sub2_mode();
    const int ds_mode       = dynamic_cast<UserInfo_B&>(b_i->userInfo()).sub3_mode();

    if ( b_mode == 0 || d_mode == 0 ) continue;

    const double Mass_B    = dynamic_cast<UserInfo_B&>(b_i->userInfo()).Mass_B();
    const double Delta_E   = dynamic_cast<UserInfo_B&>(b_i->userInfo()).Delta_E();
    const double cosT      = dynamic_cast<UserInfo_B&>(b_i->userInfo()).cos_thetaT();
    const double cosB      = dynamic_cast<UserInfo_B&>(b_i->userInfo()).cos_thetaB();
    const double chisq     = dynamic_cast<UserInfo_B&>(b_i->userInfo()).masschisq();
    const double purity    = dynamic_cast<UserInfo_B&>(b_i->userInfo()).purity();
    const int    flag_mc   = dynamic_cast<UserInfo_B&>(b_i->userInfo()).flag_mc();

   if ( purity <=0. ) continue;

   if ( m_flagHist == 0 ) continue;

   double dmass(0.), dm_dst(0.), dvchisq(0.), d_p(0.);

    if( dst_mode == 0 ){
       dmass   = dynamic_cast<UserInfo_B&>(b_i->child(0).userInfo()).mass();
       dvchisq = dynamic_cast<UserInfo_B&>(b_i->child(0).userInfo()).vtxchisq();
    }else{ 
       dm_dst  = dynamic_cast<UserInfo_B&>(b_i->child(0).userInfo()).mass() 
	         - b_i->child(0).child(0).mass();
       dmass   = dynamic_cast<UserInfo_B&>(b_i->child(0).child(0).userInfo()).mass();
       dvchisq = dynamic_cast<UserInfo_B&>(b_i->child(0).child(0).userInfo()).vtxchisq();
    }
    d_p = (frec_util::p_cm( b_i->child(0).p()) ).vect().mag();

    double dsmass(0.), dm_dsst(0.), dsvchisq(0.), ds_p(0.);
    double mphi_from_ds(0.), mkst_from_ds(0.);

    if ( abs(b_i->child(1).lund()) == Ptype_DSplus.lund() ||
	 abs(b_i->child(1).lund()) == Ptype_DSstplus.lund() ){

    if( dsst_mode == 0 ){
       dsmass   = dynamic_cast<UserInfo_B&>(b_i->child(1).userInfo()).mass();

       dsvchisq = dynamic_cast<UserInfo_B&>(b_i->child(1).userInfo()).vtxchisq();
       if ( ds_mode == 302 ){
	 mphi_from_ds = ( b_i->child(1).child(0).p() +
	                  b_i->child(1).child(1).p() ).mag();
       }else if ( ds_mode == 303 ){
	 mkst_from_ds = ( b_i->child(1).child(0).p() 
                          + b_i->child(1).child(1).p() ).mag();
       }
    }else{ 
       dm_dsst  = dynamic_cast<UserInfo_B&>(b_i->child(1).userInfo()).mass() 
	          - b_i->child(1).child(0).mass();
       dsmass   = dynamic_cast<UserInfo_B&>(b_i->child(1).child(0).userInfo()).mass();
       dsvchisq = dynamic_cast<UserInfo_B&>(b_i->child(1).child(0).userInfo()).vtxchisq();
       if ( ds_mode == 302 ){
	 mphi_from_ds = ( b_i->child(1).child(0).child(0).p() +
	                  b_i->child(1).child(0).child(1).p() ).mag();
       }else if ( ds_mode == 303 ){
           mkst_from_ds = ( b_i->child(1).child(0).child(0).p() 
                            + b_i->child(1).child(0).child(1).p() ).mag();
       }

    }

    ds_p = (frec_util::p_cm( b_i->child(1).p()) ).vect().mag();
      
    }

    double pi_p   = -1;
    double pi_pid = -1.;
    if( abs(b_i->child(1).lund()) == Ptype_PIplus.lund() ){
      pi_p = (frec_util::p_cm( b_i->child(1).p()) ).vect().mag();
      pi_pid = selkpi_noTOF.prob(b_i->child(1).mdstCharged());
     }

    double rhomass  = 0.0;
    double rho_p    = -1;
    if( abs(b_i->child(1).lund()) == Ptype_RHOplus.lund() ){
      rhomass = b_i->child(1).mass();
      rho_p = (frec_util::p_cm( b_i->child(1).p()) ).vect().mag();
      pi_pid = selkpi.prob(b_i->child(1).child(0).mdstCharged());
     }

    double a1mass   = 0.0;
    double a1_p     =  -1;
    double a1chisq  = -1.0;
    if( abs(b_i->child(1).lund()) == Ptype_A1plus.lund() ){

      a1mass   = b_i->child(1).mass();
      rhomass  = b_i->child(1).child(0).mass();
      a1chisq  = dynamic_cast<UserInfo_B&>(b_i->child(1).userInfo()).vtxchisq();

      a1_p     = (frec_util::p_cm( b_i->child(1).p() )).vect().mag();
      rho_p    = (frec_util::p_cm( b_i->child(1).child(0).p() )).vect().mag();

      pi_pid   = selkpi.prob(b_i->child(1).child(0).child(0).mdstCharged());
      double tmp_pid = selkpi.prob(b_i->child(1).child(0).child(1).mdstCharged());
      if ( tmp_pid > pi_pid ) pi_pid = tmp_pid;
      tmp_pid  = selkpi.prob(b_i->child(1).child(1).mdstCharged());
      if ( tmp_pid > pi_pid ) pi_pid = tmp_pid;

    }

    if ( mphi_from_ds > 0 ) rhomass = mphi_from_ds;
    if ( mkst_from_ds > 0 ) a1mass  = mkst_from_ds;

    // Fill Ntuple 
    nt_b->column("Mbc",      Mass_B);       // Mbc
    nt_b->column("deltaE",   Delta_E);      // Delta E
    nt_b->column("best",     flag_best);    // best flag
    nt_b->column("best2",    flag_best2);   // best flag(2)
    nt_b->column("best3",    flag_best3);   // best flag(3)
    nt_b->column("b_lund",   b_i->lund());  // lund 
    nt_b->column("b_mode",   b_mode);       // B decay ID
    nt_b->column("dst_mode", dst_mode);     // D* decay ID
    nt_b->column("d_mode",   d_mode);       // D decay ID
    nt_b->column("dst2mode", dsst_mode);    // Ds* decay ID
    nt_b->column("ds_mode",  ds_mode);      // Ds decay ID
    nt_b->column("pi_p",     pi_p );        // pi momentum
    nt_b->column("pi_pid",   pi_pid );      // pi probK
    nt_b->column("rhomass",  rhomass);      // rho mass ( phi mass for Ds case )
    nt_b->column("rho_p",    rho_p);        // rho momentum
    nt_b->column("a1mass",   a1mass);       // a1 mass ( K* mass for Ds case )
    nt_b->column("a1chisq",  a1chisq);      // a1 vetex fit chisq
    nt_b->column("a1_p",     a1_p);         // a1 momentum
    nt_b->column("dm_dst",   dm_dst);       // D* mass - D mass
    nt_b->column("dmass",    dmass);        // D mass
    nt_b->column("dm_dst2",  dm_dsst);      // Ds* mass - Ds mass
    nt_b->column("dsmass",   dsmass);       // Ds mass
    nt_b->column("dvchisq",  dvchisq );     // reduced chisq for D  vertex fit
    nt_b->column("dsvchisq", dsvchisq );    // reduced chisq for Ds vertex fit
    nt_b->column("d_p",      d_p );         // D(*)  momentum
    nt_b->column("ds_p",     ds_p );        // Ds(*) momentum
    nt_b->column("r2",       m_R2);         // R2
    nt_b->column("cosT",     cosT);         // thrust angle
    nt_b->column("cosB",     cosB);         // B flight direction angle
    nt_b->column("bnum",     bnum);         // B candidate number
    nt_b->column("mc",       flag_mc);      // flag for MC (1=on,0=off) 
    nt_b->column("chisq",    chisq);        // reduced chisq for best cand selection
    nt_b->column("purity",   purity);       // B purity
    nt_b->dumpData();

    
  } // end of loop b_i 


  return;

}

void frec_ana::FillHist_Monitor( std::vector<Particle>& B ){

  if ( m_flagHist == 0 ) return;

  const atc_pid selkpi(m_ACCq, m_TOFq, m_CDCq, KAON_CODE, PION_CODE);
  const atc_pid selkpi_noTOF(m_ACCq, -1, m_CDCq, KAON_CODE, PION_CODE);

  for( std::vector<Particle>::iterator b_i = B.begin();
       b_i != B.end(); b_i++ ) {

    const int flag_bestEach = dynamic_cast<UserInfo_B&>(b_i->userInfo()).flag_bestEach();
    if ( flag_bestEach == 0 ) continue;

    const int b_mode   = dynamic_cast<UserInfo_B&>(b_i->userInfo()).b_mode();
    const int dst_mode = dynamic_cast<UserInfo_B&>(b_i->userInfo()).sub0_mode();
    const int d_mode   = dynamic_cast<UserInfo_B&>(b_i->userInfo()).sub1_mode();
    if ( b_mode == 0 || d_mode == 0 ) continue;

    const double Mass_B  = dynamic_cast<UserInfo_B&>(b_i->userInfo()).Mass_B();
    const double Delta_E = dynamic_cast<UserInfo_B&>(b_i->userInfo()).Delta_E();
    const double cosT    = dynamic_cast<UserInfo_B&>(b_i->userInfo()).cos_thetaT();
    const double cosB    = dynamic_cast<UserInfo_B&>(b_i->userInfo()).cos_thetaB();
    const double chisq   = dynamic_cast<UserInfo_B&>(b_i->userInfo()).masschisq();
    const int    flag_mc = dynamic_cast<UserInfo_B&>(b_i->userInfo()).flag_mc();

    double dmass(-1.), dm_dst(-1.), dvchisq(-1.);
    double ms_dst(-1.), ps_dst(-1.);
    double mpi0_d(-1.), ppi0_d(-1.);
    double mks_d(-1.), pks_d(-1.);
    double kd_pid(-1.), pid_pid(-1.), tmp_pid(-1.);

    if( dst_mode == 0 ){

       dmass   = dynamic_cast<UserInfo_B&>(b_i->child(0).userInfo()).mass();
       dvchisq = dynamic_cast<UserInfo_B&>(b_i->child(0).userInfo()).vtxchisq();

       for ( unsigned i=0;i<b_i->child(0).nChildren();i++){
	 if ( abs(b_i->child(0).child(i).lund()) == Ptype_PI0.lund() ){
	   mpi0_d = b_i->child(0).child(i).mdstPi0().mass();
	   ppi0_d = (frec_util::p_cm( b_i->child(0).child(i).p())).vect().mag();
	 }
	 if ( abs(b_i->child(0).child(i).lund()) == Ptype_Ks.lund() ){
	   mks_d = b_i->child(0).child(i).p().mag();
	   pks_d = (frec_util::p_cm( b_i->child(0).child(i).p())).vect().mag();
	 }
	 if ( abs(b_i->child(0).child(i).lund()) == Ptype_Kplus.lund() ){
	   tmp_pid = selkpi.prob(b_i->child(0).child(i).mdstCharged());
	   if ( kd_pid < 0 || tmp_pid < kd_pid ) kd_pid = tmp_pid;
	 }
	 if ( abs(b_i->child(0).child(i).lund()) == Ptype_PIplus.lund() ){
	   tmp_pid = selkpi.prob(b_i->child(0).child(i).mdstCharged());
	   if ( tmp_pid > pid_pid ) pid_pid = tmp_pid;
	 }
       }

    }else{ 
       dm_dst   = dynamic_cast<UserInfo_B&>(b_i->child(0).userInfo()).mass() 
	          - b_i->child(0).child(0).mass();
       dmass    = dynamic_cast<UserInfo_B&>(b_i->child(0).child(0).userInfo()).mass();
       dvchisq  = dynamic_cast<UserInfo_B&>(b_i->child(0).child(0).userInfo()).vtxchisq();
       if ( b_i->child(0).child(1).mdstPi0() ){
	 ms_dst = b_i->child(0).child(1).mdstPi0().mass();
       }else{
	 ms_dst = b_i->child(0).child(1).p().mag();
       }
       ps_dst   = b_i->child(0).child(1).p().vect().mag();

       for ( unsigned i=0;i<b_i->child(0).child(0).nChildren();i++){
	 if ( abs(b_i->child(0).child(0).child(i).lund()) == Ptype_PI0.lund() ){
	   mpi0_d = b_i->child(0).child(0).child(i).mdstPi0().mass();
	   ppi0_d = (frec_util::p_cm( b_i->child(0).child(0).child(i).p())).vect().mag();
	 }
	 if ( abs(b_i->child(0).child(0).child(i).lund()) == Ptype_Ks.lund() ){
	   mks_d = b_i->child(0).child(0).child(i).p().mag();
	   pks_d = (frec_util::p_cm( b_i->child(0).child(0).child(i).p())).vect().mag();
	 }
	 if ( abs(b_i->child(0).child(0).child(i).lund()) == Ptype_Kplus.lund() ){
	   tmp_pid = selkpi.prob(b_i->child(0).child(0).child(i).mdstCharged());
	   if ( kd_pid < 0 || tmp_pid < kd_pid ) kd_pid = tmp_pid;
	 }
	 if ( abs(b_i->child(0).child(0).child(i).lund()) == Ptype_PIplus.lund() ){
	   tmp_pid = selkpi.prob(b_i->child(0).child(0).child(i).mdstCharged());
	   if ( tmp_pid > pid_pid ) pid_pid = tmp_pid;
	 }
       }
    }

    bool flag_dpi(false);
    bool flag_drho(false);
    bool flag_da1(false);

    int lundX = abs(b_i->child(1).pType().lund());

    if ( lundX == Ptype_PIplus.lund() )       flag_dpi  = true;
    else if ( lundX == Ptype_RHOplus.lund() ) flag_drho = true;
    else if ( lundX == Ptype_A1plus.lund() )  flag_da1  = true;

    if ( dm_dst > 0.15 ){
      flag_drho = false; flag_da1 = false;
    }

    double d_p = (frec_util::p_cm( b_i->child(0).p())).vect().mag();

    double pi_pid;

    if ( flag_dpi ){

      double pi_p   = (frec_util::p_cm( b_i->child(1).p())).vect().mag();
      pi_pid = selkpi_noTOF.prob(b_i->child(1).mdstCharged());

      double mkpi_2(0.), mkpi0_2(0.), dalitz(0.);
      if ( d_mode == 2 ){
	if (dst_mode == 0 ){
	  HepLorentzVector p_k   = b_i->child(0).child(0).p();
	  HepLorentzVector p_pi  = b_i->child(0).child(1).p();
	  HepLorentzVector p_pi0 = b_i->child(0).child(2).p();
	  mkpi_2  = (p_k + p_pi).mag2();
	  mkpi0_2 = (p_k + p_pi0 ).mag2();
	  float f_dmass(dmass), f_mkpi_2(mkpi_2), f_mkpi0_2(mkpi0_2);
	  dalitz = weight_(&f_dmass, &f_mkpi_2, &f_mkpi0_2);
	}else{
	  HepLorentzVector p_k   = b_i->child(0).child(0).child(0).p();
	  HepLorentzVector p_pi  = b_i->child(0).child(0).child(1).p();
	  HepLorentzVector p_pi0 = b_i->child(0).child(0).child(2).p();
	  mkpi_2  = (p_k + p_pi).mag2();
	  mkpi0_2 = (p_k + p_pi0 ).mag2();
	  float f_dmass(dmass), f_mkpi_2(mkpi_2), f_mkpi0_2(mkpi0_2);
	  dalitz = weight_(&f_dmass, &f_mkpi_2, &f_mkpi0_2);
	}
      }

      nt_dpi->column("Mbc",      Mass_B);          // Mbc
      nt_dpi->column("deltaE",   Delta_E);         // Delta E
      nt_dpi->column("b_lund",   b_i->lund());     // lund 
      nt_dpi->column("b_mode",   b_mode);          // B decay ID
      nt_dpi->column("dst_mode", dst_mode);        // D* decay ID
      nt_dpi->column("d_mode",   d_mode);          // D decay ID
      nt_dpi->column("dm_dst",   dm_dst);          // D* mass - D mass
      nt_dpi->column("ms_dst",   ms_dst );         // mass, p for slow pi, gamma from D*
      nt_dpi->column("ps_dst",   ps_dst );
      nt_dpi->column("pi_p",     pi_p);            // pi momentum from B
      nt_dpi->column("pi_pid",   pi_pid );         // pi probK ( from B, w/o TOF )
      nt_dpi->column("dmass",    dmass);           // D mass
      nt_dpi->column("dvchisq",  dvchisq );        // reduced chisq for D vertex fit
      nt_dpi->column("d_p",      d_p);             // D(*) momentum
      nt_dpi->column("kd_pid",   kd_pid);          // K's  minimum probK from D
      nt_dpi->column("pid_pid",  pid_pid);         // pi's maximum probK from D
      nt_dpi->column("mpi0_d",   mpi0_d );         // mass, p for pi0 from D
      nt_dpi->column("ppi0_d",   ppi0_d );          
      nt_dpi->column("mks_d",    mks_d );          // mass, p for Ks from D
      nt_dpi->column("pks_d",    pks_d );          //  ( p : lab. frame )
      nt_dpi->column("mkpi_d",   sqrt(mkpi_2)  );  // for dalitz (Kpipi0)
      nt_dpi->column("mkpi0_d",  sqrt(mkpi0_2) );
      nt_dpi->column("dalitz_d", dalitz );
      nt_dpi->column("r2",       m_R2);            // R2
      nt_dpi->column("cosT",     cosT);            // thrust angle
      nt_dpi->column("cosB",     cosB);            // B flight direction angle
      nt_dpi->column("mc",       flag_mc);         // flag for MC (1=on,0=off) 
      nt_dpi->column("chisq",    chisq);           // reduced chisq for best cand selection
      nt_dpi->dumpData();
      
    }

    if ( flag_drho ){

      const double rhomass  = b_i->child(1).mass();
      const double rho_p    = (frec_util::p_cm( b_i->child(1).p())).vect().mag();

      const double rpi_p    = (frec_util::p_cm( b_i->child(1).child(0).p())).vect().mag();
      const double rpi0_p   = (frec_util::p_cm( b_i->child(1).child(1).p())).vect().mag();
      const double rpi0_m   = b_i->child(1).child(1).mdstPi0().mass();

      pi_pid = selkpi.prob(b_i->child(1).child(0).mdstCharged());

      nt_drho->column("Mbc",      Mass_B);       // Mbc
      nt_drho->column("deltaE",   Delta_E);      // Delta E
      nt_drho->column("b_lund",   b_i->lund());  // lund 
      nt_drho->column("b_mode",   b_mode);       // B decay ID
      nt_drho->column("dst_mode", dst_mode);     // D* decay ID
      nt_drho->column("d_mode",   d_mode);       // D decay ID
      nt_drho->column("rhomass",  rhomass);      // rho mass
      nt_drho->column("rho_p",    rho_p);        // rho momentum
      nt_drho->column("rpi_p",    rpi_p);        // pi  momentum from rho
      nt_drho->column("rpi0_p",   rpi0_p);       // pi0 momentum from rho
      nt_drho->column("rpi0_m",   rpi0_m);       // pi0 mass     from rho
      nt_drho->column("pi_pid",   pi_pid);       // pi probK ( from B )
      nt_drho->column("dm_dst",   dm_dst);       // D* mass - D mass
      nt_drho->column("dmass",    dmass);        // D mass
      nt_drho->column("dvchisq",  dvchisq );     // reduced chisq for D vertex fit
      nt_drho->column("d_p",      d_p);          // D(*) momentum
      nt_drho->column("kd_pid",   kd_pid);       // K's  minimum probK from D
      nt_drho->column("pid_pid",  pid_pid);      // pi's maximum probK from D
      nt_drho->column("r2",       m_R2);         // R2
      nt_drho->column("cosT",     cosT);         // thrust angle
      nt_drho->column("cosB",     cosB);         // B flight direction angle
      nt_drho->column("mc",       flag_mc);      // flag for MC (1=on,0=off) 
      nt_drho->column("chisq",    chisq);        // reduced chisq for best cand selection
      nt_drho->dumpData();

    }

    if ( flag_da1 ){

      const double a1mass  = b_i->child(1).mass();
      const double a1_p    = (frec_util::p_cm( b_i->child(1).p())).vect().mag();

      const double rhomass = b_i->child(1).child(0).mass();
      const double rho_p   = (frec_util::p_cm( b_i->child(1).child(0).p())).vect().mag();

      const double rpi1_p  = (frec_util::p_cm( b_i->child(1).child(0).child(0).p())).vect().mag();
      const double rpi2_p  = (frec_util::p_cm( b_i->child(1).child(0).child(0).p())).vect().mag();
      const double api_p   = (frec_util::p_cm( b_i->child(1).child(1).p())).vect().mag();

      const double a1chisq = dynamic_cast<UserInfo_B&>(b_i->child(1).userInfo()).vtxchisq();

      pi_pid = selkpi.prob(b_i->child(1).child(0).child(0).mdstCharged());
      double tmp_pid = selkpi.prob(b_i->child(1).child(0).child(1).mdstCharged());
      if ( tmp_pid > pi_pid ) pi_pid = tmp_pid;
      tmp_pid = selkpi.prob(b_i->child(1).child(1).mdstCharged());
      if ( tmp_pid > pi_pid ) pi_pid = tmp_pid;

      nt_da1->column("Mbc",      Mass_B);       // Mbc
      nt_da1->column("deltaE",   Delta_E);      // Delta E
      nt_da1->column("b_lund",   b_i->lund());  // lund 
      nt_da1->column("b_mode",   b_mode);       // B decay ID
      nt_da1->column("dst_mode", dst_mode);     // D* decay ID
      nt_da1->column("d_mode",   d_mode);       // D decay ID
      nt_da1->column("a1mass",   a1mass);       // a1 mass
      nt_da1->column("a1_p",     a1_p);         // a1 momentum
      nt_da1->column("rhomass",  rhomass);      // rho mass
      nt_da1->column("rho_p",    rho_p);        // rho momentum
      nt_da1->column("rpi1_p",   rpi1_p);       // pi momentum from rho0
      nt_da1->column("rpi2_p",   rpi2_p);       
      nt_da1->column("api_p",    api_p);        // pi momentum from a1
      nt_da1->column("a1chisq",  a1chisq);      // reduced chisq for a1 vetex fit
      nt_da1->column("pi_pid",   pi_pid);       // pi's maximum probK ( from B )
      nt_da1->column("dm_dst",   dm_dst);       // D* mass - D mass
      nt_da1->column("dmass",    dmass);        // D mass
      nt_da1->column("dvchisq",  dvchisq );     // reduced chisq for D vertex fit
      nt_da1->column("d_p",      d_p);          // D(*) momentum
      nt_da1->column("kd_pid",   kd_pid);       // K's  minimum probK from D
      nt_da1->column("pid_pid",  pid_pid);      // pi's maximum probK from D
      nt_da1->column("r2",       m_R2);         // R2
      nt_da1->column("cosT",     cosT);         // thrust angle
      nt_da1->column("cosB",     cosB);         // B flight direction angle
      nt_da1->column("mc",       flag_mc);      // flag for MC (1=on, 0=off) 
      nt_da1->column("chisq",    chisq);        // reduced chisq for best cand selection
      nt_da1->dumpData();

    }

  } // end of loop b_i 

  return;

}

// -----------------------
//    count B
// -----------------------
void frec_ana::count_B( std::vector<Particle> &B_pm,
			std::vector<Particle> &B_zero ){


  int flag_b(0), flag_b2dpi(0), flag_b2drho(0), flag_b2da1(0), flag_b2dds(0);
  int flag_b_mc(0), flag_b2dpi_mc(0), flag_b2drho_mc(0), flag_b2da1_mc(0), flag_b2dds_mc(0);

  for( std::vector<Particle>::iterator i = B_pm.begin();
       i != B_pm.end(); i++){

    int lund_meson = abs(i->relation().child(1).lund());
    int flag_mc    = dynamic_cast<UserInfo_B&>(i->userInfo()).flag_mc();

    flag_b = 1;  if ( flag_mc ) flag_b_mc = 1;
    if ( lund_meson == Ptype_PIplus.lund() ) {
      flag_b2dpi  = 1; if ( flag_mc ) flag_b2dpi_mc = 1;
    }else if ( lund_meson == Ptype_RHOplus.lund() ){
      flag_b2drho = 1; if ( flag_mc ) flag_b2drho_mc = 1;
    }else if ( lund_meson == Ptype_A1plus.lund()  ){
      flag_b2da1  = 1; if ( flag_mc ) flag_b2da1_mc = 1;
    }else if ( lund_meson == Ptype_DSplus.lund() ||
              lund_meson == Ptype_DSstplus.lund() ){
      flag_b2dds  = 1; if ( flag_mc ) flag_b2dds_mc = 1;
    }
  }

  for( std::vector<Particle>::iterator i = B_zero.begin();
       i != B_zero.end(); i++){

    int lund_meson = abs(i->relation().child(1).lund());
    int flag_mc    = dynamic_cast<UserInfo_B&>(i->userInfo()).flag_mc();

    flag_b = 1;  if ( flag_mc ) flag_b_mc = 1;
    if ( lund_meson == Ptype_PIplus.lund() ){ 
      flag_b2dpi  = 1; if ( flag_mc ) flag_b2dpi_mc = 1;
    }else if ( lund_meson == Ptype_RHOplus.lund() ){
      flag_b2drho = 1; if ( flag_mc ) flag_b2drho_mc = 1;
    }else if ( lund_meson == Ptype_A1plus.lund()  ){
      flag_b2da1  = 1; if ( flag_mc ) flag_b2da1_mc = 1;
    }else if ( lund_meson == Ptype_DSplus.lund() ||
              lund_meson == Ptype_DSstplus.lund() ){
      flag_b2dds  = 1; if ( flag_mc ) flag_b2dds_mc = 1;
    }

  }

    // --- count B events
    if ( flag_b )         kid_statistics::fill_sta( m_name, m_cid_b, (float)1. );
    if ( flag_b_mc )      kid_statistics::fill_sta( m_name, m_cid_b_mc, (float)1. );
    // --- count B(D pi)
    if ( flag_b2dpi )     kid_statistics::fill_sta( m_name, m_cid_b2dpi, (float)1. );
    if ( flag_b2dpi_mc )  kid_statistics::fill_sta( m_name, m_cid_b2dpi_mc, (float)1. );
    // --- count B(D rho)
    if ( flag_b2drho )    kid_statistics::fill_sta( m_name, m_cid_b2drho, (float)1. );
    if ( flag_b2drho_mc ) kid_statistics::fill_sta( m_name, m_cid_b2drho_mc, (float)1. );
    // --- count B(D a1)
    if ( flag_b2da1 )     kid_statistics::fill_sta( m_name, m_cid_b2da1, (float)1. );  
    if ( flag_b2da1_mc )  kid_statistics::fill_sta( m_name, m_cid_b2da1_mc, (float)1. );  
    // --- count B(D Ds)
    if ( flag_b2dds )     kid_statistics::fill_sta( m_name, m_cid_b2dds, (float)1. );
    if ( flag_b2dds_mc )  kid_statistics::fill_sta( m_name, m_cid_b2dds_mc, (float)1. );  

}

// -----------------------------------------------------
//   save B candidates thorough brecon/fullrecon table
// -----------------------------------------------------
int frec_ana::save_Brecon(std::vector<Particle>& B_zero, std::vector<Particle>& B_pm  ){

  int num_B = B_zero.size() + B_pm.size();
  if ( num_B == 0 ) return num_B;

  Brecon_header_Manager & bhMgr = Brecon_header_Manager::get_manager();
  Brecon_Manager        & brMgr = Brecon_Manager::get_manager();
  Fullrecon_Manager     & frMgr = Fullrecon_Manager::get_manager();
  bhMgr.remove();
  brMgr.remove();
  frMgr.remove();

  std::vector<Particle> B_sort;
  for( std::vector<Particle>::iterator b_i = B_zero.begin();
       b_i != B_zero.end(); b_i++ ) B_sort.push_back(*b_i);
  
  for( std::vector<Particle>::iterator b_i = B_pm.begin();
       b_i != B_pm.end(); b_i++ ) B_sort.push_back(*b_i);

   std::sort( B_sort.begin(), B_sort.end(), b_lowerChisq );
   std::sort( B_sort.begin(), B_sort.end(), b_higherPurity );
   std::sort( B_sort.begin(), B_sort.end(), b_withBest );

  AnaBrecon brecon;

  int idlast(0);
  for( std::vector<Particle>::iterator b_i = B_sort.begin();
       b_i != B_sort.end(); b_i++ ) {

    int b_mode   = dynamic_cast<UserInfo_B&>(b_i->userInfo()).b_mode();
    if ( b_mode == 0 ) continue ;
    int dst_mode  = dynamic_cast<UserInfo_B&>(b_i->userInfo()).sub0_mode();
    int d_mode    = dynamic_cast<UserInfo_B&>(b_i->userInfo()).sub1_mode();
    int dsst_mode = dynamic_cast<UserInfo_B&>(b_i->userInfo()).sub2_mode();
    int ds_mode   = dynamic_cast<UserInfo_B&>(b_i->userInfo()).sub3_mode();
    if ( d_mode == 0 ) continue;

    double Mass_B     = dynamic_cast<UserInfo_B&>(b_i->userInfo()).Mass_B();
    double Delta_E    = dynamic_cast<UserInfo_B&>(b_i->userInfo()).Delta_E();

    double purity     = dynamic_cast<UserInfo_B&>(b_i->userInfo()).purity();
    double rchisq     = dynamic_cast<UserInfo_B&>(b_i->userInfo()).masschisq();
    double cosT       = dynamic_cast<UserInfo_B&>(b_i->userInfo()).cos_thetaT();
    double cosB       = dynamic_cast<UserInfo_B&>(b_i->userInfo()).cos_thetaB();

    int flag_bestEach = dynamic_cast<UserInfo_B&>(b_i->userInfo()).flag_bestEach();
    int flag_best     = dynamic_cast<UserInfo_B&>(b_i->userInfo()).flag_best();
    int flag_best3    = dynamic_cast<UserInfo_B&>(b_i->userInfo()).flag_best3();
    int flag_mc       = dynamic_cast<UserInfo_B&>(b_i->userInfo()).flag_mc();

    if ( purity <=0. ) continue;
    //    if ( flag_bestEach == 0 ) continue;

    brecon.setBrecon( *b_i );

    int id(0);
    for ( unsigned i=idlast;i<brMgr.count();i++){
      if ( brMgr[i].mother() == 0 ){
	id = i + 1; break;
       }
    }

    if ( id <= 0 ) continue;

    idlast = id - 1;
    Fullrecon& ftable = frMgr.add();
    ftable.FullBrec( brMgr[idlast] );
    idlast += 1;

    double dvchisq(0.);
    if( dst_mode == 0 ){
       dvchisq = dynamic_cast<UserInfo_B&>(b_i->child(0).userInfo()).vtxchisq();
    }else{ 
       dvchisq = dynamic_cast<UserInfo_B&>(b_i->child(0).child(0).userInfo()).vtxchisq();
    }
    double dsvchisq(0.);
    if ( abs(b_i->child(1).lund()) == Ptype_DSplus.lund() ||
	 abs(b_i->child(1).lund()) == Ptype_DSstplus.lund() ){
      if( dsst_mode == 0 ){
	dsvchisq = dynamic_cast<UserInfo_B&>(b_i->child(1).userInfo()).vtxchisq();
      }else{ 
	dsvchisq = dynamic_cast<UserInfo_B&>(b_i->child(1).child(0).userInfo()).vtxchisq();
      }
    }
    double a1chisq(0.);
    if( abs(b_i->child(1).lund()) == Ptype_A1plus.lund() ){
       a1chisq = dynamic_cast<UserInfo_B&>(b_i->child(1).userInfo()).vtxchisq();
    }

    ftable.kind(1);
    ftable.b_mode(b_mode);
    ftable.sub0_mode(dst_mode);
    ftable.sub1_mode(d_mode);

    ftable.deltE(Delta_E);
    ftable.Mbc(Mass_B);

    ftable.Aux0(purity);
    ftable.Aux1(rchisq);
    ftable.Aux2(cosT);
    ftable.Aux3(m_R2);
    ftable.Aux4(cosB);
    ftable.Aux5(dvchisq);
    ftable.Aux6(dsvchisq);
    ftable.Aux7(a1chisq);
    ftable.Aux8(0.);
    ftable.Aux9(0.);
    ftable.flag0(flag_best); 
    ftable.flag1(flag_best3); 
    ftable.flag2(flag_bestEach);
    ftable.flag3(b_i->lund());
    ftable.flag4(flag_mc);

  }

  return num_B;

}

// ---------------------------
//   Other functions
// ---------------------------
int frec_ana::get_b_mode( Particle& b ){

  int b_mode(0);

  const Particle &d = b.child(0);
  const Particle &x = b.child(1);
  const int d_lund  = abs(d.lund());
  const int x_lund  = abs(x.lund());


  if ( d_lund == Ptype_Dplus.lund() ){

    if ( x_lund == Ptype_PIplus.lund() )        b_mode =   1;
    else if ( x_lund == Ptype_RHOplus.lund()  ) b_mode =   3;
    else if ( x_lund == Ptype_A1plus.lund()   ) b_mode =   5;
    else if ( x_lund == Ptype_DSplus.lund()   ) b_mode =  11;
    else if ( x_lund == Ptype_DSstplus.lund() ) b_mode =  12;

  }else if ( d_lund == Ptype_Dstplus.lund() ){

    if ( x_lund == Ptype_PIplus.lund() )         b_mode =   2;
    else if ( x_lund == Ptype_RHOplus.lund()  )  b_mode =   4;
    else if ( x_lund == Ptype_A1plus.lund()   )  b_mode =   6;
    else if ( x_lund == Ptype_DSplus.lund()   )  b_mode =  13;
    else if ( x_lund == Ptype_DSstplus.lund() )  b_mode =  14;

  }else if ( d_lund == Ptype_D0.lund()  ) {

    if ( x_lund == Ptype_PIplus.lund() )         b_mode = 101;
    else if ( x_lund == Ptype_RHOplus.lund() )   b_mode = 103;
    else if ( x_lund == Ptype_A1plus.lund()  )   b_mode = 105;
    else if ( x_lund == Ptype_DSplus.lund()   )  b_mode = 111;
    else if ( x_lund == Ptype_DSstplus.lund() )  b_mode = 112;

  }else if ( d_lund == Ptype_Dst0.lund() ){

    if ( x_lund == Ptype_PIplus.lund() )         b_mode = 102;
    else if ( x_lund == Ptype_RHOplus.lund()  )  b_mode = 104;
    else if ( x_lund == Ptype_A1plus.lund()   )  b_mode = 106;
    else if ( x_lund == Ptype_DSplus.lund()   )  b_mode = 113;
    else if ( x_lund == Ptype_DSstplus.lund() )  b_mode = 114;

  }

 return b_mode;

}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
