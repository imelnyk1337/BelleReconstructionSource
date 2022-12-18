//
// $Id: funda_recon.cc 9932 2006-11-12 14:26:53Z katayama $
//
// $Log$
// Revision 1.8  2004/11/22 05:24:36  matumot
// updates for new kid_statistics
// module parameter, "LogFileName" is introduced
//
// Revision 1.7  2004/05/28 10:05:55  matumot
// adjusted pi0 mass cut to 0.1178 < Mpi0 < 0.1502 GeV/c^2
// ( consistent with pi0 mass cuts in current fix_mdst module )
//
// Revision 1.6  2004/05/04 21:23:20  matumot
// fixed for CC complier
// - namespace was changed : fullrecon -> fullreconNS
//   to avoid overlap for fullrecon struct used by fullrecon.tdf
// - use new operator for dynamic allocation of array
// - add std:: for sort function
//
// Revision 1.4  2004/04/19 09:27:34  matumot
// revival version
//
// =======================================================
// File Name : funda_recon.cc 
// -------------------------------------------------------
// Creation    ; 2004.04.07
// Description ; Implimentation for the frec_ana class
//               ( Fundamental particles )
// Author      ; T.Matsumoto (matumot@bmail.kek.jp)
// -------------------------------------------------------

#include "belle.h"
#include  "fullrecon/frec_ana.h"
#include  "fullrecon/frec_util.h"

#include "ip/IpProfile.h"

#include "particle/Particle.h"
#include "particle/utility.h"

#include "kfitter/kvertexfitter.h"
#include "kfitter/kmassfitter.h"

#include "kid/atc_pid.h"
#include "kid/kid_statistics.h"

#include "mdst/findKs.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


using namespace fullreconNS;

// ------------------------------
// Charged particles
// ------------------------------
void  frec_ana::recon_Charged( std::vector<Particle>& Charged_all, 
			       std::vector<Particle>& Pi_all, std::vector<Particle>& PiLoose_all,
			       std::vector<Particle>& K_all ){

  Ptype ptype;
  const atc_pid selkpi(m_ACCq, m_TOFq, m_CDCq, KAON_CODE, PION_CODE);

  Mdst_charged_Manager&  ChgMgr = Mdst_charged_Manager::get_manager();

  for(std::vector<Mdst_charged>::const_iterator i = ChgMgr.begin();
      i != ChgMgr.end(); i++) {

    bool select      = false;
    bool selectLoose = false;
    bool select_mc   = false;

    const double dr = frec_util::correct_dr(*i, IpProfile::e_position(), PION_CODE);
    const double dz = frec_util::correct_dz(*i, IpProfile::e_position(), PION_CODE);
    const Hep3Vector ch_p(i->px(), i->py(), i->pz());

   if ( m_flagHist ){
    H[12]->accumulate( dr, 1. );
    H[13]->accumulate( dz, 1. );
    H[14]->accumulate( ch_p.perp(), 1. );
   }

    if( fabs(dr) > fabs(m_Dr_cut)  || 
        fabs(dz) > fabs(m_Dz_cut)  ||   
	ch_p.perp() < m_Pt_cut ) continue;

  if ( m_flagHist ){
    H[15]->accumulate( i->charge(), 1. );
    H[16]->accumulate( selkpi.prob(*i), 1. );
  }

    if ( i->charge() > 0 )   ptype = Ptype_PIplus; 
    else                     ptype = Ptype_PIminus;
       
    Particle pi( *i, ptype ); 
    if ( m_mc ) {
      frec_util::set_pGenHepevt( pi, m_mc );
      if ( pi.genHepevt() ) select_mc = true;
    }

    bool piFlag = true;
    bool kFlag  = true;

    if ( m_dEdxSel ){
      piFlag = frec_util::with_dEdx( pi, PION_CODE, m_dEdxSigma_cut );
      kFlag  = frec_util::with_dEdx( pi, KAON_CODE, m_dEdxSigma_cut );
    }

    select      = frec_util::sel_KPI( pi, m_PID_cut_PI, selkpi );
    selectLoose = frec_util::sel_KPI( pi, m_PID_cut_PI_loose, selkpi );

    Charged_all.push_back(pi);
    if ( select && piFlag ) Pi_all.push_back(pi);
    if ( selectLoose && (piFlag || kFlag) ) PiLoose_all.push_back(pi);

    if ( m_flagHist && select ){
	const double p_pi = (frec_util::p_cm(pi.p())).vect().mag();
	H[31]->accumulate(p_pi, 1.);
	if ( select_mc ) H[41]->accumulate(p_pi, 1.);
    }

  }// end of loop i


  for(std::vector<Mdst_charged>::const_iterator i = ChgMgr.begin();
      i != ChgMgr.end(); i++) {

    bool select      = false;
    bool select_mc   = false;

    const double dr = frec_util::correct_dr(*i, IpProfile::e_position(), KAON_CODE);
    const double dz = frec_util::correct_dz(*i, IpProfile::e_position(), KAON_CODE);
    const Hep3Vector ch_p(i->px(), i->py(), i->pz());

    if( fabs(dr) > fabs(m_Dr_cut)  || 
        fabs(dz) > fabs(m_Dz_cut)  ||   
	ch_p.perp() < m_Pt_cut ) continue;

    if ( i->charge() > 0 )   ptype = Ptype_Kplus;
    else                     ptype = Ptype_Kminus;

    Particle k( *i, ptype );
    if ( m_mc ) {
      frec_util::set_pGenHepevt( k, m_mc );
      if ( k.genHepevt() ) select_mc = true;
    }

    if ( m_dEdxSel ){
      bool kFlag  = frec_util::with_dEdx( k , KAON_CODE, m_dEdxSigma_cut );
      if ( !kFlag ) continue;
    }

    select  = frec_util::sel_KPI( k,  m_PID_cut_K_loose,  selkpi );

    if ( select  ) {
      K_all.push_back(k);
      if ( m_flagHist ){
	const double p_k  = (frec_util::p_cm(k.p())).vect().mag();
	H[32]->accumulate(p_k, 1.);
	if ( select_mc ) H[42]->accumulate(p_k, 1.);
      }
    }

  }// end of loop i

  if ( m_flagHist ){
    H[21]->accumulate( Charged_all.size(), 1. );
    H[22]->accumulate( Pi_all.size(), 1. );
    H[23]->accumulate( K_all.size(), 1. );
  }

  return;

} 

// ------------------------------
//   Gamma reconstruction
// ------------------------------
void  frec_ana::recon_Gam( std::vector<Particle>& Gamma ){

  Mdst_gamma_Manager& GamMgr = Mdst_gamma_Manager::get_manager();

  for( std::vector<Mdst_gamma>::const_iterator i = GamMgr.begin();
       i != GamMgr.end(); i++ ) {

    Mdst_ecl& ecl_gam = i->ecl();
    const Hep3Vector gam_p(i->px(), i->py(), i->pz());

    if( ecl_gam.match() != 0 || ecl_gam.quality() != 0 ) continue;

    if ( m_flagHist ) H[33]->accumulate( ecl_gam.energy(), 1. );

    Particle gam(*i);
    setGammaError(gam,IpProfile::e_position(),IpProfile::e_position_err());

    frec_util::set_pGenHepevt( gam, m_mc );

    if ( m_flagHist && m_mc && gam.genHepevt() ) H[43]->accumulate( ecl_gam.energy(), 1. );


    if ( ecl_gam.energy() < m_Ene_cut ) continue;
    Gamma.push_back(gam);

  }

  if ( m_flagHist )  H[24]->accumulate(Gamma.size(), 1.);

  return;

}

// ------------------------------
//   PI0  reconstruction
// ------------------------------
void  frec_ana::recon_Pi0( std::vector<Particle>& Pi_zero, 
			   std::vector<Particle>& Slow_Pi_zero ) {

  Mdst_pi0_Manager&  pi0Mgr = Mdst_pi0_Manager::get_manager();

  for( std::vector<Mdst_pi0>::const_iterator i = pi0Mgr.begin();
       i != pi0Mgr.end(); i++ ) {

    bool select = false;
    bool select_mc = false;

    const double mass_pi0 = (*i).mass();
    if ( m_cutPi0 &&
	 ( mass_pi0 < m_Mdst_pi0_lower_limit || mass_pi0 > m_Mdst_pi0_upper_limit ) ) continue;

    const Hep3Vector gam_p1(i->gamma(0).px(), i->gamma(0).py(), i->gamma(0).pz());
    const Hep3Vector gam_p2(i->gamma(1).px(), i->gamma(1).py(), i->gamma(1).pz());
    const double e_gam1 = gam_p1.mag();
    const double e_gam2 = gam_p2.mag();

    Particle pi0(*i);
    frec_util::set_pGenHepevt( pi0, m_mc );
    if ( m_mc && pi0.genHepevt() ) select_mc = true;

    if ( m_flagHist ){
        H[34]->accumulate( e_gam1, 1. );
        H[34]->accumulate( e_gam2, 1. );
     if ( select_mc == true ){
        H[44]->accumulate( e_gam1, 1. );
        H[44]->accumulate( e_gam2, 1. );
      }
    }

    if( e_gam1 < m_Ene_cut_slowPI0 || 
        e_gam2 < m_Ene_cut_slowPI0 ) 	continue;

      // Mass constrained fit
      setGammasError(pi0, IpProfile::e_position(), IpProfile::e_position_err());
      doMassFit(pi0);

      Slow_Pi_zero.push_back(pi0);

      if( e_gam1 < m_Ene_cut || e_gam2 < m_Ene_cut ) continue;

      select = true;
      Pi_zero.push_back(pi0);

     if ( m_flagHist && select == true ){
	   double p_pi0 = (frec_util::p_cm(pi0.p())).vect().mag();
 	   H[35]->accumulate( p_pi0, 1. );
	   H[51]->accumulate( pi0.mdstPi0().mass(), 1. );
	 if ( select_mc == true ) {
	    H[45]->accumulate( p_pi0, 1. );
	    H[61]->accumulate( pi0.mdstPi0().mass(), 1. );
	 }
       }

      if ( select == true  ) {
          kid_statistics::fill_sta( m_name, m_cid_pi0, (float)1.);
      if ( select_mc  == true ) kid_statistics::fill_sta( m_name, m_cid_pi0_mc, (float)1.);
      }

   }


  if ( m_flagHist ) H[25]->accumulate( Pi_zero.size(), 1. );

  return;

}

// ------------------------------
//   Ks  reconstruction
// ------------------------------
void  frec_ana::recon_Ks( std::vector<Particle>& K_short ) {

  Mdst_vee2_Manager&  VeeMgr = Mdst_vee2_Manager::get_manager();

  for( std::vector<Mdst_vee2>::const_iterator ks_i = VeeMgr.begin();
       ks_i != VeeMgr.end(); ks_i++ ) {

    if ( ks_i->kind() != 1 ) continue;

    const HepLorentzVector vee( ks_i->px(), ks_i->py(), ks_i->pz(), ks_i->energy());
    if(ks_i->chisq()>=1.0e10)continue;

    const double ks_mass = vee.mag();

    if ( fabs(ks_mass - Ptype_Ks.mass()) > m_Ks_cut ) continue;

    bool preselect = false;
    bool select_mc = false;
    bool select    = false;

   if( m_Ks_fanfan_cut == 0 ){

     preselect = true;

   }else if( m_Ks_fanfan_cut == 1 ){

     FindKs fanfanKs;
     fanfanKs.candidates( *ks_i, IpProfile::e_position() );
     if ( fanfanKs.goodKs() )    preselect = true;

   }

   if ( preselect == false ) continue;

     Particle ks(*ks_i);

     // Vertex fit
     kvertexfitter kvf;
     kvf.initialVertex(ks.x());
     Particle  pi_p( ks_i->chgd(0), Ptype_PIplus );
     Particle  pi_m( ks_i->chgd(1), Ptype_PIminus );
     addTrack2fit(kvf, pi_p);
     addTrack2fit(kvf, pi_m);
     if(kvf.fit() != 0)  continue;

     UserInfo_B x; x.set(ks);
     dynamic_cast<UserInfo_B&>(ks.userInfo()).mass(ks.mass());
     dynamic_cast<UserInfo_B&>(ks.userInfo()).vtxchisq(kvf.chisq()/kvf.dgf());

     // mass constraint fit
     if ( m_D_massfit ){
       kmassfitter kmf;
       kmf.invariantMass(ks.pType().mass());
       kmf.vertex( ks.x() );
       addTrack2fit(kmf, pi_p );
       addTrack2fit(kmf, pi_m );
       if ( kmf.fit() != 0 ) continue;
       makeMother( kmf, ks );
     }else{
       makeMother( kvf, ks );
     }

     select = true;
     if ( frec_util::set_pGenHepevt( ks, m_mc ) == true ) select_mc = true;

     K_short.push_back(ks);

     if ( m_flagHist && select == true ){
	 const double p_ks = (frec_util::p_cm( ks.p() )).vect().mag();
	 const double m_ks = ks.p().mag();
	 H[36]->accumulate( p_ks, 1. );
	 H[52]->accumulate( m_ks, 1. );
	 if ( select_mc == true ){
	   H[46]->accumulate( p_ks, 1. );
	   H[62]->accumulate( m_ks, 1. );
	 }
     }

     if ( select == true ) {
       kid_statistics::fill_sta( m_name, m_cid_ks, (float)1. );
       if( select_mc == true ) kid_statistics::fill_sta( m_name, m_cid_ks_mc, (float) 1. );
     }

  }


  if ( m_flagHist ) H[26]->accumulate( K_short.size(), 1. );

  return;

}

// ------------------------------
//   Rho  reconstruction
// ------------------------------
void  frec_ana::recon_Rho( std::vector<Particle>& Rho, std::vector<Particle>& Pi, 
		           std::vector<Particle>& PI0 ) {

  const atc_pid selkpi(m_ACCq, m_TOFq, m_CDCq, KAON_CODE, PION_CODE);

  for( std::vector<Particle>::iterator pi_i = Pi.begin();
       pi_i != Pi.end(); pi_i++ ) {

    if( selkpi.prob(&pi_i->mdstCharged()) > m_PID_cut_PI_rho ) continue;

    for( std::vector<Particle>::iterator pi0_j = PI0.begin();
	 pi0_j != PI0.end(); pi0_j++ ) {

         bool select     = false;
	 bool select_mc  = false;
      
         const HepLorentzVector p_rho = frec_util::p_cm(pi_i->p() + pi0_j->p());

         const double pcm_rho = p_rho.vect().mag();
	 if ( pcm_rho  < m_RHO_Pstar_cut_low   || 
	      pcm_rho  > m_RHO_Pstar_cut_high   )  continue;

         const double rhomass = p_rho.mag();
	 if ( (fabs(rhomass - Ptype_RHOplus.mass()) > m_RHO_cut) ) continue;      

         select = true;

         Ptype ptype;
         if ( pi_i->charge() > 0 ) ptype = Ptype_RHOplus;
         else			   ptype = Ptype_RHOminus;

         Particle rho = frec_util::make_2body( ptype, *pi_i, *pi0_j , m_mc );
         if ( frec_util::set_pGenHepevt( rho, m_mc ) ) select_mc = true;
	 Rho.push_back(rho);

         //--- count good rho
        if ( select == true ){
	  kid_statistics::fill_sta( m_name, m_cid_rho, (float)1. );
         if ( select_mc == true ) kid_statistics::fill_sta( m_name, m_cid_rho_mc, (float)1. );        
	}

    } // end of loop pi0_j
  } // end of loop pi_i

  return;

}

// ------------------------------
//   Rho0  reconstruction
// ------------------------------
void  frec_ana::recon_Rho0( std::vector<Particle>& Rho0, 
			    std::vector<Particle>& PI ) {

  const atc_pid selkpi(m_ACCq, m_TOFq, m_CDCq, KAON_CODE, PION_CODE);

  for( std::vector<Particle>::iterator pi_i = PI.begin();
       pi_i != PI.end(); pi_i++ ) {

    if ( selkpi.prob(&pi_i->mdstCharged()) > m_PID_cut_PI_a1 ) continue;

    for( std::vector<Particle>::iterator pi_j = pi_i + 1;
	 pi_j != PI.end(); pi_j++ ) {

        if ( pi_i->charge() * pi_j->charge() > 0 ) continue;
	if ( selkpi.prob(&pi_j->mdstCharged()) > m_PID_cut_PI_a1 ) continue;

         const HepLorentzVector p_rho = frec_util::p_cm( pi_i->p() + pi_j->p() );

         const double pcm_rho0     = p_rho.vect().mag();
	 if( pcm_rho0 < m_RHO0_Pstar_cut_low   || 
	     pcm_rho0 > m_RHO0_Pstar_cut_high  )  continue;

	 const double rhomass = p_rho.mag();
	 if ( fabs ( rhomass -  Ptype_RHO0.mass() ) > m_RHO0_cut )  continue;

         Particle rho = frec_util::make_2body( Ptype_RHO0, *pi_i, *pi_j , m_mc );
	
         // Vertex fit
         kvertexfitter kvf;
	 kvf.initialVertex(IpProfile::e_position());
         addTrack2fit(kvf,*pi_i);
         addTrack2fit(kvf,*pi_j);

	 if(kvf.fit() != 0) continue;

         makeMother(kvf, rho );
	 double rchisq = kvf.chisq()/kvf.dgf();

         // Set userinfo
         UserInfo_B x; x.set(rho);
	 dynamic_cast<UserInfo_B&>(rho.userInfo()).vtxchisq(rchisq);
	 Rho0.push_back(rho);

     } // end of loop pi_j
   }// end of loop pi_i

  return;
}

// ------------------------------
//   A1  reconstruction
// ------------------------------
void  frec_ana::recon_A1(std::vector<Particle>& A1, std::vector<Particle>& RHO0, 
		      std::vector<Particle>& PI) {

  const atc_pid selkpi(m_ACCq, m_TOFq, m_CDCq, KAON_CODE, PION_CODE);

  for( std::vector<Particle>::iterator rho_i = RHO0.begin();
       rho_i != RHO0.end(); rho_i++ ) {

     const HepLorentzVector p_rho = frec_util::p_cm( rho_i->p() );

    for( std::vector<Particle>::iterator pi_j = PI.begin();
	 pi_j != PI.end(); pi_j++ ) {

      // remove duplication
      const bool chk_flag = frec_util::prohib_dupli(*rho_i, *pi_j);
      if ( chk_flag == true ) continue;

      if( selkpi.prob(&pi_j->mdstCharged()) > m_PID_cut_PI_a1 ) continue;

      bool select    = false;
      bool select_mc = false;

      const HepLorentzVector p_a1 = frec_util::p_cm( rho_i->p() + pi_j->p() );

      const double pcm_a1 = p_a1.vect().mag();
      if ( pcm_a1 < m_A1_Pstar_cut_low   || 
	   pcm_a1 > m_A1_Pstar_cut_high ) continue;

       const double a1mass = p_a1.mag();
       if ( a1mass < m_A1_cut_low || a1mass > m_A1_cut_high ) continue;
      
       Ptype ptype;
       if(pi_j->charge() > 0) ptype = Ptype_A1plus;
       else 		     ptype = Ptype_A1minus;
       Particle a1 = frec_util::make_2body( ptype, *rho_i, *pi_j , m_mc );
      
       kvertexfitter kvf;
       kvf.initialVertex(IpProfile::e_position());
       addTrack2fit(kvf,*rho_i);
       addTrack2fit(kvf,*pi_j);

       if( kvf.fit() != 0 )  continue;

       makeMother(kvf,a1);
       double rchi2 = kvf.chisq()/kvf.dgf();
       if( m_A1_Chisq_cut >0  && rchi2 > m_A1_Chisq_cut ) continue;

       select = true;
       if ( frec_util::set_pGenHepevt( a1, m_mc ) == true ) select_mc = true;
      
       UserInfo_B x; x.set(a1);
       dynamic_cast<UserInfo_B&>(a1.userInfo()).vtxchisq(rchi2);
       A1.push_back(a1);

       if ( select == true ){
	 kid_statistics::fill_sta( m_name, m_cid_a1, (float)1 );
	if ( select_mc == true ) kid_statistics::fill_sta( m_name, m_cid_a1_mc, (float)1. );
      }

    } // end of loop pi_j
  } // end of loop rho_i

  return;
}

// ------------------------------
//     All particles   
// ------------------------------
void 
frec_ana::recon_ALL( std::vector<Particle>& ALL, 
		     std::vector<Particle>& Charged_all,
		     std::vector<Particle>& Gamma ){

  for( std::vector<Particle>::iterator chg_i = Charged_all.begin();
       chg_i != Charged_all.end(); chg_i++ ) ALL.push_back(*chg_i);      

  for( std::vector<Particle>::iterator gam_i = Gamma.begin();
       gam_i != Gamma.end(); gam_i++ ) ALL.push_back(*gam_i);
    
  return;

}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
