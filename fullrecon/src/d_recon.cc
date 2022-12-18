//
// $Id: d_recon.cc 9932 2006-11-12 14:26:53Z katayama $
//
// $Log$
// Revision 1.7  2004/11/22 05:24:36  matumot
// updates for new kid_statistics
// module parameter, "LogFileName" is introduced
//
// Revision 1.6  2004/05/04 21:23:19  matumot
// fixed for CC complier
// - namespace was changed : fullrecon -> fullreconNS
//   to avoid overlap for fullrecon struct used by fullrecon.tdf
// - use new operator for dynamic allocation of array
// - add std:: for sort function
//
// Revision 1.4  2004/04/19 09:27:33  matumot
// revival version
//
// =====================================================
// File Name : d_recon.cc
// -----------------------------------------------------
// Creation    ; 2004.04.07
// Description ; Implimentation for the frec_ana class
//               ( D(*) reconstruction )
// Author      ; T.Matsumoto (matumot@bmail.kek.jp)
// ------------------------------------------------------

#include "belle.h"
#include  "fullrecon/frec_ana.h"
#include  "fullrecon/frec_util.h"

#include "ip/IpProfile.h"

#include "particle/Particle.h"
#include "particle/utility.h"
#include "kfitter/kvertexfitter.h"
#include "kfitter/kmassfitter.h"

#include "kid/kid_statistics.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


using namespace fullreconNS;

// Dalitz weight function
// ~~~~~~~~~~~~~~~~~~~~~~~~~~
extern "C" {
  extern float weight_(float*, float*, float*);
};


// -----------------------
//  D reconstruction
// -----------------------
void 
frec_ana::recon_D( std::vector<Particle>& D,   std::vector<Particle>& chgD,
		   std::vector<Particle>& DS,
	           std::vector<Particle>& K,   std::vector<Particle>& PI, 
		   std::vector<Particle>& PI0, std::vector<Particle>& Ks, 
		   std::vector<Particle>& GAM ){

  const atc_pid selkpi( m_ACCq, m_TOFq, m_CDCq, KAON_CODE, PION_CODE );
  Ptype ptype;

  int d_mode; double d_mass;

  //++++++ d__mode +++++++++++
  //   1: D0  --> Kpi
  //   7: D0  --> KK
  // 106: D+  --> K K pi
  // 302: Ds+ --> K K pi
  // 101: D+  --> K pi pi
  // 303: Ds+ --> K pi pi
  //   3: D0  --> K pi pi pi
  // 102: D+  --> K pi pi pi0 
  //   2: D0  --> K pi pi0
  // 103: D+  --> Ks pi
  // 301: Ds  --> Ks K
  //   5: D0  --> Ks pi pi
  // 105: D+  --> Ks pi pi pi 
  //   6: D0  --> Ks pi pi pi0
  // 104: D+  --> Ks pi pi0
  //   4: D0  --> ks pi0
  //+++++++++++++++++++++++++

  // Reconstruction of D0-->Kpi, K K, K K pi
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for( std::vector<Particle>::iterator k_i = K.begin();
       k_i != K.end(); k_i++ ) {

    const float chg_K  = k_i->charge();

    for( std::vector<Particle>::iterator pi_j = PI.begin();
	 pi_j != PI.end(); pi_j++ ) {

      const float chg_Pi = pi_j->charge();

      if( chg_K * chg_Pi > 0 ) continue;

      // mode 1: D0 --> K pi
      // ~~~~~~~~~~~~~~~~~~~~
      d_mode = 1;
      d_mass = ( k_i->p() + pi_j->p() ).mag();

      if ( fabs( d_mass - m_D_fitmass[d_mode] ) < m_D_cut[d_mode]  ){

	if( chg_K < 0 ) ptype = Ptype_D0; 
	else 	        ptype = Ptype_D0B;
	Particle d = frec_util::make_2body( ptype, *k_i, *pi_j, m_mc );
	if ( sel_D( d, d_mode ) == true ) D.push_back(d);

      }

    } // end of loop pi_j

    for( std::vector<Particle>::iterator k_j = k_i + 1;
	 k_j != K.end(); k_j++ ){

      if ( k_j->charge() * chg_K > 0 ) continue;

       // mode 7: D0 --> K K
       // ~~~~~~~~~~~~~~~~~~~~
	d_mode = 7;
	d_mass = ( k_i->p() + k_j->p() ).mag();

      if ( fabs( d_mass - m_D_fitmass[d_mode] ) < m_D_cut[d_mode] ){

	Particle d = frec_util::make_2body( Ptype_D0,   *k_i, *k_j, m_mc );
	if ( sel_D( d,  d_mode ) == true ) D.push_back(d);
	Particle db = frec_util::make_2body( Ptype_D0B, *k_i, *k_j, m_mc );
	if ( sel_D( db, d_mode ) == true ) D.push_back(db);

      }

    for( std::vector<Particle>::iterator pi_k = PI.begin();
	 pi_k != PI.end(); pi_k++ ){

        // avoid duplication
	if ( pi_k->mdstCharged().get_ID() == k_i->mdstCharged().get_ID() ) continue;
	if ( pi_k->mdstCharged().get_ID() == k_j->mdstCharged().get_ID() ) continue;

	// mode 106: D+ --> K K pi
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~
	d_mode = 106;
	d_mass = (k_i->p() + k_j->p() + pi_k->p() ).mag();

        if ( fabs( d_mass - m_D_fitmass[d_mode] ) < m_D_cut[d_mode]  ){

	  if ( pi_k->charge() > 0 ) ptype = Ptype_Dplus; 
	  else  	            ptype = Ptype_Dminus;

	  Particle d = frec_util::make_3body( ptype, *k_i, *k_j, *pi_k, m_mc );
	  if ( sel_D( d, d_mode ) == true ) chgD.push_back(d);

       }

       if ( m_reconDDs ){

	 // mode 302: Ds+ --> K K pi
         // ~~~~~~~~~~~~~~~~~~~~~~~~~~
	 d_mode = 302;

         if ( fabs( d_mass - m_D_fitmass[d_mode] ) < m_D_cut[d_mode]  ){

	  if ( pi_k->charge() > 0 ) ptype = Ptype_DSplus; 
	  else  	            ptype = Ptype_DSminus;

	  Particle d = frec_util::make_3body( ptype, *k_i, *k_j, *pi_k, m_mc );
	  d.usable(UNUSABLE);
	  double mphi = (k_i->p() + k_j->p()).mag();
	  if ( fabs(mphi - Ptype_Phi.mass()) < m_Phi_cut ) d.usable(USABLE);
	  if ( sel_D( d, d_mode ) == true ) DS.push_back(d);

        }

       }

      } // end of loop pi_k
    } // end of loop k_j
  } // end of loop k_i

  // Reconstruction of D-->K pi pi, K pi pi pi, 
  //                       K pi pi pi0, K pi pi0
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for( std::vector<Particle>::iterator k_i = K.begin();
       k_i != K.end(); k_i++ ) {

    const bool select_k = frec_util::sel_KPI( (*k_i), m_PID_cut_K,  selkpi );
    if ( ! select_k ) continue;

    const float chg_K  = k_i->charge();

    for( std::vector<Particle>::iterator pi_j = PI.begin();
	 pi_j != PI.end(); pi_j++ ) {

      const float chg_Pi = pi_j->charge();

      if( chg_K * chg_Pi > 0 ) continue;

       // we don't use Ds+ --> K pi pi as default for now
      if ( m_reconDDs && m_selDDs == 0 ){

      for( std::vector<Particle>::iterator pi_k = pi_j + 1;
	   pi_k != PI.end(); pi_k++ ) {

	if( chg_Pi * pi_k->charge() > 0 ) continue;

	// mode 303: Ds+ --> K pi pi
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~
	d_mode = 303;
	d_mass = ( k_i->p() + pi_j->p() + pi_k->p() ).mag();
	if ( fabs( d_mass - m_D_fitmass[d_mode] ) < m_D_cut[d_mode]  ){

	  if( chg_K > 0 ) ptype = Ptype_DSplus; 
	  else  	  ptype = Ptype_DSminus;

	  Particle d = frec_util::make_3body( ptype, *k_i, *pi_j, *pi_k, m_mc );
	  d.usable(UNUSABLE);
	  double mkst = (k_i->p() + pi_j->p()).mag();
	  if ( fabs( mkst - Ptype_Kstplus.mass() ) < m_Kst_cut ) d.usable(USABLE);
	  if ( sel_D( d, d_mode ) == true ) DS.push_back(d);

         }

      }// end of loop pi_k


     }

      for( std::vector<Particle>::iterator pi_k = pi_j + 1;
	   pi_k != PI.end(); pi_k++ ) {

	if( chg_Pi * pi_k->charge() < 0 ) continue;
	
	// mode 101: D+ --> K pi pi
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~
	d_mode = 101;
	d_mass = ( k_i->p() + pi_j->p() + pi_k->p() ).mag();

	if ( fabs( d_mass - m_D_fitmass[d_mode] ) < m_D_cut[d_mode]  ){

	  if( chg_K < 0 ) ptype = Ptype_Dplus; 
	  else  	        ptype = Ptype_Dminus;

	  Particle d = frec_util::make_3body( ptype, *k_i, *pi_j, *pi_k, m_mc );
	  if ( sel_D( d, d_mode ) == true ) chgD.push_back(d);

       }

       for( std::vector<Particle>::iterator pi_l = PI.begin();
	     pi_l != PI.end(); pi_l++ ) {

	 if( chg_K * pi_l->charge() < 0 ) continue;

	 // avoid duplication
	 if ( pi_l->mdstCharged().get_ID() == k_i->mdstCharged().get_ID()  ) continue;
	 if ( pi_l->mdstCharged().get_ID() == pi_j->mdstCharged().get_ID() ) continue;
	 if ( pi_l->mdstCharged().get_ID() == pi_k->mdstCharged().get_ID() ) continue;

         // mode 3: D0 --> K pi pi pi
	 // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  d_mode = 3;
	  d_mass = ( k_i->p() + pi_j->p() + pi_k->p() + pi_l->p() ).mag();

	  if ( fabs( d_mass - m_D_fitmass[d_mode] ) < m_D_cut[d_mode] ){

	    if( chg_K < 0 ) ptype = Ptype_D0; 
	    else 	    ptype = Ptype_D0B;

	    Particle d = frec_util::make_4body( ptype, *k_i, *pi_j, *pi_k, *pi_l, m_mc );
	    if ( sel_D( d, d_mode ) == true ) D.push_back(d);

	  }

      } // end of loop pi_l

      for( std::vector<Particle>::iterator pi0_l = PI0.begin();
	   pi0_l != PI0.end(); pi0_l++ ) {

 	 // mode 102: D+ --> K pi pi pi0 
	 // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  d_mode = 102;
	  d_mass = ( k_i->p() + pi_j->p() + pi_k->p() + pi0_l->p() ).mag();

	  if ( fabs( d_mass - m_D_fitmass[d_mode] ) < m_D_cut[d_mode]  ){

	    if ( chg_K < 0 ) ptype = Ptype_Dplus;
	    else             ptype = Ptype_Dminus;

	    Particle d = frec_util::make_4body( ptype, *k_i, *pi_j, *pi_k, *pi0_l, m_mc );
	    if ( sel_D( d, d_mode ) == true ) chgD.push_back(d);

	  }

	} // end of loop pi0_l
      } // end of loop pi_k

      for( std::vector<Particle>::iterator pi0_k = PI0.begin();
	   pi0_k != PI0.end(); pi0_k++) {

        // mode 2: D0 -> K pi pi0
        // ~~~~~~~~~~~~~~~~~~~~~~~~
	d_mode = 2;
	d_mass = ( k_i->p() + pi_j->p() + pi0_k->p() ).mag();

        if ( fabs( d_mass - m_D_fitmass[d_mode] ) < m_D_cut[d_mode]  ){

	  if( chg_K < 0 ) ptype = Ptype_D0; 
	  else 	          ptype = Ptype_D0B;

	  Particle d = frec_util::make_3body( ptype, *k_i, *pi_j, *pi0_k, m_mc );
	  if ( sel_D( d, d_mode ) == true ) D.push_back(d);

	}

       }// end of loop pi0_k
     } // end of loop pi_j
   } // end of loop k_i

  // Reconstruction of D--> Ks K, Ks pi, Ks pi pi, Ks pi pi pi
  //                        K s pi pi pi0, Ks pi pi0, Ks pi0
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for( std::vector<Particle>::iterator ks_i = Ks.begin();
       ks_i != Ks.end(); ks_i++ ) {

    for( std::vector<Particle>::iterator k_j = K.begin();
	 k_j != K.end(); k_j++ ) {

      // avoid duplication
      const bool chk_flag = frec_util::prohib_dupli( *ks_i, *k_j );
      if( chk_flag == true ) continue;

     if ( m_reconDDs ){

      // mode 301: Ds+ -> Ks K
      // ~~~~~~~~~~~~~~~~~~~~~~~~
      d_mode = 301;
      d_mass = ( ks_i->p() + k_j->p() ).mag();

      if ( fabs( d_mass - m_D_fitmass[d_mode] ) < m_D_cut[d_mode]  ){

	if( k_j->charge() > 0 ) ptype = Ptype_DSplus; 
	else 		         ptype = Ptype_DSminus;
	Particle d = frec_util::make_2body( ptype, *ks_i, *k_j, m_mc );
	if ( sel_D( d, d_mode ) == true ) DS.push_back(d);

      }

     }

    } // end of loop k_j

    for( std::vector<Particle>::iterator pi_j = PI.begin();
	 pi_j != PI.end(); pi_j++ ) {

      // avoid duplication
      const bool chk_flag = frec_util::prohib_dupli( *ks_i, *pi_j );
      if( chk_flag == true ) continue;

      // mode 103: D+ -> Ks pi
      // ~~~~~~~~~~~~~~~~~~~~~~
      d_mode = 103;
      d_mass = ( ks_i->p() + pi_j->p() ).mag();

      if ( fabs( d_mass - m_D_fitmass[d_mode] ) < m_D_cut[d_mode]  ){

	if( pi_j->charge() > 0 ) ptype = Ptype_Dplus; 
	else 		         ptype = Ptype_Dminus;
	Particle d = frec_util::make_2body( ptype, *ks_i, *pi_j, m_mc );
	if ( sel_D( d, d_mode ) == true ) chgD.push_back(d);

      }

      for( std::vector<Particle>::iterator pi_k = pi_j+ 1;
	   pi_k != PI.end(); pi_k++ ) {

	if( pi_j->charge() * pi_k->charge() > 0 ) continue;

	// avoid duplication
	const bool chk_flag2 = frec_util::prohib_dupli( *ks_i, *pi_k );
	if( chk_flag2 == true ) continue;

       // mode 5: D0 -> Ks pi pi
       // ~~~~~~~~~~~~~~~~~~~~~~
	d_mode = 5;
	d_mass = ( ks_i->p() + pi_j->p() + pi_k->p() ).mag();

       if ( fabs( d_mass - m_D_fitmass[d_mode] ) < m_D_cut[d_mode]  ){

	 Particle d = frec_util::make_3body( Ptype_D0,   *ks_i, *pi_j, *pi_k, m_mc );
	 if ( sel_D( d,  d_mode ) == true ) D.push_back(d);

	 Particle db = frec_util::make_3body( Ptype_D0B, *ks_i, *pi_j, *pi_k, m_mc );
	 if ( sel_D( db, d_mode ) == true ) D.push_back(db);

       }

      for( std::vector<Particle>::iterator pi_l = PI.begin();
	   pi_l != PI.end(); pi_l++ ){

	// avoid duplication
	const bool chk_flag3 = frec_util::prohib_dupli( *ks_i, *pi_l );
	if( chk_flag3 == true ) continue;

        // avoid duplication
	if ( pi_l->mdstCharged().get_ID() == pi_j->mdstCharged().get_ID() ) continue;
	if ( pi_l->mdstCharged().get_ID() == pi_k->mdstCharged().get_ID() ) continue;

       // mode 105: D+ --> Ks pi pi pi 
       // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	d_mode = 105;
	d_mass = ( ks_i->p() + pi_j->p() + pi_k->p() + pi_l->p() ).mag();

        if ( fabs( d_mass - m_D_fitmass[d_mode] ) < m_D_cut[d_mode]  ){

	  if ( pi_l->charge() > 0 ) ptype = Ptype_Dplus; 
	  else 		            ptype = Ptype_Dminus;
	  Particle d = frec_util::make_4body( ptype, *ks_i, *pi_j, *pi_k, *pi_l, m_mc );
	  if ( sel_D( d, d_mode ) == true ) chgD.push_back(d);

        }

      } //end of loop pi_l

      for( std::vector<Particle>::iterator pi0_l = PI0.begin();
	   pi0_l != PI0.end(); pi0_l++ ) {

	// mode 6: D0 --> Ks pi pi pi0
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	d_mode = 6;
	d_mass = ( ks_i->p() + pi_j->p() + pi_k->p() + pi0_l->p() ).mag();

        if ( fabs( d_mass - m_D_fitmass[d_mode] ) < m_D_cut[d_mode]  ){

	 Particle d  = frec_util::make_4body( Ptype_D0,  *ks_i, *pi_j, *pi_k, *pi0_l, m_mc );
	 if ( sel_D( d,  d_mode ) == true ) D.push_back(d);
	 Particle db = frec_util::make_4body( Ptype_D0B, *ks_i, *pi_j, *pi_k, *pi0_l, m_mc );
	 if ( sel_D( db, d_mode ) == true ) D.push_back(db);

       }

       } //end of loop pi0_l
      } // end of loop pi_k

     for( std::vector<Particle>::iterator pi0_k = PI0.begin();
	  pi0_k != PI0.end(); pi0_k++ ) {
      
	// mode 104: D+ --> Ks pi pi0
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	d_mode = 104;
	d_mass = ( ks_i->p() + pi_j->p() + pi0_k->p() ).mag();

        if ( fabs( d_mass - m_D_fitmass[d_mode] ) < m_D_cut[d_mode]  ){

	 if( pi_j->charge() > 0) ptype = Ptype_Dplus; 
	 else 		         ptype = Ptype_Dminus;
	 Particle d = frec_util::make_3body( ptype, *ks_i, *pi_j, *pi0_k, m_mc );
	 if ( sel_D( d, d_mode ) == true ) chgD.push_back(d);

       }

      }// end of loop pi0_k
     } // end of loop pi_j

     for( std::vector<Particle>::iterator pi0_j = PI0.begin();
	  pi0_j != PI0.end(); pi0_j++ ) {

       // mode 4: D0 -> Ks pi0
       // ~~~~~~~~~~~~~~~~~~~~~
       d_mode = 4;
       d_mass = ( ks_i->p() + pi0_j->p() ).mag();

      if ( fabs( d_mass - m_D_fitmass[d_mode] ) < m_D_cut[d_mode]  ){

	Particle d = frec_util::make_2body( Ptype_D0,   *ks_i, *pi0_j, m_mc );
	if ( sel_D( d,  d_mode ) == true ) D.push_back(d);
	Particle db = frec_util::make_2body( Ptype_D0B, *ks_i, *pi0_j, m_mc );
	if ( sel_D( db, d_mode ) == true ) D.push_back(db);

      }

     }// end of loop pi0_j
    } // end of loop ks_i

  return;
}

// -----------------------
//  D* reconstruction
// -----------------------
void 
frec_ana::recon_Dst( std::vector<Particle>& Dst, std::vector<Particle>& chgDst,
		     std::vector<Particle>& DSst,
		     std::vector<Particle>& D, std::vector<Particle>& chgD,
		     std::vector<Particle>& DS,
		     std::vector<Particle>& PI, std::vector<Particle>& PI0, 
		     std::vector<Particle>& GAM ){

  int dst_mode; double d_mass, dm;
  //++++++ dst_mode +++++++++
  // 1:  D*0  --> D0 pi0
  // 2:  D*0  -->D0 gamma
  // 3:  D*+  -->D0 pi+
  // 4:  D*+  -->D+ pi0
  // 5:  Ds*+ -->Ds+ gamma 
  //+++++++++++++++++++++++++

  Ptype ptype;

  for( std::vector<Particle>::iterator d_i = D.begin();
       d_i != D.end(); d_i++ ) {

     const int d_mode( dynamic_cast<UserInfo_B&>(d_i->userInfo()).sub1_mode());
     if ( d_mode == 0 ) continue;
     const int lund_d(d_i->lund());
    
     d_mass  = d_i->p().mag();

    for( std::vector<Particle>::iterator pi0_j = PI0.begin();
	 pi0_j != PI0.end(); pi0_j++ ) {

      // avoid duplication
      const bool chk_flag = frec_util::prohib_dupli(*d_i, *pi0_j);
      if( chk_flag == true ) continue;

      // mode 1: D*0 -> D0 pi0
      // ~~~~~~~~~~~~~~~~~~~~~~
      dst_mode = 1;
      dm = (d_i->p() + pi0_j->p()).mag() - d_mass;

       if ( fabs( dm - m_Dm[dst_mode] ) < m_Dm_cut[dst_mode] ){

        if( lund_d > 0 ) ptype = Ptype_Dst0;
        else 	         ptype = Ptype_Dst0B;

        Particle dst = frec_util::make_2body( ptype, *d_i, *pi0_j, m_mc );
        if ( sel_Dst( dst, dst_mode ) == true ) Dst.push_back(dst);

      }

    }// end of loop pi0_j

    for( std::vector<Particle>::iterator gamma_j = GAM.begin();
	 gamma_j != GAM.end(); gamma_j++ ) {

      // avoid duplication
      const bool chk_flag = frec_util::prohib_dupli( *d_i, *gamma_j );
      if( chk_flag == true ) continue;

       // mode 2: D*0 -> D0 gamma
       // ~~~~~~~~~~~~~~~~~~~~~~~~
       dst_mode = 2;
       dm = (d_i->p() + gamma_j->p()).mag() - d_mass;

       if ( fabs( dm - m_Dm[dst_mode] ) < m_Dm_cut[dst_mode] ){

         if( lund_d > 0 ) ptype = Ptype_Dst0;
	 else 	        ptype = Ptype_Dst0B;

	 Particle dst = frec_util::make_2body( ptype, *d_i, *gamma_j, m_mc );
	 if ( sel_Dst( dst, dst_mode ) == true ) Dst.push_back(dst);

       }

    } // end of loop gamma_j

    for( std::vector<Particle>::iterator pi_j = PI.begin();
	 pi_j != PI.end(); pi_j++ ) {

      if( lund_d * pi_j->charge() < 0 ) continue;

      // avoid duplication
      const bool chk_flag = frec_util::prohib_dupli( *d_i, *pi_j );
      if( chk_flag == true ) continue;

     // mode 3: D*+ -> D0 pi
     // ~~~~~~~~~~~~~~~~~~~~~
      dst_mode = 3;
      dm = (d_i->p() + pi_j->p()).mag() - d_mass;

      if ( fabs( dm - m_Dm[dst_mode] ) < m_Dm_cut[dst_mode] ){

        if( lund_d > 0 ) ptype = Ptype_Dstplus; 
	else 	         ptype = Ptype_Dstminus;

	Particle dst = frec_util::make_2body( ptype, *d_i, *pi_j, m_mc );
	if ( sel_Dst(dst, dst_mode ) == true ) chgDst.push_back(dst);

      }

    } // end of loop pi_j
  }// end of loop d_i

  for( std::vector<Particle>::iterator d_i = chgD.begin();
       d_i != chgD.end(); d_i++ ) {

    const int d_mode(dynamic_cast<UserInfo_B&>(d_i->userInfo()).sub1_mode());
    if ( d_mode == 0 ) continue;
    const int lund_d(d_i->lund());

     const HepLorentzVector p_d = d_i->p();
     double d_mass = p_d.mag();

    for( std::vector<Particle>::iterator pi0_j = PI0.begin();
	 pi0_j != PI0.end(); pi0_j++ ) {

      // mode 4: D*+ -> D+ pi0
      // ~~~~~~~~~~~~~~~~~~~~~~~
      dst_mode = 4;
      dm = (d_i->p() + pi0_j->p()).mag() - d_mass;

      if ( fabs( dm - m_Dm[dst_mode] ) < m_Dm_cut[dst_mode] ){

         if( lund_d > 0) ptype = Ptype_Dstplus; 
         else 	         ptype = Ptype_Dstminus;

         Particle dst = frec_util::make_2body( ptype, *d_i, *pi0_j, m_mc );
         if ( sel_Dst( dst, dst_mode ) == true ) chgDst.push_back( dst );

      }

    } // end of loop pi0_j
  } // end of loop d_i

 if ( m_reconDDs ){

  for( std::vector<Particle>::iterator ds_i = DS.begin();
       ds_i != DS.end(); ds_i++ ) {

     const int d_mode( dynamic_cast<UserInfo_B&>(ds_i->userInfo()).sub1_mode());
     if ( d_mode == 0 ) continue;
     const int lund_d(ds_i->lund());
    
     d_mass  = ds_i->p().mag();

    for( std::vector<Particle>::iterator gamma_j = GAM.begin();
	 gamma_j != GAM.end(); gamma_j++ ) {

      // avoid duplication
      const bool chk_flag = frec_util::prohib_dupli( *ds_i, *gamma_j );
      if( chk_flag == true ) continue;

       // mode 5: Ds*+ -> Ds+ gamma
       // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
       dst_mode = 5;
       dm = (ds_i->p() + gamma_j->p()).mag() - d_mass;

       if ( fabs( dm - m_Dm[dst_mode] ) < m_Dm_cut[dst_mode] ){

         if( lund_d > 0 ) ptype = Ptype_DSstplus;
	 else 	          ptype = Ptype_DSstminus;

	 Particle dst = frec_util::make_2body( ptype, *ds_i, *gamma_j, m_mc );
	 dst.usable(UNUSABLE);
	 if ( ds_i->usable() ) dst.usable(USABLE);
	 if ( sel_Dst( dst, dst_mode ) == true ) DSst.push_back(dst);
       }

     } // end of loop gamma_j
   }// end of loop ds_i

  }

  return;

}

// -----------------------
//  D selection
// -----------------------
bool
frec_ana::sel_D(Particle &d, const int d_mode  ){

  bool select        = false;
  bool select_mc     = false;

  double d_sigma_default(0);
  switch ( d_mode ){
    // set wider modes for pi0 including modes
    case 2: case 102: case 104: d_sigma_default = 0.015; break;
    case 4:   d_sigma_default  = 0.03;
    default : d_sigma_default = 0.006; break;
  }

  // Set MC flag
  if ( m_mc ){

      if ( d.genHepevt() ) select_mc = true;
      // check dbar also for cp states
      if ( d_mode == 4 || d_mode == 5 || d_mode == 6 || d_mode == 7 ){
        Particle dbar( d.p(), Ptype(-d.lund()));
        for (unsigned i=0; i<d.nChildren();i++)	
	  dbar.relation().append(d.child(i));
        if ( frec_util::set_pGenHepevt(dbar, m_mc) == true ) select_mc = true;
       }

   }

  // Vertex fit
  double rchisq(0.);

  if ( m_D_vertexfit ){

    if ( d_mode == 4 || d_mode == 103 || d_mode == 104 || d_mode == 301 ) {
      select = true; // for D->Kspi0, Kspi, Kspipi0 and D->KsK don't apply vertex fit
    }else {

      kvertexfitter kvf;
      kvf.initialVertex(IpProfile::e_position());
      for(unsigned i=0;i<d.nChildren();i++){
	 if ( d.child(i).mdstCharged() ) addTrack2fit(kvf,d.child(i));
      }

      if( kvf.fit() == 0 ){
	select = true;
	rchisq = kvf.chisq()/kvf.dgf();
      }

    }

    if( m_D_Chisq_cut > 0.  && rchisq > m_D_Chisq_cut ) select = false;

   }else {
     select = true;
   }

   // selection with pion PID
   const atc_pid selkpi( m_ACCq, m_TOFq, m_CDCq, KAON_CODE, PION_CODE );

   double pid_pid(-1.), tmp_pid(-1.);
   for ( unsigned i=0; i<d.nChildren(); i++ ){
      if ( abs(d.child(i).lund()) == Ptype_PIplus.lund() ){
	tmp_pid = selkpi.prob(d.child(i).mdstCharged());
	if ( tmp_pid > pid_pid ) pid_pid = tmp_pid;
      }
    }
   if ( pid_pid > m_PID_cut_PI_D ) select = false;

    // Dalitz weight cut for D0-->K-pi+pi0
    if ( d_mode == 2 ){

	Particle &k   = d.child(0);
	Particle &pi  = d.child(1);
	Particle &pi0 = d.child(2);
	const HepLorentzVector d_p    = k.p() + pi.p() + pi0.p();
	const HepLorentzVector kpi_p  = k.p() + pi.p();
	const HepLorentzVector kpi0_p = k.p() + pi0.p();

	float d_mass = (float)d.mass();
	float kpi_mag2 = (float)kpi_p.mag2();
	float kpi0_mag2 = (float)kpi0_p.mag2();
	const float dalitz = weight_(&d_mass, &kpi_mag2, &kpi0_mag2);

	if ( dalitz <= m_Dalitz_cut ) select = false;

     }

   // D candidate
   if( select == true ) {

	// mass chisq
	double d_sigma(sqrt(d.momentum().dMass()));
	if ( !(d_sigma < 100) || d_sigma == 0 ) d_sigma = d_sigma_default;
	const double masschi = (d.mass()-m_D_fitmass[d_mode])/d_sigma;

	// Set UserInfo
	UserInfo_B x; x.set(d);
	dynamic_cast<UserInfo_B&>(d.userInfo()).sub1_mode(d_mode);
	dynamic_cast<UserInfo_B&>(d.userInfo()).masschisq(masschi*masschi);
	dynamic_cast<UserInfo_B&>(d.userInfo()).vtxchisq(rchisq);
	dynamic_cast<UserInfo_B&>(d.userInfo()).mass(d.mass());

	if ( m_D_massfit ){

	  kmassfitter kmf;
	  kmf.invariantMass(d.pType().mass());
	  kmf.vertex( IpProfile::e_position() );
	  for(unsigned i=0;i<d.nChildren();i++) addTrack2fit(kmf,d.child(i));
	  if( kmf.fit() != 0 ) select = false;
	  else                 makeMother( kmf, d );

	}

    }

   if ( select == true ){
     kid_statistics::fill_sta( m_name, m_cid_d,(float)1.);
     if ( select_mc == true ) kid_statistics::fill_sta( m_name, m_cid_d_mc,(float)1.);
   }

  return select;

}

// -----------------------
//  D* selection
// -----------------------
bool
frec_ana::sel_Dst( Particle &dst, const int dst_mode ){

  bool select    = false;
  bool select_mc = false;

  const int d_mode = dynamic_cast<UserInfo_B&>(dst.child(0).userInfo()).sub1_mode();

  if ( dst_mode * d_mode == 0 ) return select;

    select = true;
    if ( m_mc && dst.genHepevt() ) select_mc = true;  

    const double dm = dst.mass() - dst.child(0).mass();
    double dm_sigma = sqrt( dst.momentum().dMass() - dst.child(0).momentum().dMass());
	    	                
    double dm_sigma_default;
    switch( dst_mode ){
      case 2: case 5: dm_sigma_default = 0.01;   break;
      case 3:  dm_sigma_default = 0.0005; break;
      default: dm_sigma_default = 0.001;  break;
    }
    if ( !(dm_sigma < 100) || dm_sigma == 0 ) dm_sigma = dm_sigma_default;

    const double masschi = (dm - m_Dm[dst_mode])/dm_sigma;

    UserInfo_B x; x.set(dst);
    dynamic_cast<UserInfo_B&>(dst.userInfo()).sub0_mode(dst_mode);
    dynamic_cast<UserInfo_B&>(dst.userInfo()).sub1_mode(d_mode);
    dynamic_cast<UserInfo_B&>(dst.userInfo()).masschisq(masschi*masschi);
    dynamic_cast<UserInfo_B&>(dst.userInfo()).mass(dst.mass());

    if ( m_D_massfit ){

      kmassfitter kmf;
      kmf.invariantMass(dst.pType().mass());
      kmf.vertex( IpProfile::e_position() );
      for(unsigned i=0;i<dst.nChildren();i++) addTrack2fit(kmf,dst.child(i));
      if ( kmf.fit() != 0 ) select = false;
      else  makeMother( kmf, dst );

    }

   if ( select == true ) {
      kid_statistics::fill_sta( m_name, m_cid_dst, (float)1. );
      if ( select_mc == true ) kid_statistics::fill_sta( m_name, m_cid_dst_mc, (float)1. );
   }

   return select;

}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
