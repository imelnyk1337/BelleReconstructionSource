// -*- C++ -*-
//
// Description: 
//   To reproduce Mdst_charged from Mdst_trk, Mdst_trk_fit, Mdst_sim_xref and 
//                                  Mdst_klm_mu_ex
//                Mdst_gamma   from Mdst_ecl
//                Mdst_pi0     from Mdst_ecl and Mdst_gamma
//                Mdst_klong   from Mdst_klm_cluster and mdst_ecl
//                Mdst_ecl_cr  from Mdst_ecl_aux and Mdst_ecl
//                (except for mass, width, and nhits)
//
// Author:      Hiroshi Hamasaki, 
// Created:     Thu Dec 21 15:27:23 JST 2000
// $Id: mdst2mdst.cc 10002 2007-02-26 06:56:17Z katayama $
//
// Revision history
//
// $Log$
// Revision 1.5  2002/02/25 01:38:00  katayama
// Bug fix(no effect)
//
// Revision 1.4  2001/12/12 02:04:10  hitoshi
// updated for mdst_ecl/trk (by Kakuno).
//
// Revision 1.3  2001/02/02 02:41:30  hitoshi
// fix by Hamasaki.
//
// Revision 1.2  2001/01/31 13:08:01  hitoshi
// commented out 2 lines in response to hkim's request.
//
// Revision 1.1  2001/01/31 11:47:54  hitoshi
// a function to reproduce mdst_charged ... tables.
//

#include "belle.h"
#include <math.h>

#include "panther/panther.h"
#include "belleCLHEP/Vector/ThreeVector.h"
#include "belleCLHEP/Vector/LorentzVector.h"
#include "ecl/Pi0Fitter.h"
#include "mdst/mdst.h"

#include MDST_H
#include ECL_H
#include ECLG_H
#include ECLP_H
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


// external functions
extern "C" {
  void recsim_mdst_propgt_( float*, float[], float[], float[],
			    float[], float[], int* );
  int rec_mdst_fill_pi0_();
};

// forward declaration
//void mdst2mdstCharged( void );
void mdst2mdstKlong( void );
void mdst2mdstGamma( void );
void mdst2mdstPi0( void );
void mdst2mdstEclCr( void );
void mdstGamma2RececlGamma( void );
void mdstEcl2RececlShower( void );
void pi0_from_2gamma( Rececl_gamma_Manager& rececl_gamma_manager,
		      Rececl_pi0_Manager& rececl_pi0_manager );

// ==========================================================================
// main function
void
mdst2mdst( bool, bool fillMdstKlong,
	   bool fillMdstGamma,   bool fillMdstPi0,
	   bool fillMdstEclCr )
{
  // Mdst_charged
  // mdst_charged should already be kept. 2001/12/10 H.Kakuno
  //if ( fillMdstCharged ) { mdst2mdstCharged(); }

  // Mdst_klong
  if ( fillMdstKlong ) { mdst2mdstKlong(); }

  // Mdst_gamma
  if ( fillMdstGamma || fillMdstPi0 ) { mdst2mdstGamma(); }

  // Mdst_pi0
  if ( fillMdstPi0 ) { mdst2mdstPi0(); }

  // Mdst_ecl_cr
  if ( fillMdstEclCr ) { mdst2mdstEclCr(); }
}

// ==========================================================================
// reproduce Mdst_charged from Mdst_trk, Mdst_trk_fit, Mdst_sim_xref and
// Mdst_klm_mu_ex
/*
void 
mdst2mdstCharged( void )
{
  // get_manager for Mdst_charged, Mdst_trk, Mdst_trk_fit, Mdst_sim_xref and
  // Mdst_klm_mu_ex
  Mdst_charged_Manager & mdst_charged_manager 
    = Mdst_charged_Manager::get_manager();
  Mdst_trk_Manager & mdst_trk_manager     
    = Mdst_trk_Manager::get_manager();
  Mdst_trk_fit_Manager & mdst_trk_fit_manager 
    = Mdst_trk_fit_Manager::get_manager();
  Mdst_sim_xref_Manager & mdst_sim_xref_manager
    = Mdst_sim_xref_Manager::get_manager();
  Mdst_klm_mu_ex_Manager & mdst_klm_mu_ex_manager
    = Mdst_klm_mu_ex_Manager::get_manager();

  // if already exist, then return
  if ( mdst_charged_manager.count() > 0 ) { return; }

  // clear Mdst_charged
  mdst_charged_manager.remove();

  // loop over the Mdst_trk
  for ( Mdst_trk_Manager::iterator mtrk = mdst_trk_manager.begin();
	mtrk != mdst_trk_manager.end(); ++mtrk ) {

    // pointer to Mdst_trk_fit with pion hypothese
    Mdst_trk_fit & mdst_trk_fit = (*mtrk).mhyp(2);

    // check quality
    if ( (*mtrk).quality() != 0 ) { continue; }

    // new entry for mdst_charged
    Mdst_charged & mdst_charged = mdst_charged_manager.add();

    // charge
    //mdst_charged.charge( sign(1., mdst_trk_fit.helix(2)) );
    if ( mdst_trk_fit.helix(2) >= 0 ) {
      mdst_charged.charge( 1. );
    } else {
      mdst_charged.charge( -1. );
    }

    // momentum at origin
    float pivot[3], helix[5], error[15], helix0[5], error0[15];
    for ( int k=0; k<15; k++ ) {
      if ( k < 3 ) { pivot[k] = mdst_trk_fit.pivot(k); }
      if ( k < 5 ) { helix[k] = mdst_trk_fit.helix(k); }
                     error[k] = mdst_trk_fit.error(k);
    }
    
    float amass = mdst_trk_fit.mass();

    // propagate helix to origin
    int iret;
    recsim_mdst_propgt_( &amass, pivot, helix,  error,
			                helix0, error0, &iret );
    if ( iret == 0 ) {
      mdst_charged.px( -sin(helix0[1])/fabs(helix0[2]) );
      mdst_charged.py(  cos(helix0[1])/fabs(helix0[2]) );
      mdst_charged.pz(      helix0[4] /fabs(helix0[2]) );
    }

    // mass
    mdst_charged.mass( amass );

    // pointer to Mdst_trk
    mdst_charged.trk( *mtrk ); 

    // pointer to Mdst_acc, Mdst_tof, Mdst_klm
    for ( Mdst_sim_xref_Manager::iterator xref = mdst_sim_xref_manager.begin();
	  xref != mdst_sim_xref_manager.end(); ++xref ) {
      if ( (*mtrk).get_ID() == (*xref).trk().get_ID() ) { 
	mdst_charged.acc( (*xref).acc() );
	mdst_charged.tof( (*xref).tof() ); 
	mdst_charged.klm( (*xref).klm() );
	break;
      }
    }

    // pointer to Mdst_muid (rec_mdst_fill_mu.cc)
    // (*) When the numbering of Mdst_muid is the same as Mdst_klm_mu_ex,
    //     this code is valid.
    for ( Mdst_klm_mu_ex_Manager::iterator mu = mdst_klm_mu_ex_manager.begin();
	  mu != mdst_klm_mu_ex_manager.end(); ++mu ) {
      if ( mdst_charged.get_ID() == (*mu).pMDST_Charged().get_ID() ) {
	mdst_charged.muid( Mdst_muid_Manager::get_manager()((*mu).get_ID()) );
	break;
      }
    }
  }
}
*/

// ==========================================================================
// reproduce Mdst_klong from Mdst_klm_cluster and Mdst_ecl
// created by M.Yamaga @ 15 Jan 2001
void
mdst2mdstKlong( void )
{
  // Code to reproduce Mdst_klong table from Mdst_klm_cluster and Mdst_ecl. 
  // 15-Jan-2001 M.Yamaga
  Mdst_klong_Manager&       klmgr  = Mdst_klong_Manager::get_manager();
  Mdst_klm_cluster_Manager& kcmgr  = Mdst_klm_cluster_Manager::get_manager();
  Mdst_ecl_Manager&         eclmgr = Mdst_ecl_Manager::get_manager();

  // if already exist, then return
  if ( klmgr.count() > 0 ) { return; }

  // clear Mdst_klong
  klmgr.remove();

  // loop over the Mdst_klm_cluster
  for( Mdst_klm_cluster_Manager::iterator kc = kcmgr.begin();
       kc != kcmgr.end(); kc++ ){
    if ( kc->trk_ID() == 0 ){ 
      // neutral cluster
      if ( kc->ecl_ID() &&       // matched to ECL cluster
           kc->layers() >= 1 ){  // #layers >= 1
        Mdst_klong&  klong = klmgr.add();
        bool  foundecl = false;
        for( Mdst_ecl_Manager::iterator ecl = eclmgr.begin();
             ecl != eclmgr.end(); ecl++ ){
          if ( ecl->get_ID() == kc->ecl_ID() ){
            klong.cos_x( sin( ecl->theta() ) * cos( ecl->phi() ) );
            klong.cos_y( sin( ecl->theta() ) * sin( ecl->phi() ) );
            klong.cos_z( cos( ecl->theta() )                     );
            klong.ecl( *ecl ); 
            klong.reset_klm();
            klong.klmc( *kc );
            foundecl = true;
            break;
          }
        }
        if ( !foundecl ){
          klong.cos_x( sin( kc->theta() ) * cos( kc->phi() ) );
          klong.cos_y( sin( kc->theta() ) * sin( kc->phi() ) );
          klong.cos_z( cos( kc->theta() )                    );
          klong.reset_ecl(); 
          klong.reset_klm();
          klong.klmc( *kc );
        }
      }
      else if ( kc->layers() >= 2 ){ // not matched to ECL cluster
        Mdst_klong&  klong = klmgr.add();
        klong.cos_x( sin( kc->theta() ) * cos( kc->phi() ) );
        klong.cos_y( sin( kc->theta() ) * sin( kc->phi() ) );
        klong.cos_z( cos( kc->theta() )                    );
        klong.reset_ecl(); 
        klong.reset_klm();
        klong.klmc( *kc );
      }
    }
  }
  // [end] Code to reproduce Mdst_klong table.
}

// ==========================================================================
// reproduce Mdst_gamma from Mdst_ecl 
// Mdst_gamma <= Rececl_gamma  (rec_mdst_fill_gamma.cc)
//            <= Rececl_shower (rec/ecl/ecl-gamma/basf_if/rececl_gamma.cc)
//            == Mdst_ecl      (rec_mdst_fill_ecl2.cc)

void
mdst2mdstGamma( void )
{
  // get manager for Mdst_ecl and Mdst_gamma
  Mdst_ecl_Manager   & mdst_ecl_manager   = Mdst_ecl_Manager::get_manager();
  Mdst_gamma_Manager & mdst_gamma_manager = Mdst_gamma_Manager::get_manager();

  // if already exist, then return
  if ( mdst_gamma_manager.count() > 0 ) { return; }

  // clear Mdst_gamma
  mdst_gamma_manager.remove();

  // loop over the Mdst_ecl to make Mdst_gamma
  for ( Mdst_ecl_Manager::iterator i = mdst_ecl_manager.begin();
	i != mdst_ecl_manager.end(); ++i ) {

    if ( (*i).quality() != 0 ) { continue; }
    if ( (*i).match()   == 1 ) { continue; }
    if ( !good_gamma( *i )   ) { continue; }

    double energy = (*i).energy();
    double theta  = (*i).theta();
    double phi    = (*i).phi();

    double px = energy * sin(theta) * cos(phi);
    double py = energy * sin(theta) * sin(phi);
    double pz = energy * cos(theta);

    // add
    Mdst_gamma & gamma = mdst_gamma_manager.add();
    gamma.px( px );
    gamma.py( py );
    gamma.pz( pz );
    gamma.ecl( *i );
  }
}

// ==========================================================================
// reproduce Mdst_pi0 from Mdst_gamma and Mdst_ecl
//  Mdst_pi0 == Rececl_pi0 (rec_mdst_fill_pi0.cc)
//           <= Rececl_gamma(==Mdst_gamma), Rececl_shower(==Mdst_ecl) 
//              (rec/ecl/ecl-pi0/basf_if/rececl_pi0.cc)
void
mdst2mdstPi0( void )
{
  // get manager for Mdst_pi0 
  Mdst_pi0_Manager & pi0Mgr = Mdst_pi0_Manager::get_manager();

  // if already exist, then return
  if ( pi0Mgr.count() > 0 ) { return; }

  // Mdst_ecl => Rececl_shower
  if ( BsCouTab( RECECL_SHOWER ) <= 0 ) {
    mdstEcl2RececlShower();
  }

  // Mdst_gamma => Rececl_gamma
  if ( BsCouTab( RECECL_GAMMA ) <= 0 ) {
    mdstGamma2RececlGamma();
  }

  // clear Rececl_pi0
  Rececl_pi0_Manager::get_manager().remove();

  // Rececl_gamma and Rececl_shower => Rececl_pi0
   Rececl_gamma_Manager& rececl_gamma_manager
      = Rececl_gamma_Manager::get_manager();
   Rececl_pi0_Manager& rececl_pi0_manager
      = Rececl_pi0_Manager::get_manager();

   pi0_from_2gamma(rececl_gamma_manager, rececl_pi0_manager);

   // Rececl_pi0 => Mdst_pi0
   Mdst_pi0_Manager::get_manager().remove();
   rec_mdst_fill_pi0_();
}

// --------------------------------------------------------------------------
// Mdst_ecl => Rececl_shower
void
mdstEcl2RececlShower( void )
{
  // get manager for Mdst_ecl and Rececl_shower
  Mdst_ecl_Manager      & mdst_ecl_manager 
    = Mdst_ecl_Manager::get_manager();
  Rececl_shower_Manager & rececl_shower_manager 
    = Rececl_shower_Manager::get_manager();

  // clear Rececl_shower
  rececl_shower_manager.remove();
  
  // loop over the Mdst_ecl to fill Rececl_shower
  for ( Mdst_ecl_Manager::iterator i = mdst_ecl_manager.begin();
	i != mdst_ecl_manager.end(); ++i ) {

    Rececl_shower & rececl_shower = rececl_shower_manager.add();

    rececl_shower.status( (*i).quality() );
    rececl_shower.energy( (*i).energy()  );
    rececl_shower.phi   ( (*i).phi()     );
    rececl_shower.theta ( (*i).theta()   );
    rececl_shower.r     ( (*i).r()       );

    for ( int j=0; j<6; ++j ) {
      rececl_shower.error( j, (*i).error(j) );
    }

  }
}

// --------------------------------------------------------------------------
// Mdst_gamma => Rececl_gamma
void
mdstGamma2RececlGamma( void )
{
  // get manager for Mdst_gamma and Rececl_gamma
  Mdst_gamma_Manager   & mdst_gamma_manager 
    = Mdst_gamma_Manager::get_manager();
  Rececl_gamma_Manager & rececl_gamma_manager 
    = Rececl_gamma_Manager::get_manager();

  // clear Rececl_gamma
  rececl_gamma_manager.remove();

  // loop over the Mdst_gamma to fill Rececl_gamma
  for ( Mdst_gamma_Manager::iterator i = mdst_gamma_manager.begin();
	i != mdst_gamma_manager.end(); ++i ) {
    
    Rececl_gamma & rececl_gamma = rececl_gamma_manager.add();

    rececl_gamma.shower( Rececl_shower_Manager::
			 get_manager()( (*i).ecl().get_ID()) );
    rececl_gamma.px( (*i).px() );
    rececl_gamma.py( (*i).py() );
    rececl_gamma.pz( (*i).pz() );

  }
}

// --------------------------------------------------------------------------
// ( == RecECL_Pi0::pi0_from_2gamma() @ b20001228_1026)
void
pi0_from_2gamma( Rececl_gamma_Manager& rececl_gamma_manager,
		 Rececl_pi0_Manager& rececl_pi0_manager )
{
  // constant
  static const float  gamma_energy_threshold(.02);
  static const float  pi0_mass_min(.134 - .0054 * 3); // 3 sigma
  static const float  pi0_mass_max(.134 + .0054 * 3); // 3 sigma
  static const int    fit_flag(1);

  const int n_gamma = rececl_gamma_manager.count();
  if (n_gamma < 2)
    return;

  for (int i = 0; i < n_gamma - 1; ++i) {
    Panther_ID ID_gamma1(i + 1);
    Rececl_gamma& gamma1 = rececl_gamma_manager(ID_gamma1);

    // cluster type should be checked here

    Hep3Vector p3_gamma1(gamma1.px(), gamma1.py(), gamma1.pz());

    // gamma energy cut
    const double E_gamma1 = p3_gamma1.mag();
    if (E_gamma1 < gamma_energy_threshold)
      continue;

    HepLorentzVector lv_gamma1(p3_gamma1, E_gamma1);

    for (int j = i + 1; j < n_gamma; ++j) {
      Panther_ID ID_gamma2(j + 1);
      Rececl_gamma& gamma2 = rececl_gamma_manager(ID_gamma2);

      // cluster type should be checked here

      Hep3Vector p3_gamma2(gamma2.px(), gamma2.py(), gamma2.pz());

      // gamma energy cut
      const double E_gamma2 = p3_gamma2.mag();
      if (E_gamma2 < gamma_energy_threshold)
	continue;

      Hep3Vector p3_rec = p3_gamma1 + p3_gamma2;

      // calculate opening angle
      // calculate minimum opening angle (acos(2 * beta^2 - 1))
      // apply opening angle cut
      // if (opening_angle < minimum_opening_angle - opening_angle_cut_margin)
      // continue;
      // ------------------  opening anglecut removed 99/5/27/skkim 

      HepLorentzVector lv_gamma2(p3_gamma2, E_gamma2);
      HepLorentzVector lv_rec = lv_gamma1 + lv_gamma2;
      const double mass = lv_rec.mag();

      if (fit_flag) {
	Pi0Fitter pi0fitter(gamma1, gamma2);
	pi0fitter.fit(Pi0Fitter::MASS_PI0);
//	if (pi0fitter.chi2() < chi2_max) {
// -----------------------  chi2 cut  removed   99/5/27/skkim
// fill RecECL_Pi0
        if (pi0_mass_min < mass && mass < pi0_mass_max) {	
	  Rececl_pi0& rececl_pi0 = rececl_pi0_manager.add();
	  rececl_pi0.gamma(0, gamma1);
	  rececl_pi0.gamma(1, gamma2);
	  rececl_pi0.px(pi0fitter.px());
	  rececl_pi0.py(pi0fitter.py());
	  rececl_pi0.pz(pi0fitter.pz());
	  rececl_pi0.energy(pi0fitter.E());
	  rececl_pi0.mass(mass);
	  rececl_pi0.massfit(pi0fitter.mass());
	  rececl_pi0.chi2(pi0fitter.chi2());
	}
      } else if (pi0_mass_min < mass && mass < pi0_mass_max) {
	// mass region cut
	// fill RecECL_Pi0
	Rececl_pi0& rececl_pi0 = rececl_pi0_manager.add();
	rececl_pi0.gamma(0, gamma1);
	rececl_pi0.gamma(1, gamma2);
	rececl_pi0.px(lv_rec.px());
	rececl_pi0.py(lv_rec.py());
	rececl_pi0.pz(lv_rec.pz());
	rececl_pi0.energy(lv_rec.e());
	rececl_pi0.mass(mass);
	rececl_pi0.massfit(0);
	rececl_pi0.chi2(0);
      }
    }
  }  
}

// ==========================================================================
// reporduce all variables in Mdst_ecl_cr except for mass, width, and nhits
// from Mdst_ecl_aux and Mdst_ecl
// created by H.W.Kim @ 2001/Jan/31
void
mdst2mdstEclCr( void )
{
   int n = 1; // this is important variable
// in order to fill mdst_ecl_cr table here
   Mdst_ecl_cr_Manager &mecr_mag = Mdst_ecl_cr_Manager::get_manager();
// to find out showers
   Mdst_ecl_Manager &me_mag = Mdst_ecl_Manager::get_manager();
   Mdst_ecl_aux_Index index( Mdst_ecl_aux_Manager::get_manager().index("cr") );
   index.update();
// to find out showers

 // if already exist, then return
   //   if ( mecr_mag.count() > 0 ) { return; }

   //   mecr_mag.remove();

   while (1) {
      Panther_ID crID( n++ ); // we convert a given integer number to panther id.
      std::vector<Mdst_ecl_aux> vPtr = point_from(crID, index); // collect showers

      if ( vPtr.size() == 0 ) break;

      Mdst_ecl_cr& mecr = mecr_mag.add();

      // initialization part
      double cr_energy = 0.0;
      double cr_theta = 0.0;
      double cr_phi = 0.0;
      double cr_r = 0.0;
      bool cr_match = false;

      Hep3Vector x_wsum;

#if 0
// this loop is for obtaining energy sum and x_mean first
      for( Mdst_ecl_aux_Manager::iterator i = vPtr.begin(); i != vPtr.end(); i++ ) {
         Mdst_ecl& me = me_mag( i->get_ID() );
         
         double px = me.r() * sin( me.theta() ) * cos( me.phi() );
         double py = me.r() * sin( me.theta() ) * sin( me.phi() );
         double pz = me.r() * cos( me.theta() );
         
         Hep3Vector xxx(px, py, pz);
         e_sum += me.energy();
         x_wsum += me.energy() * xxx;
      }

// center of gravity from shower position vectors
// but will be same as that obtained from hits.
      Hep3Vector x_min = x_wsum * (1/e_sum);
      
      double x_minus_x_min_squared = 0.0;
      double parl_axised = 0.0;
#endif

      for( Mdst_ecl_aux_Manager::iterator i = vPtr.begin(); i != vPtr.end(); i++ ) {
         Mdst_ecl& me = me_mag( i->get_ID() );
#if 0
         double px = me.r() * sin( me.theta() ) * cos( me.phi() );
         double py = me.r() * sin( me.theta() ) * sin( me.phi() );
         double pz = me.r() * cos( me.theta() );
         Hep3Vector xxx(px, py, pz);
         Hep3Vector x_minus_xmin = xxx - x_min;
         x_minus_x_min_squared += x_minus_xmin.mag2() * me.energy();

         parl_axised += ((i->width() * i->width()) + x_minus_xmin.mag2()) * me.energy();
#endif

         cr_energy += me.energy();
         //cr_nhits += i->nhits();
         cr_theta += me.energy() * me.theta();
         cr_phi   += me.energy() * me.phi();
         cr_r     += me.energy() * me.r();
         if(!cr_match && me.match() != 0 ) cr_match = true;
         
      } // end of for loop

      //double fWidth = sqrt(parl_axised / e_sum);

// variables calculated from new method
      double new_energy = cr_energy;
      int new_match = cr_match;

      double new_r = cr_r/cr_energy;
      double new_theta = cr_theta/cr_energy;
      double new_phi = cr_phi/cr_energy;

      //double new_width = (vPtr.size() == 1 ? vPtr[0].width() : fWidth);
      //double new_mass = (vPtr.size() == 1 ? vPtr[0].width() : fWidth) * new_energy / new_r;
      //int new_nhits = cr_nhits;

// filling new mdst_ecl_cr.
      mecr.energy( new_energy );
      mecr.match( new_match );
      mecr.nshowers( vPtr.size() );

      mecr.r( new_r );
      mecr.theta( new_theta );
      mecr.phi( new_phi );

      //mecr.nhits( new_nhits );
      //mecr.width( new_width );
      //mecr.mass( new_mass );

      mecr.quality(0);
      for (int i = 0; i < 10; ++i)
         mecr.property(i, 0.0);

   } // end of while loop
   
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
