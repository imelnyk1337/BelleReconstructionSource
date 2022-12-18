//
// $Id: frec_util.cc 10548 2008-07-24 13:13:11Z hitoshi $
//
// $Log$
// Revision 1.3  2004/05/04 21:23:20  matumot
// fixed for CC complier
// - namespace was changed : fullrecon -> fullreconNS
//   to avoid overlap for fullrecon struct used by fullrecon.tdf
// - use new operator for dynamic allocation of array
// - add std:: for sort function
//
// Revision 1.2  2004/04/19 11:07:20  matumot
// change $Header$ to $Id: frec_util.cc 10548 2008-07-24 13:13:11Z hitoshi $
//
// Revision 1.1  2004/04/19 09:27:34  matumot
// revival version
//
//
// ==========================================================
//  File Name :  frec_util.cc
// ----------------------------------------------------------
// Creation    ; 2004.04.07
// Description ; Implimentation for the frec_util class etc.
// Author      ; T.Matsumoto (matumot@bmail.kek.jp)
// ----------------------------------------------------------

#include  "belle.h"

#include  "fullrecon/frec_util.h"

#include  "ip/IpProfile.h"

#include "belleCLHEP/Matrix/Vector.h"
#include "belleCLHEP/Matrix/SymMatrix.h"
#include "helix/Helix.h"

#include "particle/Particle.h"
#include "particle/gammac.h"

#include "toolbox/FoxWolfr.h"
#include "toolbox/FuncPtr.h"

#include  EVTCLS_H
#include  MDST_OBS_H
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


using namespace fullreconNS;

// Define parameters
// ~~~~~~~~~~~~~~~~~~
const    int fullreconNS::ELE_CODE(0);
const    int fullreconNS::MUON_CODE(1);
const    int fullreconNS::PION_CODE(2);
const    int fullreconNS::KAON_CODE(3);
const    int fullreconNS::PROTON_CODE(4);

const  Ptype fullreconNS::Ptype_B0("B0");
const  Ptype fullreconNS::Ptype_B0B("B0B");
const  Ptype fullreconNS::Ptype_Bplus("B+");
const  Ptype fullreconNS::Ptype_Bminus("B-");
const  Ptype fullreconNS::Ptype_Dst0("D*0");
const  Ptype fullreconNS::Ptype_Dst0B("D*B");
const  Ptype fullreconNS::Ptype_Dstplus("D*+");
const  Ptype fullreconNS::Ptype_Dstminus("D*-");
const  Ptype fullreconNS::Ptype_D0("D0");
const  Ptype fullreconNS::Ptype_D0B("D0B");
const  Ptype fullreconNS::Ptype_Dplus("D+");
const  Ptype fullreconNS::Ptype_Dminus("D-");
const  Ptype fullreconNS::Ptype_DSplus("DS+");
const  Ptype fullreconNS::Ptype_DSminus("DS-");
const  Ptype fullreconNS::Ptype_DSstplus("DS*+");
const  Ptype fullreconNS::Ptype_DSstminus("DS*-");
const  Ptype fullreconNS::Ptype_Ks("K0S");
const  Ptype fullreconNS::Ptype_PI0("PI0");
const  Ptype fullreconNS::Ptype_Gamma("GAMM");
const  Ptype fullreconNS::Ptype_Kplus("K+");
const  Ptype fullreconNS::Ptype_Kminus("K-");
const  Ptype fullreconNS::Ptype_PIplus("PI+");
const  Ptype fullreconNS::Ptype_PIminus("PI-");
const  Ptype fullreconNS::Ptype_Eplus("E+");
const  Ptype fullreconNS::Ptype_Eminus("E-");
const  Ptype fullreconNS::Ptype_MUplus("MU+");
const  Ptype fullreconNS::Ptype_MUminus("MU-");
const  Ptype fullreconNS::Ptype_RHOplus("RHO+");
const  Ptype fullreconNS::Ptype_RHOminus("RHO-");
const  Ptype fullreconNS::Ptype_RHO0("RHO0");
const  Ptype fullreconNS::Ptype_A1plus("A1+");
const  Ptype fullreconNS::Ptype_A1minus("A1-");
const  Ptype fullreconNS::Ptype_Eta("ETA");
const  Ptype fullreconNS::Ptype_Etap("ETA'");
const  Ptype fullreconNS::Ptype_Kstplus("K*+");
const  Ptype fullreconNS::Ptype_Kstminus("K*-");
const  Ptype fullreconNS::Ptype_Kst0("K*0");
const  Ptype fullreconNS::Ptype_Kst0B("K*B");
const  Ptype fullreconNS::Ptype_Omega("OMEG");
const  Ptype fullreconNS::Ptype_Phi("PHI");

const double fullreconNS::E_HER(7.998213);  // GeV
const double fullreconNS::E_LER(3.499218);  // GeV
const double fullreconNS::CROSS_ANGLE(0.022);  // mrad

// ---------------------------
//  frec_util class
// ---------------------------

// initialization
// ~~~~~~~~~~~~~~~~~
const HepLorentzVector 
frec_util::m_p_beam  = HepLorentzVector(fullreconNS::E_HER * sin( fullreconNS::CROSS_ANGLE), 
					 0., 
					fullreconNS::E_HER * cos( fullreconNS::CROSS_ANGLE) 
					 -  fullreconNS::E_LER,  
					fullreconNS::E_HER + fullreconNS::E_LER );

// member functions
// ~~~~~~~~~~~~~~~~~~~

// correct dr, dz by IP position
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
const double frec_util::correct_dr( const Mdst_charged& chg, const HepPoint3D ip, 
				    const int id ){

  double co_dr = 10000.;

  Mdst_trk& trk = chg.trk();
  Mdst_trk_fit& fit = trk.mhyp(id);

  if(trk && fit){
#if defined(BELLE_DEBUG)
    try {
#endif
    // Get Helix prameter
    HepVector  a(5,0);
    a[0] = fit.helix(0);
    a[1] = fit.helix(1);
    a[2] = fit.helix(2);
    a[3] = fit.helix(3);
    a[4] = fit.helix(4);
    
    // Get the error matrix
    HepSymMatrix Ea(5,0); 
    Ea[0][0]=fit.error(0);
    Ea[1][0]=fit.error(1); Ea[1][1]= fit.error(2);
    Ea[2][0]=fit.error(3); Ea[2][1]= fit.error(4);
    Ea[2][2]=fit.error(5);
    Ea[3][0]=fit.error(6); Ea[3][1]= fit.error(7);
    Ea[3][2]=fit.error(8);
    Ea[3][3]=fit.error(9);
    Ea[4][0]=fit.error(10); Ea[4][1]=fit.error(11);
    Ea[4][2]=fit.error(12);
    Ea[4][3]=fit.error(13); Ea[4][4]=fit.error(14);
  
    const HepPoint3D pivot(fit.pivot(0), fit.pivot(1), fit.pivot(2));

    Helix crtrk(pivot, a, Ea);              // create Helix crtrk
    crtrk.pivot(ip);                        // move pivot point to vertex
  
    co_dr=crtrk.dr();                      //distance between helix and pivot
#if defined(BELLE_DEBUG)
    }
    catch(const std::string&e) {
      dout(Debugout::ERR,"frec_util") << "frec_util::correct_dr:" << e << std::endl;
    }
#endif    

  }

    return(co_dr);

}

const double frec_util::correct_dz( const Mdst_charged& chg, const HepPoint3D ip, 
				    const int id){

  double co_dz = 10000.;

  Mdst_trk& trk = chg.trk();
  Mdst_trk_fit& fit = trk.mhyp(id);

  if(trk && fit){
#if defined(BELLE_DEBUG)
    try {
#endif
    
    // Get Helix prameter
    HepVector  a(5,0);
    a[0] = fit.helix(0);
    a[1] = fit.helix(1);
    a[2] = fit.helix(2);
    a[3] = fit.helix(3);
    a[4] = fit.helix(4);

    // Get the error matrix
    HepSymMatrix Ea(5,0); 
    Ea[0][0]=fit.error(0);
    Ea[1][0]=fit.error(1); Ea[1][1]= fit.error(2);
    Ea[2][0]=fit.error(3); Ea[2][1]= fit.error(4);
    Ea[2][2]=fit.error(5);
    Ea[3][0]=fit.error(6); Ea[3][1]= fit.error(7);
    Ea[3][2]=fit.error(8);
    Ea[3][3]=fit.error(9);
    Ea[4][0]=fit.error(10); Ea[4][1]=fit.error(11);
    Ea[4][2]=fit.error(12);
    Ea[4][3]=fit.error(13); Ea[4][4]=fit.error(14);
  
    const HepPoint3D pivot(fit.pivot(0), fit.pivot(1), fit.pivot(2));

    Helix crtrk(pivot, a, Ea);            // create Helix ltrk 
    crtrk.pivot(ip);                      // move pivot point to vertex
    
    co_dz=crtrk.dz();                    // same but for Z plane
#if defined(BELLE_DEBUG)
    }
    catch(const std::string&e) {
      dout(Debugout::ERR,"frec_util") << "frec_util::correct_dr:" << e << std::endl;
    }
#endif    
  }

  return(co_dz);

}


// remove duplication
// ~~~~~~~~~~~~~~~~~~~~~~~
const bool frec_util::prohib_dupli( Particle chk1, Particle chk2 ){

  bool flag = false;

  const int chk1_nptcle = chk1.relation().nFinalStateParticles();
  const int chk2_nptcle = chk2.relation().nFinalStateParticles();

  for(int  i =0; i < chk1_nptcle; i++) {
    for(int  j =0; j < chk2_nptcle; j++) {
      Particle chk1_child = chk1.relation().finalStateParticle(i);
      Particle chk2_child = chk2.relation().finalStateParticle(j);
      const int chk1_chg = chk1_child.mdstCharged().get_ID();
      const int chk2_chg = chk2_child.mdstCharged().get_ID();
      const int chk1_gam = chk1_child.mdstGamma().get_ID();
      const int chk2_gam = chk2_child.mdstGamma().get_ID();
      if(chk1_chg*chk2_chg != 0  && (chk1_chg == chk2_chg) ) flag = true;
      if(chk1_gam*chk2_gam != 0 && (chk1_gam == chk2_gam) )  flag = true;
    }
  }

  return flag;

}

// dEdx flag
// ~~~~~~~~~~~
const bool frec_util::with_dEdx( Particle &p, int pcode, double cut_sigma ){

  bool flag = false;

  if ( pcode < 0 || pcode > 4 ) return flag;
  if ( p.mdstCharged() == NULL ) return flag;

  Mdst_trk &trk = p.mdstCharged().trk();

  if ( trk.dEdx() == 0. || trk.dEdx_exp(pcode) == 0. || trk.sigma_dEdx(pcode) == 0. ) return flag;

  double sigma = ( trk.dEdx() - trk.dEdx_exp(pcode))/trk.sigma_dEdx(pcode);

  if ( fabs(sigma) > cut_sigma ) return flag;

  return true;

}

// set generator information
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
const bool frec_util::set_pGenHepevt( Particle &p, const int flag_mc ){

 bool set = false;
 if ( flag_mc == 0 ) return set;

 switch( flag_mc ){
  case 1:
   if    ( p.mdstCharged() )  set_OldpGenHepevt_mdstCharged(p, p.pType());
   else if ( p.mdstGamma() )  set_OldpGenHepevt_mdstGamma(p);
   else if   ( p.mdstPi0() )  set_OldpGenHepevt_mdstPi0(p);
   else	if   ( p.mdstVee2() ) set_OldpGenHepevt_mdstVee2(p);
   else			      set_pGenHepevt_general( p );
   break;
  case 2:
   if    ( p.mdstCharged() )  set_pGenHepevt_mdstCharged(p, p.pType());
   else if ( p.mdstGamma() )  set_pGenHepevt_mdstGamma(p);
   else if   ( p.mdstPi0() )  set_pGenHepevt_mdstPi0(p);
   else	if  ( p.mdstVee2() )  set_pGenHepevt_mdstVee2(p);
   else			      set_pGenHepevt_general(p);
   break;
 default:
   break;

 }

 if ( p.genHepevt() ) set = true;

 return set;

}

const bool frec_util::set_pGenHepevt_mdstCharged( Particle& p, const Ptype& ptype ){

  Mdst_sim_trk_Manager& simtrkmgr = Mdst_sim_trk_Manager::get_manager();

  for( std::vector<Mdst_sim_trk>::const_iterator i = simtrkmgr.begin();
      i != simtrkmgr.end(); i++) {

    if( i->hepevt().isthep() > 0 && i->hepevt().idhep() == ptype.lund()) {
      if( p.mdstCharged().trk().get_ID() == i->trk().get_ID()) {

	const Gen_hepevt& hepevt = i ->hepevt();
	p.relation().genHepevt(hepevt);

	return true;

      }
    }
  }

  return false;

}

const bool frec_util::set_pGenHepevt_mdstGamma( Particle &p ){

  Mdst_sim_ecl_Manager& simeclmgr = Mdst_sim_ecl_Manager::get_manager();

  for( std::vector<Mdst_sim_ecl>::const_iterator i = simeclmgr.begin();
       i != simeclmgr.end(); i++) {

      if( i->hepevt().isthep() > 0  && (i->hepevt().idhep() == Ptype_Gamma.lund()) ) {

      if( p.mdstGamma().ecl().get_ID() == i->ecl().get_ID()) {

	const Gen_hepevt& hepevt = i ->hepevt();
	p.relation().genHepevt(hepevt);

	return true;

      }
    }
  }

  return false;

}

const bool frec_util::set_pGenHepevt_mdstPi0( Particle &p ){

  Mdst_sim_ecl_Manager& eclmgr = Mdst_sim_ecl_Manager::get_manager();

  Gen_hepevt* gam_evt1 = NULL;
  Gen_hepevt* gam_evt2 = NULL;

  for( std::vector<Mdst_sim_ecl>::const_iterator i = eclmgr.begin();
       i != eclmgr.end(); i++) {

    if( (i->ecl().get_ID() == p.mdstPi0().gamma(0).ecl().get_ID()) && 
        (i->hepevt().idhep() == Ptype_Gamma.lund()))
        gam_evt1 = &(i->hepevt());

    if( (i->ecl().get_ID() == p.mdstPi0().gamma(1).ecl().get_ID()) && 
        (i->hepevt().idhep() == Ptype_Gamma.lund()))
       gam_evt2 = &(i->hepevt());

  }
  
  if( gam_evt1 && gam_evt2 && gam_evt1->mother() && gam_evt2->mother()) {

    if( (gam_evt1->mother().get_ID() != gam_evt2->mother().get_ID()) ||
        (gam_evt1->mother().idhep() != Ptype_PI0.lund()) ||
        (gam_evt2->mother().idhep() != Ptype_PI0.lund())) {

        return false;

    }
  
    const Gen_hepevt& hepevt = gam_evt1->mother();
    p.relation().genHepevt(hepevt);

    return true;

  }

  return false;

}

const bool frec_util::set_pGenHepevt_mdstVee2( Particle &p ){

  Mdst_sim_trk_Manager& trkmgr = Mdst_sim_trk_Manager::get_manager();

  Gen_hepevt* chg_evt1 = NULL;
  Gen_hepevt* chg_evt2 = NULL;

  for(std::vector<Mdst_sim_trk>::const_iterator i = trkmgr.begin();
      i != trkmgr.end(); i++) {

    if( ( i->trk().get_ID()   == p.mdstVee2().chgd(0).trk().get_ID()) &&
        ( i->hepevt().idhep() == Ptype_PIplus.lund()))
      chg_evt1 = &(i->hepevt());

    if( ( i->trk().get_ID()   == p.mdstVee2().chgd(1).trk().get_ID()) &&
        ( i->hepevt().idhep() == Ptype_PIminus.lund()))
      chg_evt2 = &(i->hepevt());

  }
  
  if( chg_evt1 && chg_evt2 && chg_evt1->mother() && chg_evt2->mother()) {

    if((chg_evt1->mother().get_ID() != chg_evt2->mother().get_ID()) ||
       (chg_evt1->mother().idhep()  != Ptype_Ks.lund()) ||
       (chg_evt2->mother().idhep()  != Ptype_Ks.lund())) {

      return false;

    }

    const Gen_hepevt& hepevt = chg_evt1->mother();
    p.relation().genHepevt(hepevt);

    return true;

  }

  return false;

}    

const bool frec_util::set_pGenHepevt_general( Particle& p ){

  const Gen_hepevt* p_mother = NULL;

  const int nchildren = p.nChildren();
  if(!nchildren) return false;

  // Check that genHep references exist;
  for(int i=0; i<nchildren; ++i)

    if(!p.relation().child(i).genHepevt()) return false;

  // Check that child particles haven't same genHep reference;
  for( int i=0; i <nchildren-1; ++i )
    for( int j=i+1; j< nchildren; ++j )

      if( p.relation().child(i).genHepevt().get_ID() == 
          p.relation().child(j).genHepevt().get_ID() ) return false;

  // Seek mother by the first daughter;
  const Gen_hepevt *mother(&(p.child(0).genHepevt()));

  while( mother->mother() ){
    mother = &(mother->mother());
    if( mother->idhep() == p.pType().lund() ) break;
  }

  if( mother->idhep() != p.pType().lund() ) return false;
  
  // Check for other children that have the same mother;
  for( int i=1; i<nchildren; ++i ){

    const Gen_hepevt *tmp(&(p.child(i).genHepevt()));
    while( tmp->mother() ){
      tmp = &tmp->mother();
      if( tmp == mother ) break;
    }

    if( tmp != mother ) return false;

  }

  // Check that there are no another children from this mother;
  double e=0;
  for( int i=0; i<nchildren; ++i )
    e+=p.child(i).genHepevt().E();
  if(abs(e-mother->E())<.0001) {

    p_mother = mother;
    p.relation().genHepevt(*p_mother);

    return true;

  }
  
  return false;
}

const bool frec_util::set_OldpGenHepevt_mdstCharged( Particle& p, 
					       const Ptype& ptype ){

  Mdst_sim_xref_Manager& simxrefmgr = Mdst_sim_xref_Manager::get_manager();

  for( std::vector<Mdst_sim_xref>::const_iterator i = simxrefmgr.begin();
       i != simxrefmgr.end(); i++) {

    if( i->hepevt().isthep() > 0 && i->hepevt().idhep() == ptype.lund()) {
      if( p.mdstCharged().trk().get_ID() == i->trk().get_ID()) {

	const Gen_hepevt& hepevt = i ->hepevt();
	p.relation().genHepevt(hepevt);

	return true;

      }
    }
  }

  return false;

}

const bool frec_util::set_OldpGenHepevt_mdstGamma( Particle& p ){

  Mdst_sim_xref_Manager& simxrefmgr = Mdst_sim_xref_Manager::get_manager();

  for( std::vector<Mdst_sim_xref>::const_iterator i = simxrefmgr.begin();
       i != simxrefmgr.end(); i++) {

    if( i->hepevt().isthep() > 0  && (i->hepevt().idhep() == Ptype_Gamma.lund())) {
      if( p.mdstGamma().ecl().get_ID() == i->ecl().get_ID()) {

	const Gen_hepevt& hepevt = i ->hepevt();
	p.relation().genHepevt(hepevt);

	return true;

      }
    }
  }

  return false;

}

const bool frec_util::set_OldpGenHepevt_mdstPi0( Particle& p ){

  Mdst_sim_xref_Manager& xrefmgr = Mdst_sim_xref_Manager::get_manager();

  Gen_hepevt* gam_evt1 = NULL;
  Gen_hepevt* gam_evt2 = NULL;

  for(std::vector<Mdst_sim_xref>::const_iterator i = xrefmgr.begin();
      i != xrefmgr.end(); i++) {

    if( (i->ecl().get_ID() == p.mdstPi0().gamma(0).ecl().get_ID()) && 
       (i->hepevt().idhep() == Ptype_Gamma.lund()))
        gam_evt1 = &(i->hepevt());

    if( (i->ecl().get_ID() == p.mdstPi0().gamma(1).ecl().get_ID()) && 
        (i->hepevt().idhep() == Ptype_Gamma.lund()))
        gam_evt2 = &(i->hepevt());
  }
  
  if(gam_evt1 && gam_evt2 && gam_evt1->mother() && gam_evt2->mother()) {

    if((gam_evt1->mother().get_ID() != gam_evt2->mother().get_ID()) ||
       (gam_evt1->mother().idhep() != Ptype_PI0.lund()) ||
       (gam_evt2->mother().idhep() != Ptype_PI0.lund())) {

      return false;

    }

    const Gen_hepevt& hepevt = gam_evt1->mother();
    p.relation().genHepevt(hepevt);

    return true;

  }

  return false;

}

const bool frec_util::set_OldpGenHepevt_mdstVee2( Particle& p ){

  Mdst_sim_xref_Manager& xrefmgr = Mdst_sim_xref_Manager::get_manager();

  Gen_hepevt* chg_evt1 = NULL;
  Gen_hepevt* chg_evt2 = NULL;

  for( std::vector<Mdst_sim_xref>::const_iterator i = xrefmgr.begin();
       i != xrefmgr.end(); i++) {

    if( ( i->trk().get_ID() == p.mdstVee2().chgd(0).trk().get_ID()) &&
        ( i->hepevt().idhep() == Ptype_PIplus.lund()))
        chg_evt1 = &(i->hepevt());

    if( ( i->trk().get_ID() == p.mdstVee2().chgd(1).trk().get_ID()) &&
        ( i->hepevt().idhep() == Ptype_PIminus.lund()))
        chg_evt2 = &(i->hepevt());
  }
  
  if(chg_evt1 && chg_evt2 && chg_evt1->mother() && chg_evt2->mother()) {
    if((chg_evt1->mother().get_ID() != chg_evt2->mother().get_ID()) ||
       (chg_evt1->mother().idhep() != Ptype_Ks.lund()) ||
       (chg_evt2->mother().idhep() != Ptype_Ks.lund())) {

      return false;

    }

    const Gen_hepevt& hepevt = chg_evt1->mother();
    p.relation().genHepevt(hepevt);

    return true;

  }

  return false;

}    

// particles combinations
// ~~~~~~~~~~~~~~~~~~~~~~~~
Particle frec_util::make_2body( const Ptype ptype,
                       Particle& p1, Particle &p2, const int mc ) {

  Momentum mother_mom( p1.p() +  p2.p(),
                       p1.momentum().dp() + 
                       p2.momentum().dp());

  Particle mother(mother_mom, ptype );

  mother.relation().append( p1 );
  mother.relation().append( p2 );

  if ( mc ) frec_util::set_pGenHepevt( mother );

  return mother;

}

Particle frec_util::make_3body( const Ptype ptype,
				Particle& p1, Particle &p2,
				Particle& p3 , const int mc ) {

  Momentum mother_mom( p1.p() + p2.p() + p3.p(),
                       p1.momentum().dp() + 
                       p2.momentum().dp() + 
                       p3.momentum().dp() );

  Particle mother( mother_mom, ptype );

  mother.relation().append( p1 );
  mother.relation().append( p2 );
  mother.relation().append( p3 );

  if ( mc ) frec_util::set_pGenHepevt( mother );

  return mother;

}

Particle 
frec_util::make_4body( const Ptype ptype,
		       Particle& p1, Particle &p2, 
		       Particle& p3, Particle &p4, const int mc ) {

  Momentum mother_mom( p1.p()+ p2.p() + p3.p() + p4.p(),
                       p1.momentum().dp() + 
                       p2.momentum().dp() + 
                       p3.momentum().dp() + 
                       p4.momentum().dp() );

  Particle mother( mother_mom, ptype );

  mother.relation().append( p1 );
  mother.relation().append( p2 );
  mother.relation().append( p3 );
  mother.relation().append( p4 );

  if ( mc ) frec_util::set_pGenHepevt( mother );

  return mother;

}

// K/PI separation
// ~~~~~~~~~~~~~~~~
const bool frec_util::sel_KPI( Particle&p, const double cut_value, 
			       const atc_pid& selkpi ){

  bool select = false;

  const int h_lund = abs(p.pType().lund());

  if ( h_lund == Ptype_PIplus.lund() ){
    if ( selkpi.prob(p.mdstCharged()) <= cut_value ) select = true;
   }else if ( h_lund == Ptype_Kplus.lund() ){
    if( selkpi.prob(p.mdstCharged())  > cut_value )  select = true;
  }

  return select;

}

// beam related functions
// ~~~~~~~~~~~~~~~~~~~~~~~~
const double frec_util::Ebeam(){

   return m_p_beam.mag()/2.;

}

const HepLorentzVector frec_util::p_beam(){

   return m_p_beam;

}

const HepLorentzVector frec_util::p_cm( const HepLorentzVector& p ){

  const Hep3Vector CMboost = - m_p_beam.boostVector();

  HepLorentzVector m_p_cm = p;
  m_p_cm.boost( CMboost );
  return m_p_cm;

}

// other utilities
// ~~~~~~~~~~~~~~~~~~
const bool frec_util::HadronB(){
 
  bool flag_HadronB = false;

  Evtcls_hadronic_flag_Manager &evtclsMgr = 
             Evtcls_hadronic_flag_Manager::get_manager();

  if ( evtclsMgr.count() &&
       evtclsMgr.begin()->hadronic_flag(2) >= 10 ) 
    flag_HadronB = true;
 
 
  return flag_HadronB;

}

const double frec_util::R2( std::vector<Particle>& ALL ){

 std::vector<Hep3Vector> vlist;

  for( std::vector<Particle>::iterator i = ALL.begin();
       i != ALL.end(); i++){
      Hep3Vector pcm = frec_util::p_cm( i->p() ).vect();
      vlist.push_back( pcm );
  }

  FoxWolfram fw = foxwolfram( vlist.begin(), vlist.end(), SelfFunc(Hep3Vector()) );

  return fw.R(2);

}

void frec_util::division( Particle& B, std::vector<Particle>& ALL,
			    std::vector<HepLorentzVector>& pcm_bcand,
			    std::vector<HepLorentzVector>& pcm_other ){


    std::vector<Particle> b_child_ptcle;

    const int nchild = B.relation().nFinalStateParticles();

    for ( int i=0; i< nchild; i++ ){

      Particle child = B.relation().finalStateParticle(i);

      const HepLorentzVector child_p = frec_util::p_cm( child.p() );

      b_child_ptcle.push_back(child);
      pcm_bcand.push_back(child_p);

    }

    for( std::vector<Particle>::iterator o_i = ALL.begin();
 	 o_i != ALL.end(); o_i++ ){

      bool ptcle_flag = false;

      for( std::vector<Particle>::iterator b_j = b_child_ptcle.begin();
	   b_j != b_child_ptcle.end(); b_j++){

  	  const int chg_id2 = o_i->mdstCharged().get_ID();
	  const int chg_id1 = b_j->mdstCharged().get_ID();
	  const int gam_id1 = o_i->mdstGamma().get_ID();
	  const int gam_id2 = b_j->mdstGamma().get_ID();

	  if ( (chg_id1*chg_id2 != 0 ) && ( chg_id1 == chg_id2 ) ) 
	     ptcle_flag = true;
	  if ( (gam_id1*gam_id2 != 0 ) && ( gam_id1 == gam_id2 ) ) 
 	     ptcle_flag = true;
        
      } // loop end b_l

     if (ptcle_flag == false ){
       HepLorentzVector other_p = frec_util::p_cm( o_i->p() );
       pcm_other.push_back(other_p);
      }

    } // end loop o_i

    return;
}

void frec_util::setUserInfoB( Particle &p ){

  UserInfo_B x; x.set(p);

  return;
}

const double frec_util::Bpurity_hadronic_DX( int b_mode, int dst_mode, int d_mode ,
					     int dsst_mode, int ds_mode ){

 double val=0.;

 // B- --> D0pi-
 if ( b_mode == 101 && dst_mode == 0 ){
   if ( d_mode == 1 ){
     val = 0.8690350;
   }else if ( d_mode == 2 ){
     val = 0.6586160;
   }else if ( d_mode == 3 ){
     val = 0.6586010;
   }else if ( d_mode == 4 ){
     val = 0.7558300;
   }else if ( d_mode == 5 ){
     val = 0.6876350;
   }else if ( d_mode == 6 ){
     val = 0.3470780;
   }else if ( d_mode == 7 ){
     val = 0.7350520;
   }
 // B- --> D0rho-
 } else if ( b_mode == 103 && dst_mode == 0 ){
   if ( d_mode == 1 ){
     val = 0.5952910;
   }else if ( d_mode == 2 ){
     val = 0.3578390;
   }else if ( d_mode == 3 ){
     val = 0.3483190;
   }else if ( d_mode == 4 ){
     val = 0.4024700;
   }else if ( d_mode == 5 ){
     val = 0.3922770;
   }else if ( d_mode == 6 ){
     val = 0.1902670;
   }else if ( d_mode == 7 ){
     val = 0.3768050;
   }
 // B- --> D0a1-
 } else if ( b_mode == 105 && dst_mode == 0 ){
   if ( d_mode == 1 ){
     val = 0.4275780;
   }else if ( d_mode == 2 ){
     val = 0.000000;
   }else if ( d_mode == 3 ){
     val = 0.000000;
   }else if ( d_mode == 4 ){
     val = 0.000000;
   }else if ( d_mode == 5 ){
     val = 0.000000;
   }else if ( d_mode == 6 ){
     val = 0.000000;
   }else if ( d_mode == 7 ){
     val = 0.000000;
   }
 // B- --> D*0(D0pi0)pi-
 } else if ( b_mode == 102  && dst_mode == 1 ){
   if ( d_mode == 1 ){
     val = 0.9002910;
   }else if ( d_mode == 2 ){
     val = 0.5991010;
   }else if ( d_mode == 3 ){
     val = 0.5955250;
   }else if ( d_mode == 4 ){
     val = 0.7317170;
   }else if ( d_mode == 5 ){
     val = 0.6394860;
   }else if ( d_mode == 6 ){
     val = 0.3117250;
   }else if ( d_mode == 7 ){
     val = 0.7177770;
   }
 // B- --> D*0(D0pi0)rho-
 } else if ( b_mode == 104  && dst_mode == 1 ){
   if ( d_mode == 1 ){
     val = 0.6803420;
   }else if ( d_mode == 2 ){
     val = 0.4092550;
   }else if ( d_mode == 3 ){
     val = 0.4199640;
   }else if ( d_mode == 4 ){
     val = 0.4767610;
   }else if ( d_mode == 5 ){
     val = 0.4667280;
   }else if ( d_mode == 6 ){
     val = 0.1961840;
   }else if ( d_mode == 7 ){
     val = 0.5403430;
   }
 // B- --> D*0(D0pi0)a1-
 } else if ( b_mode == 106 && dst_mode == 1 ){
   if ( d_mode == 1 ){
     val = 0.5987170;
   }else if ( d_mode == 2 ){
     val = 0.2587050;
   }else if ( d_mode == 3 ){
     val = 0.2562690;
   }else if ( d_mode == 4 ){
     val = 0.3732290;
   }else if ( d_mode == 5 ){
     val = 0.3093630;
   }else if ( d_mode == 6 ){
     val = 0.000000;
   }else if ( d_mode == 7 ){
     val = 0.3562390;
   }
 // B- --> D*0(D0gamma)pi-
 } else if ( b_mode == 102  && dst_mode == 2 ){
   if ( d_mode == 1 ){
     val = 0.7204470;
   }else if ( d_mode == 2 ){
     val = 0.4518360;
   }else if ( d_mode == 3 ){
     val = 0.4446600;
   }else if ( d_mode == 4 ){
     val = 0.5507030;
   }else if ( d_mode == 5 ){
     val = 0.4665780;
   }else if ( d_mode == 6 ){
     val = 0.2396810;
   }else if ( d_mode == 7 ){
     val = 0.5121880;
   }
 // B- --> D*0(D0gamma)rho-
 } else if ( b_mode == 104  && dst_mode == 2 ){
   if ( d_mode == 1 ){
     val = 0.4204490;
   }else if ( d_mode == 2 ){
     val = 0.2716750;
   }else if ( d_mode == 3 ){
     val = 0.2543300;
   }else if ( d_mode == 4 ){
     val = 0.3210200;
   }else if ( d_mode == 5 ){
     val = 0.2842360;
   }else if ( d_mode == 6 ){
     val = 0.1612390;
   }else if ( d_mode == 7 ){
     val = 0.3290170;
   }
 // B- --> D*0(D0gamma)a1-
 } else if ( b_mode == 106 && dst_mode == 2 ){
   if ( d_mode == 1 ){
     val = 0.3234090;
   }else if ( d_mode == 2 ){
     val = 0.000000;
   }else if ( d_mode == 3 ){
     val = 0.000000;
   }else if ( d_mode == 4 ){
     val = 0.000000;
   }else if ( d_mode == 5 ){
     val = 0.000000;
   }else if ( d_mode == 6 ){
     val = 0.000000;
   }else if ( d_mode == 7 ){
     val = 0.000000;
   }
 // B0B --> D+pi-
 } else if ( b_mode == 1 && dst_mode == 0 ){
   if ( d_mode == 101 ){
     val = 0.7952350;
    }else if ( d_mode == 102 ){
     val = 0.2573410;
    }else if ( d_mode == 103 ){
     val = 0.8874780;
    }else if ( d_mode == 104 ){
     val = 0.4868740;
    }else if ( d_mode == 105 ){
     val = 0.4492920;
    }else if ( d_mode == 106 ){
     val = 0.4793050;
    }
 // B0B --> D+rho-
 } else if ( b_mode == 3 && dst_mode == 0 ){
   if ( d_mode == 101 ){
     val = 0.5499100;
    }else if ( d_mode == 102 ){
     val = 0.1713570;
    }else if ( d_mode == 103 ){
     val = 0.6605670;
    }else if ( d_mode == 104 ){
     val = 0.2898210;
    }else if ( d_mode == 105 ){
     val = 0.2770830;
    }else if ( d_mode == 106 ){
     val = 0.2734190;
    }
 // B0B --> D+a1-
 } else if ( b_mode == 5 && dst_mode == 0 ){
   if ( d_mode == 101 ){
     val = 0.3381760;
    }else if ( d_mode == 102 ){
     val = 0.000000;
    }else if ( d_mode == 103 ){
     val = 0.4240310;
    }else if ( d_mode == 104 ){
     val = 0.000000;
    }else if ( d_mode == 105 ){
     val = 0.000000;
    }else if ( d_mode == 106 ){
     val = 0.000000;
    }
 // B0B --> D*+(D0pi+)pi-
 } else if ( b_mode == 2 && dst_mode == 3 ){
   if ( d_mode == 1 ){
     val = 0.9559230;
   }else if ( d_mode == 2 ){
     val = 0.8037760;
   }else if ( d_mode == 3 ){
     val = 0.8096070;
   }else if ( d_mode == 4 ){
     val = 0.8963900;
   }else if ( d_mode == 5 ){
     val = 0.8605940;
   }else if ( d_mode == 6 ){
     val = 0.5391560;
   }else if ( d_mode == 7 ){
     val = 0.8860570;
   }
 // B0B --> D*+(D0pi+)rho-
 } else if ( b_mode == 4 && dst_mode == 3 ){
   if ( d_mode == 1 ){
     val = 0.8121750;
   }else if ( d_mode == 2 ){
     val = 0.6077220;
   }else if ( d_mode == 3 ){
     val = 0.6331800;
   }else if ( d_mode == 4 ){
     val = 0.7437160;
   }else if ( d_mode == 5 ){
     val = 0.7010010;
   }else if ( d_mode == 6 ){
     val = 0.3713530;
   }else if ( d_mode == 7 ){
     val = 0.7271300;
   }
 // B0B --> D*+(D0pi+)a1-
 } else if ( b_mode == 6 && dst_mode == 3 ){
   if ( d_mode == 1 ){
     val = 0.7381600;
   }else if ( d_mode == 2 ){
     val = 0.4885030;
   }else if ( d_mode == 3 ){
     val = 0.4805960;
   }else if ( d_mode == 4 ){
     val = 0.6223900;
   }else if ( d_mode == 5 ){
     val = 0.5350140;
   }else if ( d_mode == 6 ){
     val = 0.2486260;
   }else if ( d_mode == 7 ){
     val = 0.6615760;
   }
 // B0B --> D*+(D+pi0)pi-
 } else if ( b_mode == 2 && dst_mode == 4 ){
   if ( d_mode == 101 ){
     val = 0.6971820;
    }else if ( d_mode == 102 ){
     val = 0.2064640;
    }else if ( d_mode == 103 ){
     val = 0.8477080;
    }else if ( d_mode == 104 ){
     val = 0.3821190;
    }else if ( d_mode == 105 ){
     val = 0.3829890;
    }else if ( d_mode == 106 ){
     val = 0.3838650;
    }
 // B0B --> D*+(D+pi0)rho-
 } else if ( b_mode == 4 && dst_mode == 4 ){
   if ( d_mode == 101 ){
     val = 0.4474520;
    }else if ( d_mode == 102 ){
     val = 0.1492320;
    }else if ( d_mode == 103 ){
     val = 0.5956610;
    }else if ( d_mode == 104 ){
     val = 0.2499630;
    }else if ( d_mode == 105 ){
     val = 0.2001490;
    }else if ( d_mode == 106 ){
     val = 0.2452550;
    }
 // B0B --> D*+(D+pi0)a1-
 } else if ( b_mode == 6 && dst_mode == 4 ){
   if ( d_mode == 101 ){
     val = 0.3140260;
    }else if ( d_mode == 102 ){
     val = 0.000000;
    }else if ( d_mode == 103 ){
     val = 0.4243430;
    }else if ( d_mode == 104 ){
     val = 0.000000;
    }else if ( d_mode == 105 ){
     val = 0.000000;
    }else if ( d_mode == 106 ){
     val = 0.000000;
    }
 // B- --> D0Ds(KsK-)
 } else if ( b_mode == 111 && dst_mode == 0 && ds_mode == 301 && dsst_mode == 0 ){
   if ( d_mode == 1 ){
     val = 0.9405700;
   }else if ( d_mode == 2 ){
     val = 0.5635030;
   }else if ( d_mode == 3 ){
     val = 0.4841160;
   }else if ( d_mode == 4 ){
     val = 0.6808310;
   }else if ( d_mode == 5 ){
     val = 0.4894930;
   }else if ( d_mode == 6 ){
     val = 0.1855170;
   }else if ( d_mode == 7 ){
     val = 0.8742790;
   }
 // B- --> D0Ds(K+K-pi-)
 } else if ( b_mode == 111 && dst_mode == 0 && ds_mode == 302 && dsst_mode == 0 ){
   if ( d_mode == 1 ){
     val = 0.6388630;
   }else if ( d_mode == 2 ){
     val = 0.5547050;
   }else if ( d_mode == 3 ){
     val = 0.6097350;
   }else if ( d_mode == 4 ){
     val = 0.7763500;
   }else if ( d_mode == 5 ){
     val = 0.6825760;
   }else if ( d_mode == 6 ){
     val = 0.2217080;
   }else if ( d_mode == 7 ){
     val = 0.8901180;
   }
 // B- --> D*0(D0pi0)Ds(KsK-)
 } else if ( b_mode == 113  && dst_mode == 1 && ds_mode == 301 && dsst_mode == 0 ){
   if ( d_mode == 1 ){
     val = 0.6884170;
   }else if ( d_mode == 2 ){
     val = 0.5070850;
   }else if ( d_mode == 3 ){
     val = 0.5859120;
   }else if ( d_mode == 4 ){
     val = 0.6404450;
   }else if ( d_mode == 5 ){
     val = 0.5857500;
   }else if ( d_mode == 6 ){
     val = 0.2272760;
   }else if ( d_mode == 7 ){
     val = 0.8855580;
   }
 // B- --> D*0(D0pi0)Ds(K+K-pi-)
 } else if ( b_mode == 113  && dst_mode == 1 && ds_mode == 302 && dsst_mode == 0 ){
   if ( d_mode == 1 ){
     val = 0.6865980;
   }else if ( d_mode == 2 ){
     val = 0.3304040;
   }else if ( d_mode == 3 ){
     val = 0.2732160;
   }else if ( d_mode == 4 ){
     val = 0.4674510;
   }else if ( d_mode == 5 ){
     val = 0.4693790;
   }else if ( d_mode == 6 ){
     val = 0.2784570;
   }else if ( d_mode == 7 ){
     val = 0.4415400;
   }
 // B- --> D*0(D0gamma)Ds(KsK-)
 } else if ( b_mode == 113  && dst_mode == 2 && ds_mode == 301 && dsst_mode == 0 ){
   if ( d_mode == 1 ){
     val = 0.5350850;
   }else if ( d_mode == 2 ){
     val = 0.3172840;
   }else if ( d_mode == 3 ){
     val = 0.4150000;
   }else if ( d_mode == 4 ){
     val = 0.1;
   }else if ( d_mode == 5 ){
     val = 0.3237640;
   }else if ( d_mode == 6 ){
     val = 0.1;
   }else if ( d_mode == 7 ){
     val = 0.3364320;
   }
 // B- --> D*0(D0gamma)Ds(K+K-pi-)
 } else if ( b_mode == 113  && dst_mode == 2 && ds_mode == 302 && dsst_mode == 0 ){
   if ( d_mode == 1 ){
     val = 0.4476510;
   }else if ( d_mode == 2 ){
     val = 0.4032250;
   }else if ( d_mode == 3 ){
     val = 0.3659890;
   }else if ( d_mode == 4 ){
     val = 0.2982200;
   }else if ( d_mode == 5 ){
     val = 0.1055060;
   }else if ( d_mode == 6 ){
     val = 0.1342990;
   }else if ( d_mode == 7 ){
     val = 0.6244940;
   }
 // B0B --> D+Ds(KsK-)
 } else if ( b_mode == 11 && dst_mode == 0 && ds_mode == 301 && dsst_mode == 0 ){
   if ( d_mode == 101 ){
     val = 0.7972340;
    }else if ( d_mode == 102 ){
     val = 0.3015670;
    }else if ( d_mode == 103 ){
     val = 0.9628310;
    }else if ( d_mode == 104 ){
     val = 0.4915600;
    }else if ( d_mode == 105 ){
     val = 0.4012080;
    }else if ( d_mode == 106 ){
     val = 0.3483060;
    }
 // B0B --> D+Ds(K+K-pi-)
 } else if ( b_mode == 11 && dst_mode == 0 && ds_mode == 302 && dsst_mode == 0 ){
   if ( d_mode == 101 ){
     val = 0.4782100;
    }else if ( d_mode == 102 ){
     val = 0.2532120;
    }else if ( d_mode == 103 ){
     val = 0.6066900;
    }else if ( d_mode == 104 ){
     val = 0.5037810;
    }else if ( d_mode == 105 ){
     val = 0.3932110;
    }else if ( d_mode == 106 ){
     val = 0.4335030;
    }
 // B0B --> D*+(D0pi+)Ds(KsK-)
 } else if ( b_mode == 13 && dst_mode == 3 && ds_mode == 301 && dsst_mode == 0 ){
   if ( d_mode == 1 ){
     val = 0.9563140;
   }else if ( d_mode == 2 ){
     val = 0.7163940;
   }else if ( d_mode == 3 ){
     val = 0.8542810;
   }else if ( d_mode == 4 ){
     val = 0.9999850;
   }else if ( d_mode == 5 ){
     val = 0.9781170;
   }else if ( d_mode == 6 ){
     val = 0.1;
   }else if ( d_mode == 7 ){
     val = 0.9212280;
   }
 // B0B --> D*+(D0pi+)Ds(K+K-pi-)
 } else if ( b_mode == 13 && dst_mode == 3 && ds_mode == 302 && dsst_mode == 0 ){
   if ( d_mode == 1 ){
     val = 0.8741820;
   }else if ( d_mode == 2 ){
     val = 0.6522100;
   }else if ( d_mode == 3 ){
     val = 0.5578570;
   }else if ( d_mode == 4 ){
     val = 0.7509200;
   }else if ( d_mode == 5 ){
     val = 0.6098800;
   }else if ( d_mode == 6 ){
     val = 0.2693480;
   }else if ( d_mode == 7 ){
     val = 0.7979430;
   }
 // B0B --> D*+(D+pi0)Ds(KsK-)
 } else if ( b_mode == 13 && dst_mode == 4 && ds_mode == 301 && dsst_mode == 0 ){
   if ( d_mode == 101 ){
     val = 0.6321320;
    }else if ( d_mode == 102 ){
     val = 0.1983240;
    }else if ( d_mode == 103 ){
     val = 0.9622400;
    }else if ( d_mode == 104 ){
     val = 0.3959700;
    }else if ( d_mode == 105 ){
     val = 0.2847940;
    }else if ( d_mode == 106 ){
     val = 0.6034400;
    }
 // B0B --> D*+(D+pi0)Ds(K+K-pi-)
 } else if ( b_mode == 13 && dst_mode == 4 && ds_mode == 302 && dsst_mode == 0 ){
   if ( d_mode == 101 ){
     val = 0.3934260;
    }else if ( d_mode == 102 ){
     val = 0.1499030;
    }else if ( d_mode == 103 ){
     val = 0.5730740;
    }else if ( d_mode == 104 ){
     val = 0.3398020;
    }else if ( d_mode == 105 ){
     val = 0.2515470;
    }else if ( d_mode == 106 ){
     val = 0.1;
    }
 // B- --> D0Ds*(KsK-)
 } else if ( b_mode == 112 && dst_mode == 0 && ds_mode == 301 && dsst_mode == 5 ){
   if ( d_mode == 1 ){
     val = 0.5990210;
   }else if ( d_mode == 2 ){
     val = 0.4047470;
   }else if ( d_mode == 3 ){
     val = 0.5010490;
   }else if ( d_mode == 4 ){
     val = 0.3097370;
   }else if ( d_mode == 5 ){
     val = 0.3871640;
   }else if ( d_mode == 6 ){
     val = 0.2043220;
   }else if ( d_mode == 7 ){
     val = 0.4503110;
   }
 // B- --> D0Ds*(K+K-pi-)
 } else if ( b_mode == 112 && dst_mode == 0 && ds_mode == 302 && dsst_mode == 5 ){
   if ( d_mode == 1 ){
     val = 0.5425770;
   }else if ( d_mode == 2 ){
     val = 0.4022400;
   }else if ( d_mode == 3 ){
     val = 0.4015190;
   }else if ( d_mode == 4 ){
     val = 0.5158320;
   }else if ( d_mode == 5 ){
     val = 0.4715510;
   }else if ( d_mode == 6 ){
     val = 0.2203000;
   }else if ( d_mode == 7 ){
     val = 0.3966250;
   }
 // B- --> D*0(D0pi0)Ds*(KsK-)
 } else if ( b_mode == 114  && dst_mode == 1 && ds_mode == 301 && dsst_mode == 5 ){
   if ( d_mode == 1 ){
     val = 0.7361600;
   }else if ( d_mode == 2 ){
     val = 0.4953760;
   }else if ( d_mode == 3 ){
     val = 0.4975020;
   }else if ( d_mode == 4 ){
     val = 0.5063600;
   }else if ( d_mode == 5 ){
     val = 0.5200950;
   }else if ( d_mode == 6 ){
     val = 0.1749270;
   }else if ( d_mode == 7 ){
     val = 0.7288820;
   }
 // B- --> D*0(D0pi0)Ds*(K+K-pi-)
 } else if ( b_mode == 114  && dst_mode == 1 && ds_mode == 302 && dsst_mode == 5 ){
   if ( d_mode == 1 ){
     val = 0.6865980;
   }else if ( d_mode == 2 ){
     val = 0.3304040;
   }else if ( d_mode == 3 ){
     val = 0.2732160;
   }else if ( d_mode == 4 ){
     val = 0.4674510;
   }else if ( d_mode == 5 ){
     val = 0.4693790;
   }else if ( d_mode == 6 ){
     val = 0.2784570;
   }else if ( d_mode == 7 ){
     val = 0.4415400;
   }
 // B- --> D*0(D0gamma)Ds*(KsK-)
 } else if ( b_mode == 114  && dst_mode == 2 && ds_mode == 301 && dsst_mode == 5 ){
   if ( d_mode == 1 ){
     val = 0.4024090;
   }else if ( d_mode == 2 ){
     val = 0.3054990;
   }else if ( d_mode == 3 ){
     val = 0.3352910;
   }else if ( d_mode == 4 ){
     val = 0.2836480;
   }else if ( d_mode == 5 ){
     val = 0.4127720;
   }else if ( d_mode == 6 ){
     val = 0.1487900;
   }else if ( d_mode == 7 ){
     val = 0.1;
   }
 // B- --> D*0(D0gamma)Ds*(K+K-pi-)
 } else if ( b_mode == 114  && dst_mode == 2 && ds_mode == 302 && dsst_mode == 5 ){
   if ( d_mode == 1 ){
     val = 0.4191170;
   }else if ( d_mode == 2 ){
     val = 0.3279280;
   }else if ( d_mode == 3 ){
     val = 0.3031690;
   }else if ( d_mode == 4 ){
     val = 0.4053630;
   }else if ( d_mode == 5 ){
     val = 0.3062710;
   }else if ( d_mode == 6 ){
     val = 0.1764860;
   }else if ( d_mode == 7 ){
     val = 0.5440060;
   }
 // B0B --> D+Ds*(KsK-)
 } else if ( b_mode == 12 && dst_mode == 0 && ds_mode == 301 && dsst_mode == 5 ){
   if ( d_mode == 101 ){
     val = 0.5381500;
    }else if ( d_mode == 102 ){
     val = 0.1783790;
    }else if ( d_mode == 103 ){
     val = 0.3958610;
    }else if ( d_mode == 104 ){
     val = 0.3921310;
    }else if ( d_mode == 105 ){
     val = 0.3824790;
    }else if ( d_mode == 106 ){
     val = 0.3207780;
    }
 // B0B --> D+Ds*(K+K-pi-)
 } else if ( b_mode == 12 && dst_mode == 0 && ds_mode == 302 && dsst_mode == 5 ){
   if ( d_mode == 101 ){
     val = 0.3269810;
    }else if ( d_mode == 102 ){
     val = 0.1827030;
    }else if ( d_mode == 103 ){
     val = 0.4295950;
    }else if ( d_mode == 104 ){
     val = 0.3626490;
    }else if ( d_mode == 105 ){
     val = 0.3395500;
    }else if ( d_mode == 106 ){
     val = 0.3326640;
    }
 // B0B --> D*+(D0pi+)Ds*(KsK-)
 } else if ( b_mode == 14 && dst_mode == 3 && ds_mode == 301 && dsst_mode == 5 ){
   if ( d_mode == 1 ){
     val = 0.9456130;
   }else if ( d_mode == 2 ){
     val = 0.7461420;
   }else if ( d_mode == 3 ){
     val = 0.8097560;
   }else if ( d_mode == 4 ){
     val = 0.9261800;
   }else if ( d_mode == 5 ){
     val = 0.9938560;
   }else if ( d_mode == 6 ){
     val = 0.1;
   }else if ( d_mode == 7 ){
     val = 0.6394380;
   }
 // B0B --> D*+(D0pi+)Ds*(K+K-pi-)
 } else if ( b_mode == 14 && dst_mode == 3 && ds_mode == 302 && dsst_mode == 5 ){
   if ( d_mode == 1 ){
     val = 0.7764520;
   }else if ( d_mode == 2 ){
     val = 0.5747060;
   }else if ( d_mode == 3 ){
     val = 0.5514500;
   }else if ( d_mode == 4 ){
     val = 0.7561710;
   }else if ( d_mode == 5 ){
     val = 0.7377220;
   }else if ( d_mode == 6 ){
     val = 0.2757920;
   }else if ( d_mode == 7 ){
     val = 0.8292380;
   }
 // B0B --> D*+(D+pi0)Ds*(KsK-)
 } else if ( b_mode == 14 && dst_mode == 4 && ds_mode == 301 && dsst_mode == 5 ){
   if ( d_mode == 101 ){
     val = 0.8647090;
    }else if ( d_mode == 102 ){
     val = 0.2114160;
    }else if ( d_mode == 103 ){
     val = 0.9389900;
    }else if ( d_mode == 104 ){
     val = 0.3462740;
    }else if ( d_mode == 105 ){
     val = 0.1390550;
    }else if ( d_mode == 106 ){
     val = 0.3652850;
    }
 // B0B --> D*+(D+pi0)Ds*(K+K-pi-)
 } else if ( b_mode == 14 && dst_mode == 4 && ds_mode == 302 && dsst_mode == 5 ){
   if ( d_mode == 101 ){
     val = 0.3714510;
    }else if ( d_mode == 102 ){
     val = 0.1718110;
    }else if ( d_mode == 103 ){
     val = 0.3947190;
    }else if ( d_mode == 104 ){
     val = 0.3056260;
    }else if ( d_mode == 105 ){
     val = 0.5800770;
    }else if ( d_mode == 106 ){
     val = 0.3503170;
    }
}

 return val;

}

// ------------------------
//     thrust  class
// ------------------------

// constructors
// ~~~~~~~~~~~~~~~
thrust::thrust( HepAList<HepLorentzVector>& plist ){

  m_status = 0;

  int np   = 0;
  int len  = 4;
  int ist0 = 0;
  float pl[100][4];

  int ntrk = plist.length();
  if( ntrk < 2 ) {
    if(ntrk==0){
      m_status=0;
      m_thr=1;
      m_thr_axis=Hep3Vector(0,0,0);

      // I'm not sure following values are proper
      // they are not used, anyway
      m_sph_axis=Hep3Vector(0,0,0);
      m_obl=0;
      m_sph=0;
      m_apl=0;
    }
    if(ntrk==1){
      m_status=1;
      m_thr=1;
      m_thr_axis=plist[0]->vect().unit();

      // I'm not sure following values are proper
      // they are not used, anyway
      m_sph_axis=plist[0]->vect().unit();
      m_obl=0;
      m_sph=0;
      m_apl=0;
    }
    
    return;
  }

  for(int i = 0; i < ntrk; i++){
    pl[np][0]=plist[i]->x();
    pl[np][1]=plist[i]->y();
    pl[np][2]=plist[i]->z();
    pl[np][3]=plist[i]->t();
    np++;
  }

  float thr0, obl0;
  float eig1[3][4];
  FORTRAN_ROUTINE(thrust)(&np,pl,&len,&ist0,&thr0,&obl0,eig1);
  if( thr0 == -1.0 || obl0 == -1.0 ){ return; }
  Hep3Vector thr_axis0(eig1[0][1],eig1[0][2],eig1[0][3]);
  float sph0, apl0;
  float eig2[3][4];
  FORTRAN_ROUTINE(spher)(&np,pl,&len,&ist0,&sph0,&apl0,eig2);
  Hep3Vector sph_axis0(eig2[0][1],eig2[0][2],eig2[0][3]);

  m_thr_axis = thr_axis0;
  m_sph_axis = sph_axis0;
  m_thr  = thr0;
  m_obl  = obl0;
  m_sph  = sph0;
  m_apl  = apl0;
  m_status = 1;
 
}

thrust::thrust(std::vector<HepLorentzVector>& ptcle_p){

   m_status = 0;

  int np   = 0;
  int len  = 4;
  int ist0 = 0;
  int ntrk = ptcle_p.size();
  float pl[1000][4];

  if(ntrk > 100) dout(Debugout::WARN,"frec_util") << "ptcle over>> " << ntrk << std::endl;

  if(ntrk >= 1000) {
    dout(Debugout::ERR,"frec_util") << "ptcle over returning>> " << ntrk << std::endl;
    return;
  }
  
  if( ntrk < 2 ) {
    if(ntrk==0){
      m_status=0;
      m_thr=1;
      m_thr_axis=Hep3Vector(0,0,0);

      // I'm not sure following values are proper
      // they are not used, anyway
      m_sph_axis=Hep3Vector(0,0,0);
      m_obl=0;
      m_sph=0;
      m_apl=0;
    }
    if(ntrk==1){
      m_status=1;
      m_thr=1;
      m_thr_axis=ptcle_p[0].vect().unit();

      // I'm not sure following values are proper
      // they are not used, anyway
      m_sph_axis=ptcle_p[0].vect().unit();
      m_obl=0;
      m_sph=0;
      m_apl=0;
    }

    return;
  }

  for(std::vector<HepLorentzVector>::const_iterator i = ptcle_p.begin(); 
      i != ptcle_p.end(); i++) {
    pl[np][0] = i->px();
    pl[np][1] = i->py();
    pl[np][2] = i->pz();
    pl[np][3] = i->e();
    np++;
  }

  float thr0, obl0;
  float eig1[3][4];
  FORTRAN_ROUTINE(thrust)(&np,pl,&len,&ist0,&thr0,&obl0,eig1);
  if( thr0 == -1.0 || obl0 == -1.0 ){ return; }
  Hep3Vector thr_axis0(eig1[0][1],eig1[0][2],eig1[0][3]);
  float sph0, apl0;
  float eig2[3][4];
  FORTRAN_ROUTINE(spher)(&np,pl,&len,&ist0,&sph0,&apl0,eig2);
  Hep3Vector sph_axis0(eig2[0][1],eig2[0][2],eig2[0][3]);

  m_thr_axis = thr_axis0;
  m_sph_axis = sph_axis0;
  m_thr  = thr0;
  m_obl  = obl0;
  m_sph  = sph0;
  m_apl  = apl0;
  m_status = 1;
 
} 


// ------------------------
//    UserInfo_B class
// ------------------------

// default constructor
// ~~~~~~~~~~~~~~~~~~~~~~~
UserInfo_B::UserInfo_B() {

  m_b_mode        = 0;
  m_sub0_mode     = 0;
  m_sub1_mode     = 0;
  m_sub2_mode     = 0;
  m_sub3_mode     = 0;

  m_Mass_B        = -5.;
  m_Delta_E       = -5.;
  m_cosT          = -5.;
  m_cosB          = -5.;
  m_mass          = -5.;
  m_masschisq     = -5.;
  m_vtxchisq      = -5.;
  m_purity        = -5.;

  m_flag_bestEach = 0;
  m_flag_best     = 0;
  m_flag_best2    = 0;
  m_flag_best3    = 0;
  m_flag_mc       = 0;

}

// copy constructor
// ~~~~~~~~~~~~~~~~~~~
UserInfo_B::UserInfo_B(const UserInfo_B& x) {

  m_b_mode        = x.m_b_mode;
  m_sub0_mode     = x.m_sub0_mode;
  m_sub1_mode     = x.m_sub1_mode;
  m_sub2_mode     = x.m_sub2_mode;
  m_sub3_mode     = x.m_sub3_mode;

  m_Mass_B        = x.m_Mass_B;
  m_Delta_E       = x.m_Delta_E;
  m_cosT          = x.m_cosT;
  m_cosB          = x.m_cosB;
  m_mass          = x.m_mass;
  m_masschisq     = x.m_masschisq;
  m_vtxchisq      = x.m_vtxchisq;
  m_purity        = x.m_purity;

  m_flag_bestEach = x.m_flag_bestEach;
  m_flag_best     = x.m_flag_best;
  m_flag_best2    = x.m_flag_best2;
  m_flag_best3    = x.m_flag_best3;
  m_flag_mc       = x.m_flag_mc;

}

// construct self object
// ~~~~~~~~~~~~~~~~~~~~~~~~
UserInfo_B* UserInfo_B::clone(void) const {

  UserInfo_B *x = new UserInfo_B(*this);
  return x;

}

// operator
// ~~~~~~~~~~~~
UserInfo_B& UserInfo_B::operator = (const UserInfo_B& x){

  m_b_mode        = x.m_b_mode;
  m_sub0_mode     = x.m_sub0_mode;
  m_sub1_mode     = x.m_sub1_mode;
  m_sub2_mode     = x.m_sub2_mode;
  m_sub3_mode     = x.m_sub3_mode;

  m_Mass_B        = x.m_Mass_B;
  m_Delta_E       = x.m_Delta_E;
  m_cosT          = x.m_cosT;
  m_cosB          = x.m_cosB;
  m_mass          = x.m_mass;
  m_masschisq     = x.m_masschisq;
  m_vtxchisq      = x.m_vtxchisq;
  m_purity        = x.m_purity;

  m_flag_bestEach = x.m_flag_bestEach;
  m_flag_best     = x.m_flag_best;
  m_flag_best2    = x.m_flag_best2;
  m_flag_best3    = x.m_flag_best3;
  m_flag_mc       = x.m_flag_mc;

  return *this;

}

// set Particle UserInfo
// ~~~~~~~~~~~~~~~~~~~~~~~~~
void UserInfo_B::set(std::vector<Particle>& B){

  for(std::vector<Particle>::iterator i = B.begin(); i != B.end(); i++)  set(*i);

  return;

}

void UserInfo_B::set(Particle& B){

  const UserInfo_B info;
  B.userInfo(info);

  return;

}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
