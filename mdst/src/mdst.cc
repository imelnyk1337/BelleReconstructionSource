// 1999/5/28 skkim
// added good_gamma
//
// $Id: mdst.cc 10712 2009-01-11 16:18:37Z hitoshi $
//
// $Log$
// Revision 1.43  2003/04/15 02:48:16  katayama
// pi0resol added
//
// Revision 1.42  2002/03/14 01:19:32  katayama
// correction functions are moved to fix_mdst
//
// Revision 1.41  2002/02/27 01:39:26  katayama
// uninitialized variables/removal of unused declarations.
//
// Revision 1.40  2002/02/25 10:01:16  hitoshi
// added an option for MC 20011214_0817 (by J.Tanaka).
//
// Revision 1.39  2002/02/25 01:38:00  katayama
// Bug fix(no effect)
//
// Revision 1.38  2002/02/22 03:26:10  hitoshi
// updated Benergy further using higher statstics (by Sanjay).
//
// Revision 1.37  2002/02/21 00:47:10  hitoshi
// updated correct_ecl (by Miyabayashi).
//
// Revision 1.36  2002/02/20 01:32:19  katayama
// std::
//
// Revision 1.35  2002/02/16 02:48:18  hitoshi
// updated Benergy for e15 (by Sanjay).
//
// Revision 1.34  2002/02/15 03:22:33  hitoshi
// added scale error for e15 (by J.Tanaka).
//
// Revision 1.33  2002/02/14 09:28:29  hitoshi
// added factors for exp15 (mode=2).
//
// Revision 1.32  2001/12/25 18:55:02  hitoshi
// added mdst_good_event function.
//
// Revision 1.31  2001/12/25 11:12:30  katayama
// more cleaning up
//
// Revision 1.30  2001/12/15 06:39:22  hitoshi
// updated tracing back to hepevt (by Kakuno).
//
// Revision 1.29  2001/12/12 02:04:10  hitoshi
// updated for mdst_ecl/trk (by Kakuno).
//
// Revision 1.28  2001/12/05 05:51:24  katayama
//  for CC
//
// Revision 1.27  2001/12/05 04:56:56  hitoshi
// commented out redundant ifndef's in scale_error (requested by J.Tanaka).
//
// Revision 1.26  2001/12/04 12:07:26  hitoshi
// added scale_error, scale gamma energy, benergy.
//
// Revision 1.25  2001/11/29 18:34:24  katayama
// CC compatibility
//
// Revision 1.24  2001/11/05 13:44:19  hitoshi
// updated factor for e13 for mmode=1
//
// Revision 1.23  2001/11/05 12:39:44  hitoshi
// updated so that it can handle factor per 100 runs.
//
// Revision 1.22  2000/12/10 03:19:16  hitoshi
// added a function remove_duplicates.
//
// Revision 1.21  2000/11/17 15:32:21  hitoshi
// Fixed a bug. So far this couldn't handle old MDST w/o mdst_event_add table.
//
// Revision 1.20  2000/10/21 14:33:33  hitoshi
// turn on scale flag when scale_momenta is called with scale !=1. Save scale
// factor.
//
// Revision 1.19  2000/10/11 18:44:45  hitoshi
// added protection not to scale momentum multiplly.
//
// Revision 1.18  2000/05/26 09:07:00  hitoshi
// added check of existence of belle_event.
//
// Revision 1.17  2000/05/19 09:26:52  hitoshi
// added a function to copy vee to vee2.
//
// Revision 1.16  2000/05/16 12:06:54  hitoshi
// just cosmetic change.
//
// Revision 1.15  2000/05/02 10:09:59  hitoshi
// added momentum scaling func.
//
// Revision 1.14  2000/03/07 11:10:25  katayama
// compatibility with CC5.0
//
// Revision 1.13  2000/02/08 07:50:12  katayama
// Added Fang Fang's Ks finder for now
//
// Revision 1.12  1999/07/13 03:16:29  katayama
// added rectrk->mdst_charged (temporary)
//
// Revision 1.11  1999/07/10 12:17:45  katayama
// added trk->charged
//
// Revision 1.10  1999/05/29 05:05:42  katayama
// from skim san
//
// Revision 1.9  1999/04/27 23:45:15  katayama
// Use const keyword
//
// Revision 1.8  1999/03/17 07:33:52  hitoshi
// do not calculate prob when ndf=0.
//
// Revision 1.7  1999/01/16 10:31:12  katayama
// clean up includes
//
// Revision 1.6  1999/01/13 15:25:41  hitoshi
// implemented good trk selection logic.
//
// Revision 1.5  1999/01/13 00:15:59  katayama
// added good_charged
//
// Revision 1.4  1998/10/01 06:44:14  katayama
// Keep up with rec2mdst changes
//
// Revision 1.3  1997/09/19 00:49:43  katayama
// Added Id and Log
//
//
#include "belleCLHEP/Vector/ThreeVector.h"
#include "belleCLHEP/Matrix/Vector.h"
#include "helix/Helix.h"

#include "mdst/mdst.h"
#include "fix_mdst/fix_mdst.h"
#include <cmath>
#include <cfloat>

#include "belle.h"
#include MDST_H
#include MDST_OBS_H
#include HEPEVT_H
#include TRK_H

#include "event/BelleEvent.h"
#include BELLETDF_H

#include "belleCLHEP/Matrix/Matrix.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  static const int PDG_PI0(111);
  static const int PDG_K0S(310);
  static const int PDG_LAM(3122);
  static const int PDG_ALAM(-3122);
  static const int PDG_GAMM(22);
  

extern "C" {
  float prob_(const float&, const int&);
}

// trace back to isthep>=0 or isthep=-10(decay-in-flight)
static const Gen_hepevt &gen_level_step1(const Gen_hepevt & gen){
  return (gen && gen.mother() && gen.isthep()<0 && gen.isthep()!=-10) 
    ? gen_level_step1(gen.mother()) : gen;
}

// trace back to isthep>0 level
static const Gen_hepevt &gen_level_step2(const Gen_hepevt & gen){
  return (gen && gen.mother() && gen.isthep()<0)
    ? gen_level_step2(gen.mother()) : gen;
}

const Gen_hepevt &gen_level(const Gen_hepevt & gen){
  const Gen_hepevt &level1(gen_level_step1(gen));
  if (!level1 || level1.isthep() >= 0){
    return level1;
  }else if (level1.isthep() == -10){
    const Gen_hepevt &level2(gen_level_step2(gen));
    // If "gen" is made by Ks(which is made by generator) daughter, returns daughter
    if (level2.idhep() == 310 || abs(level2.idhep()) == 3122) return level1;
    // Otherwize, trace back to generator level
    else return level2;
  }
  return level1;
}

const Gen_hepevt &get_hepevt(const Mdst_trk&trk, int i) {
  Mdst_sim_trk_Manager & xrefmgr = Mdst_sim_trk_Manager::get_manager();
  Mdst_sim_trk_Index index(xrefmgr.index( "trk" ));
  index.update();
  Panther_ID trackID(trk.get_ID());
  std::vector<Mdst_sim_trk> point = point_from ( trackID, index );
  if(i>=0 && i<(int)point.size())  return point[i].hepevt();
  else return Gen_hepevt_Manager::get_manager().get_NULL();
}

const Gen_hepevt &get_hepevt(const Mdst_charged&ch, int i) {
  if (!Mdst_sim_trk_Manager::get_manager().count()){
    // for backward compatibility
    Mdst_sim_xref_Manager & xrefmgr = Mdst_sim_xref_Manager::get_manager();
    Mdst_sim_xref_Index index(xrefmgr.index( "charged" ));
    index.update();
    Panther_ID trackID(ch.get_ID());
    std::vector<Mdst_sim_xref> point = point_from ( trackID, index );
    if(i>=0 && i<(int)point.size())  return point[i].hepevt();
    else return Gen_hepevt_Manager::get_manager().get_NULL();
  }
  return get_hepevt(ch.trk(),i);
}

const Gen_hepevt &get_hepevt(const Mdst_ecl&ecl, int i) {
  Mdst_sim_ecl_Manager & xrefmgr = Mdst_sim_ecl_Manager::get_manager();
  Mdst_sim_ecl_Index index(xrefmgr.index( "ecl" ));
  index.update();
  Panther_ID eclID(ecl.get_ID());
  std::vector<Mdst_sim_ecl> point = point_from ( eclID, index );
  if(i>=0 && i<(int)point.size())  return point[i].hepevt();
  else return Gen_hepevt_Manager::get_manager().get_NULL();
}

const Gen_hepevt &get_hepevt(const Mdst_gamma&gamma, int i) {
  if (!Mdst_sim_ecl_Manager::get_manager().count()){
    // do NOT support backward compatibility
    // because the returned "Gen_hepevt & " may not be correct
    return Gen_hepevt_Manager::get_manager().get_NULL();
  }
  return get_hepevt(gamma.ecl(),i);
}

const Gen_hepevt &get_hepevt(const Mdst_pi0 &pi0, int) {
  const Gen_hepevt &hepevt1 = gen_level(get_hepevt(pi0.gamma(0)));
  const Gen_hepevt &hepevt2 = gen_level(get_hepevt(pi0.gamma(1)));
  if( hepevt1 && hepevt1.mother() && hepevt2 && hepevt2.mother() ) {
    Gen_hepevt &mother1 = hepevt1.mother();
    if( mother1.get_ID() == hepevt2.mother().get_ID()
          && mother1.idhep() == PDG_PI0 ) return mother1;
  }
  return Gen_hepevt_Manager::get_manager().get_NULL();
}

const Gen_hepevt &get_hepevt(const Mdst_vee2& vee, int) {
  const Gen_hepevt &hepevt1 = get_hepevt(vee.chgd(0));
  const Gen_hepevt &hepevt2 = get_hepevt(vee.chgd(1));
  if( hepevt1 && hepevt1.mother() && hepevt2 && hepevt2.mother() ) {
    const Gen_hepevt &mother1 = hepevt1.mother();
    if( mother1.get_ID() == hepevt2.mother().get_ID() ) {
        int idhep = mother1.idhep();
        int kind = vee.kind();
        if( (kind==1 && idhep==PDG_K0S) || (kind==2 && idhep==PDG_LAM) ||
            (kind==3 && idhep==PDG_ALAM) || (kind==4 && idhep==PDG_GAMM) )
          return mother1;
    }
  }
  return Gen_hepevt_Manager::get_manager().get_NULL();
}


Mdst_charged &mdst_charged(const Mdst_trk&trk, int ith) {
  Mdst_charged_Manager & chmgr = Mdst_charged_Manager::get_manager();
  Mdst_charged_Index index(chmgr.index( "trk" ));
  index.update();
  Panther_ID trkID(trk.get_ID());
  std::vector<Mdst_charged> point = point_from ( trkID, index );
  if(ith>=0 && ith<(int)point.size())  return point[ith];
  else return Mdst_charged_Manager::get_manager().get_NULL();
}
  

Mdst_charged &mdst_charged(const Rectrk&trk, int ith) {
  Mdst_charged_Manager & chmgr = Mdst_charged_Manager::get_manager();
  Mdst_charged_Index index(chmgr.index( "trk" ));
  index.update();
  Panther_ID trkID(trk.get_ID());
  std::vector<Mdst_charged> point = point_from ( trkID, index );
//   Panther_ID HepID(point[i].hepevt_ID());
//   Gen_hepevt_Manager &hepmgr = Gen_hepevt_Manager::get_manager();
//   return hepmgr ( HepID );
  if(ith>=0 && ith<(int)point.size())  return point[ith];
  else return Mdst_charged_Manager::get_manager().get_NULL();
  //  return point[ith].charged();
}
  

bool good_charged(const Mdst_charged& charged, float cl_cut, float dz_cut, 
float dr_cut) {

  Mdst_trk & trk = charged.trk();

  Mdst_trk_fit & trkfit = trk.mhyp(2); 

  bool ret = true;

  // confidence level cut
  // dout(Debugout::DDEBUG,"mdst") << "prob_=" << prob_(trkfit.chisq(), trkfit.ndf()) << std::endl;
  //  if ( 0 == trkfit.ndf()) {
  if ( trkfit.ndf() <= 0 || trkfit.chisq() < 0.) {
    ret = false;
  } else {  
    if ( prob_(trkfit.chisq(), trkfit.ndf()) <= cl_cut) ret = false;
  }
  // dz cut
  // dout(Debugout::DDEBUG,"mdst") << "z=" << trkfit.helix(3) << std::endl;
  if ( fabs(trkfit.helix(3)) > dz_cut ) ret = false;
  // dr cut
  // dout(Debugout::DDEBUG,"mdst") << "z=" << trkfit.helix(0) << std::endl;
  if ( fabs(trkfit.helix(0)) > dr_cut ) ret = false;

  return ret;
}

bool good_gamma(const Mdst_gamma& gamma, float ecut, float e925cut, 
float widcut, int nhcut, float) {

  bool ret = true;

  // check shower status first

  Mdst_ecl & shower = gamma.ecl();
  Mdst_ecl_aux_Manager& aux = Mdst_ecl_aux_Manager::get_manager();
  Mdst_ecl_aux & shower_aux(aux(Panther_ID(shower.get_ID())));

  int match = shower.match();
  double energy = shower.energy();
  double theta = shower.theta();
  float ftheta = theta*180./M_PI;
  float fe9oe25 = shower_aux.e9oe25();
  float fwidth = shower_aux.width();
  int nhit = shower_aux.nhits();

  if(ftheta <  17.0) ret=false;
  if(ftheta > 150.0) ret=false;
  if(shower.quality() != 0) ret=false;
  if(match == 1) ret=false;
  if(energy < ecut) ret=false;

  if(energy < 0.5) {
     if(nhit <= nhcut) ret=false;
     if(fe9oe25 < e925cut) ret=false;
     if(fwidth > widcut) ret=false;
  }
  return ret;
}

bool good_gamma(const Mdst_ecl& shower, float ecut, float e925cut, 
float widcut, int nhcut, float) {

  bool ret = true;

  Mdst_ecl_aux_Manager& aux = Mdst_ecl_aux_Manager::get_manager();
  Mdst_ecl_aux & shower_aux(aux(Panther_ID(shower.get_ID())));

  int match = shower.match();
  double energy = shower.energy();
  double theta = shower.theta();
  float ftheta = theta*180./M_PI;
  float fe9oe25 = shower_aux.e9oe25();
  float fwidth = shower_aux.width();
  int nhit = shower_aux.nhits();

  if(ftheta <  17.0) ret=false;
  if(ftheta > 150.0) ret=false;
  if(shower.quality() != 0) ret=false;
  if(match == 1) ret=false;
  if(energy < ecut) ret=false;

  if(energy < 0.5) {
     if(nhit <= nhcut) ret=false;
     if(fe9oe25 < e925cut) ret=false;
     if(fwidth > widcut) ret=false;
  }
  return ret;
}

//======================================================
void copy_vee_to_vee2 () {
//======================================================
// Copy MDST_Vee contents to MDST_Vee2.
//======================================================

  // dout(Debugout::DDEBUG,"mdst") << BsCouTab(MDST_VEE2) << std::endl;
  // dout(Debugout::DDEBUG,"mdst") << BsCouTab(MDST_VEE ) << std::endl;

  if (!BsCouTab(MDST_VEE2) && BsCouTab(MDST_VEE)) {
  } else {
    // dout(Debugout::DDEBUG,"mdst") << "return" << std::endl; 
      return;
  }

  Mdst_vee_Manager& mdst_vee_manager = Mdst_vee_Manager::get_manager();
  Mdst_vee2_Manager& mdst_vee2_manager = Mdst_vee2_Manager::get_manager();

  for (Mdst_vee_Manager::iterator i = mdst_vee_manager.begin();
       i != mdst_vee_manager.end(); i++) {
    Mdst_vee& vee = *i;
    Mdst_vee2& vee2 = mdst_vee2_manager.add();
    //
    vee2.chgd(0,vee.dau(0));
    vee2.chgd(1,vee.dau(1));
    //
    vee2.kind(vee.kind());
    //
    vee2.px(vee.px());
    vee2.py(vee.py());
    vee2.pz(vee.pz());
    vee2.energy(vee.energy());
    //
    vee2.vx(vee.vx());
    vee2.vy(vee.vy());
    vee2.vz(vee.vz());
    //
    vee2.z_dist(vee.z_dist());
    vee2.chisq(1.e10);
    //
    vee2.type(vee.type());
    //
  }
  return;
}

//==========================================================
void remove_duplicates () {
//==========================================================
// Removes duplicated entities in exp.7 MDSTs.
// Should be called every event at the beginnig of user_ana.
//==========================================================

  // dout(Debugout::DDEBUG,"mdst") << "remove_duplicates is called" << std::endl;

//Check existence of Belle_Event
  //  Belle_event_Manager& evtmgr = Belle_event_Manager::get_manager();
  //  if(0 == evtmgr.count()) return; //do nothing if not exist

//Check exp. no., data or MC
  //  if(evtmgr[0].ExpNo()!=7) return;
  //  if(evtmgr[0].ExpMC()==2) return;

  //temp
  //  return;
  //temp

  Mdst_event_add_Manager& mevtmgr = Mdst_event_add_Manager::get_manager();
  int nent = mevtmgr.count();
  if( nent != 2) {
    return;
  } else {
    //            mevtmgr.dump();
    std::vector<Mdst_event_add>::iterator i=mevtmgr.end();
    mevtmgr.remove(--i);
    //    dout(Debugout::INFO,"mdst") << mevtmgr.size() << std::endl;
    //            mevtmgr.dump();
  }

  //temp
  //  return;
  //temp


  Mdst_vee2_Manager& vee2mgr = Mdst_vee2_Manager::get_manager();
  nent = vee2mgr.count();
  //  dout(Debugout::INFO,"mdst") << nent << " " << nent%2 << std::endl;
  int ient;
  if( nent <= 0) {
  } else if(nent > 200) {
    dout(Debugout::ERR,"mdst") << "remove_duplicates error: size(mdst_vee2)>200" << std::endl;    
  } else if((nent%2)==1) {
    dout(Debugout::ERR,"mdst") << "remove_duplicates error: mdst_event_add and mdst_vee2 are inconsistent" << std::endl;    
  } else { 
    nent /= 2;
    //    dout(Debugout::INFO,"mdst") << "max=" << nent << std::endl;
    //    dout(Debugout::INFO,"mdst") << "begin=" << int(vee2mgr.begin()) << std::endl;
    //    dout(Debugout::INFO,"mdst") << "end  =" << int(vee2mgr.end()  ) << std::endl;
    //        vee2mgr.dump();
    ient = 1; 
    for(std::vector<Mdst_vee2>::iterator i=vee2mgr.end(); ient<=nent; ient++) {
      vee2mgr.remove(--i);
    }
    //    dout(Debugout::INFO,"mdst") << "count=" << vee2mgr.count() << std::endl;
    //    dout(Debugout::INFO,"mdst") << "begin=" << int(vee2mgr.begin()) << std::endl;
    //    dout(Debugout::INFO,"mdst") << "end  =" << int(vee2mgr.end()  ) << std::endl;
    //        vee2mgr.dump();
  }
  //  dout(Debugout::INFO,"mdst") << vee2mgr.count() << std::endl;


  Mdst_vee_daughters_Manager& veedmgr = Mdst_vee_daughters_Manager::get_manager();
  nent = veedmgr.count();
  if( nent <= 0) {
  } else if(nent > 100) {
    dout(Debugout::ERR,"mdst") << "remove_duplicates error: size(mdst_vee_daughters)>100" << std::endl;    
  } else if((nent%2)==1) {
    dout(Debugout::ERR,"mdst") << "remove_duplicates error: mdst_event_add and mdst_vee_daughters are inconsistent" << std::endl;    
  } else {
    //        veedmgr.dump();
    nent /= 2;
    ient = 1; 
    for(std::vector<Mdst_vee_daughters>::iterator i= veedmgr.end(); 
	ient<=nent; ient++) {
      veedmgr.remove(--i);
    }
    //        veedmgr.dump();
  }
  //  dout(Debugout::INFO,"mdst") << veedmgr.count() << std::endl;


  Mdst_klm_mu_ex_Manager& klmmgr = Mdst_klm_mu_ex_Manager::get_manager();
  nent = klmmgr.count();
  if( nent <= 0) {
  } else if(nent > 100) {
    dout(Debugout::ERR,"mdst") << "remove_duplicates error: size(mdst_klm_mu_ex)>100" << std::endl;    
  } else if((nent%2)==1) {
    dout(Debugout::ERR,"mdst") << "remove_duplicates error: mdst_event_add and mdst_klm_mu_ex are inconsistent" << std::endl;    
  } else {
    //        klmmgr.dump();
    nent /= 2;
    ient = 1; 
    for(std::vector<Mdst_klm_mu_ex>::iterator i= klmmgr.end(); 
	ient<=nent; ient++) {
      klmmgr.remove(--i);
    }
    //        klmmgr.dump();
  }
  //  dout(Debugout::INFO,"mdst") << klmmgr.count() << std::endl;

  return;
}

double Benergy(void) {
  return fix_mdst::get_Benergy();
}

// pi0 mass resolution function.
double Pi0resol(double p, double theta, char* side, bool mcdata, 
		int exp, int option ) {
   return fix_mdst::get_pi0resol( p, theta, side, mcdata, exp, option );
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
