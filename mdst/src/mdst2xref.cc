// -*- C++ -*-
//
// Package:     mdst
// Module:      mdst2xref
// 
// Description: Function to generate mdst_sim_xref table from other mdst tables
//
// Implimentation:
//     <Notes on implimentation>
//
// Author:      Hidekazu Kakuno
// Created:     Fri Dec  7 19:17:46 JST 2001
// $Id: mdst2xref.cc 9932 2006-11-12 14:26:53Z katayama $
//
// Revision history
//
// $Log$
// Revision 1.1  2001/12/12 02:05:03  hitoshi
// added for mdst_ecl/trk (by Kakuno).
//

// system include files
#include "belle.h"
#include "panther/panther.h"

// user include files
#include "mdst/mdst.h"

#include BELLETDF_H
#include HEPEVT_H
#include MDST_H
#include MDST_OBS_H
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

void
mdst2xref ( bool fillTrk, bool fillEcl, bool fillKlong )
{

  Mdst_sim_xref_Manager & xrefMgr = Mdst_sim_xref_Manager::get_manager();

  // skip if the mdst_sim_xref table is already filled
  if (xrefMgr.size()) return;

  belle_event *b = 
    (belle_event *)BsGetEnt(BELLE_EVENT, BsCouTab(BELLE_EVENT),BBS_No_Index);
  const bool real_data = (b && b->m_ExpMC < 2) ? 1 : 0;

  // managers
  //Mdst_charged_Manager & chMgr = Mdst_charged_Manager::get_manager();
  Mdst_trk_Manager & trkMgr = Mdst_trk_Manager::get_manager();
  Mdst_ecl_Manager & eclMgr = Mdst_ecl_Manager::get_manager();

  Mdst_charged_Index chgdIndx =
    Mdst_charged_Manager::get_manager().index("trk");
  Mdst_sim_trk_Index simtIndx =
    Mdst_sim_trk_Manager::get_manager().index("trk");
  Mdst_ecl_trk_Index ecltrkIndx = 
    Mdst_ecl_trk_Manager::get_manager().index("trk");
  Mdst_gamma_Index gammaIndx = 
    Mdst_gamma_Manager::get_manager().index("ecl");
  Mdst_sim_ecl_Index simeIndx =
    Mdst_sim_ecl_Manager::get_manager().index("ecl");
  Mdst_klong_Index klongIndx =
    Mdst_klong_Manager::get_manager().index("ecl");

  chgdIndx.update();
  gammaIndx.update();
  ecltrkIndx.update();
  simtIndx.update();
  simeIndx.update();
  klongIndx.update();

  // fill for ecl 
  // relations between ecl and trk won't be saved to aviod the
  // confliction of relation to hepevt(ecl/trk)
  if (fillEcl){
    for (std::vector<Mdst_ecl>::iterator i = eclMgr.begin();
	 i != eclMgr.end(); i++){
      bool filled(false);

      std::vector<Mdst_sim_ecl> vsime = point_from((*i).get_ID(), simeIndx);
      std::vector<Mdst_gamma> vgam = point_from((*i).get_ID(), gammaIndx);

      if (fillKlong){
	std::vector<Mdst_klong> vkl = point_from((*i).get_ID(), klongIndx);
	for (std::vector<Mdst_klong>::iterator j = vkl.begin();
	     j != vkl.end(); j++){
	  
	  // fill Mdst_sim_xref table (ecl+klm)
	  Mdst_sim_xref & xref = xrefMgr.add();
	  xref.ecl(*i);
	  if (vsime.size()){
	    xref.hepevt(vsime[0].hepevt());
	    xref.frac(vsime[0].f2(0));
	  }
	  if (vgam.size()) xref.gamma(vgam[0]);
	  xref.klong(*j);
	  xref.klm((*j).klm());
	  filled = true;
	} // end of klong loop
	if (filled) continue;
      }
      
      // fill Mdst_sim_xref table (ecl)
      // do not skip for MC even if there no relations to avoid
      // taking the wrong hepevt by get_hepevt function
      if (real_data && !vsime.size() && !vgam.size()) continue;
      Mdst_sim_xref & xref = xrefMgr.add();
      xref.ecl(*i);
      if (vsime.size()){
	xref.hepevt(vsime[0].hepevt());
	xref.frac(vsime[0].f2(0));
      }
      if (vgam.size()) xref.gamma(vgam[0]);
    }
  }

  // fill for trk
  if (fillTrk){
    for (std::vector<Mdst_trk>::iterator i = trkMgr.begin();
	 i != trkMgr.end(); i++){
      bool filled(false);
      std::vector<Mdst_sim_trk> vsimt = point_from((*i).get_ID(), simtIndx);
      std::vector<Mdst_charged> vchgd = point_from((*i).get_ID(), chgdIndx);

      if (fillEcl){
	std::vector<Mdst_ecl_trk> vecltrk = 
	  point_from((*i).get_ID(), ecltrkIndx);
	for (std::vector<Mdst_ecl_trk>::iterator j = vecltrk.begin();
	     j != vecltrk.end(); j++){
	  bool filled2(false);
	  std::vector<Mdst_gamma> vgam = 
	    point_from((*j).ecl().get_ID(), gammaIndx);

	  if (fillKlong){
	    std::vector<Mdst_klong> vkl =
	      point_from((*j).ecl().get_ID(), klongIndx);
	    for (std::vector<Mdst_klong>::iterator k = vkl.begin();
		 k != vkl.end(); k++){

	      // fill Mdst_sim_xref table (trk+ecl+klm)
	      Mdst_sim_xref & xref = xrefMgr.add();
	      xref.trk(*i);
	      if (vchgd.size()){
		xref.charged(vchgd[0]);
		xref.acc(vchgd[0].acc());
		xref.tof(vchgd[0].tof());
	      }
	      if (vsimt.size()){
		xref.hepevt(vsimt[0].hepevt());
		xref.frac(vsimt[0].frac());
	      }
	      xref.ecl((*j).ecl());
	      if (vgam.size()) xref.gamma(vgam[0]);
	      xref.klong(*k);
	      xref.klm((*k).klm());
	      filled = true;
	      filled2 = true;
	    }  // end of klong loop
	    if (filled2) continue;
	  }

	  // fill Mdst_sim_xref table (trk+ecl)
	  Mdst_sim_xref & xref = xrefMgr.add();
	  xref.trk(*i);
	  if (vchgd.size()){
	    xref.charged(vchgd[0]);
	    xref.acc(vchgd[0].acc());
	    xref.tof(vchgd[0].tof());
	  }
	  if (vsimt.size()){
	    xref.hepevt(vsimt[0].hepevt());
	    xref.frac(vsimt[0].frac());
	  }
	  xref.ecl((*j).ecl());
	  if (vgam.size()) xref.gamma(vgam[0]);
	  filled = true;

	} // end of ecltlk loop
	if (filled) continue;
      }

      // fill Mdst_sim_xref table (trk)
      if (real_data && !vsimt.size() && !vchgd.size()) continue;
      Mdst_sim_xref & xref = xrefMgr.add();
      xref.trk(*i);
      if (vchgd.size()){
	xref.charged(vchgd[0]);
	xref.acc(vchgd[0].acc());
	xref.tof(vchgd[0].tof());
      }
      if (vsimt.size()){
	xref.hepevt(vsimt[0].hepevt());
	xref.frac(vsimt[0].frac());
      }
      filled = true;
    } // end of trk loop
  }
}      
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
