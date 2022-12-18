//
// $Id: IpProfile.cc 10202 2007-06-25 10:35:15Z hitoshi $
//
// $Log$
// Revision 1.14  2002/03/29 08:32:29  tomura
// To avoid compile error.
//
// Revision 1.13  2001/12/24 12:03:26  katayama
// gcc 3.0 and headers are cleaned up
//
// Revision 1.12  2001/12/23 09:58:21  katayama
// removed Strings.h
//
// Revision 1.11  2001/06/18 11:34:23  tomura
// ip_tube is added. (by Sumisawa-san.)
// Update for new archive database "beam_profile_arc2".
// Modify usable() to do check_beginrun().
//
// Revision 1.10  2001/04/27 12:29:16  tomura
// Update for event#-dependent IP. (Default is run-by-run.)
//
// Revision 1.9  2000/11/10 09:38:33  tomura
// Update for archive database.
//
// Revision 1.8  2000/06/24 14:20:45  nakadair
// Bug fix of b_life_smaserd() in MC case.
//
// Revision 1.7  2000/04/26 09:41:42  katayama
// added std::
//
// Revision 1.6  2000/04/18 03:40:08  nakadair
// add crossing angle for data analysis.
// set beam energy from Belle_Runhead and turn HER by m_xing_angle
//
// Revision 1.5  2000/04/13 12:34:57  katayama
// Added std:: to cout,cerr,endl etc.
//
// Revision 1.4  2000/03/30 15:09:36  nakadair
//  Correct flight length of B meson, Suppress waning massages, Modify dump()
//
// Revision 1.3  2000/03/07 11:08:32  katayama
// compatibility with CC5.0
//
// Revision 1.2  2000/02/15 07:03:27  tomura
// Updated by Nakadaira-san.
//
//

#include "belle.h"
#include <iostream>
#include "event/BelleEvent.h"

#include "basf/module.h"
#include "basf/module_descr.h"

#include "panther/panther.h"
#include BELLETDF_H
#include RUN_INFO_H

#include "pntdb/TPntFDDB.h"
#include "pntdb/TPntDB.h"

#include "belleCLHEP/Vector/LorentzVector.h"
#include "belleCLHEP/Vector/ThreeVector.h"
#include "belleCLHEP/Vector/Rotation.h"
#include "belleCLHEP/Matrix/Vector.h"
#include "belleCLHEP/Matrix/SymMatrix.h"
#include "belleCLHEP/Matrix/Matrix.h"
#include "belleCLHEP/Geometry/Point3D.h"

#include "benergy/BeamEnergy.h"

#include "ip/IpProfile.h"
#include "belleutil/belutil.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif



bool         IpProfile::m_MC = false;
int          IpProfile::m_exp = -1;
int          IpProfile::m_run = -1;
int          IpProfile::m_evtbin = -1;
int          IpProfile::m_version = 0;
bool         IpProfile::m_valid = false;
bool         IpProfile::m_usable = false;
bool         IpProfile::m_b_life_smeared_usable = false;
HepSymMatrix IpProfile::m_b_life_smear_matrix(3, 0);

HepPoint3D   IpProfile::m_position;
HepSymMatrix IpProfile::m_position_err(3, 0);
HepSymMatrix IpProfile::m_position_err_b_life_smeared(3, 0);

HepPoint3D   IpProfile::m_e_position;
HepSymMatrix IpProfile::m_e_position_err(3, 0);
HepSymMatrix IpProfile::m_e_position_err_b_life_smeared(3, 0);

double       IpProfile::m_blife_smear_xy = 2.1e-3;
double       IpProfile::m_blife_smear_z  = 0.0e-2;
double       IpProfile::m_xing_angle = 2.2e-2;
double       IpProfile::m_position_offset[3] = {0.0, 0.0, 0.0};
double       IpProfile::m_position_err_offset[3] = {0.0, 0.0, 0.0};
double       IpProfile::m_err_min = 1.0e-5;

HepLorentzVector IpProfile::m_her;
HepLorentzVector IpProfile::m_ler;

double       IpProfile::m_theta_x = 0;
double       IpProfile::m_theta_y = m_xing_angle/2.;
double       IpProfile::m_theta_z = 0;


void
IpProfile::define_global(Module_descr* dscr){
  dscr->define_param("ip_profile_version",
		     "ip_profile_version", &IpProfile::m_version);
  dscr->define_param("ip_err_smear_xy",
		     "ip_err_smear_xy", &IpProfile::m_blife_smear_xy);
  dscr->define_param("ip_err_smear_z",
		     "ip_err_smear_z", &IpProfile::m_blife_smear_z);
  dscr->define_param("xing_angle",
		     "xing_angle", &IpProfile::m_xing_angle);
  dscr->define_param("ip_position_offset_x",
		     "Offset for IP position X", &IpProfile::m_position_offset[0]);
  dscr->define_param("ip_position_offset_y",
		     "Offset for IP position Y", &IpProfile::m_position_offset[1]);
  dscr->define_param("ip_position_offset_z",
		     "Offset for IP position Z", &IpProfile::m_position_offset[2]);
  dscr->define_param("ip_err_offset_x",
		     "Offset for IP error X", &IpProfile::m_position_err_offset[0]);
  dscr->define_param("ip_err_offset_y",
		     "Offset for IP error Y", &IpProfile::m_position_err_offset[1]);
  dscr->define_param("ip_err_offset_z",
		     "Offset for IP error Z", &IpProfile::m_position_err_offset[2]);
  dscr->define_param("ip_err_min",
		     "Minimum value of IP size", &IpProfile::m_err_min);

  return;
};


int
IpProfile::begin_run(const int ver){
  m_evtbin = -1;

  BeamEnergy::begin_run();

  Belle_runhead_Manager& evtmgr = Belle_runhead_Manager::get_manager();
  Belle_runhead_Manager::const_iterator belleevt = evtmgr.begin();
  if(belleevt == evtmgr.end()){
    dout(Debugout::ERR,"IpProfile") << " Error: Cannot read \"Belle_RunHead\"." << std::endl;
    m_usable = false;
    return 3;
  }

  if((*belleevt)&&(belleevt->ExpMC()==1)){
    // Experimental Data
    int exp = belleevt->ExpNo(); 
    int run = belleevt->RunNo(); 

    if(m_exp==exp && m_run==run && m_version == ver && m_valid){
      m_usable = true;
      return 0;
    }

    m_MC = false;
    m_exp = exp;
    m_run = run;
    m_version = ver;
    m_valid = false;
    m_usable = false;
    m_b_life_smeared_usable = false;

    if(int status = GetIPprofile(ver)){
      dout(Debugout::ERR,"IpProfile") << " There is no IP profile data." << std::endl;
      return status;
    }
    else{
      Ip_profile_general_Manager& ipgmgr = Ip_profile_general_Manager::get_manager();
      Ip_profile_general_Manager::const_iterator ipg = ipgmgr.begin();
      if(ipg == ipgmgr.end()){
	dout(Debugout::ERR,"IpProfile") << " There is no IP profile general data." << std::endl;
	return 4;
      }

      m_position.setX(ipg->IPx() + m_position_offset[0]);
      m_position.setY(ipg->IPy() + m_position_offset[1]);
      m_position.setZ(ipg->IPz() + m_position_offset[2]);

      if(m_position_err_offset[0]){
	double size_tmp = sqrt(ipg->IPerror(0)) + m_position_err_offset[0];
	m_position_err(1, 1) = (size_tmp > m_err_min) ?
	  size_tmp*size_tmp : m_err_min*m_err_min;
      }
      else
	m_position_err(1, 1) = ipg->IPerror(0);
      m_position_err(2, 1) = ipg->IPerror(1);
      if(m_position_err_offset[1]){
	double size_tmp = sqrt(ipg->IPerror(2)) + m_position_err_offset[1];
	m_position_err(2, 2) = (size_tmp > m_err_min) ?
	  size_tmp*size_tmp : m_err_min*m_err_min;
      }
      else
	m_position_err(2, 2) = ipg->IPerror(2);
      m_position_err(3, 1) = ipg->IPerror(3);
      m_position_err(3, 2) = ipg->IPerror(4);
      if(m_position_err_offset[2]){
	double size_tmp = sqrt(ipg->IPerror(5)) + m_position_err_offset[2];
	m_position_err(3, 3) = (size_tmp > m_err_min) ?
	  size_tmp*size_tmp : m_err_min*m_err_min;
      }
      else
	m_position_err(3, 3) = ipg->IPerror(5);

      Ip_profile_hadron_Manager& iphmgr = Ip_profile_hadron_Manager::get_manager();
      Ip_profile_hadron_Manager::const_iterator iph = iphmgr.begin();
      if(iph == iphmgr.end()){
	dout(Debugout::ERR,"IpProfile") << " There is no IP profile hadron data." << std::endl;
	return 4;
      }
      m_theta_x = iph->Theta_x();
      m_theta_y = iph->Theta_y();
      //m_theta_z = iph->Theta_z();
      m_theta_z = 0.;

      m_valid  = true;
      m_usable = true;

      if(calc_b_life_smeared())
	return 5;

      m_position_err_b_life_smeared = m_position_err + m_b_life_smear_matrix;
      m_b_life_smeared_usable = true;

      return 0; 
    }
  }
  else{
    // Monte Calro
    m_MC = true;

    m_valid = false;
    m_usable = false;
    m_b_life_smeared_usable = false;

    Belle_nominal_beam_Manager& bnbmgr = Belle_nominal_beam_Manager::get_manager();
    Belle_nominal_beam_Manager::const_iterator bnb = bnbmgr.begin();
    if(bnb == bnbmgr.end()){
      dout(Debugout::ERR,"IpProfile") << " Can not access to table [Belle_Nominal_Beam]." << std::endl;
      return 5;
    }

    m_position.setX(bnb->ip_x());
    m_position.setY(bnb->ip_y());
    m_position.setZ(bnb->ip_z());

    m_e_position = m_position;

    HepSymMatrix ip_sigma(3, 0);
    ip_sigma(1, 1) = pow(bnb->sigma_ip_x(), 2);
    ip_sigma(2, 2) = pow(bnb->sigma_ip_y(), 2);
    ip_sigma(3, 3) = pow(bnb->sigma_ip_z(), 2);
    HepRotation rot_y;
    rot_y.rotateY(bnb->angle_ip_zx());
    HepMatrix rotM(3, 3, 0);
    rotM = rot_y;
    m_position_err = ip_sigma.similarity(rotM);

    m_e_position_err = m_position_err;

    m_valid = true;
    m_usable = true;

    m_theta_x = 0.;
    m_theta_y = bnb->angle_ip_zx();
    m_theta_z = 0.;

    if(calc_b_life_smeared())
      return 5;

    m_position_err_b_life_smeared = m_position_err + m_b_life_smear_matrix;
    m_b_life_smeared_usable = true;

    m_e_position_err_b_life_smeared = m_position_err_b_life_smeared;

    return 1;
  }
};


int
IpProfile::calc_b_life_smeared(){
  Belle_nominal_beam_Manager& bnbmgr = Belle_nominal_beam_Manager::get_manager();
  Belle_nominal_beam_Manager::const_iterator bnb = bnbmgr.begin();

  //double e_her = 7.9965;
  //double e_ler = 3.5000;
  double e_her = BeamEnergy::E_HER();
  double e_ler = BeamEnergy::E_LER();

  //Belle_runhead_Manager& rhdmgr = Belle_runhead_Manager::get_manager();
  //Belle_runhead_Manager::const_iterator rhd = rhdmgr.begin();
  //if(rhd != rhdmgr.end()){
  //  e_her = rhd->EHER();
  //  e_ler = rhd->ELER();
  //}
  HepLorentzVector p4_her(0.0, 0.0,  e_her, e_her);
  p4_her.rotateY(m_xing_angle);
  HepLorentzVector p4_ler(0.0, 0.0, -e_ler, e_ler);

  if(bnb == bnbmgr.end()){
// dout(Debugout::DDEBUG,"IpProfile") << " Can not access to table [Belle_Nominal_Beam]." << std::endl;
// dout(Debugout::DDEBUG,"IpProfile") << " regarded as TDF values" << std::endl;
  }
  else{
    Hep3Vector p3_her(bnb->px_high(), bnb->py_high(), bnb->pz_high());
    Hep3Vector p3_ler(bnb->px_low(), bnb->py_low(), bnb->pz_low());
    p4_ler = HepLorentzVector(p3_ler, p3_ler.mag());
    p4_her = HepLorentzVector(p3_her, p3_her.mag());
  }
// dout(Debugout::DDEBUG,"IpProfile") << p4_her << std::endl;
// dout(Debugout::DDEBUG,"IpProfile") << p4_ler << std::endl;

  m_her = p4_her;
  m_ler = p4_ler;

  HepSymMatrix blife_cor(3, 0);
  blife_cor(1, 1) = pow(m_blife_smear_xy, 2);
  blife_cor(2, 2) = pow(m_blife_smear_xy, 2);
  blife_cor(3, 3) = pow(m_blife_smear_z, 2);

  Hep3Vector boostaxis((p4_ler+p4_her).boostVector());
  Hep3Vector z_axis(0.0, 0.0, 1.0);
  HepRotation rot;

  rot.rotate(boostaxis.angle(z_axis), z_axis.cross(boostaxis));
  dout(Debugout::INFO,"IpProfile") << boostaxis.angle(z_axis) << std::endl;

  HepMatrix rotM(3, 3, 0);
  rotM = rot;
  m_b_life_smear_matrix = blife_cor.similarity(rotM);

  return 0;
}


void
IpProfile::dump(std::ostream& sout){
  (void *)sout;
  dout(Debugout::ERR,"IpProfile")
    << " dump() function does not accept an argument any more." << std::endl;
}

void IpProfile::dump() {
  if( m_usable ) {
    dout(Debugout::INFO,"IpProfile")
      << " IpProfile is available." << std::endl;
  } else {
    dout(Debugout::WARN,"IpProfile")
      << " IpProfile is !NOT! available." << std::endl;
  }

  dout(Debugout::INFO,"IpProfile") 
    << " Exp No." << m_exp << ", Run No." << m_run << std::endl;
  dout(Debugout::INFO,"IpProfile")
    << " Version No. " << m_version << " (0 means latest version.)"
    << std::endl;
  dout(Debugout::INFO,"IpProfile")
    << " IP position : " << m_position << std::endl;
  dout(Debugout::INFO,"IpProfile")
    << " IP position error: " << m_position_err << std::flush;
  dout(Debugout::INFO,"IpProfile")
    << " Flight length of B meson : x,y = " << m_blife_smear_xy << " (cm)"
    << std::endl
    << "                                             : z   = "
    << m_blife_smear_z << "(cm)"<< std::endl;
  dout(Debugout::INFO,"IpProfile")
    << " IP position error(smeard): " << m_position_err_b_life_smeared
    << std::flush;
  return;
}


void
IpProfile::e_dump(std::ostream& sout){
  (void *)sout;
  dout(Debugout::ERR,"IpProfile")
    << " e_dump() function does not accept an argument any more." << std::endl;
}

void IpProfile::e_dump(){
  if( m_usable ) {
    dout(Debugout::INFO,"IpProfile")
      << " IpProfile of this event is available." << std::endl;
  } else {
    dout(Debugout::WARN,"IpProfile")
      << " IpProfile of this event is !NOT! available." << std::endl;
  }
  dout(Debugout::INFO,"IpProfile")
    << " Exp No." << m_exp << ", Run No." << m_run 
    << ", Evt BIN No." << m_evtbin << std::endl;
  dout(Debugout::INFO,"IpProfile")
    << " Version No. " << m_version << " (0 means latest version.)"
    << std::endl;
  dout(Debugout::INFO,"IpProfile")
    << " IP position : " << m_e_position << std::endl;
  dout(Debugout::INFO,"IpProfile")
    << " IP position error: " << m_e_position_err << std::flush;
  dout(Debugout::INFO,"IpProfile")
    << " Flight length of B meson : x,y = " << m_blife_smear_xy << " (cm)"
    << std::endl
    << "                                             : z   = "
    << m_blife_smear_z << "(cm)"<< std::endl;
  dout(Debugout::INFO,"IpProfile")
    << " IP position error(smeard): " << m_e_position_err_b_life_smeared
    << std::flush;
  return;
}


int
IpProfile::GetIPprofile(const int ver){
  Ip_profile_general_Manager& General_mgr = Ip_profile_general_Manager::get_manager();
  Ip_profile_hadron_Manager& Hadron_mgr = Ip_profile_hadron_Manager::get_manager();
  Ip_profile_bhabha_Manager& Bhabha_mgr = Ip_profile_bhabha_Manager::get_manager();
  Accelerator_info_Manager& Ainfo_mgr = Accelerator_info_Manager::get_manager();

  Ainfo_mgr.remove();
  Bhabha_mgr.remove();
  Hadron_mgr.remove();
  General_mgr.remove();

  TPntFDDB fddb("rif");
  TPntDB db(fddb, "beam_profile");
  if(!db.IsOK()){
    dout(Debugout::ERR,"IpProfile") << "Error : Cannot access rif::beam_profile.\n";
    return 2;
  }

  if(db.Get(m_exp, m_run, ver) <= 0){
    db.Close();
    db.Open(fddb, "beam_profile_arc1");
    if(!db.IsOK()){
      dout(Debugout::ERR,"IpProfile") << "Error : Cannot access rif::beam_profile_arc1.\n";
      return 3;
    }

    if(db.Get(m_exp, m_run, ver) <= 0){
      db.Close();
      db.Open(fddb, "beam_profile_arc2");
      if(!db.IsOK()){
	dout(Debugout::ERR,"IpProfile") << "Error : Cannot access rif::beam_profile_arc2.\n";
	return 4;
      }

      if(db.Get(m_exp, m_run, ver) <= 0){
	dout(Debugout::ERR,"IpProfile") << "Error : Cannot get data from database correctly. (Exp"
		  << m_exp << ",Run" << m_run << ")\n";
	return 1;
      }
    }
  }

  db.non_clear(IP_PROFILE_GENERAL);
  db.non_clear(IP_PROFILE_HADRON);
  db.non_clear(IP_PROFILE_BHABHA);
  db.non_clear(ACCELERATOR_INFO);

  return 0;
}


void
IpProfile::check_runno(){
  if(check_beginrun() == -1){
    dout(Debugout::ERR,"IpProfile") << " Warning : Calling IpProfile::begin_run() again."
	      << std::endl;

    begin_run();
  }
}


int
IpProfile::set_evtbin_number(){
  // If MC, do nothing.
  if(m_MC) return 0;

  Belle_event_Manager& evt_mgr = Belle_event_Manager::get_manager();
  Belle_event_Manager::const_iterator itevt = evt_mgr.begin();
  if(itevt == evt_mgr.end()){
    dout(Debugout::ERR,"IpProfile") << " Error : There is no Belle_event table."
	      << std::endl;

    return 2;
  }

  int evtno = itevt->EvtNo()&0x0FFFFFFF;

  return set_evtbin_number(evtno);
}


int
IpProfile::set_evtbin_number(const int event){
  // If MC, do nothing
  if(m_MC) return 0;

  if(event/BIN_EVENTS != m_evtbin){
    m_evtbin = event/BIN_EVENTS;

    Ip_profile_general_Manager& ipgmgr = Ip_profile_general_Manager::get_manager();

    if(ipgmgr.count() <= 0){
      dout(Debugout::ERR,"IpProfile") << " Error : Ip_profile_general table is empty."
		<< std::endl;

      return 1;
    }

    Ip_profile_general *ipg(NULL);

    if(m_evtbin+2 > ipgmgr.count())
      ipg = &(ipgmgr.back());
    else{
      Panther_ID ID(m_evtbin+2);
      ipg = &ipgmgr(ID);
    }

    m_e_position.setX(ipg->IPx() + m_position_offset[0]);
    m_e_position.setY(ipg->IPy() + m_position_offset[1]);
    m_e_position.setZ(ipg->IPz() + m_position_offset[2]);

    if(m_position_err_offset[0]){
      double size_tmp = sqrt(ipg->IPerror(0)) + m_position_err_offset[0];
      m_e_position_err(1, 1) = (size_tmp > m_err_min) ?
	size_tmp*size_tmp : m_err_min*m_err_min;
    }
    else
      m_e_position_err(1, 1) = ipg->IPerror(0);
    m_e_position_err(2, 1) = ipg->IPerror(1);
    if(m_position_err_offset[1]){
      double size_tmp = sqrt(ipg->IPerror(2)) + m_position_err_offset[1];
      m_e_position_err(2, 2) = (size_tmp > m_err_min) ?
	size_tmp*size_tmp : m_err_min*m_err_min;
    }
    else
      m_e_position_err(2, 2) = ipg->IPerror(2);
    m_e_position_err(3, 1) = ipg->IPerror(3);
    m_e_position_err(3, 2) = ipg->IPerror(4);
    if(m_position_err_offset[2]){
      double size_tmp = sqrt(ipg->IPerror(5)) + m_position_err_offset[2];
      m_e_position_err(3, 3) = (size_tmp > m_err_min) ?
	size_tmp*size_tmp : m_err_min*m_err_min;
    }
    else
      m_e_position_err(3, 3) = ipg->IPerror(5);

    m_e_position_err_b_life_smeared = m_e_position_err + m_b_life_smear_matrix;
  }

  return 0;
}


const HepPoint3D&
IpProfile::e_position(){
  set_evtbin_number();

  return m_e_position;
}


const HepSymMatrix&
IpProfile::e_position_err(){
  set_evtbin_number();

  return m_e_position_err;
}


const HepSymMatrix&
IpProfile::e_position_err_b_life_smeared(){
  set_evtbin_number();

  return m_e_position_err_b_life_smeared;
}


Particle
IpProfile::ip_tube(const int flag){
  Ptype gamm( "GAMM" );
  if( m_b_life_smear_matrix.trace() == 0 ){
    dout(Debugout::ERR,"IpProfile")<<" Error: Please call IpProfile::begin_run()"
	     <<" in your BASF begin_run function."<<std::endl;
    return Particle();
  }

	// modified by T.H 2006/0509
  // HepLorentzVector tmp_mom( 0., 0., 1., 1. );
  HepLorentzVector tmp_mom( 0., 0., 1e10, 1e10 );
  tmp_mom.rotateX(m_theta_x);
  tmp_mom.rotateY(m_theta_y);

  Particle tmp_prt;
  tmp_prt.pType( gamm );
  tmp_prt.momentum().momentum(tmp_mom);
  tmp_prt.momentum().position(position(flag),
			      position_err_b_life_smeared(flag));

  return tmp_prt;
}


kfitterparticle
IpProfile::ip_tube_kfitterparticle(const int flag){
  Ptype gamm( "GAMM" );
  if( m_b_life_smear_matrix.trace() == 0 ){
    dout(Debugout::ERR,"IpProfile")<<" Error: Please call IpProfile::begin_run()"
	     <<" in your BASF begin_run function."<<std::endl;
    return kfitterparticle( HepLorentzVector(0.,0.,0.,0.),
			    HepPoint3D(0.,0.,0.),
			    HepSymMatrix(7,0), 0.0, 0.0 );
  }
  
  HepSymMatrix tmp_E( 7, 0 );

	// modified by T.H 2006/0509
  // HepLorentzVector tmp_mom( 0., 0., 1., 1. );
  HepLorentzVector tmp_mom( 0., 0., 1e10, 1e10 );
  tmp_mom.rotateX(m_theta_x);
  tmp_mom.rotateY(m_theta_y);

  tmp_E.sub( 5, position_err_b_life_smeared(flag) );

  return kfitterparticle(tmp_mom, position(flag), tmp_E, 0, gamm.mass() );
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
