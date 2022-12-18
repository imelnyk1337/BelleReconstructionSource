//
// $Id: IpProfile.h 10143 2007-05-25 23:12:57Z katayama $
//
// $Log$
// Revision 1.7  2002/01/03 11:04:36  katayama
// Point3D and other header files are cleaned
//
// Revision 1.6  2001/06/18 11:34:24  tomura
// ip_tube is added. (by Sumisawa-san.)
// Update for new archive database "beam_profile_arc2".
// Modify usable() to do check_beginrun().
//
// Revision 1.5  2001/04/27 12:30:46  tomura
// Update for event#-dependent IP. (Default is run-by-run.)
//
// Revision 1.4  2000/04/18 03:38:25  nakadair
// add crossing angle for data analysis.
//
// Revision 1.3  2000/04/13 12:34:56  katayama
// Added std:: to cout,cerr,endl etc.
//
// Revision 1.2  2000/02/15 07:05:03  tomura
// Updated by Nakadaira-san.
//
//
#ifndef __IPPROFILE_H__
#define __IPPROFILE_H__

#include "belle.h"
#include "basf/module.h"


#include "belleCLHEP/Vector/LorentzVector.h"
#include "belleCLHEP/Matrix/Vector.h"
#include "belleCLHEP/Matrix/SymMatrix.h"
#include "belleCLHEP/Matrix/Matrix.h"
#ifndef CLHEP_POINT3D_H
#include "belleCLHEP/Geometry/Point3D.h"
#endif

#include "particle/Particle.h"
#include "kfitter/kfitterparticle.h"

#include <iostream>
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

#define BIN_EVENTS 10000


class IpProfile{
public:
  // constructor
  IpProfile(){};

  // destructor
  virtual ~IpProfile(){};

  // register variables as module parameter
  static void define_global(Module_descr* dscr);

  // Set IP profile
  static int begin_run(const int version = m_version);

  // Calculate smearing matrix of B lifetime
  static int calc_b_life_smeared();

  // Set event-depend IP profile
  static int set_evtbin_number();
  static int set_evtbin_number(const int event);

  // IP position
  static const HepPoint3D& position() {
    return m_position;
  };
  static const HepPoint3D& e_position();
  static const HepPoint3D& position(const int flag) {
    return flag ? e_position() : position();
  };

  // Error matrix of IP position
  static const HepSymMatrix& position_err() {
    return m_position_err;
  };
  static const HepSymMatrix& e_position_err();
  static const HepSymMatrix& position_err(const int flag) {
    return flag ? e_position_err() : position_err();
  };

  // Error matrix of IP position smeared by B lifetime
  static const HepSymMatrix& position_err_b_life_smeared() {
    return m_position_err_b_life_smeared;
  };
  static const HepSymMatrix& e_position_err_b_life_smeared();
  static const HepSymMatrix& position_err_b_life_smeared(const int flag) {
    return flag ?
      e_position_err_b_life_smeared() : position_err_b_life_smeared();
  };


  // IP tube
  static Particle ip_tube(const int flag);
  static Particle ip_tube(){
    return ip_tube(0);
  }
  static Particle e_ip_tube(){
    return ip_tube(1);
  }
  static kfitterparticle ip_tube_kfitterparticle(const int flag);
  static kfitterparticle ip_tube_kfitterparticle(){
    return ip_tube_kfitterparticle(0);
  }
  static kfitterparticle e_ip_tube_kfitterparticle(){
    return ip_tube_kfitterparticle(1);
  }


  // Return Experimental number
  static int ExpNo() {
    return m_exp;
  };
  // Return Run number
  static int RunNo() {
    return m_run;
  };
  // Version of IP profile
  static int version() {
    return m_version;
  };
  // Return Event bin number
  static int EvtBinNo() {
    return m_evtbin;
  };
  // Return B life smearing matrix
  static HepSymMatrix blife_smear_matrix() {
    return m_b_life_smear_matrix;
  }

  // Optical dump of data
  static void dump();
  static void e_dump();
  // followings are obsolete
  static void dump(std::ostream& sout);
  static void e_dump(std::ostream& sout);

  // Flight length of B along the X/Y axis.
  static double blife_smear_xy() {
    return m_blife_smear_xy;
  };
  // Flight length of B along the Z axis.
  static double blife_smear_z() {
    return m_blife_smear_z;
  };

  // Crossing angle of HER
  static double xing_angle() {
    return m_xing_angle;
  };

  // IP profile is valid for run #m_run
  static bool valid() {
    return m_valid;
  };
  // IP profile is usable for analysis
  static bool usable() {
    check_runno();
    return m_usable;
  };

  // IP profile  smeared by B lifetime is usable for analysis
  static bool b_life_smeared_usable() {
    return m_b_life_smeared_usable;
  };

public:
  // Modifiers

  // set Exp# and run#
  static void ExpRunNo(const int exp, const int  run) {
    m_exp = exp;
    m_run = run;
    return;
  };
  // set version
  static int version(const int ver) {
    return m_version = ver;
  };
  static bool valid(const bool vld) {
    return m_valid = vld;
  };
  static bool usable(const bool usbl) {
    return m_usable = usbl;
  };
  static bool b_life_smeared_usable(const bool b_usbl) {
    return m_b_life_smeared_usable = b_usbl;
  };

  // set IP position
  static void position(const HepPoint3D& x) {
    m_position = x;
    return;
  };
  // set Error matrix of IP opsition
  static void position_err(const HepSymMatrix& x) {
    m_position_err = x;
    return;
  };
  // set Error matrix of IP opsition smeared by B lifetime
  static void position_err_b_life_smeared(const HepSymMatrix& x) {
    m_position_err_b_life_smeared = x;
    return;
  };

  // set Flight length of B along the X/Y axis.
  static double blife_smear_xy(const double x) {
    return m_blife_smear_xy = x;
  };
  // set Flight length of B along the Z axis.
  static double blife_smear_z(const double x) {
    return m_blife_smear_z = x;
  };
  // set B life smearing matrix
  static HepSymMatrix blife_smear_matrix(const HepSymMatrix x) {
    return m_b_life_smear_matrix = x;
  }

  // set crossing angle of HER
  static double xing_angle(const double x) {
    return m_xing_angle = x;
  };

public:
  static int m_version;
  static double m_blife_smear_xy;
  static double m_blife_smear_z;
  static double m_xing_angle;
  static double m_position_offset[3];
  static double m_position_err_offset[3];
  static double m_err_min;

private:
  static int GetIPprofile(const int version = 0);
  static void check_runno();

  static bool         m_MC;
  static int          m_exp;
  static int          m_run;
  static bool         m_valid;
  static bool         m_usable;
  static bool         m_b_life_smeared_usable;
  static HepPoint3D   m_position;
  static HepSymMatrix m_position_err;
  static HepSymMatrix m_position_err_b_life_smeared;

  static int          m_evtbin;
  static HepPoint3D   m_e_position;
  static HepSymMatrix m_e_position_err;
  static HepSymMatrix m_e_position_err_b_life_smeared;

  static HepSymMatrix m_b_life_smear_matrix;

  static HepLorentzVector m_her;
  static HepLorentzVector m_ler;

  static double m_theta_x;
  static double m_theta_y;
  static double m_theta_z;
};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* __IPPROFILE_H__ */
