//
// $Id: frec_util.h 10013 2007-02-27 06:55:37Z katayama $
//
// $Log$
// Revision 1.3  2004/05/04 21:22:37  matumot
// fixed for CC complier
// - namespace was changed : fullrecon -> fullreconNS
//   to avoid overlap for fullrecon struct used by fullrecon.tdf
// - use new operator for dynamic allocation of array
// - add std:: for sort function
//
// Revision 1.2  2004/04/19 11:07:20  matumot
// change $Header$ to $Id: frec_util.h 10013 2007-02-27 06:55:37Z katayama $
//
// Revision 1.1  2004/04/19 09:27:32  matumot
// revival version
//
//
// ======================================================
// File Name : frec_util.h
// ------------------------------------------------------
// Creation    ; 2004.04.07
// Description ; Definition for the frec_util class etc.
// Author      ; T.Matsumoto (matumot@bmail.kek.jp)
// ------------------------------------------------------

#ifndef  __FREC_UTIL__
#define  __FREC_UTIL__

#include <vector>
#include <iosfwd>

#include "belle.h"

#include "particle/Particle.h"
#include "particle/ParticleUserInfo.h"

#include "kid/atc_pid.h"

// Use CLHEP
// ~~~~~~~~~~
#include "belleCLHEP/Vector/LorentzVector.h"
#include "belleCLHEP/Vector/ThreeVector.h"
#include <belleCLHEP/Utilities/fortran.h>
#include "belleCLHEP/Alist/AList.h"

#include HEPEVT_H
#include MDST_H

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// Define parameters
// ~~~~~~~~~~~~~~~~~~
namespace fullreconNS{

extern const int ELE_CODE;
extern const int MUON_CODE;
extern const int PION_CODE;
extern const int KAON_CODE;
extern const int PROTON_CODE;

extern const Ptype Ptype_B0;
extern const Ptype Ptype_B0B;
extern const Ptype Ptype_Bplus;
extern const Ptype Ptype_Bminus;
extern const Ptype Ptype_Dst0;
extern const Ptype Ptype_Dst0B;
extern const Ptype Ptype_Dstplus;
extern const Ptype Ptype_Dstminus;
extern const Ptype Ptype_D0;
extern const Ptype Ptype_D0B;
extern const Ptype Ptype_Dplus;
extern const Ptype Ptype_Dminus;
extern const Ptype Ptype_DSplus;
extern const Ptype Ptype_DSminus;
extern const Ptype Ptype_DSstplus;
extern const Ptype Ptype_DSstminus;
extern const Ptype Ptype_Ks;
extern const Ptype Ptype_PI0;
extern const Ptype Ptype_Gamma;
extern const Ptype Ptype_Kplus;
extern const Ptype Ptype_Kminus;
extern const Ptype Ptype_PIplus;
extern const Ptype Ptype_PIminus;
extern const Ptype Ptype_Eplus;
extern const Ptype Ptype_Eminus;
extern const Ptype Ptype_MUplus;
extern const Ptype Ptype_MUminus;
extern const Ptype Ptype_RHOplus;
extern const Ptype Ptype_RHOminus;
extern const Ptype Ptype_RHO0;
extern const Ptype Ptype_A1plus;
extern const Ptype Ptype_A1minus;
extern const Ptype Ptype_Eta;
extern const Ptype Ptype_Etap;
extern const Ptype Ptype_Kstplus;
extern const Ptype Ptype_Kstminus;
extern const Ptype Ptype_Kst0;
extern const Ptype Ptype_Kst0B;
extern const Ptype Ptype_Omega;
extern const Ptype Ptype_Phi;

extern const double E_HER;
extern const double E_LER;
extern const double CROSS_ANGLE;

}

// --------------------------
//  frec_util class
// --------------------------
class frec_util {

 public:

  // member functions
  // ~~~~~~~~~~~~~~~~~
  static const double correct_dr(const Mdst_charged&, const HepPoint3D, const int);
  static const double correct_dz(const Mdst_charged&, const HepPoint3D, const int);

  static const bool prohib_dupli(Particle, Particle);
  static const bool with_dEdx( Particle&, const int, const double );

  static const bool  set_pGenHepevt( Particle&, const int = 2 );
  static const bool  set_pGenHepevt_mdstCharged( Particle&, const Ptype& );
  static const bool  set_pGenHepevt_mdstGamma( Particle& );
  static const bool  set_pGenHepevt_mdstPi0( Particle& );
  static const bool  set_pGenHepevt_mdstVee2( Particle& );
  static const bool  set_pGenHepevt_general( Particle& );

  static const bool  set_OldpGenHepevt_mdstCharged( Particle&, const Ptype& );
  static const bool  set_OldpGenHepevt_mdstGamma( Particle& );
  static const bool  set_OldpGenHepevt_mdstPi0( Particle& );
  static const bool  set_OldpGenHepevt_mdstVee2( Particle& );

  static Particle make_2body( const Ptype, Particle&, Particle& , const int = 0 );
  static Particle make_3body( const Ptype, Particle&, Particle&, Particle&, const int = 0 );
  static Particle make_4body( const Ptype, Particle&, Particle&, 
                                    Particle&, Particle&, const int = 0 );

  static const bool sel_KPI( Particle&, const double, const atc_pid& );

  static const bool HadronB();

  static const double Ebeam();
  static const HepLorentzVector p_beam();

  static const HepLorentzVector p_cm( const HepLorentzVector& );

  static const double R2( std::vector<Particle>& );

  static void division( Particle&, std::vector<Particle>&,
			std::vector<HepLorentzVector>&,
			std::vector<HepLorentzVector>& );

  static void setUserInfoB( Particle& );

  static const double Bpurity_hadronic_DX( const int, const int, const int ,
					   const int, const int );

 private :

  // variables
  // ~~~~~~~~~~~~
  static const HepLorentzVector m_p_beam;

};


// -------------------
//  thrust class
// -------------------

// forward declarations
// ~~~~~~~~~~~~~~~~~~~~~~
extern "C" {
  extern void FORTRAN_ROUTINE(thrust)(int*,float[100][4],int*,int*
                                      ,float*,float*,float[3][4]);
  extern void FORTRAN_ROUTINE(spher)(int*,float[100][4],int*,int*
                                     ,float*,float*,float[3][4]);
}

class thrust
{
  
public:

  // constructors and destructor
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  thrust();
  thrust( HepAList<HepLorentzVector>& );
  thrust( std::vector<HepLorentzVector>& );
  virtual ~thrust(){}

  // member functions
  // ~~~~~~~~~~~~~~~~
  inline Hep3Vector thr_axis(){ return m_thr_axis; }
  inline Hep3Vector sph_axis(){ return m_sph_axis; }
  inline float thr(){ return m_thr; }
  inline float obl(){ return m_obl; }
  inline float sph(){ return m_sph; }
  inline float apl(){ return m_apl; }
  inline int status(){ return m_status; }

private:

  // constructors and destructor
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  thrust( const thrust& );

  // assignment operator(s)
  // ~~~~~~~~~~~~~~~~~~~~~~
  const thrust& operator=( const thrust& );

  // comparison operators
  // ~~~~~~~~~~~~~~~~~~~~
  bool operator==( const thrust& ) const;
  bool operator!=( const thrust& ) const;

  // data members
  // ~~~~~~~~~~~~
  int m_status;
  float m_thr;
  float m_obl;
  float m_sph;
  float m_apl;
  Hep3Vector m_thr_axis;
  Hep3Vector m_sph_axis;

};

// -----------------------------------------
//  UserInfo_B class (Particle UserInfo)
// -----------------------------------------
class UserInfo_B : public ParticleUserInfo{

public :

  // default constructor & destructor
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     UserInfo_B();
    ~UserInfo_B(){}

  // copy constructor
  // ~~~~~~~~~~~~~~~~
  UserInfo_B( const UserInfo_B& );

  // constructs self object
  // ~~~~~~~~~~~~~~~~~~~~~~
  UserInfo_B* clone(void) const;

  // copy operator
  // ~~~~~~~~~~~~~
  UserInfo_B& operator = (const UserInfo_B&);

public :

  //  set ParticleUserInfo
  // ~~~~~~~~~~~~~~~~~~~~~~~
  void set( std::vector<Particle>& );
  void set( Particle& );

  // set variables
  // ~~~~~~~~~~~~~~~
  void b_mode( const int& x ){ m_b_mode = x; };
  void sub0_mode( const int& x ){ m_sub0_mode = x; };
  void sub1_mode( const int& x ){ m_sub1_mode = x; };
  void sub2_mode( const int& x ){ m_sub2_mode = x; };
  void sub3_mode( const int& x ){ m_sub3_mode = x; };

  void Mass_B( const double& x ){ m_Mass_B = x; };
  void Delta_E( const double& x ){ m_Delta_E = x; };
  void cos_thetaT( const double& x ){ m_cosT = x; };
  void cos_thetaB( const double& x ){ m_cosB = x; };
  void mass( const double& x ){ m_mass = x; };
  void masschisq( const double& x ){ m_masschisq = x;};
  void vtxchisq( const double& x ){ m_vtxchisq = x;};
  void purity( const double&x ){ m_purity = x;}

  void flag_bestEach( const int& x ){ m_flag_bestEach = x; } ;
  void flag_best( const int& x ){ m_flag_best = x; } ;
  void flag_best2( const int& x ){ m_flag_best2 = x; } ;
  void flag_best3( const int& x ){ m_flag_best3 = x; } ;
  void flag_mc( const int& x ){ m_flag_mc = x; } ;

  //  get variables
  // ~~~~~~~~~~~~~~~
  const int& b_mode( void )             const { return m_b_mode; };
  const int& sub0_mode( void )          const { return m_sub0_mode; };
  const int& sub1_mode( void )          const { return m_sub1_mode; };
  const int& sub2_mode( void )          const { return m_sub2_mode; };
  const int& sub3_mode( void )          const { return m_sub3_mode; };

  const double& Mass_B( void )          const { return m_Mass_B; };
  const double& Delta_E( void )         const { return m_Delta_E; };
  const double& cos_thetaT( void )      const { return m_cosT; };
  const double& cos_thetaB( void )      const { return m_cosB; };
  const double& mass( void )            const { return m_mass; };
  const double& masschisq( void )       const { return m_masschisq; };
  const double& vtxchisq( void )        const { return m_vtxchisq; };
  const double& purity( void )          const { return m_purity; };

  const int& flag_bestEach( void )      const { return m_flag_bestEach; };
  const int& flag_best( void )          const { return m_flag_best; };
  const int& flag_best2( void )         const { return m_flag_best2; };
  const int& flag_best3( void )         const { return m_flag_best3; };
  const int& flag_mc( void )            const { return m_flag_mc; };

private :

  int    m_b_mode;
  int    m_sub0_mode;
  int    m_sub1_mode;
  int    m_sub2_mode;
  int    m_sub3_mode;

  double m_Mass_B;
  double m_Delta_E;
  double m_cosT;
  double m_cosB;
  double m_mass;
  double m_masschisq;
  double m_vtxchisq;
  double m_purity;

  int    m_flag_bestEach;
  int    m_flag_best;
  int    m_flag_best2;
  int    m_flag_best3;
  int    m_flag_mc;

};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif 
