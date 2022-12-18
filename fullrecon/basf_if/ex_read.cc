// ===================================================================
//  File Name : ex_read.cc
//  ------------------------------------------------------------------
//  Creation    ; 2004.04.07
//  Description ; example program for analysis with fullrecon table
//  Author      ; T.Matsumoto (matumot@bmail.kek.jp)
// -------------------------------------------------------------------

#include <vector>
#include <iosfwd>

#include "belle.h"

#include "basf/module.h"
#include "basf/module_descr.h"
#include "tuple/BelleTupleManager.h"

#include "particle/Particle.h"
#include "hamlet/AnaBrecon.h"
#include "fullrecon/frec_util.h"
#include "ip/IpProfile.h"

#include BRECON_H
#include FULLRECON_H

#include "panther/panther.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


using namespace fullreconNS;

// Class definition
class ex_read : public Module {
public:

  ex_read ( void );
  ~ex_read ( void ){};
  void init ( int *status ){};
  void term ( void ){};
  void disp_stat ( const char* ){};
  void hist_def ( void );
  void event ( BelleEvent*, int* );
  void begin_run ( BelleEvent*, int *status );
  void end_run ( BelleEvent*, int *status ){};
  void other ( int*, BelleEvent*, int* ){};

public:
  
  int m_debug;
  int m_version;

private:

  int m_ACCq;
  int m_TOFq;
  int m_CDCq;

  double m_PID_cut_K;
  double m_PID_cut_PI;

  double m_Dr_cut;
  double m_Dz_cut;

  double m_cut_dE_low_ver1;
  double m_cut_dE_high_ver1;
  double m_cut_dE_low;
  double m_cut_dE_high;
  double m_cut_Mbc_low;
  double m_cut_Mbc_high;

  BelleHistogram *H[200];

};

// BASF Interface Function 
extern "C" Module_descr *mdcl_ex_read ()
{

  /* main */
  
  ex_read *module = new ex_read;
  Module_descr *dscr = new Module_descr ( "ex_read", module );

  dscr->define_param ( "debug", "debug_parameter (0:off, 1:on)", 
		       &module->m_debug );
  dscr->define_param ( "version", "version number", 
		       &module->m_version );

  return dscr;
  
}

// Member Functions

// Constructor
ex_read::ex_read ( void ) {

  m_debug   = 0;
  m_version = 3;

  m_ACCq = 3;
  m_TOFq = 1;
  m_CDCq = 5;

  m_PID_cut_K  = 0.7;
  m_PID_cut_PI = 0.7;

  m_Dr_cut = 2.0;   // cm
  m_Dz_cut = 5.0;   // cm

  m_cut_dE_low_ver1   = -0.05;
  m_cut_dE_high_ver1  =  0.05;

  m_cut_dE_low   = -0.08;
  m_cut_dE_high  =  0.06;

  m_cut_Mbc_low  =  5.27;
  m_cut_Mbc_high =  5.29;

}

void ex_read::begin_run( BelleEvent*, int * ){

  IpProfile::begin_run();
  dout(Debugout::INFO,"ex_read") << "corrected IP -->> " << IpProfile::position() << std::endl;

  return;

}

// hist_def function
void ex_read::hist_def ( void ) {

  extern BelleTupleManager* BASF_Histogram;
  BelleTupleManager* tm = BASF_Histogram;

  int    id;
  int    bin,  bin_Mbc(50),   bin_dE(60); 
  double low,  Mbc_low(5.2),  dE_low(-0.3);
  double high, Mbc_high(5.3), dE_high(0.3);

  // Counter
  id = 1;
  H[id] =tm->histogram("Event Counter (All)",             6, -0.5, 5.5, id );
  id = 2;
  H[id] =tm->histogram("Ncand B    ",  101, -0.5, 100.5, id );
  id = 3;
  H[id] =tm->histogram("Nentry Brecon",  200, 0, 1000, id );

  id = 4;
  H[id] =tm->histogram("Event Counter (bestB(chisq), neutral)",  6, -0.5, 5.5, id );
  id = 5;
  H[id] =tm->histogram("Event Counter (bestB(chisq), charged)",  6, -0.5, 5.5, id );
  id = 6;
  H[id] =tm->histogram("Event Counter (bestB(purity), neutral)",  6, -0.5, 5.5, id );
  id = 7;
  H[id] =tm->histogram("Event Counter (bestB(purity), charged)",  6, -0.5, 5.5, id );

  // Mbc, deltaE
  id = 10;
  H[id]=tm->histogram("Mbc vs. deltaE (neutral B) ", 
		      bin_Mbc, Mbc_low, Mbc_high, bin_dE, dE_low, dE_high, id);
  id = 11;
  H[id]=tm->histogram("deltaE (neutral B) ", bin_dE, dE_low, dE_high, id);
  id = 12;
  H[id]=tm->histogram("Mbc with best chisq  (neutral B) ", bin_Mbc, Mbc_low, Mbc_high, id);
  id = 13;
  H[id]=tm->histogram("Mbc with best purity (neutral B) ", bin_Mbc, Mbc_low, Mbc_high, id);
  id = 14;
  H[id]=tm->histogram("Mbc(sideband) with best chisq  (neutral B) ", bin_Mbc, Mbc_low, Mbc_high, id);

  id = 20;
  H[id]=tm->histogram("Mbc vs. deltaE (charge B) ", 
		      bin_Mbc, Mbc_low, Mbc_high, bin_dE, dE_low, dE_high, id);
  id = 21;
  H[id]=tm->histogram("deltaE (charged B) ", bin_dE, dE_low, dE_high, id);
  id = 22;
  H[id]=tm->histogram("Mbc with best chisq  (charged B) ", bin_Mbc, Mbc_low, Mbc_high, id);
  id = 23;
  H[id]=tm->histogram("Mbc with best purity (charged B) ", bin_Mbc, Mbc_low, Mbc_high, id);
  id = 24;
  H[id]=tm->histogram("Mbc(sideband) with best chisq (charged B) ", bin_Mbc, Mbc_low, Mbc_high, id);


  // Information for recoil B side
  bin = 16, low = -0.5; high = 15.5;
  id = 31;
  H[id]=tm->histogram("pi- multiplicity in B0B", bin, low, high, id );
  id = 32;
  H[id]=tm->histogram("pi+ multiplicity in B0B", bin, low, high, id );
  id = 33;
  H[id]=tm->histogram("pi- multiplicity in B-",  bin, low, high, id );
  id = 34;
  H[id]=tm->histogram("pi+ multiplicity in B-",  bin, low, high, id );
  id = 41;
  H[id]=tm->histogram("K- multiplicity in B0B",  bin, low, high, id );
  id = 42;
  H[id]=tm->histogram("K+ multiplicity in B0B",  bin, low, high, id );
  id = 43;
  H[id]=tm->histogram("K- multiplicity in B-",  bin, low, high, id );
  id = 44;
  H[id]=tm->histogram("K+ multiplicity in B-",  bin, low, high, id );

  bin = 150, low = 0.; high = 3.0;
  id = 51;
  H[id]=tm->histogram("pi- Pcm in B0B",  bin, low, high, id );
  id = 52;
  H[id]=tm->histogram("pi+ Pcm in B0B",  bin, low, high, id );
  id = 53;
  H[id]=tm->histogram("pi- Pcm in B-",  bin, low, high, id );
  id = 54;
  H[id]=tm->histogram("pi+ Pcm in B-",  bin, low, high, id );
  id = 61;
  H[id]=tm->histogram("K- Pcm in B0B",  bin, low, high, id );
  id = 62;
  H[id]=tm->histogram("K+ Pcm in B0B",  bin, low, high, id );
  id = 63;
  H[id]=tm->histogram("K- Pcm in B-",  bin, low, high, id );
  id = 64;
  H[id]=tm->histogram("K+ Pcm in B-",  bin, low, high, id );

  bin = 100, low = Ptype_D0.mass() - 0.1; high = Ptype_D0.mass() + 0.1;
  id = 71;
  H[id]=tm->histogram("D0(K-pi+)  in B0B",  bin, low, high, id );
  id = 72;
  H[id]=tm->histogram("D0B(K+pi-) in B0B",  bin, low, high, id );
  id = 73;
  H[id]=tm->histogram("D0(K-pi+)  in B-",  bin, low, high, id );
  id = 74;
  H[id]=tm->histogram("D0B(K+pi-) in B-",  bin, low, high, id );

  bin = 100, low = Ptype_Dplus.mass() - 0.1; high = Ptype_Dplus.mass() + 0.1;
  id = 81;
  H[id]=tm->histogram("D+(K-pi+pi+) in B0B",  bin, low, high, id );
  id = 82;
  H[id]=tm->histogram("D-(K+pi-pi-) in B0B",  bin, low, high, id );
  id = 83;
  H[id]=tm->histogram("D+(K-pi+pi+) in B-",  bin, low, high, id );
  id = 84;
  H[id]=tm->histogram("D-(K+pi-pi-) in B-",  bin, low, high, id );

  bin = 225, low = 0.; high = 4.5;
  id = 91;
  H[id]=tm->histogram("M(X) in B0B -> Xpi-",  bin, low, high, id );
  id = 92;
  H[id]=tm->histogram("M(X) in B0B -> Xpi+",  bin, low, high, id );
  id = 93;
  H[id]=tm->histogram("M(X) in B-  -> Xpi-",  bin, low, high, id );
  id = 94;
  H[id]=tm->histogram("M(X) in B-  -> Xpi+",  bin, low, high, id );

  id = 101;
  H[id]=tm->histogram("M(X) in B0B -> Xpi- (Ebeam constraint)",  bin, low, high, id );
  id = 102;
  H[id]=tm->histogram("M(X) in B0B -> Xpi+ (Ebeam constraint)",  bin, low, high, id );
  id = 103;
  H[id]=tm->histogram("M(X) in B-  -> Xpi- (Ebeam constraint)",  bin, low, high, id );
  id = 104;
  H[id]=tm->histogram("M(X) in B-  -> Xpi+ (Ebeam constraint)",  bin, low, high, id );

}

// event function
void ex_read::event ( BelleEvent* evptr, int* status ) {
  
  Fullrecon_Manager     & frMgr  = Fullrecon_Manager::get_manager();
  Brecon_header_Manager & bhMgr  = Brecon_header_Manager::get_manager();
  Brecon_Manager        & brMgr  = Brecon_Manager::get_manager();

  Mdst_charged_Manager  & chgMgr = Mdst_charged_Manager::get_manager();
  if ( frMgr.count() <= 0 || bhMgr.count() <= 0 ) return;

  // Counter
  // ~~~~~~~
  H[1]->accumulate(  1, (float) 1. );
  int nB = bhMgr[0].nBrecon();
  H[2]->accumulate( nB, (float) 1. );
  int nBrecon = brMgr.count();
  H[3]->accumulate( nBrecon, (float) 1. );

  // Select best B candidates from Fullrecon table
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  std::vector<Particle> bestB_chisq_neutral;
  std::vector<Particle> bestB_chisq_charged;
  std::vector<Particle> bestB_purity_neutral;
  std::vector<Particle> bestB_purity_charged;
  std::vector<Particle> bestB;

  AnaBrecon brecon;

  for( std::vector<Fullrecon>::iterator f_i = frMgr.begin();
       f_i != frMgr.end(); f_i++ ){

    Fullrecon& frecB = *f_i;

    int kind      = frecB.kind();
    // select B to charm decays
    if ( kind != 1 ) continue;

    double deltaE = frecB.deltE();
    double Mbc    = frecB.Mbc();

    bool   flag_select_deltaE    = false;
    bool   flag_select_Mbc       = false;
    bool   flag_select           = false;
    bool   flag_sideband_deltaE  = false;

    if ( m_version > 1 ){
      if ( deltaE > m_cut_dE_low  && deltaE < m_cut_dE_high  ) flag_select_deltaE   = true;
    }else{
      if ( deltaE > m_cut_dE_low_ver1  && deltaE < m_cut_dE_high_ver1  ) flag_select_deltaE = true;
    }

    if ( deltaE > 0.1           && deltaE < 0.3  )           flag_sideband_deltaE = true;
    if ( Mbc    > m_cut_Mbc_low && Mbc    < m_cut_Mbc_high ) flag_select_Mbc      = true;
    if ( flag_select_deltaE  && flag_select_Mbc )            flag_select          = true;

    int flag_bestPurity(0);
    int flag_bestChisq(0);

    if ( m_version > 1 ){
      flag_bestPurity = frecB.flag0();
      flag_bestChisq  = frecB.flag1();
    }else if ( m_version == 1 ){
      flag_bestChisq = frecB.flag0();
    }

    if ( !(flag_bestChisq || flag_bestChisq) ) continue;

    bool flagB_charged = false;
    bool flagB_neutral = false;

    const Particle & bcand = brecon.getParticle( (int)frecB.get_ID() );

    if ( abs(bcand.lund()) == Ptype_Bplus.lund() )   flagB_charged = true;
    else if ( abs(bcand.lund()) == Ptype_B0.lund() ) flagB_neutral = true;

    if ( flag_bestPurity ){

      if ( flag_select_deltaE ){
	if ( flagB_neutral ) H[13]->accumulate( Mbc, (float)1. );
	if ( flagB_charged ) H[23]->accumulate( Mbc, (float)1. );
      }

      if ( flag_select ){
	if ( flagB_neutral )      bestB_purity_neutral.push_back(bcand);
	else if ( flagB_charged ) bestB_purity_charged.push_back(bcand);
      }

    }

    if ( flag_bestChisq != 1 ) continue;

      if ( flagB_neutral ){
        H[10]->accumulate( Mbc, deltaE, (float)1. );
	if ( flag_select_Mbc )      H[11]->accumulate( deltaE, (float)1. );
	if ( flag_select_deltaE )   H[12]->accumulate( Mbc, (float)1. );
	if ( flag_sideband_deltaE ) H[14]->accumulate( Mbc, (float)1. );
      }

      if ( flagB_charged ){
        H[20]->accumulate( Mbc, deltaE, (float)1. );
	if ( flag_select_Mbc )       H[21]->accumulate( deltaE, (float)1. );
	if ( flag_select_deltaE )    H[22]->accumulate( Mbc, (float)1. );
	if ( flag_sideband_deltaE )  H[24]->accumulate( Mbc, (float)1. );
      }

      if ( flag_select ){
	if ( flagB_neutral )       bestB_chisq_neutral.push_back(bcand);
	else if ( flagB_charged )  bestB_chisq_charged.push_back(bcand);
      }

  } // end of loop f_i

  if ( bestB_chisq_neutral.size() )  H[4]->accumulate( 1., (float)1. );
  if ( bestB_chisq_charged.size() )  H[5]->accumulate( 1., (float)1. );
  if ( bestB_purity_neutral.size() ) H[6]->accumulate( 1., (float)1. );
  if ( bestB_purity_charged.size() ) H[7]->accumulate( 1., (float)1. );

  // use bestB candidates with best purity or best chisq
  if ( m_version > 1 ){
    //    if ( bestB_purity_neutral.size() ) bestB.push_back( bestB_purity_neutral[0] );
    //    if ( bestB_purity_charged.size() ) bestB.push_back( bestB_purity_charged[0] );
    if ( bestB_chisq_neutral.size() ) bestB.push_back( bestB_chisq_neutral[0] );
    if ( bestB_chisq_charged.size() ) bestB.push_back( bestB_chisq_charged[0] );
  }else if ( m_version == 1 ){
    if ( bestB_chisq_neutral.size() ) bestB.push_back( bestB_chisq_neutral[0] );
    if ( bestB_chisq_charged.size() ) bestB.push_back( bestB_chisq_charged[0] );
  }

  // Analyze recoil B side [ select pi, K and 
  // check inclusive pi,k and D0(Kpi), D+(Kpipi) and missing mass in Xpi ]
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Ptype ptype;
  const atc_pid selkpi( m_ACCq, m_TOFq, m_CDCq, KAON_CODE, PION_CODE );

  for( std::vector<Particle>::iterator b_i = bestB.begin();
       b_i != bestB.end(); b_i++ ){

    if ( m_debug )  b_i->dump( "recursive" );

    bool flagB_neutral = false;
    bool flagB_charged = false;

    Particle & bcand = *b_i;

    Particle brecoil( frec_util::p_beam() - bcand.p(), Ptype(-bcand.lund()) );

    if ( abs(brecoil.lund()) == Ptype_B0.lund() )    flagB_neutral = true;
    if ( abs(brecoil.lund()) == Ptype_Bplus.lund() ) flagB_charged = true;

    std::vector<Particle> PI_all;
    std::vector<Particle> K_all;

    // reconstruct pion
    for ( std::vector<Mdst_charged>::const_iterator i = chgMgr.begin();
	  i != chgMgr.end(); i++ ){

      const double dr = frec_util::correct_dr(*i, IpProfile::position(), PION_CODE);
      const double dz = frec_util::correct_dr(*i, IpProfile::position(), PION_CODE);

      if ( fabs(dr) > m_Dr_cut || fabs(dz) > m_Dz_cut ) continue;

      if ( (*i).charge() > 0 ) ptype = Ptype_PIplus;
      else                     ptype = Ptype_PIminus;
      Particle pi(*i, ptype );
      if ( frec_util::sel_KPI( pi, m_PID_cut_PI, selkpi ) == false ) continue;
      if ( frec_util::prohib_dupli( pi, bcand ) == true ) continue;

      PI_all.push_back(pi);

    }

    // reconstruct kaon
    for ( std::vector<Mdst_charged>::const_iterator i = chgMgr.begin();
	  i != chgMgr.end(); i++ ){

      const double dr = frec_util::correct_dr(*i, IpProfile::position(), KAON_CODE);
      const double dz = frec_util::correct_dr(*i, IpProfile::position(), KAON_CODE);

      if ( fabs(dr) > m_Dr_cut || fabs(dz) > m_Dz_cut ) continue;

      if ( (*i).charge() > 0 ) ptype = Ptype_Kplus;
      else                     ptype = Ptype_Kminus;
      Particle k(*i, ptype );
      if ( frec_util::sel_KPI( k, m_PID_cut_K, selkpi ) == false ) continue;
      if ( frec_util::prohib_dupli( k, bcand ) == true ) continue;

      K_all.push_back(k);

    }

   int count_b2piMinus(0);
   int count_b2piPlus(0);
   for ( std::vector<Particle>::iterator pi_i = PI_all.begin();
	 pi_i != PI_all.end(); pi_i++ ){

     double pcm = frec_util::p_cm( pi_i->p() ).vect().mag();

     bool flag_b2piMinus = false;
     bool flag_b2piPlus  = false;
     if ( brecoil.lund() * pi_i->lund() > 0 ) flag_b2piMinus = true;
     else                                     flag_b2piPlus  = true;

     if ( flag_b2piMinus ){
       count_b2piMinus += 1;
       if ( flagB_neutral ) H[51]->accumulate( pcm, (float)1. );
       if ( flagB_charged ) H[53]->accumulate( pcm, (float)1. );
     }

     if ( flag_b2piPlus ){
       count_b2piPlus += 1;
       if ( flagB_neutral ) H[52]->accumulate( pcm, (float)1. );
       if ( flagB_charged ) H[54]->accumulate( pcm, (float)1. );
     }

   }

   if ( flagB_neutral ){
     H[31]->accumulate( count_b2piMinus, (float)1. );
     H[32]->accumulate( count_b2piPlus,  (float)1. );
   }
   if ( flagB_charged ){
     H[33]->accumulate( count_b2piMinus, (float)1. );
     H[34]->accumulate( count_b2piPlus,  (float)1. );
   }

   int count_b2kMinus(0);
   int count_b2kPlus(0);
   for ( std::vector<Particle>::iterator k_i = K_all.begin();
	 k_i != K_all.end(); k_i++ ){

     double pcm = frec_util::p_cm( k_i->p() ).vect().mag();

     bool flag_b2kMinus = false;
     bool flag_b2kPlus  = false;
     if ( brecoil.lund() * k_i->lund() > 0 )  flag_b2kMinus = true;
     else                                     flag_b2kPlus  = true;

     if ( flag_b2kMinus ){
       count_b2kMinus += 1;
       if ( flagB_neutral ) H[61]->accumulate( pcm, (float)1. );
       if ( flagB_charged ) H[63]->accumulate( pcm, (float)1. );
     }

     if ( flag_b2kPlus ){
       count_b2kPlus += 1;
       if ( flagB_neutral ) H[62]->accumulate( pcm, (float)1. );
       if ( flagB_charged ) H[64]->accumulate( pcm, (float)1. );
     }

   }

   if ( flagB_neutral ){
     H[41]->accumulate( count_b2kMinus, (float)1. );
     H[42]->accumulate( count_b2kPlus,  (float)1. );
   }
   if ( flagB_charged ){
     H[43]->accumulate( count_b2kMinus, (float)1. );
     H[44]->accumulate( count_b2kPlus,  (float)1. );
   }

    // reconstruct D-->Kpi and Kpipi
    double dmass;
    for( std::vector<Particle>::iterator k_i = K_all.begin();
	 k_i != K_all.end(); k_i++ ){

      const float chg_K = k_i->charge();

      for ( std::vector<Particle>::iterator pi_j = PI_all.begin();
	    pi_j != PI_all.end(); pi_j++ ){
   
        const float chg_PI = pi_j->charge();

	if ( chg_K < 0 ) ptype = Ptype_D0;
	else             ptype = Ptype_D0B;
	dmass = (k_i->p() + pi_j->p()).mag();

	bool flag_b2d0  = false;
	bool flag_b2d0b = false;
	if ( brecoil.lund() * ptype.lund() < 0 ) flag_b2d0  = true;
	else                                     flag_b2d0b = true;

	if ( flag_b2d0 ){
	  if ( flagB_neutral ) H[71]->accumulate( dmass, (float)1. );
          if ( flagB_charged ) H[73]->accumulate( dmass, (float)1. );
        }
	if ( flag_b2d0b ){
	  if ( flagB_neutral ) H[72]->accumulate( dmass, (float)1. );
          if ( flagB_charged ) H[74]->accumulate( dmass, (float)1. );
        }

	for ( std::vector<Particle>::iterator pi_k = pi_j+1;
	      pi_k != PI_all.end(); pi_k++ ){

	  if ( chg_PI * pi_k->charge() < 0 ) continue;

	  if ( chg_K < 0 ) ptype = Ptype_Dplus;
	  else             ptype = Ptype_Dminus;
	  dmass = (k_i->p() + pi_j->p() + pi_k->p()).mag();

  	  bool flag_b2dminus  = false;
	  bool flag_b2dplus   = false;
	  if ( brecoil.lund() * ptype.lund() > 0 ) flag_b2dminus = true;
	  else                                     flag_b2dplus  = true;

	if ( flag_b2dplus ){
	  if ( flagB_neutral ) H[81]->accumulate( dmass, (float)1. );
          if ( flagB_charged ) H[83]->accumulate( dmass, (float)1. );
        }

	if ( flag_b2dminus ){
	  if ( flagB_neutral ) H[82]->accumulate( dmass, (float)1. );
          if ( flagB_charged ) H[84]->accumulate( dmass, (float)1. );
        }

      } // end of loop pi_k
    } // end of loop pi_j
   } // end of loop pi_i

   /// reconstruct M(X) in B-->Xpi

   // pcm for brecoil with Ebeam constraint
   HepLorentzVector pcm_brecoil( - frec_util::p_cm(bcand.p()).vect(), frec_util::Ebeam() );

   for ( std::vector<Particle>::iterator pi_i = PI_all.begin();
	      pi_i != PI_all.end(); pi_i++ ){

     int nsvdHit_rphi = pi_i->mdstCharged().trk().mhyp(PION_CODE).nhits(3);
     int nsvdHit_z    = pi_i->mdstCharged().trk().mhyp(PION_CODE).nhits(4);

     bool flag_svdHit = false;
     if ( nsvdHit_rphi >= 2 && nsvdHit_z >= 2 ) flag_svdHit = true;

     double Mx  = (brecoil.p() - pi_i->p()).mag();
     double Mx2 = (pcm_brecoil - frec_util::p_cm(pi_i->p())).mag();

     bool flag_b2piMinus = false;
     bool flag_b2piPlus  = false;
     
     if ( brecoil.lund() * pi_i->lund() > 0 )  flag_b2piMinus = true;
     else                                      flag_b2piPlus  = true;

     if ( flag_svdHit ){

     if ( flag_b2piMinus ){
	if ( flagB_neutral ) {
	  H[91]->accumulate( Mx, (float)1. );
	  H[101]->accumulate( Mx2, (float)1. );
	}
        if ( flagB_charged ){
	  H[93]->accumulate( Mx, (float)1. );
	  H[103]->accumulate( Mx2, (float)1. );
	}
     }
     
     if ( flag_b2piPlus ){
	if ( flagB_neutral ){
	  H[92]->accumulate( Mx, (float)1. );
	  H[102]->accumulate( Mx2, (float)1. );
	}
        if ( flagB_charged ){
	  H[94]->accumulate( Mx, (float)1. );
	  H[104]->accumulate( Mx2, (float)1. );
	}
     }

     }

   } // end of loop pi_i
  } // end of loop b_i


  return;

}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
