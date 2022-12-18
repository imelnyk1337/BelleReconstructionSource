//
//
//
#include "belle.h"
#include "particle/Particle.h"
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"

#include "panther/panther.h"
#include MDST_H
#include HEPEVT_H

#include "exUserInfo.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


#define DEBUG_JPSIKS 1

//utilities
double
beam_energy_constraint(const Particle &b, double e = 8., double p = 3.5)
{
  //b -- generally B0 or B0B
  //e -- electron beam energy
  //p -- positron beam enegry
  HepLorentzVector boost_vector(0., 0., e-p, e+p);
  HepLorentzVector b0(b.p());
  b0.boost(-(boost_vector.boostVector()));
  double mass2 = e*p - b0.vect().mag2();
  double mass  = (mass2 > 0.) ? sqrt(mass2) :  -sqrt(-mass2);
  return mass;
}

bool
checkB0(const Particle &jpsi, const Particle ks)
{
  //if(jpsi.relation().child(0).relation().isIdenticalWith(jpsi.relation().child(1).relation(),PC_CHARGED))return false;
  if(jpsi.relation().child(0).relation().isIdenticalWith(ks.relation().child(0).relation(),PC_CHARGED))return false;
  if(jpsi.relation().child(0).relation().isIdenticalWith(ks.relation().child(1).relation(),PC_CHARGED))return false;
  if(jpsi.relation().child(1).relation().isIdenticalWith(ks.relation().child(0).relation(),PC_CHARGED))return false;
  if(jpsi.relation().child(1).relation().isIdenticalWith(ks.relation().child(1).relation(),PC_CHARGED))return false;
  //if(ks.relation().child(0).relation().isIdenticalWith(ks.relation().child(1).relation(),PC_CHARGED))return false;
  return true;
}

void
setUserInfo(Particle &p)
{
  p.userInfo(*(new UserInfo));
}

void 
setUserInfo(vector<Particle> &p)
{
  for(int i=0;i<p.size();++i)setUserInfo(p[i]);
}


// Module class
class test_module : public Module
{
public:
  test_module(void){};
  ~test_module(void){};
  void init(int *status){};
  void term(void){};
  void disp_stat(const char*){};
  void hist_def(void);
  void event(BelleEvent*, int*);
  void begin_run(BelleEvent*, int *status){};
  void end_run(BelleEvent*, int *status){};
  void other(int*, BelleEvent*, int*){};

private:
  BelleHistogram *m_jpsi, *m_vee, *m_b, *m_b_beam;
};


extern "C" Module_descr *mdcl_exjpsiks()
{
  test_module *module = new test_module;
  Module_descr *dscr = new Module_descr ( "exjpsiks", module );
  return dscr;
}


void 
test_module::hist_def(void)
{
  extern BelleTupleManager *BASF_Histogram;
  m_jpsi   = BASF_Histogram->histogram("Mass(J/psi - GeV)",100,2.5,3.5);
  m_vee    = BASF_Histogram->histogram("Mass(Vee - GeV)",100,0.,1.);
  m_b      = BASF_Histogram->histogram("Mass(B - GeV)",50,5.0,5.5);
  m_b_beam = BASF_Histogram->histogram("Mass(B.Beam - GeV)",50,5.0,5.5);
}


void 
test_module::event(BelleEvent *evptr, int *status)
{
  Mdst_charged_Manager &charged_mag = Mdst_charged_Manager::get_manager();
  Mdst_vee_Manager     &vee_mag     = Mdst_vee_Manager::get_manager();

  //list of Particle
  vector<Particle> electron;
  vector<Particle> positron;
  vector<Particle> muon_minus;
  vector<Particle> muon_plus;
  vector<Particle> pi_minus;
  vector<Particle> pi_plus;
  vector<Particle> kshort;
  vector<Particle> jpsi;
  vector<Particle> b0;

  //Particle Type
  Ptype ptype_pi_plus("PI+");
  Ptype ptype_pi_minus("PI-");
  Ptype ptype_muon_plus("MU+");
  Ptype ptype_muon_minus("MU-");
  Ptype ptype_positron("E+");
  Ptype ptype_electron("E-");
  Ptype ptype_psi("PSI");
  Ptype ptype_b0("B0");

  //fill muon and electron lists with MDST_Charged Data Base
  for(vector<Mdst_charged>::iterator i = charged_mag.begin();
      i != charged_mag.end(); ++i){
    if(i->muid() && (*i).muid().muon() >= 2){
      if((*i).charge() > 0.){
	  Particle tmp(*i, ptype_muon_plus);
	  muon_plus.push_back(tmp);
	}else{
	  Particle tmp(*i, ptype_muon_minus);
	  muon_minus.push_back(tmp);
	}
    }else{
      if((*i).charge() > 0.){
	positron.push_back(Particle(*i, ptype_positron));
	pi_plus.push_back(Particle(*i, ptype_pi_plus));
      }else{
	electron.push_back(Particle(*i, ptype_electron));
	pi_minus.push_back(Particle(*i, ptype_pi_minus));
      }
    }
  }

#if 1
  //set UserInfo if we want to use. 
  //(in this program, it is not important. ONLY example!)
  setUserInfo(muon_plus);
  for(int i=0;i<muon_plus.size();++i){
    dout(Debugout::INFO,"exjpsiks") << dynamic_cast<UserInfo&>(muon_plus[i].userInfo()).cl() << " --> ";
    dynamic_cast<UserInfo&>(muon_plus[i].userInfo()).cl(1.0);
    dout(Debugout::INFO,"exjpsiks") << dynamic_cast<UserInfo&>(muon_plus[i].userInfo()).cl() << std::endl;
  }
#endif

  //fill vee lists with MDST_Vee Data Base
  for(vector<Mdst_vee>::iterator i = vee_mag.begin();
      i != vee_mag.end(); ++i){
    m_vee->accumulate((*i).mass(),1.);
    if((*i).kind() == 1)kshort.push_back(Particle(*i));
  }

  //make J/Psi from electron
  for(vector<Particle>::iterator i = electron.begin();
      i != electron.end(); ++i){
    for(vector<Particle>::iterator j = positron.begin();
	j != positron.end(); ++j){
      double mass = ((*i).momentum().p() + (*j).momentum().p()).mag();
      m_jpsi->accumulate(mass,1.);
      if(3.09693 - 0.05 <= mass && mass <= 3.09693 + 0.05){ // rough
	Particle jpsi_cand((*i).momentum().p() + (*j).momentum().p(),ptype_psi);
	jpsi_cand.relation().append(*i);
	jpsi_cand.relation().append(*j);
	jpsi.push_back(jpsi_cand);
      }
    }
  }

  //make J/Psi from muon
  for(vector<Particle>::iterator i = muon_minus.begin();
      i != muon_minus.end(); ++i){
    for(vector<Particle>::iterator j = muon_plus.begin();
	j != muon_plus.end(); ++j){
      double mass = ((*i).momentum().p() + (*j).momentum().p()).mag();
      m_jpsi->accumulate(mass,1.);
      if(3.09693 - 0.05 <= mass && mass <= 3.09693 + 0.05){ // rough
	Particle jpsi_cand((*i).momentum().p() + (*j).momentum().p(),ptype_psi);
	jpsi_cand.relation().append(*i);
	jpsi_cand.relation().append(*j);
	jpsi.push_back(jpsi_cand);
      }
    }
  }

  //make B0 from J/Psi and Vee
  for(vector<Particle>::iterator i = kshort.begin();
      i != kshort.end(); ++i){
    for(vector<Particle>::iterator j = jpsi.begin();
	j != jpsi.end(); ++j){
      if(checkB0(*j,*i)){
	double mass = ((*i).momentum().p() + (*j).momentum().p()).mag();
	m_b->accumulate(mass,1.);
	Particle b0_cand((*i).momentum().p() + (*j).momentum().p(),ptype_b0);
	b0_cand.relation().append(*i);
	b0_cand.relation().append(*j);
	b0.push_back(b0_cand);
      }
    }
  }

#if DEBUG_JPSIKS
  if(b0.size() != 0){
    dout(Debugout::INFO,"exjpsiks") << "B0 particle # = " << b0.size() << std::endl;
    for(vector<Particle>::iterator i = b0.begin();
	i != b0.end(); ++i){
      (*i).dump("mass momentum recursive");
      (*i).relation().dump("recursive");
    }
  }
#endif

  //beam constraint
  for(vector<Particle>::iterator i = b0.begin();
      i != b0.end(); ++i){
    double mass = beam_energy_constraint(*i);
    if(mass >= 0.)m_b_beam->accumulate(mass,1.);
  }
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
