#include "particle/Particle.h"
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"

#define FOR_FSIM 1

#include "belle.h"
#include MDST_H
#if FOR_FSIM
#include MDST_OBS_H
#endif
#include HEPEVT_H

#include "particle/combination.h"
#include "particle/utility.h"
#include "tables/belletdf.h"
#include "tables/evtcls.h"
#include "mdst/mdst.h"
#include "tables/evtvtx.h"

#include "kfitter/kmakemother.h"

#include "userinfo.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


const double BF = 1.5;

#if FOR_FSIM
void 
setGenHepInfoGfsim(std::vector<Particle> &list){
  //...gamma
  if(list.size() == 0)return;

  Gen_hepevt_Manager    &genMgr  = Gen_hepevt_Manager::get_manager();
  Mdst_sim_xref_Manager &xrefMgr = Mdst_sim_xref_Manager::get_manager();
  for(std::vector<Particle>::iterator i = list.begin();
      i != list.end(); ++i){
    //...mdst_gamma
    if(i->mdstGamma()){
      for(Mdst_sim_xref_Manager::iterator j = xrefMgr.begin();
          j != xrefMgr.end(); ++j){
        if(j->gamma() && 
           i->mdstGamma().get_ID() == j->gamma().get_ID() && 
           j->hepevt() &&
           i->pType().lund() == j->hepevt().idhep()){
          i->relation().genHepevt(genMgr[(int)(j->hepevt().get_ID())-1]);
          break;
        }
      }
    }
  }
}

void 
setGenHepInfoPfsim(Particle &p){
  Gen_hepevt_Manager    &genMgr  = Gen_hepevt_Manager::get_manager();
  Mdst_sim_xref_Manager &xrefMgr = Mdst_sim_xref_Manager::get_manager();
  if( !p.child(0) || !p.child(1) ) return;
  if( !p.child(0).mdstGamma() || !p.child(1).mdstGamma() ) return;
  if(  p.child(0).mdstGamma() ==  p.child(1).mdstGamma() ) return;

  for(Mdst_sim_xref_Manager::iterator i = xrefMgr.begin();
      i != xrefMgr.end(); ++i){
    if( !i->gamma() || !i->hepevt() || (i->hepevt().idhep() != 22) ) continue;
    if( p.child(0).mdstGamma() == i->gamma() )
      p.child(0).relation().genHepevt(i->hepevt());
    if( p.child(1).mdstGamma() == i->gamma() )
      p.child(1).relation().genHepevt(i->hepevt());
    if( p.child(0).genHepevt() && p.child(1).genHepevt() )
      if( p.child(0).genHepevt() !=  p.child(1).genHepevt() &&
          p.child(0).genHepevt().mother() == 
          p.child(1).genHepevt().mother() &&
          p.child(0).genHepevt().mother().idhep() == 111 )
        p.relation().genHepevt(p.child(0).genHepevt().mother());
  }
}

void 
setGenHepInfoPfsim(std::vector<Particle> &p_list){
  for(int i=0; i<p_list.size(); ++i)
    setGenHepInfoPfsim(p_list[i]);
}
#endif

class d0kpipi0_mv_module : public Module
{
public:
  d0kpipi0_mv_module(void);
  ~d0kpipi0_mv_module(void){};
  void init(int*);
  void term(void){};
  void disp_stat(const char*){};
  void hist_def(void);
  void event(BelleEvent*, int*);
  void begin_run(BelleEvent*, int*);
  void end_run(BelleEvent*, int*){};
  void other(int*, BelleEvent*, int*){};

public:

private:
  void endEvent(void);
  void fitD0Event(vector<Particle>&);
  void writeHisto(vector<Particle>&, BelleTuple*);
  void fitPI0Event(vector<Particle>&, BelleTuple*);

private:
  // PType
  Ptype m_ptypeD0;
  Ptype m_ptypePI0;

  vector<Particle> m_kaonP, m_kaonM;
  vector<Particle> m_pionP, m_pionM;
  vector<Particle> m_gamma;

private:
  // Histgram and NTuple
  BelleHistogram *m_event_counter;

  BelleTuple     *m_hist;
  BelleTuple     *m_hist2;

  int m_expNo, m_runNo, m_evtNo, m_farmNo;
  int m_mc;
  double m_r2;
};

extern "C" Module_descr *mdcl_d0kpipi0_mv()
{
  d0kpipi0_mv_module *module = new d0kpipi0_mv_module;
  Module_descr *dscr = new Module_descr ( "d0kpipi0_mv", module );
  return dscr;
}

d0kpipi0_mv_module::d0kpipi0_mv_module(void)
  : m_ptypeD0("D0"),
    m_ptypePI0("PI0"),
    m_mc(1)
{
}

void 
d0kpipi0_mv_module::hist_def(void)
{
  extern BelleTupleManager *BASF_Histogram;

  m_event_counter = BASF_Histogram->histogram("Event Counter", 11, -5.5, 5.5);

  m_hist = BASF_Histogram->ntuple("Event Info", 
				  // 7
				  "mc lund expN runN evtN farmN R2 "
				  // 12
				  "ga0px0 ga0py0 ga0pz0 "
				  "ga0epx02 ga0epy02 ga0epz02 "
				  "ga0px1 ga0py1 ga0pz1 "
				  "ga0epx2 ga0epy2 ga0epz2 "
				  // 12
				  "ga1px0 ga1py0 ga1pz0 "
				  "ga1epx02 ga1epy02 ga1epz02 "
				  "ga1px1 ga1py1 ga1pz1 "
				  "ga1epx2 ga1epy2 ga1epz2 "
				  // 17
				  "pi0cl pi0chisq pi0ndf "
				  "pi0m0 pi0px0 pi0py0 pi0pz0 "
				  "pi0epx02 pi0epy02 pi0epz02 "
				  "pi0m1 pi0px1 pi0py1 pi0pz1 "
				  "pi0epx2 pi0epy2 pi0epz2 "
				  // 20
				  "cl chisq ndf "
				  "m0 px0 py0 pz0 "
				  "m1 px1 py1 pz1 "
				  "epx2 epy2 epz2 "
				  "vx vy vz "
				  "evx2 evy2 evz2 "
				  // 6
				  "vx0 vy0 vz0 "
				  "vx1 vy1 vz1 "
				  // 5
				  "gga0 "
				  "gga0px gga0py gga0pz gga0m "
				  // 5
				  "gga1 "
				  "gga1px gga1py gga1pz gga1m "
				  // 11
				  "gpi0 "
				  "gpi0px gpi0py gpi0pz gpi0m "
				  "gpi0dvx gpi0dvy gpi0dvz "
				  "gpi0pvx gpi0pvy gpi0pvz "
				  // 11
				  "gd0 "
				  "gpx gpy gpz gm "
				  "gdvx gdvy gdvz "
				  "gpvx gpvy gpvz "
				  );
  m_hist2= BASF_Histogram->ntuple("Event Info2", 
				  // 7
				  "mc lund expN runN evtN farmN R2 "
				  // 12
				  "ga0px0 ga0py0 ga0pz0 "
				  "ga0epx02 ga0epy02 ga0epz02 "
				  "ga0px1 ga0py1 ga0pz1 "
				  "ga0epx2 ga0epy2 ga0epz2 "
				  // 12
				  "ga1px0 ga1py0 ga1pz0 "
				  "ga1epx02 ga1epy02 ga1epz02 "
				  "ga1px1 ga1py1 ga1pz1 "
				  "ga1epx2 ga1epy2 ga1epz2 "
				  // 17
				  "pi0cl pi0chisq pi0ndf "
				  "pi0m0 pi0px0 pi0py0 pi0pz0 "
				  "pi0epx02 pi0epy02 pi0epz02 "
				  "pi0m1 pi0px1 pi0py1 pi0pz1 "
				  "pi0epx2 pi0epy2 pi0epz2 "
				  // 5
				  "gga0 "
				  "gga0px gga0py gga0pz gga0m "
				  // 5
				  "gga1 "
				  "gga1px gga1py gga1pz gga1m "
				  // 11
				  "gpi0 "
				  "gpi0px gpi0py gpi0pz gpi0m "
				  );
}

void
d0kpipi0_mv_module::begin_run(BelleEvent *dummy,
			      int *dummyStatus)
{
  //eid::init_data();
  return;
}

void
d0kpipi0_mv_module::init(int *dummyStatus)
{
  //dummy
  Ptype dummy("E-");
}

void 
d0kpipi0_mv_module::endEvent(void)
{
  eraseVector(m_kaonP);
  eraseVector(m_kaonM);
  eraseVector(m_pionP);
  eraseVector(m_pionM);
  eraseVector(m_gamma);
}

void
d0kpipi0_mv_module::event(BelleEvent *evptr, int *status)
{
  *status = -1;
  m_event_counter->accumulate(0.,1.);

  // MC or Data
  if(m_mc){
    Belle_event_Manager &evtMgr = Belle_event_Manager::get_manager();
    if(evtMgr.count() != 0){
      if(evtMgr[0].ExpMC() != 2){
        m_mc = 0; // not MC
        dout(Debugout::INFO,"d0kpipi0_mv") << "(module) This analysis is in REAL DATA not MC DATA." << std::endl;
      }
    }
  }

  // Event Information
  m_expNo = m_runNo = m_evtNo = m_farmNo = 0;
  Belle_event_Manager &evtMgr = Belle_event_Manager::get_manager();
  if(evtMgr.count()){
    const int MASK28BIT = 0x0FFFFFFF;
    m_expNo  = evtMgr[0].ExpNo();
    m_runNo  = evtMgr[0].RunNo();
    m_evtNo  = (int)(evtMgr[0].EvtNo() & MASK28BIT);
    m_farmNo = evtMgr[0].EvtNo() >> 28;
  }

  // Event Shape
  m_r2 = -1.;
  Evtcls_hadron_info_Manager &clsMgr = Evtcls_hadron_info_Manager::get_manager();
  if(clsMgr.count())m_r2 = clsMgr[0].R2();

  // Make Charged Particles (with No Cut)
  makeKPi(m_kaonP, m_kaonM, m_pionP, m_pionM, 0);
  setUserInfo(m_kaonP);
  setUserInfo(m_kaonM);
  setUserInfo(m_pionP);
  setUserInfo(m_pionM);
  if(m_mc){
    setGenHepInfoF(m_kaonP);setGenHepInfoF(m_kaonM);
    setGenHepInfoF(m_pionP);setGenHepInfoF(m_pionM);
  }

#if FOR_FSIM
  // src/sim/fsim/src/feclset.F
  Mdst_ecl_Manager &eclMgr = Mdst_ecl_Manager::get_manager();
  const double theta1 = atan2(125.0, 192.5);
  const double theta2 = atan2(125.0, -94.5);
  //dout(Debugout::INFO,"d0kpipi0_mv") << "Theta : " << theta1 << ", " << theta2 << std::endl;
  for(Mdst_ecl_Manager::iterator j = eclMgr.begin();
      j != eclMgr.end(); ++j){
    double theta = j->theta();
    double Rt = 0.0; // cm
    if(theta1 < theta && theta < theta2){
      Rt = 125.0;
      //j->r(125.0); // cm
    }else if(theta <= theta1){
      Rt = 196.2*tan(theta);
      //j->r(fabs(196.2*tan(theta)));
    }else{
      Rt = -102.2*tan(theta);
      //j->r(fabs(-102.2*tan(theta)));
    }
    j->r(fabs(Rt)*sqrt(1.+1./tan(theta)/tan(theta)));
    //dout(Debugout::INFO,"d0kpipi0_mv") << j->r() << "cm" << std::endl;
    if(theta < 0.){
      dout(Debugout::INFO,"d0kpipi0_mv") << "Error ecl theta = " << theta << std::endl;
    }
  }
#endif

  // Make Gamma
  makeGamma(m_gamma);
  setUserInfo(m_gamma);
  if(m_mc){
#if FOR_FSIM
    setGenHepInfoGfsim(m_gamma);
#else
    setGenHepInfoG(m_gamma);
#endif
    /*
    for(int i=0;i<m_gamma.size();++i){
      if(m_gamma[i].genHepevt()){
	dout(Debugout::INFO,"d0kpipi0_mv") << "g" << i << " " << m_gamma[i].genHepevt().get_ID() << std::endl;
      }
    }
    */
  }

  // Make Pi0
  vector<Particle> PI0;
  combination(PI0, m_ptypePI0, m_gamma, m_gamma, 0.05);
  setUserInfo(PI0);
  if(m_mc){
    setGenHepInfoPfsim(PI0);
    /*
    for(int i=0;i<PI0.size();++i){
      if(PI0[i].genHepevt()){
	dout(Debugout::INFO,"d0kpipi0_mv") << "p" << i << " " << PI0[i].genHepevt().get_ID() << std::endl;
      }
    }
    */
  }

  fitPI0Event(PI0, m_hist2);

  // Make D0 (with No Cut)
  vector<Particle> D0;
  combination(D0, m_ptypeD0, m_kaonM, m_pionP, PI0, 0.15);
  setUserInfo(D0);
  if(m_mc){
    setGenHepInfoR(D0);
  }

  fitD0Event(D0);

  writeHisto(D0, m_hist);

  endEvent();
}

void
d0kpipi0_mv_module::fitD0Event(vector<Particle> &plist)
{
  for(int i=0;i<plist.size();++i){
    // find a vertex
    kvertexfitter kvf;
    kvf.magneticField(BF);
    addTrack2fit(kvf, plist[i].child(0));
    addTrack2fit(kvf, plist[i].child(1));
    unsigned err = kvf.fit();
    if(err != 0)continue;
    //dout(Debugout::INFO,"d0kpipi0_mv") << "0 " << kvf.vertex() << ", m =  " << plist[i].child(2).mass() << std::endl;

    // mass fit of PI0
    Particle g0 = plist[i].child(2).child(0);
    Particle g1 = plist[i].child(2).child(1);
    Particle tmpg0 = plist[i].child(2).child(0);
    Particle tmpg1 = plist[i].child(2).child(1);
    HepSymMatrix tmpErr(3,0);
    HepPoint3D org(0.,0.,0.);
    setGammaError(tmpg0, org, tmpErr);
    setGammaError(tmpg1, org, tmpErr);

    //tmpErr[0][0] = 1.; // 1cm^2
    //tmpErr[1][1] = 1.;
    //tmpErr[2][2] = 1.; --> wrong!!!
    /*
    dout(Debugout::INFO,"d0kpipi0_mv") << "Before" << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g0.momentum().p() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g0.momentum().x() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g0.momentum().dpx() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g1.momentum().p() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g1.momentum().x() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g1.momentum().dpx() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g0.mdstGamma().ecl().energy() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g0.mdstGamma().ecl().phi() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g0.mdstGamma().ecl().theta() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g0.mdstGamma().ecl().r() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g1.mdstGamma().ecl().energy() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g1.mdstGamma().ecl().phi() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g1.mdstGamma().ecl().theta() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g1.mdstGamma().ecl().r() << std::endl;
    */
    //if(0.575919 > g0.mdstGamma().ecl().theta())dout(Debugout::INFO,"d0kpipi0_mv") << "LLLLL" << std::endl;
    
    setGammaError(g0, kvf.vertex(), kvf.errVertex());
    setGammaError(g1, kvf.vertex(), kvf.errVertex());
    //setGammaError(g0, kvf.vertex(), tmpErr);
    //setGammaError(g1, kvf.vertex(), tmpErr);

    /*
    dout(Debugout::INFO,"d0kpipi0_mv") << "After" << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g0.momentum().p() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g0.momentum().x() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g0.momentum().dpx() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g1.momentum().p() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g1.momentum().x() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g1.momentum().dpx() << std::endl;
    */

    if(plist[i].genHepevt()){
      dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).px(tmpg0.momentum().p().x(),0);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).py(tmpg0.momentum().p().y(),0);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).pz(tmpg0.momentum().p().z(),0);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).epx(tmpg0.momentum().dp()[0][0],0);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).epy(tmpg0.momentum().dp()[1][1],0);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).epz(tmpg0.momentum().dp()[2][2],0);

      dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).px(tmpg1.momentum().p().x(),0);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).py(tmpg1.momentum().p().y(),0);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).pz(tmpg1.momentum().p().z(),0);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).epx(tmpg1.momentum().dp()[0][0],0);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).epy(tmpg1.momentum().dp()[1][1],0);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).epz(tmpg1.momentum().dp()[2][2],0);

      dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).px(g0.momentum().p().x(),1);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).py(g0.momentum().p().y(),1);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).pz(g0.momentum().p().z(),1);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).epx(g0.momentum().dp()[0][0],1);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).epy(g0.momentum().dp()[1][1],1);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).epz(g0.momentum().dp()[2][2],1);

      dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).px(g1.momentum().p().x(),1);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).py(g1.momentum().p().y(),1);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).pz(g1.momentum().p().z(),1);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).epx(g1.momentum().dp()[0][0],1);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).epy(g1.momentum().dp()[1][1],1);
      dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).epz(g1.momentum().dp()[2][2],1);
    }

    //setGammaError(g0, kvf.vertex(), tmpErr);
    //setGammaError(g1, kvf.vertex(), tmpErr);
    kmassfitter kmf;
    kmf.magneticField(BF);
    addTrack2fit(kmf, g0);
    addTrack2fit(kmf, g1);
    kmf.invariantMass(m_ptypePI0.mass());
    kmf.vertex(kvf.vertex());
    kmf.errVertex(kvf.errVertex());
    //kmf.atDecayPoint();
    err = kmf.fit();
    /*
    dout(Debugout::INFO,"d0kpipi0_mv") << "err = " << err << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g0.momentum().dpx() << std::endl;
    dout(Debugout::INFO,"d0kpipi0_mv") << g1.momentum().dpx() << std::endl;
    */
    if(err != 0)continue;
    //dout(Debugout::INFO,"d0kpipi0_mv") << "1 " << kvf.vertex() << " " << kmf.vertex() << std::endl;
    Particle pi0 = plist[i].child(2);
    kmakemother kmm2;
    kmm2.magneticField(BF);
    if(plist[i].genHepevt()){
      dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).mass(plist[i].child(2).mass(),0);
      dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).px(plist[i].child(2).p().x(),0);
      dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).py(plist[i].child(2).p().y(),0);
      dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).pz(plist[i].child(2).p().z(),0);
      dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).epx(plist[i].child(2).momentum().dp()[0][0],0);
      dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).epy(plist[i].child(2).momentum().dp()[1][1],0);
      dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).epz(plist[i].child(2).momentum().dp()[2][2],0);
    }
    makeMother(kmm2,kmf,pi0,0);
    if(plist[i].genHepevt()){
      dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).cl(kmf.cl());
      dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).chisq(kmf.chisq());
      dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).ndf(kmf.dgf());
      dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).mass(pi0.mass(),1);
      dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).px(pi0.p().x(),1);
      dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).py(pi0.p().y(),1);
      dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).pz(pi0.p().z(),1);
      dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).epx(pi0.momentum().dp()[0][0],1);
      dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).epy(pi0.momentum().dp()[1][1],1);
      dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).epz(pi0.momentum().dp()[2][2],1);
    }

    // mass and vertex fit
#if 0
    kmassvertexfitter kf;
    kf.magneticField(BF);
    addTrack2fit(kf, plist[i].child(0));
    addTrack2fit(kf, plist[i].child(1));
    addTrack2fit(kf, pi0);
    kf.invariantMass(m_ptypeD0.mass());
#else
    kmassfitter kf;
    kf.magneticField(BF);
    addTrack2fit(kf, plist[i].child(0));
    addTrack2fit(kf, plist[i].child(1));
    addTrack2fit(kf, pi0);
    kf.invariantMass(m_ptypeD0.mass());
    kf.vertex(kvf.vertex());
    kf.errVertex(kvf.errVertex());
#endif
    err = kf.fit();
    if(err == 0){
      //dout(Debugout::INFO,"d0kpipi0_mv") << "2 " << kvf.vertex() << " " << kmf.vertex() << " " << kf.vertex() << std::endl;
      kmakemother kmm;
      kmm.magneticField(BF);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).mass(plist[i].mass(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).px(plist[i].p().x(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).py(plist[i].p().y(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).pz(plist[i].p().z(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).vx(kvf.vertex().x(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).vy(kvf.vertex().y(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).vz(kvf.vertex().z(),0);
      makeMother(kmm, kf, plist[i], 0);// change Momentum Info.
      dynamic_cast<UserInfo&>(plist[i].userInfo()).cl(kf.cl());
      dynamic_cast<UserInfo&>(plist[i].userInfo()).chisq(kf.chisq());
      dynamic_cast<UserInfo&>(plist[i].userInfo()).ndf(kf.dgf());
      dynamic_cast<UserInfo&>(plist[i].userInfo()).mass(plist[i].mass(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).px(plist[i].p().x(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).py(plist[i].p().y(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).pz(plist[i].p().z(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).vx(kf.vertex().x(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).vy(kf.vertex().y(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).vz(kf.vertex().z(),1);
    }else{
      HepPoint3D vtx(999.,999.,999.);
      HepSymMatrix errVtx(3,0);
      plist[i].momentum().decayVertex(vtx,errVtx);
    }
  }
}

void
d0kpipi0_mv_module::writeHisto(vector<Particle> &plist,
			       BelleTuple *hist)
{
  for(int i=0;i<plist.size();++i){
    // lund(charge) info.
    hist->column("lund", plist[i].lund());

    // event info.
    hist->column("expN",  m_expNo); // Exp #
    hist->column("runN",  m_runNo); // Run #
    hist->column("evtN",  m_evtNo); // Event #
    hist->column("farmN", m_farmNo);// Farm #

    // event shape
    hist->column("R2", m_r2);

    // gamma
    hist->column("ga0px0",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).px(0));
    hist->column("ga0px1",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).px(1));
    hist->column("ga0py0",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).py(0));
    hist->column("ga0py1",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).py(1));
    hist->column("ga0pz0",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).pz(0));
    hist->column("ga0pz1",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).pz(1));

    hist->column("ga0epx02",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).epx(0));
    hist->column("ga0epy02",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).epy(0));
    hist->column("ga0epz02",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).epz(0));
    hist->column("ga0epx2",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).epx(1));
    hist->column("ga0epy2",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).epy(1));
    hist->column("ga0epz2",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).epz(1));

    hist->column("ga1px0",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).px(0));
    hist->column("ga1px1",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).px(1));
    hist->column("ga1py0",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).py(0));
    hist->column("ga1py1",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).py(1));
    hist->column("ga1pz0",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).pz(0));
    hist->column("ga1pz1",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).pz(1));

    hist->column("ga1epx02",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).epx(0));
    hist->column("ga1epy02",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).epy(0));
    hist->column("ga1epz02",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).epz(0));
    hist->column("ga1epx2",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).epx(1));
    hist->column("ga1epy2",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).epy(1));
    hist->column("ga1epz2",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).epz(1));

    // pi0
    hist->column("pi0cl",   dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).cl());
    hist->column("pi0chisq",dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).chisq());
    hist->column("pi0ndf",  dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).ndf());

    hist->column("pi0m0",dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).mass(0));
    hist->column("pi0m1",dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).mass(1));
    hist->column("pi0px0",dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).px(0));
    hist->column("pi0px1",dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).px(1));
    hist->column("pi0py0",dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).py(0));
    hist->column("pi0py1",dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).py(1));
    hist->column("pi0pz0",dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).pz(0));
    hist->column("pi0pz1",dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).pz(1));

    hist->column("pi0epx02",dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).epx(0));
    hist->column("pi0epy02",dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).epy(0));
    hist->column("pi0epz02",dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).epz(0));
    hist->column("pi0epx2",dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).epx(1));
    hist->column("pi0epy2",dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).epy(1));
    hist->column("pi0epz2",dynamic_cast<UserInfo&>(plist[i].child(2).userInfo()).epz(1));

    // D0
    hist->column("cl",   dynamic_cast<UserInfo&>(plist[i].userInfo()).cl());
    hist->column("chisq",dynamic_cast<UserInfo&>(plist[i].userInfo()).chisq());
    hist->column("ndf",  dynamic_cast<UserInfo&>(plist[i].userInfo()).ndf());

    hist->column("m0",dynamic_cast<UserInfo&>(plist[i].userInfo()).mass(0));
    hist->column("m1",dynamic_cast<UserInfo&>(plist[i].userInfo()).mass(1));
    hist->column("px0",dynamic_cast<UserInfo&>(plist[i].userInfo()).px(0));
    hist->column("px1",dynamic_cast<UserInfo&>(plist[i].userInfo()).px(1));
    hist->column("py0",dynamic_cast<UserInfo&>(plist[i].userInfo()).py(0));
    hist->column("py1",dynamic_cast<UserInfo&>(plist[i].userInfo()).py(1));
    hist->column("pz0",dynamic_cast<UserInfo&>(plist[i].userInfo()).pz(0));
    hist->column("pz1",dynamic_cast<UserInfo&>(plist[i].userInfo()).pz(1));

    hist->column("epx2",plist[i].momentum().dp()[0][0]);
    hist->column("epy2",plist[i].momentum().dp()[1][1]);
    hist->column("epz2",plist[i].momentum().dp()[2][2]);

    hist->column("vx",plist[i].momentum().decayVertex().x());
    hist->column("vy",plist[i].momentum().decayVertex().y());
    hist->column("vz",plist[i].momentum().decayVertex().z());

    hist->column("evx2",plist[i].momentum().dDecayVertex()[0][0]);
    hist->column("evy2",plist[i].momentum().dDecayVertex()[1][1]);
    hist->column("evz2",plist[i].momentum().dDecayVertex()[2][2]);

    // 0: Kpi, 1: Kpipi0
    hist->column("vx0",dynamic_cast<UserInfo&>(plist[i].userInfo()).vx(0));
    hist->column("vx1",dynamic_cast<UserInfo&>(plist[i].userInfo()).vx(1));
    hist->column("vy0",dynamic_cast<UserInfo&>(plist[i].userInfo()).vy(0));
    hist->column("vy1",dynamic_cast<UserInfo&>(plist[i].userInfo()).vy(1));
    hist->column("vz0",dynamic_cast<UserInfo&>(plist[i].userInfo()).vz(0));
    hist->column("vz1",dynamic_cast<UserInfo&>(plist[i].userInfo()).vz(1));

    // MC info.
    if(m_mc){
      hist->column("mc",1);
      if(plist[i].child(2).child(0).genHepevt()){
	hist->column("gga0",  1);
	hist->column("gga0px",  plist[i].child(2).child(0).genHepevt().PX());
	hist->column("gga0py",  plist[i].child(2).child(0).genHepevt().PY());
	hist->column("gga0pz",  plist[i].child(2).child(0).genHepevt().PZ());
	hist->column("gga0m",   plist[i].child(2).child(0).genHepevt().M());
      }else{
	hist->column("gga0",  0);
      }
      if(plist[i].child(2).child(1).genHepevt()){
	hist->column("gga1",  1);
	hist->column("gga1px",  plist[i].child(2).child(1).genHepevt().PX());
	hist->column("gga1py",  plist[i].child(2).child(1).genHepevt().PY());
	hist->column("gga1pz",  plist[i].child(2).child(1).genHepevt().PZ());
	hist->column("gga1m",   plist[i].child(2).child(1).genHepevt().M());
      }else{
	hist->column("gga1",  0);
      }
      if(plist[i].child(2).genHepevt()){
	hist->column("gpi0",  1);
	hist->column("gpi0px",  plist[i].child(2).genHepevt().PX());
	hist->column("gpi0py",  plist[i].child(2).genHepevt().PY());
	hist->column("gpi0pz",  plist[i].child(2).genHepevt().PZ());
	hist->column("gpi0m",   plist[i].child(2).genHepevt().M());
	hist->column("gpi0dvx", plist[i].child(2).child(0).genHepevt().VX()*0.1);
	hist->column("gpi0dvy", plist[i].child(2).child(0).genHepevt().VY()*0.1);
	hist->column("gpi0dvz", plist[i].child(2).child(0).genHepevt().VZ()*0.1);
	hist->column("gpi0pvx", plist[i].child(2).genHepevt().VX()*0.1);
	hist->column("gpi0pvy", plist[i].child(2).genHepevt().VY()*0.1);
	hist->column("gpi0pvz", plist[i].child(2).genHepevt().VZ()*0.1);
      }else{
	hist->column("gpi0",  0);
      }
      if(plist[i].genHepevt()){
	hist->column("gd0",  1);
	hist->column("gpx",  plist[i].genHepevt().PX());
	hist->column("gpy",  plist[i].genHepevt().PY());
	hist->column("gpz",  plist[i].genHepevt().PZ());
	hist->column("gm",   plist[i].genHepevt().M());
	hist->column("gdvx", plist[i].child(0).genHepevt().VX()*0.1);
	hist->column("gdvy", plist[i].child(0).genHepevt().VY()*0.1);
	hist->column("gdvz", plist[i].child(0).genHepevt().VZ()*0.1);
	hist->column("gpvx", plist[i].genHepevt().VX()*0.1);
	hist->column("gpvy", plist[i].genHepevt().VY()*0.1);
	hist->column("gpvz", plist[i].genHepevt().VZ()*0.1);
      }else{
	hist->column("gd0",0);
      }
    }else{
      // Real Data
      hist->column("mc",-1);
    }
    hist->dumpData();
  }
}


void
d0kpipi0_mv_module::fitPI0Event(vector<Particle> &plist,
				BelleTuple *hist)
{
  for(int i=0;i<plist.size();++i){
    // mass fit of PI0
    Particle g0 = plist[i].child(0);
    Particle g1 = plist[i].child(1);
    HepSymMatrix tmpErr(3,0);
    HepPoint3D org(0.,0.,0.);
    setGammaError(g0, org, tmpErr);
    setGammaError(g1, org, tmpErr);

    kmassfitter kmf;
    kmf.magneticField(BF);
    addTrack2fit(kmf, g0);
    addTrack2fit(kmf, g1);
    kmf.invariantMass(m_ptypePI0.mass());
    kmf.vertex(org);
    kmf.atDecayPoint();
    int err = kmf.fit();

    if(err != 0)continue;
    Particle pi0 = plist[i];
    kmakemother kmm2;
    kmm2.magneticField(BF);
    if(plist[i].genHepevt()){
      hist->column("pi0m0",plist[i].mass());
      hist->column("pi0px0",plist[i].p().x());
      hist->column("pi0py0",plist[i].p().y());
      hist->column("pi0pz0",plist[i].p().z());

      hist->column("pi0epx02",plist[i].momentum().dp()[0][0]);
      hist->column("pi0epy02",plist[i].momentum().dp()[1][1]);
      hist->column("pi0epz02",plist[i].momentum().dp()[2][2]);
    }
    makeMother(kmm2,kmf,pi0,0);
    hist->column("pi0cl",   kmf.cl());
    hist->column("pi0chisq",kmf.chisq());
    hist->column("pi0ndf",  kmf.dgf());
    if(plist[i].genHepevt()){
      hist->column("pi0m1",pi0.mass());
      hist->column("pi0px1",pi0.p().x());
      hist->column("pi0py1",pi0.p().y());
      hist->column("pi0pz1",pi0.p().z());

      hist->column("pi0epx2",pi0.momentum().dp()[0][0]);
      hist->column("pi0epy2",pi0.momentum().dp()[1][1]);
      hist->column("pi0epz2",pi0.momentum().dp()[2][2]);
    }

    if(plist[i].genHepevt()){
      // lund(charge) info.
      hist->column("lund", plist[i].lund());

      // event info.
      hist->column("expN",  m_expNo); // Exp #
      hist->column("runN",  m_runNo); // Run #
      hist->column("evtN",  m_evtNo); // Event #
      hist->column("farmN", m_farmNo);// Farm #
      
      // event shape
      hist->column("R2", m_r2);

      /*
      // gamma
      hist->column("ga0px0",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).px(0));
      hist->column("ga0px1",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).px(1));
      hist->column("ga0py0",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).py(0));
      hist->column("ga0py1",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).py(1));
      hist->column("ga0pz0",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).pz(0));
      hist->column("ga0pz1",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).pz(1));
      
      hist->column("ga0epx02",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).epx(0));
      hist->column("ga0epy02",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).epy(0));
      hist->column("ga0epz02",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).epz(0));
      hist->column("ga0epx2",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).epx(1));
      hist->column("ga0epy2",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).epy(1));
      hist->column("ga0epz2",dynamic_cast<UserInfo&>(plist[i].child(2).child(0).userInfo()).epz(1));
      
      hist->column("ga1px0",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).px(0));
      hist->column("ga1px1",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).px(1));
      hist->column("ga1py0",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).py(0));
      hist->column("ga1py1",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).py(1));
      hist->column("ga1pz0",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).pz(0));
      hist->column("ga1pz1",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).pz(1));
      
      hist->column("ga1epx02",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).epx(0));
      hist->column("ga1epy02",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).epy(0));
      hist->column("ga1epz02",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).epz(0));
      hist->column("ga1epx2",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).epx(1));
      hist->column("ga1epy2",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).epy(1));
      hist->column("ga1epz2",dynamic_cast<UserInfo&>(plist[i].child(2).child(1).userInfo()).epz(1));
      */

      if(m_mc){
	hist->column("mc",1);
	if(plist[i].child(0).genHepevt()){
	  hist->column("gga0",  1);
	  hist->column("gga0px",  plist[i].child(0).genHepevt().PX());
	  hist->column("gga0py",  plist[i].child(0).genHepevt().PY());
	  hist->column("gga0pz",  plist[i].child(0).genHepevt().PZ());
	  hist->column("gga0m",   plist[i].child(0).genHepevt().M());
	}else{
	  hist->column("gga0",  0);
	}
	if(plist[i].child(1).genHepevt()){
	  hist->column("gga1",  1);
	  hist->column("gga1px",  plist[i].child(1).genHepevt().PX());
	  hist->column("gga1py",  plist[i].child(1).genHepevt().PY());
	  hist->column("gga1pz",  plist[i].child(1).genHepevt().PZ());
	  hist->column("gga1m",   plist[i].child(1).genHepevt().M());
	}else{
	  hist->column("gga1",  0);
	}
	if(plist[i].genHepevt()){
	  hist->column("gpi0",  1);
	  hist->column("gpi0px",  plist[i].genHepevt().PX());
	  hist->column("gpi0py",  plist[i].genHepevt().PY());
	  hist->column("gpi0pz",  plist[i].genHepevt().PZ());
	  hist->column("gpi0m",   plist[i].genHepevt().M());
	}else{
	  hist->column("gpi0",  0);
	}
      }else{
	// Real Data
	hist->column("mc",-1);
      }
    }
    hist->dumpData();
  }
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
