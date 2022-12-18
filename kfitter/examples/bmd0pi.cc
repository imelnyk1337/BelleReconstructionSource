#include "particle/Particle.h"
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"

#include "belle.h"
#include MDST_H
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

class bmd0pi_module : public Module
{
public:
  bmd0pi_module(void);
  ~bmd0pi_module(void){};
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
  void fitBEvent(vector<Particle>&);
  void writeHisto(vector<Particle>&, BelleTuple*);

private:
  // PType
  Ptype m_ptypeD0;
  Ptype m_ptypeBM;

  vector<Particle> m_kaonP, m_kaonM;
  vector<Particle> m_pionP, m_pionM;

private:
  // Histgram and NTuple
  BelleHistogram *m_event_counter;

  BelleTuple     *m_hist;

  int m_expNo, m_runNo, m_evtNo, m_farmNo;
  int m_mc;
  double m_r2;
};

extern "C" Module_descr *mdcl_bmd0pi()
{
  bmd0pi_module *module = new bmd0pi_module;
  Module_descr *dscr = new Module_descr ( "bmd0pi", module );
  return dscr;
}

bmd0pi_module::bmd0pi_module(void)
  : m_ptypeD0("D0"),
    m_ptypeBM("B-"),
    m_mc(1)
{
}

void 
bmd0pi_module::hist_def(void)
{
  extern BelleTupleManager *BASF_Histogram;

  m_event_counter = BASF_Histogram->histogram("Event Counter", 11, -5.5, 5.5);

  m_hist = BASF_Histogram->ntuple("Event Info", 
				  // 7
				  "mc lund expN runN evtN farmN R2 "
				  // 20
				  "dcl dchisq dndf "
				  "dm0 dpx0 dpy0 dpz0 "
				  "dm1 dpx1 dpy1 dpz1 "
				  "edpx2 edpy2 edpz2 "
				  "dvx dvy dvz "
				  "edvx2 edvy2 edvz2 "
				  // 20
				  "cl chisq ndf "
				  "m0 px0 py0 pz0 "
				  "m1 px1 py1 pz1 "
				  "epx2 epy2 epz2 "
				  "vx vy vz "
				  "evx2 evy2 evz2 "
				  // 11
				  "gbm "
				  "gpx gpy gpz gm "
				  "gbdvx gbdvy gbdvz "
				  "gbpvx gbpvy gbpvz "
				  // 11
				  "gd0 "
				  "gdpx gdpy gdpz gdm "
				  "gddvx gddvy gddvz "
				  "gdpvx gdpvy gdpvz "
				  );
}

void
bmd0pi_module::begin_run(BelleEvent *dummy,
			  int *dummyStatus)
{
  //eid::init_data();
  return;
}

void
bmd0pi_module::init(int *dummyStatus)
{
  //dummy
  Ptype dummy("E-");
}

void 
bmd0pi_module::endEvent(void)
{
  eraseVector(m_kaonP);
  eraseVector(m_kaonM);
  eraseVector(m_pionP);
  eraseVector(m_pionM);
}

void
bmd0pi_module::event(BelleEvent *evptr, int *status)
{
  *status = -1;
  m_event_counter->accumulate(0.,1.);

  // MC or Data
  if(m_mc){
    Belle_event_Manager &evtMgr = Belle_event_Manager::get_manager();
    if(evtMgr.count() != 0){
      if(evtMgr[0].ExpMC() != 2){
        m_mc = 0; // not MC
        dout(Debugout::INFO,"bmd0pi") << "(module) This analysis is in REAL DATA not MC DATA." << std::endl;
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

  // Make D0 (with No Cut)
  vector<Particle> D0;
  combination(D0, m_ptypeD0, m_kaonM, m_pionP);
  setUserInfo(D0);
  if(m_mc){
    setGenHepInfoR(D0);
  }
  fitD0Event(D0);

  // Make B- (with No Cut)
  vector<Particle> BM;
  combination(BM, m_ptypeBM, D0, m_pionM);
  setUserInfo(BM);
  if(m_mc){
    setGenHepInfoR(BM);
  }
  fitBEvent(BM);

  writeHisto(BM, m_hist);

  endEvent();
}

void
bmd0pi_module::fitD0Event(vector<Particle> &plist)
{
  for(int i=0;i<plist.size();++i){
    // mass and vertex fit
    kmassvertexfitter kf;
    kf.magneticField(BF);
    addTrack2fit(kf, plist[i].child(0));
    addTrack2fit(kf, plist[i].child(1));
    kf.invariantMass(m_ptypeD0.mass());
    unsigned err = kf.fit();
    if(err == 0){
      kmakemother kmm;
      kmm.magneticField(BF);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).mass(plist[i].mass(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).px(plist[i].p().x(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).py(plist[i].p().y(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).pz(plist[i].p().z(),0);
      makeMother(kmm, kf, plist[i], 0);// change Momentum Info.
      dynamic_cast<UserInfo&>(plist[i].userInfo()).cl(kf.cl());
      dynamic_cast<UserInfo&>(plist[i].userInfo()).chisq(kf.chisq());
      dynamic_cast<UserInfo&>(plist[i].userInfo()).ndf(kf.dgf());
      dynamic_cast<UserInfo&>(plist[i].userInfo()).mass(plist[i].mass(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).px(plist[i].p().x(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).py(plist[i].p().y(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).pz(plist[i].p().z(),1);
    }else{
      HepPoint3D vtx(999.,999.,999.);
      HepSymMatrix errVtx(3,0);
      plist[i].momentum().decayVertex(vtx,errVtx);
    }
  }
}

void
bmd0pi_module::fitBEvent(vector<Particle> &plist)
{
  for(int i=0;i<plist.size();++i){
    // vertex fit
    kvertexfitter kf;
    kf.magneticField(BF);
    addTrack2fit(kf, plist[i].child(0));
    addTrack2fit(kf, plist[i].child(1));
    unsigned err = kf.fit();
    if(err == 0){
      kmakemother kmm;
      kmm.magneticField(BF);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).mass(plist[i].mass(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).px(plist[i].p().x(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).py(plist[i].p().y(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).pz(plist[i].p().z(),0);
      makeMother(kmm, kf, plist[i], 0);// change Momentum Info.
      dynamic_cast<UserInfo&>(plist[i].userInfo()).cl(kf.cl());
      dynamic_cast<UserInfo&>(plist[i].userInfo()).chisq(kf.chisq());
      dynamic_cast<UserInfo&>(plist[i].userInfo()).ndf(kf.dgf());
      dynamic_cast<UserInfo&>(plist[i].userInfo()).mass(plist[i].mass(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).px(plist[i].p().x(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).py(plist[i].p().y(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).pz(plist[i].p().z(),1);
    }else{
      HepPoint3D vtx(999.,999.,999.);
      HepSymMatrix errVtx(3,0);
      plist[i].momentum().decayVertex(vtx,errVtx);
    }
  }
}

void
bmd0pi_module::writeHisto(vector<Particle> &plist,
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

    // D0
    m_hist->column("dcl",   dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).cl());
    m_hist->column("dchisq",dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).chisq());
    m_hist->column("dndf",  dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).ndf());

    m_hist->column("dm0",dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).mass(0));
    m_hist->column("dm1",dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).mass(1));
    m_hist->column("dpx0",dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).px(0));
    m_hist->column("dpx1",dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).px(1));
    m_hist->column("dpy0",dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).py(0));
    m_hist->column("dpy1",dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).py(1));
    m_hist->column("dpz0",dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).pz(0));
    m_hist->column("dpz1",dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).pz(1));

    m_hist->column("edpx2",plist[i].child(0).momentum().dp()[0][0]);
    m_hist->column("edpy2",plist[i].child(0).momentum().dp()[1][1]);
    m_hist->column("edpz2",plist[i].child(0).momentum().dp()[2][2]);

    m_hist->column("dvx",plist[i].child(0).momentum().decayVertex().x());
    m_hist->column("dvy",plist[i].child(0).momentum().decayVertex().y());
    m_hist->column("dvz",plist[i].child(0).momentum().decayVertex().z());

    m_hist->column("edvx2",plist[i].child(0).momentum().dDecayVertex()[0][0]);
    m_hist->column("edvy2",plist[i].child(0).momentum().dDecayVertex()[1][1]);
    m_hist->column("edvz2",plist[i].child(0).momentum().dDecayVertex()[2][2]);

    // B-
    m_hist->column("cl",   dynamic_cast<UserInfo&>(plist[i].userInfo()).cl());
    m_hist->column("chisq",dynamic_cast<UserInfo&>(plist[i].userInfo()).chisq());
    m_hist->column("ndf",  dynamic_cast<UserInfo&>(plist[i].userInfo()).ndf());

    m_hist->column("m0",dynamic_cast<UserInfo&>(plist[i].userInfo()).mass(0));
    m_hist->column("m1",dynamic_cast<UserInfo&>(plist[i].userInfo()).mass(1));
    m_hist->column("px0",dynamic_cast<UserInfo&>(plist[i].userInfo()).px(0));
    m_hist->column("px1",dynamic_cast<UserInfo&>(plist[i].userInfo()).px(1));
    m_hist->column("py0",dynamic_cast<UserInfo&>(plist[i].userInfo()).py(0));
    m_hist->column("py1",dynamic_cast<UserInfo&>(plist[i].userInfo()).py(1));
    m_hist->column("pz0",dynamic_cast<UserInfo&>(plist[i].userInfo()).pz(0));
    m_hist->column("pz1",dynamic_cast<UserInfo&>(plist[i].userInfo()).pz(1));

    m_hist->column("epx2",plist[i].momentum().dp()[0][0]);
    m_hist->column("epy2",plist[i].momentum().dp()[1][1]);
    m_hist->column("epz2",plist[i].momentum().dp()[2][2]);

    m_hist->column("vx",plist[i].momentum().decayVertex().x());
    m_hist->column("vy",plist[i].momentum().decayVertex().y());
    m_hist->column("vz",plist[i].momentum().decayVertex().z());

    m_hist->column("evx2",plist[i].momentum().dDecayVertex()[0][0]);
    m_hist->column("evy2",plist[i].momentum().dDecayVertex()[1][1]);
    m_hist->column("evz2",plist[i].momentum().dDecayVertex()[2][2]);

    // MC info.
    if(m_mc){
      hist->column("mc",1);
      if(plist[i].child(0).genHepevt()){
	m_hist->column("gd0",   1);
	m_hist->column("gdpx",  plist[i].child(0).genHepevt().PX());
	m_hist->column("gdpy",  plist[i].child(0).genHepevt().PY());
	m_hist->column("gdpz",  plist[i].child(0).genHepevt().PZ());
	m_hist->column("gdm",   plist[i].child(0).genHepevt().M());
	m_hist->column("gddvx", plist[i].child(0).child(0).genHepevt().VX()*0.1);
	m_hist->column("gddvy", plist[i].child(0).child(0).genHepevt().VY()*0.1);
	m_hist->column("gddvz", plist[i].child(0).child(0).genHepevt().VZ()*0.1);
	m_hist->column("gdpvx", plist[i].child(0).genHepevt().VX()*0.1);
	m_hist->column("gdpvy", plist[i].child(0).genHepevt().VY()*0.1);
	m_hist->column("gdpvz", plist[i].child(0).genHepevt().VZ()*0.1);
      }else{
	m_hist->column("gd0",   0);
      }
      if(plist[i].genHepevt()){
	m_hist->column("gbm",  1);
	m_hist->column("gpx",  plist[i].genHepevt().PX());
	m_hist->column("gpy",  plist[i].genHepevt().PY());
	m_hist->column("gpz",  plist[i].genHepevt().PZ());
	m_hist->column("gm",   plist[i].genHepevt().M());
	m_hist->column("gbdvx", plist[i].child(0).genHepevt().VX()*0.1);
	m_hist->column("gbdvy", plist[i].child(0).genHepevt().VY()*0.1);
	m_hist->column("gbdvz", plist[i].child(0).genHepevt().VZ()*0.1);
	m_hist->column("gbpvx", plist[i].genHepevt().VX()*0.1);
	m_hist->column("gbpvy", plist[i].genHepevt().VY()*0.1);
	m_hist->column("gbpvz", plist[i].genHepevt().VZ()*0.1);
      }else{
	m_hist->column("gbm",0);
      }
    }else{
      // Real Data
      hist->column("mc",-1);
    }
    hist->dumpData();
  }
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
