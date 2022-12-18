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

class dprhopi_module : public Module
{
public:
  dprhopi_module(void);
  ~dprhopi_module(void){};
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
  void fitRHO0Event(vector<Particle>&);
  void fitDPEvent(vector<Particle>&);
  void writeHisto(vector<Particle>&, BelleTuple*);

private:
  // PType
  Ptype m_ptypeRHO0;
  Ptype m_ptypeDP;

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

extern "C" Module_descr *mdcl_dprhopi()
{
  dprhopi_module *module = new dprhopi_module;
  Module_descr *dscr = new Module_descr ( "dprhopi", module );
  return dscr;
}

dprhopi_module::dprhopi_module(void)
  : m_ptypeRHO0("RHO0"),
    m_ptypeDP("D+"),
    m_mc(1)
{
}

void 
dprhopi_module::hist_def(void)
{
  extern BelleTupleManager *BASF_Histogram;

  m_event_counter = BASF_Histogram->histogram("Event Counter", 11, -5.5, 5.5);

  m_hist = BASF_Histogram->ntuple("Event Info", 
				  // 7
				  "mc lund expN runN evtN farmN R2 "
				  // 20
				  "rcl rchisq rndf "
				  "rm0 rpx0 rpy0 rpz0 "
				  "rm1 rpx1 rpy1 rpz1 "
				  "erpx2 erpy2 erpz2 "
				  "rvx rvy rvz "
				  "ervx2 ervy2 ervz2 "
				  // 20
				  "cl chisq ndf "
				  "m0 px0 py0 pz0 "
				  "m1 px1 py1 pz1 "
				  "epx2 epy2 epz2 "
				  "vx vy vz "
				  "evx2 evy2 evz2 "
				  // 11
				  "gr "
				  "grpx grpy grpz grm "
				  "grdvx grdvy grdvz "
				  "grpvx grpvy grpvz "
				  // 11
				  "gdstar "
				  "gpx gpy gpz gm "
				  "gdvx gdvy gdvz "
				  "gpvx gpvy gpvz "
				  );
}

void
dprhopi_module::begin_run(BelleEvent *dummy,
			  int *dummyStatus)
{
  //eid::init_data();
  return;
}

void
dprhopi_module::init(int *dummyStatus)
{
  //dummy
  Ptype dummy("E-");
}

void 
dprhopi_module::endEvent(void)
{
  eraseVector(m_kaonP);
  eraseVector(m_kaonM);
  eraseVector(m_pionP);
  eraseVector(m_pionM);
}

void
dprhopi_module::event(BelleEvent *evptr, int *status)
{
  *status = -1;
  m_event_counter->accumulate(0.,1.);

  // MC or Data
  if(m_mc){
    Belle_event_Manager &evtMgr = Belle_event_Manager::get_manager();
    if(evtMgr.count() != 0){
      if(evtMgr[0].ExpMC() != 2){
        m_mc = 0; // not MC
        dout(Debugout::INFO,"dprhopi") << "(module) This analysis is in REAL DATA not MC DATA." << std::endl;
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
  vector<Particle> RHO0;
  combination(RHO0, m_ptypeRHO0, m_pionM, m_pionP);
  setUserInfo(RHO0);
  if(m_mc){
    setGenHepInfoR(RHO0);
  }
  fitRHO0Event(RHO0);

  // Make D*+ (with No Cut)
  vector<Particle> DP;
  combination(DP, m_ptypeDP, RHO0, m_pionP);
  setUserInfo(DP);
  if(m_mc){
    setGenHepInfoR(DP);
  }

  fitDPEvent(DP);

  writeHisto(DP, m_hist);

  endEvent();
}

void
dprhopi_module::fitRHO0Event(vector<Particle> &plist)
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
dprhopi_module::fitDPEvent(vector<Particle> &plist)
{
  for(int i=0;i<plist.size();++i){
    // mass fit
    kmassfitter kf;
    kf.magneticField(BF);
    addTrack2fit(kf, plist[i].child(0)); // RHO0
    addTrack2fit(kf, plist[i].child(1)); // PI+
    kf.invariantMass(m_ptypeDP.mass());
    kf.atDecayPoint();
    kf.vertex(plist[i].momentum().decayVertex());
#if 1
    // new
    kf.unfixMass();
    kf.fixMass();
#endif
    unsigned err = kf.fit();

    if(err == 0){
      kmakemother kmm;
      kmm.magneticField(BF);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).mass(plist[i].mass(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).px(plist[i].p().x(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).py(plist[i].p().y(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).pz(plist[i].p().z(),0);
      makeMother(kmm, kf, plist[i], 0); // change Momentum Info.
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
dprhopi_module::writeHisto(vector<Particle> &plist,
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

    // RHO
    m_hist->column("rcl",   dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).cl());
    m_hist->column("rchisq",dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).chisq());
    m_hist->column("rndf",  dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).ndf());

    m_hist->column("rm0",dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).mass(0));
    m_hist->column("rm1",dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).mass(1));
    m_hist->column("rpx0",dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).px(0));
    m_hist->column("rpx1",dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).px(1));
    m_hist->column("rpy0",dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).py(0));
    m_hist->column("rpy1",dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).py(1));
    m_hist->column("rpz0",dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).pz(0));
    m_hist->column("rpz1",dynamic_cast<UserInfo&>(plist[i].child(0).userInfo()).pz(1));

    m_hist->column("erpx2",plist[i].child(0).momentum().dp()[0][0]);
    m_hist->column("erpy2",plist[i].child(0).momentum().dp()[1][1]);
    m_hist->column("erpz2",plist[i].child(0).momentum().dp()[2][2]);

    m_hist->column("rvx",plist[i].child(0).momentum().decayVertex().x());
    m_hist->column("rvy",plist[i].child(0).momentum().decayVertex().y());
    m_hist->column("rvz",plist[i].child(0).momentum().decayVertex().z());

    m_hist->column("ervx2",plist[i].child(0).momentum().dDecayVertex()[0][0]);
    m_hist->column("ervy2",plist[i].child(0).momentum().dDecayVertex()[1][1]);
    m_hist->column("ervz2",plist[i].child(0).momentum().dDecayVertex()[2][2]);

    // D+
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
	m_hist->column("gr",   1);
	m_hist->column("grpx", plist[i].child(0).genHepevt().PX());
	m_hist->column("grpy", plist[i].child(0).genHepevt().PY());
	m_hist->column("grpz", plist[i].child(0).genHepevt().PZ());
	m_hist->column("grm", plist[i].child(0).genHepevt().M());
	m_hist->column("grdvx", plist[i].child(0).child(0).genHepevt().VX()*0.1);
	m_hist->column("grdvy", plist[i].child(0).child(0).genHepevt().VY()*0.1);
	m_hist->column("grdvz", plist[i].child(0).child(0).genHepevt().VZ()*0.1);
	m_hist->column("grpvx", plist[i].child(0).genHepevt().VX()*0.1);
	m_hist->column("grpvy", plist[i].child(0).genHepevt().VY()*0.1);
	m_hist->column("grpvz", plist[i].child(0).genHepevt().VZ()*0.1);
      }else{
	m_hist->column("gr",   0);
      }
      if(plist[i].genHepevt()){
	m_hist->column("gdstar",  1);
	m_hist->column("gpx",  plist[i].genHepevt().PX());
	m_hist->column("gpy",  plist[i].genHepevt().PY());
	m_hist->column("gpz",  plist[i].genHepevt().PZ());
	m_hist->column("gm",   plist[i].genHepevt().M());
	m_hist->column("gdvx", plist[i].child(0).genHepevt().VX()*0.1);
	m_hist->column("gdvy", plist[i].child(0).genHepevt().VY()*0.1);
	m_hist->column("gdvz", plist[i].child(0).genHepevt().VZ()*0.1);
	m_hist->column("gpvx", plist[i].genHepevt().VX()*0.1);
	m_hist->column("gpvy", plist[i].genHepevt().VY()*0.1);
	m_hist->column("gpvz", plist[i].genHepevt().VZ()*0.1);
      }else{
	m_hist->column("gdstar",0);
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
