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

class d0kpipipi_v_mv_module : public Module
{
public:
  d0kpipipi_v_mv_module(void);
  ~d0kpipipi_v_mv_module(void){};
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
  void massD0Event(vector<Particle>&);
  void writeHisto(vector<Particle>&, BelleTuple*);

private:
  // PType
  Ptype m_ptypeD0;

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

extern "C" Module_descr *mdcl_d0kpipipi_v_mv()
{
  d0kpipipi_v_mv_module *module = new d0kpipipi_v_mv_module;
  Module_descr *dscr = new Module_descr ( "d0kpipipi_v_mv", module );
  return dscr;
}

d0kpipipi_v_mv_module::d0kpipipi_v_mv_module(void)
  : m_ptypeD0("D0"),
    m_mc(1)
{
}

void 
d0kpipipi_v_mv_module::hist_def(void)
{
  extern BelleTupleManager *BASF_Histogram;

  m_event_counter = BASF_Histogram->histogram("Event Counter", 11, -5.5, 5.5);

  m_hist = BASF_Histogram->ntuple("Event Info", 
				  // 7
				  "mc lund expN runN evtN farmN R2 "
				  // 20
				  "cl chisq ndf "
				  "m0 px0 py0 pz0 "
				  "m1 px1 py1 pz1 "
				  "epx2 epy2 epz2 "
				  "vx vy vz "				  
				  "evx2 evy2 evz2 "
				  // 21
				  "dvx2 dvy2 dvz2 devx2 devy2 devz2 dcl2 "
				  "dvx3 dvy3 dvz3 devx3 devy3 devz3 dcl3 "
				  "dvx4 dvy4 dvz4 devx4 devy4 devz4 dcl4 "
				  // 11
				  "gd0 "
				  "gpx gpy gpz gm "
				  "gdvx gdvy gdvz "
				  "gpvx gpvy gpvz "
				  );
}

void
d0kpipipi_v_mv_module::begin_run(BelleEvent *dummy,
			  int *dummyStatus)
{
  //eid::init_data();
  return;
}

void
d0kpipipi_v_mv_module::init(int *dummyStatus)
{
  //dummy
  Ptype dummy("E-");
}

void 
d0kpipipi_v_mv_module::endEvent(void)
{
  eraseVector(m_kaonP);
  eraseVector(m_kaonM);
  eraseVector(m_pionP);
  eraseVector(m_pionM);
}

void
d0kpipipi_v_mv_module::event(BelleEvent *evptr, int *status)
{
  *status = -1;
  m_event_counter->accumulate(0.,1.);

  // MC or Data
  if(m_mc){
    Belle_event_Manager &evtMgr = Belle_event_Manager::get_manager();
    if(evtMgr.count() != 0){
      if(evtMgr[0].ExpMC() != 2){
        m_mc = 0; // not MC
        dout(Debugout::INFO,"d0kpipipi_v_mv") << "(module) This analysis is in REAL DATA not MC DATA." << std::endl;
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
  combination(D0,  m_ptypeD0, m_kaonM, m_pionM, m_pionP, m_pionP, 0.08);
  setUserInfo(D0);

  if(m_mc){
    setGenHepInfoR(D0);
  }

  massD0Event(D0);

  writeHisto(D0, m_hist);

  endEvent();
}

void
d0kpipipi_v_mv_module::massD0Event(vector<Particle> &plist)
{
  for(int i=0;i<plist.size();++i){
    // vertex fit (2prog)
    kvertexfitter kf2;
    kf2.magneticField(BF);
    addTrack2fit(kf2, plist[i].child(0));
    addTrack2fit(kf2, plist[i].child(1));
    unsigned test_err2 = kf2.fit();

    // vertex fit (4prog)
    kvertexfitter kf4;
    kf4.magneticField(BF);
    addTrack2fit(kf4, plist[i].child(0));
    addTrack2fit(kf4, plist[i].child(1));
    addTrack2fit(kf4, plist[i].child(2));
    addTrack2fit(kf4, plist[i].child(3));
    unsigned test_err4 = kf4.fit();

    // vertex fit (3prog)
    kvertexfitter kf;
    kf.magneticField(BF);
    addTrack2fit(kf, plist[i].child(0));
    addTrack2fit(kf, plist[i].child(1));
    addTrack2fit(kf, plist[i].child(2));
    unsigned err = kf.fit();
    if(err != 0)continue;

    // mass fit
    kmassfitter kmf;
    kmf.magneticField(BF);
    addTrack2fit(kmf, plist[i].child(0));
    addTrack2fit(kmf, plist[i].child(1));
    addTrack2fit(kmf, plist[i].child(2));
    addTrack2fit(kmf, plist[i].child(3));
    kmf.invariantMass(m_ptypeD0.mass());
    //kmf.atDecayPoint();
    kmf.vertex(kf.vertex());
    kmf.errVertex(kf.errVertex());
    err = kmf.fit();

    if(err == 0){
      kmakemother kmm;
      kmm.magneticField(BF);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).mass(plist[i].mass(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).px(plist[i].p().x(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).py(plist[i].p().y(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).pz(plist[i].p().z(),0);
      makeMother(kmm, kmf, plist[i], 0);// change Momentum Info.
      dynamic_cast<UserInfo&>(plist[i].userInfo()).cl(kmf.cl());
      dynamic_cast<UserInfo&>(plist[i].userInfo()).chisq(kmf.chisq());
      dynamic_cast<UserInfo&>(plist[i].userInfo()).ndf(kmf.dgf());
      dynamic_cast<UserInfo&>(plist[i].userInfo()).mass(plist[i].mass(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).px(plist[i].p().x(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).py(plist[i].p().y(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).pz(plist[i].p().z(),1);

      dynamic_cast<UserInfo&>(plist[i].userInfo()).vx(kf2.vertex().x(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).vy(kf2.vertex().y(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).vz(kf2.vertex().z(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).vx(kf.vertex().x(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).vy(kf.vertex().y(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).vz(kf.vertex().z(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).vx(kf4.vertex().x(),2);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).vy(kf4.vertex().y(),2);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).vz(kf4.vertex().z(),2);

      dynamic_cast<UserInfo&>(plist[i].userInfo()).evx(kf2.errVertex()[0][0],0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).evy(kf2.errVertex()[1][1],0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).evz(kf2.errVertex()[2][2],0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).evx(kf.errVertex()[0][0],1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).evy(kf.errVertex()[1][1],1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).evz(kf.errVertex()[2][2],1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).evx(kf4.errVertex()[0][0],2);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).evy(kf4.errVertex()[1][1],2);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).evz(kf4.errVertex()[2][2],2);

      dynamic_cast<UserInfo&>(plist[i].userInfo()).CL(kf2.cl(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).CL(kf.cl(),1);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).CL(kf4.cl(),2);
    }else{
      HepPoint3D vtx(999.,999.,999.);
      HepSymMatrix errVtx(3,0);
      plist[i].momentum().decayVertex(vtx,errVtx);
    }
  }
}

void
d0kpipipi_v_mv_module::writeHisto(vector<Particle> &plist,
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

    // D0 v+m
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

    // 2 prong
    m_hist->column("dvx2", dynamic_cast<UserInfo&>(plist[i].userInfo()).vx(0));
    m_hist->column("dvy2", dynamic_cast<UserInfo&>(plist[i].userInfo()).vy(0));
    m_hist->column("dvz2", dynamic_cast<UserInfo&>(plist[i].userInfo()).vz(0));
    m_hist->column("devx2",dynamic_cast<UserInfo&>(plist[i].userInfo()).evx(0));
    m_hist->column("devy2",dynamic_cast<UserInfo&>(plist[i].userInfo()).evy(0));
    m_hist->column("devz2",dynamic_cast<UserInfo&>(plist[i].userInfo()).evz(0));
    m_hist->column("dcl2", dynamic_cast<UserInfo&>(plist[i].userInfo()).CL(0));

    // 3 prong
    m_hist->column("dvx3", dynamic_cast<UserInfo&>(plist[i].userInfo()).vx(1));
    m_hist->column("dvy3", dynamic_cast<UserInfo&>(plist[i].userInfo()).vy(1));
    m_hist->column("dvz3", dynamic_cast<UserInfo&>(plist[i].userInfo()).vz(1));
    m_hist->column("devx3",dynamic_cast<UserInfo&>(plist[i].userInfo()).evx(1));
    m_hist->column("devy3",dynamic_cast<UserInfo&>(plist[i].userInfo()).evy(1));
    m_hist->column("devz3",dynamic_cast<UserInfo&>(plist[i].userInfo()).evz(1));
    m_hist->column("dcl3", dynamic_cast<UserInfo&>(plist[i].userInfo()).CL(1));

    // 4 prong
    m_hist->column("dvx4", dynamic_cast<UserInfo&>(plist[i].userInfo()).vx(2));
    m_hist->column("dvy4", dynamic_cast<UserInfo&>(plist[i].userInfo()).vy(2));
    m_hist->column("dvz4", dynamic_cast<UserInfo&>(plist[i].userInfo()).vz(2));
    m_hist->column("devx4",dynamic_cast<UserInfo&>(plist[i].userInfo()).evx(2));
    m_hist->column("devy4",dynamic_cast<UserInfo&>(plist[i].userInfo()).evy(2));
    m_hist->column("devz4",dynamic_cast<UserInfo&>(plist[i].userInfo()).evz(2));
    m_hist->column("dcl4", dynamic_cast<UserInfo&>(plist[i].userInfo()).CL(2));

    // MC info.
    if(m_mc){
      hist->column("mc",1);
      if(plist[i].genHepevt()){
	m_hist->column("gd0",  1);
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
	m_hist->column("gd0",0);
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
