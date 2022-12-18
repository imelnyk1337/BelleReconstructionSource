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
#include "ip/IpProfile.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


const double BF = 1.5;

class d0kpi_v_module : public Module
{
public:
  d0kpi_v_module(void);
  ~d0kpi_v_module(void){};
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
  void vertexD0Event(vector<Particle>&);
  void writeHisto(vector<Particle>&, BelleTuple*);
  HepPoint3D beamPosition(Particle &p,
			  const HepPoint3D &ipPosition, const HepSymMatrix &ipError,
			  HepSymMatrix &errBp,
			  double &bcl, double &bchisq, unsigned &bndf,
			  HepLorentzVector &newP, HepPoint3D &newX,
			  HepSymMatrix &newErr, HepMatrix &newErrV);

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

extern "C" Module_descr *mdcl_d0kpi_v()
{
  d0kpi_v_module *module = new d0kpi_v_module;
  Module_descr *dscr = new Module_descr ( "d0kpi_v", module );
  IpProfile::define_global(dscr);
  return dscr;
}

d0kpi_v_module::d0kpi_v_module(void)
  : m_ptypeD0("D0"),
    m_mc(1)
{
}

void 
d0kpi_v_module::hist_def(void)
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
				  // 15
				  "pvx pvy pvz "
				  "epvx2 epvy2 epvz2 "
				  "pcl pchisq pndf "
				  "px2 py2 pz2 "
				  "epx22 epy22 epz22 "
				  // 12
				  "gip gd0 "
				  "gpx gpy gpz gm "
				  "gdvx gdvy gdvz "
				  "gpvx gpvy gpvz "
				  );
}

void
d0kpi_v_module::begin_run(BelleEvent *dummy,
			  int *dummyStatus)
{
  //eid::init_data();

  // Get IP profile data from $BELLE_POSTGRES_SERVER 
  IpProfile::begin_run();
  // Dump IP profile data to STDOUT (optional) 
  IpProfile::dump();
  return;
}

void
d0kpi_v_module::init(int *dummyStatus)
{
  //dummy
  Ptype dummy("E-");
}

void 
d0kpi_v_module::endEvent(void)
{
  eraseVector(m_kaonP);
  eraseVector(m_kaonM);
  eraseVector(m_pionP);
  eraseVector(m_pionM);
}

void
d0kpi_v_module::event(BelleEvent *evptr, int *status)
{
  *status = -1;
  m_event_counter->accumulate(0.,1.);

  // MC or Data
  if(m_mc){
    Belle_event_Manager &evtMgr = Belle_event_Manager::get_manager();
    if(evtMgr.count() != 0){
      if(evtMgr[0].ExpMC() != 2){
        m_mc = 0; // not MC
        dout(Debugout::INFO,"d0kpi_v") << "(module) This analysis is in REAL DATA not MC DATA." << std::endl;
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
  vector<Particle> D0, D0B;
  combination(D0,  m_ptypeD0,  m_kaonM, m_pionP);
  setUserInfo(D0);

  if(m_mc){
    setGenHepInfoR(D0);
  }

  vertexD0Event(D0);

  writeHisto(D0, m_hist);

  endEvent();
}

void
d0kpi_v_module::vertexD0Event(vector<Particle> &plist)
{
  for(int i=0;i<plist.size();++i){
    kvertexfitter kvf;
    kvf.magneticField(BF);
    addTrack2fit(kvf, plist[i].child(0));
    addTrack2fit(kvf, plist[i].child(1));
    unsigned err = kvf.fit();
    if(err == 0){
      kmakemother kmm;
      kmm.magneticField(BF);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).mass(plist[i].mass(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).px(plist[i].p().x(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).py(plist[i].p().y(),0);
      dynamic_cast<UserInfo&>(plist[i].userInfo()).pz(plist[i].p().z(),0);
      makeMother(kmm, kvf, plist[i], 0);// change Momentum Info.
      dynamic_cast<UserInfo&>(plist[i].userInfo()).cl(kvf.cl());
      dynamic_cast<UserInfo&>(plist[i].userInfo()).chisq(kvf.chisq());
      dynamic_cast<UserInfo&>(plist[i].userInfo()).ndf(kvf.dgf());
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
d0kpi_v_module::writeHisto(vector<Particle> &plist,
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

    if(IpProfile::usable()){
      HepSymMatrix errBp(3,0);
      double bcl,bchisq;
      unsigned bndf;
      HepLorentzVector newP;
      HepPoint3D newX;
      HepSymMatrix newErr(7,0);
      HepMatrix newErrV(3,7,0);
      HepPoint3D bp = beamPosition(plist[i],
				   IpProfile::position(1), IpProfile::position_err(1),
				   errBp,
				   bcl,bchisq,bndf,
				   newP,newX,
				   newErr,newErrV);

      m_hist->column("pvx",bp.x());
      m_hist->column("pvy",bp.y());
      m_hist->column("pvz",bp.z());
      m_hist->column("epvx2",errBp[0][0]);
      m_hist->column("epvy2",errBp[1][1]);
      m_hist->column("epvz2",errBp[2][2]);
      m_hist->column("pcl",   bcl);
      m_hist->column("pchisq",bchisq);
      m_hist->column("pndf",  bndf);
      m_hist->column("px2", newP.x());
      m_hist->column("py2", newP.y());
      m_hist->column("pz2", newP.z());
      m_hist->column("epx22", newErr[0][0]);
      m_hist->column("epy22", newErr[1][1]);
      m_hist->column("epz22", newErr[2][2]);
    }

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
	if(plist[i].genHepevt().mother()){
	  m_hist->column("gip",  0);
	}else{
	  m_hist->column("gip",  1);
	}
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

HepPoint3D
d0kpi_v_module::beamPosition(Particle &p, 
			     const HepPoint3D &ipPosition, const HepSymMatrix &ipError,
			     HepSymMatrix &errBp,
			     double &bcl, double &bchisq, unsigned &bndf,
			     HepLorentzVector &newP, HepPoint3D &newX, 
			     HepSymMatrix &newErr, HepMatrix &newErrV)			     
{
  HepPoint3D bp(999.,999.,999.);
  bcl = -1.0;
  kvertexfitter kvf;
  kvf.magneticField(BF);
  kvf.addTrack(p.momentum().p(),
               p.momentum().x(),
               p.momentum().dpx(),
               p.pType().charge(),
               p.mass());
  kvf.initialVertex(ipPosition);
  kvf.vertexProfile(ipError);
  int err = kvf.fit();
  if(err == 0){
    bp      = kvf.vertex();
    errBp   = kvf.errVertex();
    bcl     = kvf.cl();
    bchisq  = kvf.chisq();
    bndf    = kvf.dgf();
    newP    = kvf.momentum(0);
    newX    = kvf.position(0);
    newErr  = kvf.error(0);
    newErrV = kvf.errVertexTrack(0);
  }
  return bp;
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
