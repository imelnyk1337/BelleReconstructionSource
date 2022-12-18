//
// $Id: fill_ip.cc 10258 2007-09-21 06:02:59Z hitoshi $
//
// $Log$
// Revision 1.8  2004/09/16 02:59:30  katayama
// New version from Kakuno san
//
// Revision 1.7  2001/12/24 12:03:26  katayama
// gcc 3.0 and headers are cleaned up
//
// Revision 1.6  2000/11/15 12:00:50  tomura
// Check whether BelleTuple* is NULL or not.
//
// Revision 1.5  2000/11/13 09:52:23  tomura
// Add new variables to ntuples.
//
// Revision 1.4  2000/05/16 09:15:18  katayama
// Use std:: etc.
//
// Revision 1.3  2000/05/14 05:22:47  katayama
// Added std::
//
// Revision 1.2  2000/04/25 07:13:09  tomura
// Some parameters are added.
//
// Revision 1.1  2000/04/04 11:58:31  tomura
// New basf module and new Panther table are added.
//

#include "belle.h"
#include "belleCLHEP/Geometry/Point3D.h"
#include "belleCLHEP/Matrix/SymMatrix.h"

#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"

#include "basf/module.h"
#include "basf/module_descr.h"

#include "mdst/mdst.h"
#include "particle/Particle.h"
#include "particle/Ptype.h"
#include "particle/utility.h"
#include "kfitter/kvertexfitter.h"

#include "panther/panther.h"
#include BELLETDF_H
#include MDST_H
#include EVTCLS_H
#include EVTVTX_H
#include IP_H
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif



#define NTUPLE_ID 9501

#define NSVD_RPHI 2
#define NSVD_Z 2

#define CHISQ_THRESHOLD 40.0
#define PSTAR_CUT 5.0


// Class definition
class fill_ip : public Module {
public:
  fill_ip(void);
  ~fill_ip(void){};
  void init(int*);
  void term(void);
  void disp_stat(const char*){};
  void hist_def(void);
  void event(BelleEvent*, int*);
  void begin_run(BelleEvent*, int*);
  void end_run(BelleEvent*, int*){};
  void other(int*, BelleEvent*, int*){};

public:
  int nhit_svd_rphi, nhit_svd_z;
  int m_book_flag;
  int m_use_histo;
  double m_thr_chisq, m_thr_ps;

private:
  int hadron_ip(void);
  int bhabha_ip(void);
  int mumu_ip(void);
  int other_ip(void);
  void save_table(void);

private:
  int m_exp, m_run, m_event, m_starttime;
  int m_flag;
  HepPoint3D m_position;
  HepSymMatrix m_error;
  int m_ndf;
  double m_chisq, m_cl;

  BelleTuple* T_hadron;
  BelleTuple* T_bhabha;
  BelleTuple* T_mumu;

  BelleHistogram* H_hadron[10];
  BelleHistogram* H_bhabha[10];
  BelleHistogram* H_mumu[10];
  BelleHistogram* H_all[10];
};

// BASF Interface Function 
extern "C" Module_descr *mdcl_fill_ip()
{
  fill_ip *module = new fill_ip;
  Module_descr *dscr = new Module_descr("fill_ip", module);

  dscr->define_param("book_ntuple", "Book Ntuple or not", &module->m_book_flag);
  dscr->define_param("use_histo", "Use Histogram (for IP monitor)", &module->m_use_histo);
  dscr->define_param("nhit_SVD_rphi", "No. of hits in SVD-rphi", &module->nhit_svd_rphi);
  dscr->define_param("nhit_SVD_z", "No. of hits in SVD-z", &module->nhit_svd_z);
  dscr->define_param("Chisq_Threshold", "Threshold for chisq cut", &module->m_thr_chisq);
  dscr->define_param("Pstar_Threshold", "Threshold for Pstar cut", &module->m_thr_ps);

  return dscr;
}

// Member Functions

// Constructor
fill_ip::fill_ip(void)
{
  m_book_flag = 1;
  m_use_histo = 1;

  nhit_svd_rphi = NSVD_RPHI;
  nhit_svd_z = NSVD_Z;

  m_thr_chisq = CHISQ_THRESHOLD;
  m_thr_ps = PSTAR_CUT;

  m_starttime = 0;
}

// init function
void fill_ip::init(int* status) 
{
  // dummy construction
  Ptype ptype_dummy("VPHO");
}

// term function
void fill_ip::term (void) 
{
}

// begin_run function
void fill_ip::begin_run(BelleEvent *evptr, int *status)
{
  Belle_runhead_Manager& runhead_mgr = Belle_runhead_Manager::get_manager();

  Belle_runhead_Manager::iterator itrh = runhead_mgr.begin();
  if(itrh == runhead_mgr.end()){
    dout(Debugout::ERR,"fill_ip") << "[fill_ip] Error : There is no Belle_RunHead table." << std::endl;
    return;
  }

  m_starttime = itrh->Time();

  return;
}

// hist_def function
void fill_ip::hist_def(void)
{
  extern BelleTupleManager* BASF_Histogram;
  BelleTupleManager* tm = BASF_Histogram;
  
  T_hadron = tm->ntuple("IP from Hadron ", "expno run event time Vx Vy Vz Exx Eyx Eyy Ezx Ezy Ezz ndf chisq cl", NTUPLE_ID);
  T_bhabha = tm->ntuple("IP from Bhabha", "expno run event time Vx Vy Vz Exx Eyx Eyy Ezx Ezy Ezz ndf chisq cl theta phi E P Px Py Pz");
  T_mumu = tm->ntuple("IP from Mu-pair", "expno run event time Vx Vy Vz Exx Eyx Eyy Ezx Ezy Ezz ndf chisq cl theta phi E P Px Py Pz");
  if (m_use_histo){
    H_all[0] = tm->histogram("Vx", 200, -0.2, 0.2);
    H_all[1] = tm->histogram("Vy", 200, -0.2, 0.2);
    H_all[2] = tm->histogram("Vz", 200, -5., 5.);
    H_all[3] = tm->histogram("Exx", 200, 0., 0.02);
    H_all[4] = tm->histogram("Eyy", 200, 0., 0.02);
    H_all[5] = tm->histogram("Ezz", 200, 0., 0.02);
    H_all[6] = tm->histogram("VxvsVz", 200, -3., 3., 200, -0.2, 0.2);
    H_all[7] = tm->histogram("VyvsVz", 200, -3., 3., 200, -0.2, 0.2);
    H_all[8] = tm->histogram("VyvsVx", 200, -0.2, 0.2, 200, -0.2, 0.2);
  }
}

// event function
void fill_ip::event(BelleEvent* evptr, int* status)
{
  Fill_ip_vertex_Manager& ipv_mgr = Fill_ip_vertex_Manager::get_manager();
  ipv_mgr.remove();

  Belle_event_Manager& event_mgr = Belle_event_Manager::get_manager();

  Belle_event_Manager::iterator itevt = event_mgr.begin();
  if(itevt == event_mgr.end()){
    dout(Debugout::ERR,"fill_ip") << "[fill_ip] Error : There is no Belle_event table." << std::endl;
    return;
  }

  m_exp = itevt->ExpNo();
  m_run = itevt->RunNo();
  m_event = itevt->EvtNo()&0x0FFFFFFF;

  Evtcls_flag_Manager& clsflg_mgr = Evtcls_flag_Manager::get_manager();

  Evtcls_flag_Manager::iterator itflg = clsflg_mgr.begin();
  if(itflg == clsflg_mgr.end()){
    dout(Debugout::ERR,"fill_ip") << "[fill_ip] Error : There is no Evtcls_flag table." << std::endl;
    return;
  }

  m_flag = 0;

  if(itflg->flag(0) >= 30){
    int err = hadron_ip();
    if(err) return;
  }
  else if(itflg->flag(1) == 10 || itflg->flag(1) == 20 || itflg->flag(1) == 30
	  || itflg->flag(1) == 40 || itflg->flag(1) == 110 ||
	  itflg->flag(1) == 120 || itflg->flag(1) == 130 ||
	  itflg->flag(1) == 140){
    int err = bhabha_ip();
    if(err) return;
  }
  else if(itflg->flag(3) == 10 || itflg->flag(3) == 11 || itflg->flag(4) == 20
	  || itflg->flag(3) == 30 || itflg->flag(3) == 40){
    int err = mumu_ip();
    if(err) return;
  }
  else{
    int err = other_ip();
    if(err) return;
  }

  save_table();

  return;
}

int fill_ip::hadron_ip(void)
{
  Evtvtx_primary_vertex_Manager& evtvtx_mgr
    = Evtvtx_primary_vertex_Manager::get_manager();

  Evtvtx_trk_Manager& evttrk_mgr = Evtvtx_trk_Manager::get_manager();

  Evtvtx_primary_vertex_Manager::iterator itvtx = evtvtx_mgr.begin();
  if(itvtx == evtvtx_mgr.end()) return 4;

  if(!(itvtx->quality() == 2 || itvtx->quality() ==3)) return 3;

  kvertexfitter kv;

  HepPoint3D pvtx(itvtx->PV_x(), itvtx->PV_y(), itvtx->PV_z());
  kv.initialVertex(pvtx);

  std::vector<Particle> pion;

  Ptype ptype_pi_plus("PI+");
  Ptype ptype_pi_minus("PI-");

  for(Evtvtx_trk_Manager::iterator i = evttrk_mgr.begin();
	i != evttrk_mgr.end(); ++i){
    //      if(!good_charged(i->charged())) continue;

    if(i->charged().trk().mhyp(2).nhits(3) < nhit_svd_rphi) continue;
    if(i->charged().trk().mhyp(2).nhits(4) < nhit_svd_z) continue;

    Particle tmp(i->charged(),
		 (i->charged().charge()>0.0 ? ptype_pi_plus : ptype_pi_minus));
    pion.push_back(tmp);
  }

  if(pion.size() < 2) return 2;

  addTrack2fit(kv, pion);

  unsigned int err = kv.fit();
  if(err) return 1;

  m_flag = 0x0001;

  m_position = kv.vertex();
  m_error = kv.errVertex();
  m_ndf = kv.dgf();
  m_chisq = kv.chisq();
  m_cl = kv.cl();

  if(m_ndf > 0 && m_chisq/m_ndf < m_thr_chisq){
    m_flag |= 0x0002;

    if(m_book_flag && T_hadron){
      T_hadron->column("expno", m_exp);
      T_hadron->column("run", m_run);
      T_hadron->column("event", m_event);

      Belle_event_Manager& event_mgr = Belle_event_Manager::get_manager();

      Belle_event_Manager::iterator itevt = event_mgr.begin();
      if(itevt != event_mgr.end())
	T_hadron->column("time", itevt->Time() - m_starttime);
      else
	T_hadron->column("time", -999);

      T_hadron->column("Vx", kv.vertex().x());
      T_hadron->column("Vy", kv.vertex().y());
      T_hadron->column("Vz", kv.vertex().z());

      T_hadron->column("Exx", kv.errVertex()[0][0]);
      T_hadron->column("Eyx", kv.errVertex()[1][0]);
      T_hadron->column("Eyy", kv.errVertex()[1][1]);
      T_hadron->column("Ezx", kv.errVertex()[2][0]);
      T_hadron->column("Ezy", kv.errVertex()[2][1]);
      T_hadron->column("Ezz", kv.errVertex()[2][2]);

      T_hadron->column("ndf", kv.dgf());
      T_hadron->column("chisq", kv.chisq());
      T_hadron->column("cl", kv.cl());

      T_hadron->dumpData();

      if (m_use_histo){
	H_all[0]->accumulate(kv.vertex().x());
	H_all[1]->accumulate(kv.vertex().y());
	H_all[2]->accumulate(kv.vertex().z());
	H_all[3]->accumulate(std::sqrt(kv.errVertex()[0][0]));
	H_all[4]->accumulate(std::sqrt(kv.errVertex()[1][1]));
	H_all[5]->accumulate(std::sqrt(kv.errVertex()[2][2]));
	double vx(kv.vertex().x()), vy(kv.vertex().y()), vz(kv.vertex().z());
	H_all[6]->accumulate(vz, vx, (float)1.);
	H_all[7]->accumulate(vz, vy, (float)1.);
	H_all[8]->accumulate(vx, vy, (float)1.);
      }
    }
  }

  return 0;
}

int fill_ip::bhabha_ip(void)
{
  Mdst_charged_Manager& charged_mgr = Mdst_charged_Manager::get_manager();
  if(charged_mgr.count() != 2) return 3;

  std::vector<Particle> electron;

  Ptype ptype_electron("E-");
  Ptype ptype_positron("E+");
  Ptype ptype_upsilon4s("UPS4");

  float Qtot = 0.0;

  double theta, phi;
  bool ps_cut = true;

  for(Mdst_charged_Manager::iterator i = charged_mgr.begin();
      i != charged_mgr.end(); ++i) {
    //if (!good_charged(*i)) continue;

    Qtot += i->charge();

    if(i->trk().mhyp(0).nhits(3) < nhit_svd_rphi) continue;
    if(i->trk().mhyp(0).nhits(4) < nhit_svd_z) continue;

    Particle tmp(*i, (i->charge() > 0.0 ? ptype_positron : ptype_electron));
    if(i->charge() > 0.0){
      theta = tmp.p().theta();
      phi = tmp.p().phi();
    }
    if(pStar(tmp).rho() <= m_thr_ps) ps_cut = false;

    electron.push_back(tmp);
  }

  if(Qtot != 0.0 || electron.size() != 2) return 2;

  Particle Y(electron[0].p()+electron[1].p(), ptype_upsilon4s);
  Y.relation().append(electron[0]);
  Y.relation().append(electron[1]);

  kvertexfitter kv;
  addTrack2fit(kv, electron);

  unsigned int err = kv.fit();
  if(err) return 1;

  makeMother(kv, Y);

  m_flag = 0x0010;

  m_position = kv.vertex();
  m_error = kv.errVertex();
  m_ndf = kv.dgf();
  m_chisq = kv.chisq();
  m_cl = kv.cl();

  if(m_ndf > 0 && m_chisq/m_ndf < m_thr_chisq && ps_cut){
    m_flag |= 0x0020;

    if(m_book_flag && T_bhabha){
      T_bhabha->column("expno", m_exp);
      T_bhabha->column("run", m_run);
      T_bhabha->column("event", m_event);

      Belle_event_Manager& event_mgr = Belle_event_Manager::get_manager();

      Belle_event_Manager::iterator itevt = event_mgr.begin();
      if(itevt != event_mgr.end())
	T_bhabha->column("time", itevt->Time() - m_starttime);
      else
	T_bhabha->column("time", -999);

      T_bhabha->column("Vx", kv.vertex().x());
      T_bhabha->column("Vy", kv.vertex().y());
      T_bhabha->column("Vz", kv.vertex().z());

      T_bhabha->column("Exx", kv.errVertex()[0][0]);
      T_bhabha->column("Eyx", kv.errVertex()[1][0]);
      T_bhabha->column("Eyy", kv.errVertex()[1][1]);
      T_bhabha->column("Ezx", kv.errVertex()[2][0]);
      T_bhabha->column("Ezy", kv.errVertex()[2][1]);
      T_bhabha->column("Ezz", kv.errVertex()[2][2]);

      T_bhabha->column("ndf", kv.dgf());
      T_bhabha->column("chisq", kv.chisq());
      T_bhabha->column("cl", kv.cl());

      T_bhabha->column("theta", theta);
      T_bhabha->column("phi", phi);

      T_bhabha->column("E", Y.e());
      T_bhabha->column("P", Y.ptot());
      T_bhabha->column("Px", Y.px());
      T_bhabha->column("Py", Y.py());
      T_bhabha->column("Pz", Y.pz());

      T_bhabha->dumpData();

      if (m_use_histo){
	H_all[0]->accumulate(kv.vertex().x());
	H_all[1]->accumulate(kv.vertex().y());
	H_all[2]->accumulate(kv.vertex().z());
	H_all[3]->accumulate(std::sqrt(kv.errVertex()[0][0]));
	H_all[4]->accumulate(std::sqrt(kv.errVertex()[1][1]));
	H_all[5]->accumulate(std::sqrt(kv.errVertex()[2][2]));
	double vx(kv.vertex().x()), vy(kv.vertex().y()), vz(kv.vertex().z());
	H_all[6]->accumulate(vz, vx, (float)1.);
	H_all[7]->accumulate(vz, vy, (float)1.);
	H_all[8]->accumulate(vx, vy, (float)1.);
      }
    }
  }

  return 0;
}

int fill_ip::mumu_ip(void)
{
  Mdst_charged_Manager& charged_mgr = Mdst_charged_Manager::get_manager();
  if(charged_mgr.count() != 2) return 3;

  std::vector<Particle> muon;

  Ptype ptype_mu_plus("MU+");
  Ptype ptype_mu_minus("MU-");
  Ptype ptype_upsilon4s("UPS4");

  float Qtot = 0.0;

  double theta, phi;
  bool ps_cut = true;

  for(Mdst_charged_Manager::iterator i = charged_mgr.begin();
      i != charged_mgr.end(); ++i) {
    //if (!good_charged(*i)) continue;

    Qtot += i->charge();

    if(i->trk().mhyp(1).nhits(3) < nhit_svd_rphi) continue;
    if(i->trk().mhyp(1).nhits(4) < nhit_svd_z) continue;

    Particle tmp(*i, (i->charge() > 0.0 ? ptype_mu_plus : ptype_mu_minus));
    if(i->charge() > 0.0){
      theta = tmp.p().theta();
      phi = tmp.p().phi();
    }
    if(pStar(tmp).rho() <= m_thr_ps) ps_cut = false;

    muon.push_back(tmp);
  }

  if(Qtot != 0.0 || muon.size() != 2) return 2;

  Particle Y(muon[0].p()+muon[1].p(), ptype_upsilon4s);
  Y.relation().append(muon[0]);
  Y.relation().append(muon[1]);

  kvertexfitter kv;
  addTrack2fit(kv, muon);

  unsigned int err = kv.fit();
  if(err) return 1;

  makeMother(kv, Y);

  m_flag = 0x0100;

  m_position = kv.vertex();
  m_error = kv.errVertex();
  m_ndf = kv.dgf();
  m_chisq = kv.chisq();
  m_cl = kv.cl();

  if(m_ndf > 0 && m_chisq/m_ndf < m_thr_chisq && ps_cut){
    m_flag |= 0x0200;

    if(m_book_flag && T_mumu){
      T_mumu->column("expno", m_exp);
      T_mumu->column("run", m_run);
      T_mumu->column("event", m_event);

      Belle_event_Manager& event_mgr = Belle_event_Manager::get_manager();

      Belle_event_Manager::iterator itevt = event_mgr.begin();
      if(itevt != event_mgr.end())
	T_mumu->column("time", itevt->Time() - m_starttime);
      else
	T_mumu->column("time", -999);

      T_mumu->column("Vx", kv.vertex().x());
      T_mumu->column("Vy", kv.vertex().y());
      T_mumu->column("Vz", kv.vertex().z());

      T_mumu->column("Exx", kv.errVertex()[0][0]);
      T_mumu->column("Eyx", kv.errVertex()[1][0]);
      T_mumu->column("Eyy", kv.errVertex()[1][1]);
      T_mumu->column("Ezx", kv.errVertex()[2][0]);
      T_mumu->column("Ezy", kv.errVertex()[2][1]);
      T_mumu->column("Ezz", kv.errVertex()[2][2]);

      T_mumu->column("ndf", kv.dgf());
      T_mumu->column("chisq", kv.chisq());
      T_mumu->column("cl", kv.cl());

      T_mumu->column("theta", theta);
      T_mumu->column("phi", phi);

      T_mumu->column("E", Y.e());
      T_mumu->column("P", Y.ptot());
      T_mumu->column("Px", Y.px());
      T_mumu->column("Py", Y.py());
      T_mumu->column("Pz", Y.pz());

      T_mumu->dumpData();

      if (m_use_histo){
	H_all[0]->accumulate(kv.vertex().x());
	H_all[1]->accumulate(kv.vertex().y());
	H_all[2]->accumulate(kv.vertex().z());
	H_all[3]->accumulate(std::sqrt(kv.errVertex()[0][0]));
	H_all[4]->accumulate(std::sqrt(kv.errVertex()[1][1]));
	H_all[5]->accumulate(std::sqrt(kv.errVertex()[2][2]));
	double vx(kv.vertex().x()), vy(kv.vertex().y()), vz(kv.vertex().z());
	H_all[6]->accumulate(vz, vx, (float)1.);
	H_all[7]->accumulate(vz, vy, (float)1.);
	H_all[8]->accumulate(vx, vy, (float)1.);
      }
    }
  }

  return 0;
}

int fill_ip::other_ip(void)
{
  Evtvtx_primary_vertex_Manager& evtvtx_mgr
    = Evtvtx_primary_vertex_Manager::get_manager();

  Evtvtx_trk_Manager& evttrk_mgr = Evtvtx_trk_Manager::get_manager();

  Evtvtx_primary_vertex_Manager::iterator itvtx = evtvtx_mgr.begin();
  if(itvtx == evtvtx_mgr.end()) return 4;

  if(!(itvtx->quality() == 2 || itvtx->quality() ==3)) return 3;

  kvertexfitter kv;

  HepPoint3D pvtx(itvtx->PV_x(), itvtx->PV_y(), itvtx->PV_z());
  kv.initialVertex(pvtx);

  std::vector<Particle> pion;

  Ptype ptype_pi_plus("PI+");
  Ptype ptype_pi_minus("PI-");

  for(Evtvtx_trk_Manager::iterator i = evttrk_mgr.begin();
	i != evttrk_mgr.end(); ++i){
    //      if(!good_charged(i->charged())) continue;

    if(i->charged().trk().mhyp(2).nhits(3) < nhit_svd_rphi) continue;
    if(i->charged().trk().mhyp(2).nhits(4) < nhit_svd_z) continue;

    Particle tmp(i->charged(),
		 (i->charged().charge()>0.0 ? ptype_pi_plus : ptype_pi_minus));
    pion.push_back(tmp);
  }

  if(pion.size() < 2) return 2;

  addTrack2fit(kv, pion);

  unsigned int err = kv.fit();
  if(err) return 1;

  m_flag = 0x1000;

  m_position = kv.vertex();
  m_error = kv.errVertex();
  m_ndf = kv.dgf();
  m_chisq = kv.chisq();
  m_cl = kv.cl();

  return 0;
}

void fill_ip::save_table(void)
{
  Fill_ip_vertex_Manager& ipv_mgr = Fill_ip_vertex_Manager::get_manager();

  Fill_ip_vertex& entity = ipv_mgr.add();

  entity.flag(m_flag);
  for(int i=0; i<3; i++) entity.v(i, m_position(i));
  for(int i=0; i<3; i++)
    for(int j=0; j<=i; j++)
      entity.error(i*(i+1)/2+j , m_error[i][j]);
  entity.ndf(m_ndf);
  entity.chisq(m_chisq);
  entity.cl(m_cl);

  return;
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
