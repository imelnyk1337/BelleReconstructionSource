#include "belle.h"
#include "particle/Particle.h"
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"

#include "combination.h"
#include "utility.h"

#define EXKFIT exkfit
//#define EXKFIT std

#include HEPEVT_H
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


class exkfitter : public Module {
public:
  exkfitter(void);
  ~exkfitter(void){};
  void init(int *status);
  void term(void){};
  void disp_stat(const char*){};
  void hist_def(void);
  void event(BelleEvent*, int*);
  void begin_run(BelleEvent*, int *status){};
  void end_run(BelleEvent*, int *status){};
  void other(int*, BelleEvent*, int*){};

public:
  int m_sim;
  int m_debug;

private:
  // PType
  Ptype m_ptypeDM; // D-
  Ptype m_ptypeB0; // D-

private:
  void doKmvFit(vector<Particle> &plist);
  unsigned doKmvFit(Particle &p);
  void doKvFit(vector<Particle> &plist);
  unsigned doKvFit(Particle &p);

private:
  // Histgram and NTuple
  BelleHistogram *m_event_counter;
  BelleHistogram *m_cl_dm_s, *m_cl_dm_n;
  BelleHistogram *m_vtx_dm_x, *m_vtx_dm_y, *m_vtx_dm_z;
  BelleHistogram *m_pullvtx_dm_x, *m_pullvtx_dm_y, *m_pullvtx_dm_z;
  BelleHistogram *m_cl_b0_s, *m_cl_b0_n;
  BelleHistogram *m_vtx_b0_x, *m_vtx_b0_y, *m_vtx_b0_z;
  BelleHistogram *m_pullvtx_b0_x, *m_pullvtx_b0_y, *m_pullvtx_b0_z;
};

extern "C" Module_descr *mdcl_exkfitter() {
  exkfitter *module = new exkfitter;
  Module_descr *dscr   = new Module_descr ( "exkfitter", module );
  // 0 = fsim, 1 = gsim
  dscr->define_param("SIM_TYPE", "SIM_TYPE", &module->m_sim);
  // 0 = no print, 1 or more = print
  dscr->define_param("DEBUG", "DEBUG", &module->m_debug);
  return dscr;
}

exkfitter::exkfitter(void)
  : m_sim(1), m_debug(0), m_ptypeDM("D-"), m_ptypeB0("B0") {
}

void
exkfitter::init(void){
  dout(Debugout::INFO,"exkfitter") << std::endl;
  dout(Debugout::INFO,"exkfitter") << "Example of Kfitter(exkfitter)" << std::endl;
  dout(Debugout::INFO,"exkfitter") << "  B0 -> D-pi-pi+pi+          " << std::endl;
  dout(Debugout::INFO,"exkfitter") << "        D- -> K+pi-pi-       " << std::endl;
  dout(Debugout::INFO,"exkfitter") << "  Simulator(0:fsim, 1:gsim) = " << m_sim << std::endl;
  dout(Debugout::INFO,"exkfitter") << "  Debug(0:no print, 1 or more:print) = " << m_debug << std::endl;
  dout(Debugout::INFO,"exkfitter") << std::endl;
}

void 
exkfitter::hist_def(void) {
  extern BelleTupleManager *BASF_Histogram;

  // 8 characters
  m_event_counter = BASF_Histogram->histogram("Event Counter", 11, -5.5, 5.5);
  m_cl_dm_s  = BASF_Histogram->histogram("CL(D- : signal)", 100, 0.0, 1.0);
  m_cl_dm_n  = BASF_Histogram->histogram("CL(D- : noise)", 100, 0.0, 1.0);
  m_vtx_dm_x = BASF_Histogram->histogram("x of VTX(D-)", 100, -0.1, 0.1);
  m_vtx_dm_y = BASF_Histogram->histogram("y of VTX(D-)", 100, -0.1, 0.1);
  m_vtx_dm_z = BASF_Histogram->histogram("z of VTX(D-)", 100, -0.1, 0.1);
  m_pullvtx_dm_x = BASF_Histogram->histogram("pull x of VTX(D-)", 100, -5, 5);
  m_pullvtx_dm_y = BASF_Histogram->histogram("pull y of VTX(D-)", 100, -5, 5);
  m_pullvtx_dm_z = BASF_Histogram->histogram("pull z of VTX(D-)", 100, -5, 5);
  m_cl_b0_s  = BASF_Histogram->histogram("CL(B0 : signal)", 100, 0.0, 1.0);
  m_cl_b0_n  = BASF_Histogram->histogram("CL(B0 : noise)", 100, 0.0, 1.0);
  m_vtx_b0_x = BASF_Histogram->histogram("x of VTX(B0)", 100, -0.1, 0.1);
  m_vtx_b0_y = BASF_Histogram->histogram("y of VTX(B0)", 100, -0.1, 0.1);
  m_vtx_b0_z = BASF_Histogram->histogram("z of VTX(B0)", 100, -0.1, 0.1);
  m_pullvtx_b0_x = BASF_Histogram->histogram("pull x of VTX(B0)", 100, -5, 5);
  m_pullvtx_b0_y = BASF_Histogram->histogram("pull y of VTX(B0)", 100, -5, 5);
  m_pullvtx_b0_z = BASF_Histogram->histogram("pull z of VTX(B0)", 100, -5, 5);
}

void 
exkfitter::event(BelleEvent *evptr, int *status) {
  // Event Counter
  m_event_counter->accumulate(0.,1.);

  vector<Particle> kaonP, kaonM;
  vector<Particle> pionP, pionM;

  // Make Charged Particles (with No Cut)
  EXKFIT::makeKPi(kaonP, kaonM, pionP, pionM);

  // Set Generator Level Information to Kaon and Pion
  EXKFIT::setGenHepInfoF(kaonP);
  //EXKFIT::setGenHepInfoF(kaonM);
  EXKFIT::setGenHepInfoF(pionP);
  EXKFIT::setGenHepInfoF(pionM);
  if(m_sim){
    EXKFIT::setUniqueGenHepInfoFBySvdAndDeltaP(kaonP);
    //EXKFIT::setUniqueGenHepInfoFBySvdAndDeltaP(kaonM);
    EXKFIT::setUniqueGenHepInfoFBySvdAndDeltaP(pionP);
    EXKFIT::setUniqueGenHepInfoFBySvdAndDeltaP(pionM);
  }

  // Make D-
  vector<Particle> DM, FITDM;
  EXKFIT::combination(DM, m_ptypeDM, kaonP, pionM, pionM);

  // Set Generator Level Information to D-
  EXKFIT::setGenHepInfoR(DM);

  // Mass and Vertex Fit (D-)
  EXKFIT::deepCopy(DM,FITDM);
  doKmvFit(FITDM);

  vector<Particle> B0, FITB0;

  // Make B0
  EXKFIT::combination(B0, m_ptypeB0, FITDM, pionM, pionP, pionP);

  // Set Generator Level Information to B0
  EXKFIT::setGenHepInfoR(B0);

  // Vertex Fit (B0)
  EXKFIT::deepCopy(B0,FITB0);
  doKvFit(FITB0);

  EXKFIT::deleteDeepCopiedObjects(FITDM);
  EXKFIT::deleteDeepCopiedObjects(FITB0);
}

void
exkfitter::doKmvFit(vector<Particle> &plist) {
  for(unsigned i=0;i<plist.size();++i)
    doKmvFit(plist[i]);
}

unsigned
exkfitter::doKmvFit(Particle &p) {
  kmassvertexfitter kmv;
  kmv.invariantMass(p.pType().mass());
  for(unsigned j=0;j<p.relation().nChildren();++j){
    EXKFIT::addTrack2fit(kmv,p.relation().child(j));
  }
  unsigned err = kmv.fit();
  if(err)return 0;
  if(p.genHepevt()){
    m_cl_dm_s->accumulate(kmv.cl(),1.);
    m_vtx_dm_x->accumulate(kmv.vertex().x()-p.relation().child(0).genHepevt().VX()*0.1,1.);
    m_vtx_dm_y->accumulate(kmv.vertex().y()-p.relation().child(0).genHepevt().VY()*0.1,1.);
    m_vtx_dm_z->accumulate(kmv.vertex().z()-p.relation().child(0).genHepevt().VZ()*0.1,1.);
    m_pullvtx_dm_x->accumulate((kmv.vertex().x()-p.relation().child(0).genHepevt().VX()*0.1)
			       /sqrt(kmv.errVertex()[0][0]),1.);
    m_pullvtx_dm_y->accumulate((kmv.vertex().y()-p.relation().child(0).genHepevt().VY()*0.1)
			       /sqrt(kmv.errVertex()[1][1]),1.);
    m_pullvtx_dm_z->accumulate((kmv.vertex().z()-p.relation().child(0).genHepevt().VZ()*0.1)
			       /sqrt(kmv.errVertex()[2][2]),1.);
  }else{
    m_cl_dm_n->accumulate(kmv.cl(),1.);
  }
  if(m_debug){
    if(p.genHepevt())dout(Debugout::INFO,"exkfitter") << "kmv: signal" << std::endl;
    else dout(Debugout::INFO,"exkfitter") << "kmv: noise" << std::endl;
    dout(Debugout::INFO,"exkfitter") << "chisq = " << kmv.chisq() << std::endl;
    double total(0.);
    for(unsigned j=0;j<p.relation().nChildren();++j){
      dout(Debugout::INFO,"exkfitter") << "  " << j << " : " << kmv.chisq(j) << std::endl;
      total += kmv.chisq(j);
    }
    dout(Debugout::INFO,"exkfitter") << "sum of chisq = " << total << std::endl;
  }
  return EXKFIT::makeMother(kmv,p);
}

void
exkfitter::doKvFit(vector<Particle> &plist) {
  for(unsigned i=0;i<plist.size();++i)
    doKvFit(plist[i]);
}

unsigned
exkfitter::doKvFit(Particle &p) {
  kvertexfitter kv;
  for(unsigned j=0;j<p.relation().nChildren();++j){
    EXKFIT::addTrack2fit(kv,p.relation().child(j));
  }
  unsigned err = kv.fit();
  if(err)return 0;
  if(p.genHepevt()){
    m_cl_b0_s->accumulate(kv.cl(),1.);
    m_vtx_b0_x->accumulate(kv.vertex().x()-p.relation().child(0).genHepevt().VX()*0.1,1.);
    m_vtx_b0_y->accumulate(kv.vertex().y()-p.relation().child(0).genHepevt().VY()*0.1,1.);
    m_vtx_b0_z->accumulate(kv.vertex().z()-p.relation().child(0).genHepevt().VZ()*0.1,1.);
    m_pullvtx_b0_x->accumulate((kv.vertex().x()-p.relation().child(0).genHepevt().VX()*0.1)
			       /sqrt(kv.errVertex()[0][0]),1.);
    m_pullvtx_b0_y->accumulate((kv.vertex().y()-p.relation().child(0).genHepevt().VY()*0.1)
			       /sqrt(kv.errVertex()[1][1]),1.);
    m_pullvtx_b0_z->accumulate((kv.vertex().z()-p.relation().child(0).genHepevt().VZ()*0.1)
			       /sqrt(kv.errVertex()[2][2]),1.);
  }else{
    m_cl_b0_n->accumulate(kv.cl(),1.);
  }
  if(m_debug){
    if(p.genHepevt())dout(Debugout::INFO,"exkfitter") << "kv: signal" << std::endl;
    else dout(Debugout::INFO,"exkfitter") << "kv: noise" << std::endl;
    dout(Debugout::INFO,"exkfitter") << "chisq = " << kv.chisq() << std::endl;
    double total(0.);
    for(unsigned j=0;j<p.relation().nChildren();++j){
      dout(Debugout::INFO,"exkfitter") << "  " << j << " : " << kv.chisq(j) << std::endl;
      total += kv.chisq(j);
    }
    dout(Debugout::INFO,"exkfitter") << "sum of chisq = " << total << std::endl;
  }
  return EXKFIT::makeMother(kv,p);
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
