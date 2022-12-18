#include "belle.h"
#if 0
#include "particle/Particle.h"
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"

#include "particle/combination.h"
#include "particle/utility.h"
//#include "fbtg/fbtag_main.h"
#include "kfitter/kmakemother.h"

#include HEPEVT_H
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


class jpdln_module : public Module {
public:
  jpdln_module(void);
  ~jpdln_module(void){};
  void init(int *status);
  void term(void){};
  void disp_stat(const char*){};
  void hist_def(void);
  void event(BelleEvent*, int*);
  void begin_run(BelleEvent*, int *){};
  void end_run(BelleEvent*, int *){};
  void other(int*, BelleEvent*, int*){};

public:

private:
  void genhepAnalysis(void);
  unsigned dlnEvent(unsigned*, unsigned*, unsigned*);
  double missMass2(Particle&, Particle&, double&);
  void directionCosTheta(Particle&, Particle&, double*);
  int flavor(Particle&, Particle&, unsigned&);
  unsigned vertexing(kvertexfitter&, Particle&, Particle&, vector<Particle>&, vector<Particle>&);
  void mainAnalysis(Particle&, Particle&, vector<Particle>&, vector<Particle>&, unsigned);
  void dStarLeptonCheck(Particle&, Particle&);

  // For D* and lepton
  unsigned doKVFit(kvertexfitter&, Particle&, Particle&);
  void addTrack2fit(kvertexfitter&, Particle&, Particle&);
  unsigned makeDStar(kvertexfitter&, Particle&);
  void fittedLepton(kvertexfitter&, Particle&);
  void genCheck(Particle&, Particle&);
  
  // For D0
  void doKMFit(vector<Particle> &);
  unsigned doKMFit(Particle &);

private:
  // PType
  Ptype m_ptypeD0;
  Ptype m_ptypeD0B;
  Ptype m_ptypeDStarP;
  Ptype m_ptypeDStarM;

private:
  // Supported Class by BELLE(not me)
  //  Fbtag_main m_fbtag;

private:
  // Histgram and NTuple
  BelleHistogram *m_event_counter;
  BelleHistogram *m_inclusive_d0_event_counter;
  BelleHistogram *m_inclusive_lepton_event_counter;
  BelleHistogram *m_dln_event_counter;
  BelleHistogram *m_dln_b0_event_counter;
  BelleHistogram *m_dln_b0b_event_counter;
  BelleHistogram *m_dln_kpi_event_counter;
  BelleTuple     *m_tuple;
  BelleTuple     *m_tuple2;
};

extern "C" Module_descr *mdcl_exdln() {
  jpdln_module *module = new jpdln_module;
  Module_descr *dscr   = new Module_descr ( "exdln", module );
  return dscr;
}

jpdln_module::jpdln_module(void)
  : m_ptypeD0("D0"),
    m_ptypeD0B("D0B"),
    m_ptypeDStarP("D*+"),
    m_ptypeDStarM("D*-") {
}

void 
jpdln_module::hist_def(void) {
  extern BelleTupleManager *BASF_Histogram;

  // 8 characters
  m_tuple         = BASF_Histogram->ntuple("Dln Event Analysis", "expnum runnum evtnum trig date time ltype lmom lmoms kmom kmoms pmom pmoms spmom spmoms d0mass dsmass massdif dlnvtxx dlnvtxy dlnvtxz dlndvtxx dlndvtxy dlndvtxz dlncl dlnchisq dlndof d0vtxx d0vtxy d0vtxz d0dvtxx d0dvtxy d0dvtxz d0cl d0chisq d0dof hemisp0 hemisp1 missmass misspp lepchg dschg flvl flvk flvo tagvtxx tagvtxy tagvtxz tagdvtxx tagdvtxy tagdvtxz gdschk gdsbf gbvtxx gbvxty gbvtxz gd0vtxx gd0vxty gd0vtxz evtmon snmon ldschk glepchk glepbf");
  m_tuple2        = BASF_Histogram->ntuple("GenHep Dln Event Analysis", 
					   "lmom kmom pmom lmomt kmomt pmomt");

  m_event_counter = BASF_Histogram->histogram("Event Counter", 11, -5.5, 5.5);
  m_inclusive_d0_event_counter = BASF_Histogram->histogram("Inc D0 Event Counter", 11, -5.5, 5.5);
  m_inclusive_lepton_event_counter = BASF_Histogram->histogram("Inc Lepton Event Counter", 11, -5.5, 5.5);
  m_dln_event_counter = BASF_Histogram->histogram("Dln Event Counter", 11, -5.5, 5.5);
  m_dln_b0_event_counter = BASF_Histogram->histogram("Dln(B0) Event Counter", 11, -5.5, 5.5);
  m_dln_b0b_event_counter = BASF_Histogram->histogram("Dln(B0B) Event Counter", 11, -5.5, 5.5);
  m_dln_kpi_event_counter = BASF_Histogram->histogram("Dln(KPi) Event Counter", 11, -5.5, 5.5);
}

void 
jpdln_module::init(int *status) {
  *status = 0;
  // Flavor Tagging Routine
  double kid_thr    = 0.9;
  double proton_thr = 0.9;
  double elid_thr   = 0.998;
  double muid_thr   = 1.0;
  double pstar_el   = 1.1;
  double pstar_mu   = 1.2;

  // Initilize Fbtag_main
//    m_fbtag.f_init( (float)  kid_thr,      // KID threshold
//  		  (float)  proton_thr,   // proton ID threshold
//  		  (float)  elid_thr,     // ELID threshold
//  		  (int)    muid_thr,     // MUID threshold
//  		  (double) pstar_el,     // p* threshold for electorons
//  		  (double) pstar_mu   ); // p* threshold for muons
}

void 
jpdln_module::event(BelleEvent *evptr, int *status) {
  // Event Counter
  m_event_counter->accumulate(0.,1.);

  //static int events = 1;
  //dout(Debugout::INFO,"exdln") << "DEBUG DUMP MODE events# = " << events++ << std::endl; 

  // Check D*ln Event in Generator Level.
  genhepAnalysis();

#define DOANA 0

#if DOANA

  vector<Particle> kaonP, kaonM;
  vector<Particle> pionP, pionM;
  vector<Particle> elecP, elecM;
  vector<Particle> muonP, muonM;

  // Make Charged Particles (with No Cut)
  makeKPi(kaonP, kaonM, pionP, pionM);
  makeLepton(elecP, elecM, muonP, muonM);

  // Event Shape : FoxWolfram
  //double foxWolfram[3];
  //vector<Particle> pList(pionP);
  //for(unsigned i=0;i<pionM.size();++i)pList.push_back(pionM[i]);
  //...gamma etc??
  //calcuFoxWolfram(pList,foxWolfram);

  // Muon ID Routine
  withMuonId(muonP);withMuonId(muonM);
  
  // Electron ID Routine
  withEId(elecP);withEId(elecM);

  // Momentum* Cut to Muon and Electron
  withPSCut(muonP,0.7,3.5);  withPSCut(muonM,0.7,3.5);
  withPSCut(elecP,0.7,3.5);  withPSCut(elecM,0.7,3.5);

  // Requre SVD Association to Muon and Electron
  withSVD(muonP,2);  withSVD(muonM,2);
  withSVD(elecP,2);  withSVD(elecM,2);

  // Inclusive Lepton Counter
  if(muonP.size()+muonM.size()+elecP.size()+elecM.size())
    m_inclusive_lepton_event_counter->accumulate(0.,1.);

  // Set Generator Level Information to Lepton
  setGenHepInfoF(muonP);
  setGenHepInfoF(muonM);
  setGenHepInfoF(elecP);
  setGenHepInfoF(elecM);
  setUniqueGenHepInfoFBySvdAndDeltaP(muonP);
  setUniqueGenHepInfoFBySvdAndDeltaP(muonM);
  setUniqueGenHepInfoFBySvdAndDeltaP(elecP);
  setUniqueGenHepInfoFBySvdAndDeltaP(elecM);

  // Requre SVD Association to Kaon
  withSVD(kaonP,2);
  withSVD(kaonM,2);

  vector<Particle> D0, D0B;

  // Set Generator Level Information to Kaon and Pion
  setGenHepInfoF(kaonP);
  setGenHepInfoF(kaonM);
  setGenHepInfoF(pionP);
  setGenHepInfoF(pionM);
  setUniqueGenHepInfoFBySvdAndDeltaP(kaonP);
  setUniqueGenHepInfoFBySvdAndDeltaP(kaonM);
  setUniqueGenHepInfoFBySvdAndDeltaP(pionP);
  setUniqueGenHepInfoFBySvdAndDeltaP(pionM);

  // Make D0 and D0B (with No Cut)
  combination(D0, m_ptypeD0, kaonM, pionP);
  combination(D0B, m_ptypeD0B, kaonP, pionM);

  // Inclusive D0 Counter
  if(D0.size() != 0 || D0B.size() != 0)
    m_inclusive_d0_event_counter->accumulate(0.,1.);

  // Mass Cut to D0 and D0B (very large)
  withMassCut(D0, 1.77, 1.98);
  withMassCut(D0B, 1.77, 1.98);

  // Requre SVD Association to Pion of D0 and D0B
  withSVD(D0, 2, 1);
  withSVD(D0B, 2, 1);

  // Set Generator Level Information to D0 and D0B
  setGenHepInfoR(D0);
  setGenHepInfoR(D0B);

  vector<Particle> N_D0B;
  vector<Particle> N_D0;

  // Mass Cut to D0(B) (large : about +-30MeV)
  // D0(B)(1.83~1.91MeV), N_D0(B)(others)
  withMassCut(D0B,N_D0B,1.80,1.95);
  withMassCut(D0,N_D0,1.80,1.95);

  vector<Particle> FITD0B;
  vector<Particle> FITN_D0B;
  vector<Particle> FITD0;
  vector<Particle> FITN_D0;

  // Mass and Vertex Fit to Signal Region D0(B)
  deepCopy(D0B,FITD0B);
  deepCopy(D0,FITD0);
  doKMFit(FITD0B);
  doKMFit(FITD0);

  // Mass and Vertex Fit to Noise Region D0(B)
  deepCopy(N_D0B,FITN_D0B);
  deepCopy(N_D0,FITN_D0);
  doKMFit(FITN_D0B);
  doKMFit(FITN_D0);
  
  // Mass Cut to Fitted D0B (large : +-25MeV)
  withMassCut(FITD0B, 
	      m_ptypeD0B.mass()-35.0, 
	      m_ptypeD0B.mass()+35.0);
  withMassCut(FITD0, 
	      m_ptypeD0.mass()-35.0, 
	      m_ptypeD0.mass()+35.0);

  // Momentum* Cut to Fitted D0(B)
  withPSCut(FITD0B,0.0,2.5);
  withPSCut(FITD0,0.0,2.5);

  vector<Particle> DSTM;
  vector<Particle> DSTP;
  vector<Particle> MISSDSTM;
  vector<Particle> MISSDSTP;

  // Make D*+- (with No Cut)
  combination(DSTM, m_ptypeDStarM, FITD0B, pionM);
  combination(DSTP, m_ptypeDStarP, FITD0, pionP);
  combination(MISSDSTM, m_ptypeDStarM, FITD0, pionM);
  combination(MISSDSTP, m_ptypeDStarP, FITD0B, pionP);

  // Set Generator Level Information to D*+-
  setGenHepInfoR(DSTM);
  setGenHepInfoR(DSTP);
  setGenHepInfoR(MISSDSTM);
  setGenHepInfoR(MISSDSTP);

  // D*+- and Lepton Found!!
  // Let's Start Analysis
  // Like Signal
  if(DSTM.size() > 0 && elecP.size() > 0){
    for(unsigned i=0;i<DSTM.size();++i){
      for(unsigned j=0;j<elecP.size();++j){
	mainAnalysis(DSTM[i],elecP[j],pionP,pionM,0);
      }
    }
  }
  if(DSTM.size() > 0 && muonP.size() > 0){
    for(unsigned i=0;i<DSTM.size();++i){
      for(unsigned j=0;j<muonP.size();++j){
	mainAnalysis(DSTM[i],muonP[j],pionP,pionM,0);
      }
    }
  }
  if(DSTP.size() > 0 && elecM.size() > 0){
    for(unsigned i=0;i<DSTP.size();++i){
      for(unsigned j=0;j<elecM.size();++j){
	mainAnalysis(DSTP[i],elecM[j],pionP,pionM,0);
      }
    }
  }
  if(DSTP.size() > 0 && muonM.size() > 0){
    for(unsigned i=0;i<DSTP.size();++i){
      for(unsigned j=0;j<muonM.size();++j){
	mainAnalysis(DSTP[i],muonM[j],pionP,pionM,0);
      }
    }
  }

  // Like Noise
  // Lepton and D* charge miss-match
  if(DSTM.size() > 0 && elecM.size() > 0){
    for(unsigned i=0;i<DSTM.size();++i){
      for(unsigned j=0;j<elecM.size();++j){
	mainAnalysis(DSTM[i],elecM[j],pionP,pionM,1);
      }
    }
  }
  if(DSTM.size() > 0 && muonM.size() > 0){
    for(unsigned i=0;i<DSTM.size();++i){
      for(unsigned j=0;j<muonM.size();++j){
	mainAnalysis(DSTM[i],muonM[j],pionP,pionM,1);
      }
    }
  }
  if(DSTP.size() > 0 && elecP.size() > 0){
    for(unsigned i=0;i<DSTP.size();++i){
      for(unsigned j=0;j<elecP.size();++j){
	mainAnalysis(DSTP[i],elecP[j],pionP,pionM,1);
      }
    }
  }
  if(DSTP.size() > 0 && muonP.size() > 0){
    for(unsigned i=0;i<DSTP.size();++i){
      for(unsigned j=0;j<muonP.size();++j){
	mainAnalysis(DSTP[i],muonP[j],pionP,pionM,1);
      }
    }
  }
  // Lepton and slow pion charge miss-match
  if(MISSDSTM.size() > 0 && elecM.size() > 0){
    for(unsigned i=0;i<MISSDSTM.size();++i){
      for(unsigned j=0;j<elecM.size();++j){
	mainAnalysis(MISSDSTM[i],elecM[j],pionP,pionM,2);
      }
    }
  }
  if(MISSDSTM.size() > 0 && muonM.size() > 0){
    for(unsigned i=0;i<MISSDSTM.size();++i){
      for(unsigned j=0;j<muonM.size();++j){
	mainAnalysis(MISSDSTM[i],muonM[j],pionP,pionM,2);
      }
    }
  }
  if(MISSDSTP.size() > 0 && elecP.size() > 0){
    for(unsigned i=0;i<MISSDSTP.size();++i){
      for(unsigned j=0;j<elecP.size();++j){
	mainAnalysis(MISSDSTP[i],elecP[j],pionP,pionM,2);
      }
    }
  }
  if(MISSDSTP.size() > 0 && muonP.size() > 0){
    for(unsigned i=0;i<MISSDSTP.size();++i){
      for(unsigned j=0;j<muonP.size();++j){
	mainAnalysis(MISSDSTP[i],muonP[j],pionP,pionM,2);
      }
    }
  }

  deleteDeepCopiedObjects(FITD0B);
  deleteDeepCopiedObjects(FITD0);
  deleteDeepCopiedObjects(FITN_D0B);
  deleteDeepCopiedObjects(FITN_D0);
#endif /* DOANA */
}
 
void 
jpdln_module::directionCosTheta(Particle &dStar, Particle &lepton, double *cosTheta){
  // cosTheta[0] ... in 3D
  // cosTheta[1] ... in 2D
  
  HepLorentzVector LVecDStar(pStar(dStar));
  HepLorentzVector LVecLepton(pStar(lepton));

  Hep3Vector VecDStar(LVecDStar.vect());
  Hep3Vector VecLepton(LVecLepton.vect());

  cosTheta[0] = VecDStar.unit()*VecLepton.unit();
  
  VecDStar.setZ(0.0);
  VecLepton.setZ(0.0);

  cosTheta[1] = VecDStar.unit()*VecLepton.unit();
#if 0
  if(cosTheta[1] > 0.95){
    dout(Debugout::INFO,"exdln") << cosTheta[0] << ", " << cosTheta[1] << std::endl;
    dout(Debugout::INFO,"exdln") << LVecDStar << std::endl;
    dout(Debugout::INFO,"exdln") << LVecLepton << std::endl;
    dout(Debugout::INFO,"exdln") << lepton.p() << std::endl;
    dout(Debugout::INFO,"exdln") << std::endl;
  }
#endif
}

double 
jpdln_module::missMass2(Particle &dStar, Particle &lepton, double &ppTerm){
  double e = 7.996;
  double p = 3.5;
  double upsilon_2 = sqrt(e*p); // upsilon mass*0.5
  double pB = sqrt(upsilon_2*upsilon_2-5.28*5.28);
  HepLorentzVector boost_vector(0., 0., e-p, e+p);
  HepLorentzVector dStarLepton(dStar.momentum().p()+lepton.momentum().p());
  dStarLepton.boost(-boost_vector.boostVector());
  double mass = dStarLepton.mag();
  double energy = dStarLepton.t();
  double mom = dStarLepton.vect().mag();

  ppTerm = 2.0*pB*mom;
  
  return 5.28*5.28+mass*mass-2.0*energy*upsilon_2;
}

int
jpdln_module::flavor(Particle &dStar, Particle &lepton, unsigned &method){
  int usedID[4];
  usedID[0] = lepton.relation().mdstCharged().get_ID();
  for(unsigned k=0;k<dStar.relation().nFinalStateParticles();++k){
    if(k >= 3){dout(Debugout::WARN,"exdln") << "Really?? BUG in D*ln or Particle Class!" << std::endl;break;}
    usedID[1+k] = dStar.relation().finalStateParticle(k).relation().mdstCharged().get_ID();
  }
  int gamma;
  //  m_fbtag.f_accum_list( 4, usedID, 0, &gamma );
  // Do flavor-tagging
  //  m_fbtag.f_tag();
  // Lepton or Kaon?
  char* chlepton = "lepton   ";
  char* chkaon   = "kaon     ";
  method = 2;
//    if(strcmp(m_fbtag.mode(),chlepton)==0){
//      method = 0;
//    }else if(strcmp(m_fbtag.mode(),chkaon)==0){
//      method = 1;  
//    }
//    return m_fbtag.flavor();
  return 0;
}

unsigned
jpdln_module::vertexing(kvertexfitter &kv,
			Particle &dStar, 
			Particle &lepton,
			vector<Particle> &pionP,
			vector<Particle> &pionM){
  vector<Particle> pionPList(pionP);
  vector<Particle> pionMList(pionM);
  removeParticle(pionPList,lepton);
  for(unsigned k=0;k<dStar.relation().nFinalStateParticles();++k){
    if(dStar.relation().finalStateParticle(k).pType().charge() > 0.)
      removeParticle(pionPList,dStar.relation().finalStateParticle(k));
    else
      removeParticle(pionMList,dStar.relation().finalStateParticle(k));
  }
  for(unsigned k=0;k<pionPList.size();++k)
    ::addTrack2fit(kv, pionPList[k]);
  for(unsigned k=0;k<pionMList.size();++k)
    ::addTrack2fit(kv, pionMList[k]);
  unsigned err = kv.fit();
  if(err != 0)return 0;
  return 1;
}

void 
jpdln_module::mainAnalysis(Particle &dStar, 
			   Particle &lepton,
			   vector<Particle> &pionP,
			   vector<Particle> &pionM,
			   unsigned sn){
#define JPDLN_DUMP_INFO 0

  // Boost Vector
  double e = 7.996;
  double p = 3.5;
  HepLorentzVector boost_vector(0., 0., e-p, e+p);

  // Vertex Fit to D* and Lepton, and Make D*
  kvertexfitter kv;
  Particle fittedDStar = dStar.deepCopy();
  unsigned err = doKVFit(kv,fittedDStar,lepton);
  if(err == 0){
    fittedDStar.deepDelete();
    m_tuple->clearData();
    return;
  }

  // Mass Dif Cut to D* (large : +-8MeV)
  if(!withMassDifCut(fittedDStar, 
		     m_ptypeDStarM.mass()-
		     m_ptypeD0B.mass()-0.020,
		     m_ptypeDStarM.mass()-
		     m_ptypeD0B.mass()+0.100, 0)){
    fittedDStar.deepDelete();
    m_tuple->clearData();
    return;
  }
#if JPDLN_DUMP_INFO
  dout(Debugout::INFO,"exdln") << "Fitted D*  , Mom : " << fittedDStar.momentum().p() << std::endl;
  dout(Debugout::INFO,"exdln") << "Fitted D0  , Mom : " << fittedDStar.relation().child(0).fittedMomentum().p() << std::endl;
  dout(Debugout::INFO,"exdln") << "Slow PI    , Mom : " << fittedDStar.relation().child(1).momentum().p() << std::endl;
  dout(Debugout::INFO,"exdln") << "Lepton     , Mom : " << lepton.fittedMomentum().p() << std::endl;
  dout(Debugout::INFO,"exdln") << "D*ln Side Vertex : " << kv.vertex() << std::endl;
  dout(Debugout::INFO,"exdln") << "D*ln Side dVertex: " << kv.errVertex() << std::endl;
#endif
  Mdst_event_Manager &evtMgr  = Mdst_event_Manager::get_manager();
  m_tuple->column("expnum", evtMgr.begin()->exp_no());
  m_tuple->column("runnum", evtMgr.begin()->run_no());
  m_tuple->column("evtnum", evtMgr.begin()->event_no());
  m_tuple->column("trig", evtMgr.begin()->trigger());
  m_tuple->column("date", evtMgr.begin()->date());
  m_tuple->column("time", evtMgr.begin()->time());
  m_tuple->column("snmon", sn);
  HepLorentzVector tmp;
  tmp = lepton.fittedMomentum().p();
  tmp.boost(-boost_vector.boostVector());
  m_tuple->column("ltype", lepton.pType().lund());
  m_tuple->column("lmom", lepton.fittedMomentum().p().vect().mag());
  m_tuple->column("lmoms", tmp.vect().mag());
  tmp = fittedDStar.relation().child(0).relation().child(0).fittedMomentum().p();
  tmp.boost(-boost_vector.boostVector());
  m_tuple->column("kmom", fittedDStar.relation().child(0).relation().child(0).fittedMomentum().p().vect().mag());
  m_tuple->column("kmoms", tmp.vect().mag());
  tmp = fittedDStar.relation().child(0).relation().child(1).fittedMomentum().p();
  tmp.boost(-boost_vector.boostVector());
  m_tuple->column("pmom", fittedDStar.relation().child(0).relation().child(1).fittedMomentum().p().vect().mag());
  m_tuple->column("pmoms", tmp.vect().mag());
  tmp = fittedDStar.relation().child(1).momentum().p();
  tmp.boost(-boost_vector.boostVector());
  m_tuple->column("spmom", fittedDStar.relation().child(1).momentum().p().vect().mag());
  m_tuple->column("spmoms", tmp.vect().mag());
  m_tuple->column("d0mass", fittedDStar.relation().child(0).momentum().mass());
  m_tuple->column("dsmass", fittedDStar.momentum().mass());
  m_tuple->column("massdif", fittedDStar.momentum().mass()-fittedDStar.relation().child(0).momentum().mass());
  m_tuple->column("dlnvtxx", kv.vertex().x());
  m_tuple->column("dlnvtxy", kv.vertex().y());
  m_tuple->column("dlnvtxz", kv.vertex().z());
  m_tuple->column("dlndvtxx", kv.errVertex()[0][0]);
  m_tuple->column("dlndvtxy", kv.errVertex()[1][1]);
  m_tuple->column("dlndvtxz", kv.errVertex()[2][2]);
  m_tuple->column("dlncl", kv.cl());
  m_tuple->column("dlnchisq", kv.chisq());
  m_tuple->column("dlndof", kv.dgf());
  m_tuple->column("d0vtxx", fittedDStar.relation().child(0).momentum().decayVertex().x());
  m_tuple->column("d0vtxy", fittedDStar.relation().child(0).momentum().decayVertex().y());
  m_tuple->column("d0vtxz", fittedDStar.relation().child(0).momentum().decayVertex().z());
  m_tuple->column("d0dvtxx", fittedDStar.relation().child(0).momentum().dDecayVertex()[0][0]);
  m_tuple->column("d0dvtxy", fittedDStar.relation().child(0).momentum().dDecayVertex()[1][1]);
  m_tuple->column("d0dvtxz", fittedDStar.relation().child(0).momentum().dDecayVertex()[2][2]);
  m_tuple->column("d0cl", fittedDStar.relation().child(0).relation().child(0).fittedMomentum().cl());
  m_tuple->column("d0chisq", fittedDStar.relation().child(0).relation().child(0).fittedMomentum().chisq());
  m_tuple->column("d0dof", fittedDStar.relation().child(0).relation().child(0).fittedMomentum().dof());

  dStarLeptonCheck(fittedDStar,lepton);

  genCheck(fittedDStar,lepton);

  // Calcu Info. of Direction between D* and Lepton
  double cosTheta[2];
  directionCosTheta(fittedDStar,lepton,cosTheta);
#if JPDLN_DUMP_INFO
  dout(Debugout::INFO,"exdln") << "Hemisphere       : " << cosTheta[0] << ", " << cosTheta[1] << std::endl;
#endif  
  m_tuple->column("hemisp0", cosTheta[0]);
  m_tuple->column("hemisp1", cosTheta[1]);

  // Calcu Info. of Miss Nuetrino...No Cut in this Routine
  double ppTerm;
  double missMassSqr = missMass2(fittedDStar,lepton,ppTerm);
#if JPDLN_DUMP_INFO
  dout(Debugout::INFO,"exdln") << "Miss Mass        : " << missMassSqr << ", " << ppTerm << std::endl;
#endif
  m_tuple->column("missmass", missMassSqr);
  m_tuple->column("misspp", ppTerm);

  // Flavor Tagging
  unsigned methodFlavor;
  int dlnFlavor = flavor(fittedDStar,lepton,methodFlavor);
#if JPDLN_DUMP_INFO
  dout(Debugout::INFO,"exdln") << "Flavor of D*ln   : " << lepton.pType().charge() << std::endl;
  dout(Debugout::INFO,"exdln") << "Flavor by Tag    : " << dlnFlavor << std::endl;
#endif
  m_tuple->column("lepchg", lepton.pType().charge());
  m_tuple->column("dschg", fittedDStar.pType().charge());
  if(methodFlavor == 0){
    //lepton
    m_tuple->column("flvl", dlnFlavor);
    m_tuple->column("flvk", 0.5);
    m_tuple->column("flvo", 0.5);
  }else if(methodFlavor == 1){
    //kaon
    m_tuple->column("flvl", 0.5);
    m_tuple->column("flvk", dlnFlavor);
    m_tuple->column("flvo", 0.5);
  }else{
    m_tuple->column("flvl", 0.5);
    m_tuple->column("flvk", 0.5);
    m_tuple->column("flvo", dlnFlavor);
  }
  // TagSide Vertexing
  kvertexfitter tagKv;
  err = vertexing(tagKv,fittedDStar,lepton,pionP,pionM);
#if JPDLN_DUMP_INFO
  if(err == 1){
    dout(Debugout::INFO,"exdln") << "Tag Side Vertex  : " << tagKv.vertex() << std::endl;
    dout(Debugout::INFO,"exdln") << "Tag Side dVertex : " << tagKv.errVertex() << std::endl;
  }
  dout(Debugout::INFO,"exdln") << std::endl;
#endif
  if(err == 1){
    m_tuple->column("tagvtxx", tagKv.vertex().x());
    m_tuple->column("tagvtxy", tagKv.vertex().y());
    m_tuple->column("tagvtxz", tagKv.vertex().z());
    m_tuple->column("tagdvtxx", tagKv.errVertex()[0][0]);
    m_tuple->column("tagdvtxy", tagKv.errVertex()[1][1]);
    m_tuple->column("tagdvtxz", tagKv.errVertex()[2][2]);
  }

  m_tuple->dumpData();

  // Delete Deep Copied Objects
  fittedDStar.deepDelete();
}

void
jpdln_module::doKMFit(vector<Particle> &plist){
  for(unsigned i=0;i<plist.size();++i){
    doKMFit(plist[i]);
  }
}

unsigned
jpdln_module::doKMFit(Particle &p){
  kmassvertexfitter kmv;
  kmv.invariantMass(p.pType().mass());
  for(unsigned j=0;j<p.relation().nChildren();++j){
    ::addTrack2fit(kmv,p.relation().child(j));
  }
  unsigned err = kmv.fit();
  if(err)return 0;
  return ::makeMother(kmv,p);
}

unsigned
jpdln_module::doKVFit(kvertexfitter &kv,
		      Particle &dStar,
		      Particle &lepton) {
  //...[doKVFit, addTrackKVFit, makeDStar fittedLepton]
  addTrack2fit(kv,dStar,lepton);
  unsigned err = kv.fit();
  if(err)return 0;
  fittedLepton(kv,lepton);
  return makeDStar(kv,dStar);
}

void 
jpdln_module::addTrack2fit(kvertexfitter &kv,
			   Particle &dStar,
			   Particle &lepton) {
  //...[doKVFit, addTrackKVFit, makeDStar fittedLepton]
  ::addTrack2fit(kv,dStar.relation().child(0));
  //::addTrackKVFit(kv,dStar.relation().child(1));
  ::addTrack2fit(kv,lepton);
}

unsigned
jpdln_module::makeDStar(kvertexfitter &kv,
			Particle &dStar) {
  //...[doKVFit, addTrackKVFit, makeDStar fittedLepton]
  kmakemother kmm;
  // Fitted D0
  kmm.addTrack(kv.momentum(0),
	       kv.position(0),
	       kv.error(0),
	       dStar.relation().child(0).pType().charge());
  kmm.errVertexTrack(kv.errVertexTrack(0));
  // Not Fitted Slow PI
  kmm.addTrack(dStar.relation().child(1).momentum().p(),
	       dStar.relation().child(1).momentum().x(),
	       dStar.relation().child(1).momentum().dpx(),
	       dStar.relation().child(1).pType().charge());
  kmm.errVertexTrack();
  // Vertex Info.
  kmm.vertex(kv.vertex());
  kmm.errVertex(kv.errVertex());
  unsigned err = kmm.make();
  if(err != 0)return 0;
  dStar.momentum().momentumPosition(kmm.momentum(),
				    kmm.position(),
				    kmm.error());
  dStar.momentum().decayVertex(kv.vertex(), kv.errVertex());
 
  //...fitted children of D*, that is, D0
  dStar.relation().child(0).fittedMomentum().momentumPosition(kv.momentum(0),
							      kv.position(0),
							      kv.error(0));
  dStar.relation().child(0).fittedMomentum().vertex(kv.vertex(), kv.errVertex());
  HepMatrix tmp1(kv.errVertexTrack(0));
  dStar.relation().child(0).fittedMomentum().coVertex(tmp1);

  return 1;
}

void
jpdln_module::fittedLepton(kvertexfitter &kv,
			   Particle &lepton) {
  //...[doKVFit, addTrackKVFit, makeDStar fittedLepton]
  //...fitted lepton
  lepton.fittedMomentum().momentumPosition(kv.momentum(1),
					   kv.position(1),
					   kv.error(1));
  lepton.fittedMomentum().vertex(kv.vertex(), kv.errVertex());
  HepMatrix tmp1(kv.errVertexTrack(1));
  lepton.fittedMomentum().coVertex(tmp1);
}

void
jpdln_module::genCheck(Particle &dStar, Particle &lepton) {
  // lepton
  if(lepton.relation().genHepevt() && lepton.relation().genHepevt().mother() &&
     abs(lepton.relation().genHepevt().mother().idhep()) == 511){
    m_tuple->column("glepchk", 1.0);
    if(lepton.relation().genHepevt().mother().idhep() == 511)
      m_tuple->column("glepbf", 1.0);
    else
      m_tuple->column("glepbf", -1.0);
  }else{
    m_tuple->column("glepchk", 0.0);
    m_tuple->column("glepbf", 0.0);
  }
  // D*
  if(dStar.relation().genHepevt() && dStar.relation().genHepevt().mother() &&
     abs(dStar.relation().genHepevt().mother().idhep()) == 511){
    m_tuple->column("gdschk", 1.0);
    if(dStar.relation().genHepevt().mother().idhep() == 511)
      m_tuple->column("gdsbf", 1.0);
    else
      m_tuple->column("gdsbf", -1.0);
    // Decay Point of B meson
    m_tuple->column("gbvtxx", dStar.relation().genHepevt().VX()*0.1);
    m_tuple->column("gbvtxy", dStar.relation().genHepevt().VY()*0.1);
    m_tuple->column("gbvtxz", dStar.relation().genHepevt().VZ()*0.1);
    // Decay Point of D0 meson
    m_tuple->column("gd0vtxx", dStar.relation().child(0).child(0).genHepevt().VX()*0.1);
    m_tuple->column("gd0vtxy", dStar.relation().child(0).child(0).genHepevt().VY()*0.1);
    m_tuple->column("gd0vtxz", dStar.relation().child(0).child(0).genHepevt().VZ()*0.1);
  }else{
    m_tuple->column("gdschk", 0.0);
    m_tuple->column("gdsbf", 0.0);
  }
}

void 
jpdln_module::dStarLeptonCheck(Particle &dStar, Particle &lepton) {
  if(lepton.mdstCharged().get_ID() == dStar.relation().finalStateParticle(0).mdstCharged().get_ID() ||
     lepton.mdstCharged().get_ID() == dStar.relation().finalStateParticle(1).mdstCharged().get_ID() ||
     lepton.mdstCharged().get_ID() == dStar.relation().finalStateParticle(2).mdstCharged().get_ID()){
    // bad
    m_tuple->column("ldschk", 0.);
  }else{
    // good
    m_tuple->column("ldschk", 1.);
  }
}

void
jpdln_module::genhepAnalysis(void) {
  unsigned sameFlavor, dlnFlag;
  unsigned id[6];// (511)l+, K+, PI-, (-511)l-, K-, PI+

  unsigned flag = dlnEvent(&sameFlavor, id, &dlnFlag);
  if(dlnFlag)m_dln_event_counter->accumulate(0.,1.);
  if(flag == 1)m_dln_b0_event_counter->accumulate(0.,1.);
  if(flag == 2)m_dln_b0b_event_counter->accumulate(0.,1.);
  m_tuple->column("evtmon",flag);

  //dout(Debugout::INFO,"exdln") << flag << std::endl;
  if(!flag)return;

  m_dln_kpi_event_counter->accumulate(0.,1.);

  Gen_hepevt_Manager &genMgr  = Gen_hepevt_Manager::get_manager();
  if(flag == 1){
    Hep3Vector l(genMgr[id[0]-1].PX(),genMgr[id[0]-1].PY(),genMgr[id[0]-1].PZ());
    m_tuple2->column("lmom", l.mag());
    m_tuple2->column("lmomt", l.perp());
    Hep3Vector k(genMgr[id[1]-1].PX(),genMgr[id[1]-1].PY(),genMgr[id[1]-1].PZ());
    m_tuple2->column("kmom", k.mag());
    m_tuple2->column("kmomt", k.perp());
    Hep3Vector p(genMgr[id[2]-1].PX(),genMgr[id[2]-1].PY(),genMgr[id[2]-1].PZ());
    m_tuple2->column("pmom", p.mag());
    m_tuple2->column("pmomt", p.perp());
    m_tuple2->dumpData();
  }else if(flag == 2){
    Hep3Vector l(genMgr[id[3]-1].PX(),genMgr[id[3]-1].PY(),genMgr[id[3]-1].PZ());
    m_tuple2->column("lmom", l.mag());
    m_tuple2->column("lmomt", l.perp());
    Hep3Vector k(genMgr[id[4]-1].PX(),genMgr[id[4]-1].PY(),genMgr[id[4]-1].PZ());
    m_tuple2->column("kmom", k.mag());
    m_tuple2->column("kmomt", k.perp());
    Hep3Vector p(genMgr[id[5]-1].PX(),genMgr[id[5]-1].PY(),genMgr[id[5]-1].PZ());
    m_tuple2->column("pmom", p.mag());
    m_tuple2->column("pmomt", p.perp());
    m_tuple2->dumpData();
  }else if(flag == 3 || flag == 4 || flag == 5){
    {
      Hep3Vector l(genMgr[id[0]-1].PX(),genMgr[id[0]-1].PY(),genMgr[id[0]-1].PZ());
      m_tuple2->column("lmom", l.mag());
      m_tuple2->column("lmomt", l.perp());
      Hep3Vector k(genMgr[id[1]-1].PX(),genMgr[id[1]-1].PY(),genMgr[id[1]-1].PZ());
      m_tuple2->column("kmom", k.mag());
      m_tuple2->column("kmomt", k.perp());
      Hep3Vector p(genMgr[id[2]-1].PX(),genMgr[id[2]-1].PY(),genMgr[id[2]-1].PZ());
      m_tuple2->column("pmom", p.mag());
      m_tuple2->column("pmomt", p.perp());
      m_tuple2->dumpData();
    }
    {      
      Hep3Vector l(genMgr[id[3]-1].PX(),genMgr[id[3]-1].PY(),genMgr[id[3]-1].PZ());
      m_tuple2->column("lmom", l.mag());
      m_tuple2->column("lmomt", l.perp());
      Hep3Vector k(genMgr[id[4]-1].PX(),genMgr[id[4]-1].PY(),genMgr[id[4]-1].PZ());
      m_tuple2->column("kmom", k.mag());
      m_tuple2->column("kmomt", k.perp());
      Hep3Vector p(genMgr[id[5]-1].PX(),genMgr[id[5]-1].PY(),genMgr[id[5]-1].PZ());
      m_tuple2->column("pmom", p.mag());
      m_tuple2->column("pmomt", p.perp());
      m_tuple2->dumpData();
    }
  }
}

unsigned 
jpdln_module::dlnEvent(unsigned *sameFlavor, unsigned *id, unsigned *dlnFlag) {
  //return 0 ... not dln event
  //       1 ... dln event (B0)
  //       2 ... dln event (B0B)
  //       3 ... dln event (B0,  B0B)
  //       4 ... dln event (B0,  B0)
  //       5 ... dln event (B0B, B0B)
  Gen_hepevt_Manager &genMgr  = Gen_hepevt_Manager::get_manager();

  *sameFlavor = 0;
  *dlnFlag    = 0;
  unsigned flagB0(0);
  unsigned flagB0B(0);

  const int elecID = 11;  //11  -
  const int muonID = 13;  //13  -
  const int nelecID = 12; //12  +
  const int nmuonID = 14; //14  +
  const int dStarID = 413; //413 +

  for(vector<Gen_hepevt>::iterator i = genMgr.begin();
      i != genMgr.end(); ++i){
    unsigned tmpId[3];
    // B0
    if(i->idhep() == 511){
      if(abs(i->daLast()-i->daFirst()+1) == 3){
	unsigned countDln(0);
	unsigned dStarId(0);
	if(((genMgr[i->daFirst()-1].idhep() ==  nelecID &&
	     genMgr[i->daFirst()].idhep()  == -elecID) ||
	    (genMgr[i->daFirst()-1].idhep() ==  nmuonID &&
	     genMgr[i->daFirst()].idhep()  == -muonID)) &&
	   genMgr[i->daFirst()+1].idhep() == -dStarID){
	  countDln = 2;
	  tmpId[0] = genMgr[i->daFirst()].get_ID();
	  dStarId = genMgr[i->daFirst()+1].get_ID();
	}else if(((genMgr[i->daFirst()].idhep() ==  nelecID &&
		   genMgr[i->daFirst()-1].idhep()  == -elecID) ||
		  (genMgr[i->daFirst()].idhep() ==  nmuonID &&
		   genMgr[i->daFirst()-1].idhep()  == -muonID)) &&
		 genMgr[i->daFirst()+1].idhep()  == -dStarID){
	  countDln = 2;
	  tmpId[0] = genMgr[i->daFirst()-1].get_ID();
	  dStarId = genMgr[i->daFirst()+1].get_ID();
	}else if(((genMgr[i->daFirst()-1].idhep() ==  nelecID &&
		   genMgr[i->daFirst()+1].idhep()  == -elecID) ||
		  (genMgr[i->daFirst()-1].idhep() ==  nmuonID &&
		   genMgr[i->daFirst()+1].idhep()  == -muonID)) &&
		 genMgr[i->daFirst()].idhep()    == -dStarID){
	  countDln = 2;
	  tmpId[0] = genMgr[i->daFirst()+1].get_ID();
	  dStarId = genMgr[i->daFirst()].get_ID();
	}else if(((genMgr[i->daFirst()+1].idhep() ==  nelecID &&
		   genMgr[i->daFirst()-1].idhep()  == -elecID) ||
		  (genMgr[i->daFirst()+1].idhep() ==  nmuonID &&
		   genMgr[i->daFirst()-1].idhep()  == -muonID)) &&
		 genMgr[i->daFirst()].idhep()    == -dStarID){
	  countDln = 2;
	  tmpId[0] = genMgr[i->daFirst()-1].get_ID();
	  dStarId = genMgr[i->daFirst()].get_ID();
	}else if(((genMgr[i->daFirst()+1].idhep() ==  nelecID &&
		   genMgr[i->daFirst()].idhep()  == -elecID) ||
		  (genMgr[i->daFirst()+1].idhep() ==  nmuonID &&
		   genMgr[i->daFirst()].idhep()  == -muonID)) &&
		 genMgr[i->daFirst()-1].idhep() == -dStarID){
	  countDln = 2;
	  tmpId[0] = genMgr[i->daFirst()].get_ID();
	  dStarId = genMgr[i->daFirst()-1].get_ID();
	}else if(((genMgr[i->daFirst()].idhep() ==  nelecID &&
		   genMgr[i->daFirst()+1].idhep()  == -elecID) ||
		  (genMgr[i->daFirst()].idhep() ==  nmuonID &&
		   genMgr[i->daFirst()+1].idhep()  == -muonID)) &&
		 genMgr[i->daFirst()-1].idhep()  == -dStarID){
	  countDln = 2;
	  tmpId[0] = genMgr[i->daFirst()+1].get_ID();
	  dStarId = genMgr[i->daFirst()-1].get_ID();
	}
	if(countDln == 2){
	  *dlnFlag    = 1;
	  // dstar(-)
	  if(abs(genMgr[dStarId-1].daLast()-genMgr[dStarId-1].daFirst()+1) == 2){
	    if(genMgr[genMgr[dStarId-1].daFirst()-1].idhep() == -421 && 
	       genMgr[genMgr[dStarId-1].daLast()-1].idhep()  == -211){
	      // d0bar
	      unsigned tmp = genMgr[dStarId-1].daFirst();
	      if(abs(genMgr[tmp-1].daLast()-genMgr[tmp-1].daFirst()+1) == 2){
		unsigned id0 = genMgr[tmp-1].daFirst();
		unsigned id1 = id0+1;
		// k+, pi-
		if(genMgr[id0-1].idhep() ==  321 &&
		   genMgr[id1-1].idhep() == -211){
		  tmpId[1] = genMgr[id0-1].get_ID();
		  tmpId[2] = genMgr[id1-1].get_ID();
		  ++countDln;
		}else if(genMgr[id0-1].idhep() == -211 &&
			 genMgr[id1-1].idhep() ==  321){
		  tmpId[1] = genMgr[id1-1].get_ID();
		  tmpId[2] = genMgr[id0-1].get_ID();
		  ++countDln;
		}
	      }
	    }else if(genMgr[genMgr[dStarId-1].daFirst()-1].idhep() == -211 && 
		     genMgr[genMgr[dStarId-1].daLast()-1].idhep()  == -421){
	      // d0bar
	      unsigned tmp = genMgr[dStarId-1].daLast();
	      if(abs(genMgr[tmp-1].daLast()-genMgr[tmp-1].daFirst()+1) == 2){
		unsigned id0 = genMgr[tmp-1].daFirst();
		unsigned id1 = id0+1;
		// k+, pi-
		if(genMgr[id0-1].idhep() ==  321 &&
		   genMgr[id1-1].idhep() == -211){
		  tmpId[1] = genMgr[id0-1].get_ID();
		  tmpId[2] = genMgr[id1-1].get_ID();
		  ++countDln;
		}else if(genMgr[id0-1].idhep() == -211 &&
			 genMgr[id1-1].idhep() ==  321){
		  tmpId[1] = genMgr[id1-1].get_ID();
		  tmpId[2] = genMgr[id0-1].get_ID();
		  ++countDln;
		}
	      }
	    }
	  }
	}
	if(countDln == 3)++flagB0;
      }
      if(flagB0 == 1){
	id[0] = tmpId[0];
	id[1] = tmpId[1];
	id[2] = tmpId[2];
      }else if(flagB0 == 2){
	*sameFlavor = 1;
	id[3] = tmpId[0];
	id[4] = tmpId[1];
	id[5] = tmpId[2];
      }else if(flagB0 >= 3){
	dout(Debugout::WARN,"exdln") << "3 or more B->D*ln decay.....Strange....." << std::endl;
      }
    }
    // B0B
    if(i->idhep() == -511){
      if(abs(i->daLast()-i->daFirst()+1) == 3){
	unsigned countDln(0);
	unsigned dStarId(0);
	if(((genMgr[i->daFirst()-1].idhep() == -nelecID &&
	     genMgr[i->daFirst()].idhep()  == elecID) ||
	    (genMgr[i->daFirst()-1].idhep() == -nmuonID &&
	     genMgr[i->daFirst()].idhep()  == muonID)) &&
	   genMgr[i->daFirst()+1].idhep() ==  dStarID){
	  countDln = 2;
	  tmpId[0] = genMgr[i->daFirst()].get_ID();
	  dStarId = genMgr[i->daFirst()+1].get_ID();
	}else if(((genMgr[i->daFirst()].idhep() == -nelecID &&
		   genMgr[i->daFirst()-1].idhep()  == elecID) ||
		  (genMgr[i->daFirst()].idhep() == -nmuonID &&
		   genMgr[i->daFirst()-1].idhep()  == muonID)) &&
		 genMgr[i->daFirst()+1].idhep()  ==  dStarID){
	  countDln = 2;
	  tmpId[0] = genMgr[i->daFirst()-1].get_ID();
	  dStarId = genMgr[i->daFirst()+1].get_ID();
	}else if(((genMgr[i->daFirst()-1].idhep() == -nelecID &&
		   genMgr[i->daFirst()+1].idhep()  == elecID) ||
		  (genMgr[i->daFirst()-1].idhep() == -nmuonID &&
		   genMgr[i->daFirst()+1].idhep()  == muonID)) &&
		 genMgr[i->daFirst()].idhep()    ==  dStarID){
	  countDln = 2;
	  tmpId[0] = genMgr[i->daFirst()+1].get_ID();
	  dStarId = genMgr[i->daFirst()].get_ID();
	}else if(((genMgr[i->daFirst()+1].idhep() == -nelecID &&
		   genMgr[i->daFirst()-1].idhep()  == elecID) ||
		  (genMgr[i->daFirst()+1].idhep() == -nmuonID &&
		   genMgr[i->daFirst()-1].idhep()  == muonID)) &&
		 genMgr[i->daFirst()].idhep()    ==  dStarID){
	  countDln = 2;
	  tmpId[0] = genMgr[i->daFirst()-1].get_ID();
	  dStarId = genMgr[i->daFirst()].get_ID();
	}else if(((genMgr[i->daFirst()+1].idhep() == -nelecID &&
		   genMgr[i->daFirst()].idhep()  == elecID) ||
		  (genMgr[i->daFirst()+1].idhep() == -nmuonID &&
		   genMgr[i->daFirst()].idhep()  == muonID)) &&
		 genMgr[i->daFirst()-1].idhep() ==  dStarID){
	  countDln = 2;
	  tmpId[0] = genMgr[i->daFirst()].get_ID();
	  dStarId = genMgr[i->daFirst()-1].get_ID();
	}else if(((genMgr[i->daFirst()].idhep() == -nelecID &&
		   genMgr[i->daFirst()+1].idhep()  == elecID) ||
		  (genMgr[i->daFirst()].idhep() == -nmuonID &&
		   genMgr[i->daFirst()+1].idhep()  == muonID)) &&
		 genMgr[i->daFirst()-1].idhep()  ==  dStarID){
	  countDln = 2;
	  tmpId[0] = genMgr[i->daFirst()+1].get_ID();
	  dStarId = genMgr[i->daFirst()-1].get_ID();
	}
	if(countDln == 2){
	  *dlnFlag    = 1;
	  // dstar(+)
	  if(abs(genMgr[dStarId-1].daLast()-genMgr[dStarId-1].daFirst()+1) == 2){
	    if(genMgr[genMgr[dStarId-1].daFirst()-1].idhep() == 421 && 
	       genMgr[genMgr[dStarId-1].daLast()-1].idhep()  == 211){
	      // d0
	      unsigned tmp = genMgr[dStarId-1].daFirst();
	      if(abs(genMgr[tmp-1].daLast()-genMgr[tmp-1].daFirst()+1) == 2){
		unsigned id0 = genMgr[tmp-1].daFirst();
		unsigned id1 = id0+1;
		// k+, pi-
		if(genMgr[id0-1].idhep() == -321 &&
		   genMgr[id1-1].idhep() ==  211){
		  tmpId[1] = genMgr[id0-1].get_ID();
		  tmpId[2] = genMgr[id1-1].get_ID();
		  ++countDln;
		}else if(genMgr[id0-1].idhep() ==  211 &&
			 genMgr[id1-1].idhep() == -321){
		  tmpId[1] = genMgr[id1-1].get_ID();
		  tmpId[2] = genMgr[id0-1].get_ID();
		  ++countDln;
		}
	      }
	    }else if(genMgr[genMgr[dStarId-1].daFirst()-1].idhep() == 211 && 
		     genMgr[genMgr[dStarId-1].daLast()-1].idhep()  == 421){
	      // d0bar
	      unsigned tmp = genMgr[dStarId-1].daLast();
	      if(abs(genMgr[tmp-1].daLast()-genMgr[tmp-1].daFirst()+1) == 2){
		unsigned id0 = genMgr[tmp-1].daFirst();
		unsigned id1 = id0+1;
		// k+, pi-
		if(genMgr[id0-1].idhep() == -321 &&
		   genMgr[id1-1].idhep() ==  211){
		  tmpId[1] = genMgr[id0-1].get_ID();
		  tmpId[2] = genMgr[id1-1].get_ID();
		  ++countDln;
		}else if(genMgr[id0-1].idhep() ==  211 &&
			 genMgr[id1-1].idhep() == -321){
		  tmpId[1] = genMgr[id1-1].get_ID();
		  tmpId[2] = genMgr[id0-1].get_ID();
		  ++countDln;
		}
	      }
	    }
	  }
	}
	if(countDln == 3)++flagB0B;
      }
      if(flagB0B == 1){
	id[3] = tmpId[0];
	id[4] = tmpId[1];
	id[5] = tmpId[2];
      }else if(flagB0B == 2){
	*sameFlavor = 1;
	id[0] = tmpId[0];
	id[1] = tmpId[1];
	id[2] = tmpId[2];
      }else if(flagB0B >= 3){
	dout(Debugout::WARN,"exdln") << "3 or more B->D*ln decay.....Strange....." << std::endl;
      }
    }
  }
  if(flagB0 == 0 && flagB0B == 0)return 0;
  if(flagB0 == 1 && flagB0B == 0)return 1;
  if(flagB0 == 0 && flagB0B == 1)return 2;
  if(flagB0 == 1 && flagB0B == 1)return 3;
  if(flagB0 == 2 && flagB0B == 0)return 4;
  if(flagB0 == 0 && flagB0B == 2)return 5;
  return 0;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
