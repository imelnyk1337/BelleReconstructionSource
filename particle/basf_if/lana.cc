//-----------------------------------------------------------------------------
// $Id: lana.cc 9944 2006-11-29 07:36:07Z katayama $
//-----------------------------------------------------------------------------
// Filename : lepton.cc
// Section  : particle analysis
// Owner    : Yoshi Iwasaki
// Email    : yoshihito.iwasaki@kek.jp
//-----------------------------------------------------------------------------
// Description : Lepton analysis module.
//-----------------------------------------------------------------------------
// $Log$
// Revision 1.7  2001/12/13 15:31:53  katayama
// MDST_OBS
//
// Revision 1.6  2000/04/13 12:41:54  katayama
// Added std:: to cout,cerr,endl etc.
//
// Revision 1.5  2000/03/18 07:40:51  katayama
// fix the infinite loop problem of disp_stat
//
// Revision 1.4  2000/03/17 05:05:00  katayama
// Module functions modified
//
// Revision 1.3  2000/03/08 12:41:28  katayama
// compatibility with CC
//
// Revision 1.2  1999/01/16 10:31:26  katayama
// clean up includes
//
// Revision 1.1  1998/09/22 15:38:30  yiwasaki
// lana.cc(example of ParticleManager) added
//
//-----------------------------------------------------------------------------

#include "belle.h"
#include "particle/Particle.h"
#include "particle/ParticleManager.h"
#include "event/BelleEvent.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "tuple/BelleTupleManager.h"

#include "panther/panther.h"
#include MDST_H
#include HEPEVT_H
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


// Module class
class LeptonAnalysis : public Module {

  public:
    LeptonAnalysis(void);
    ~LeptonAnalysis(void);
    void init(int *status);
    void term(void);
    void disp_stat(const char *);
    void disp_stat ( char* ) { disp_stat((const char *)""); }
    void event(BelleEvent *, int *);
  void begin_run(BelleEvent *, int *){}
  void end_run(BelleEvent *, int *){}
  

  private:
    unsigned _nEvents;
};

extern "C" Module_descr * mdcl_lana() {
  LeptonAnalysis * m = new LeptonAnalysis();
  Module_descr * d = new Module_descr ("lana", m);
  return d;
}

LeptonAnalysis::LeptonAnalysis(void) : _nEvents(0) {

    //...Do nothing...
}

LeptonAnalysis::~LeptonAnalysis(void) {

    //...Do nothing...
}

void 
LeptonAnalysis::init(int *status) {
  *status = 0;
    //...Initialize histograms...
//     extern BelleTupleManager *BASF_Histogram;
//     m_jpsi   = BASF_Histogram->histogram("Mass(J/psi - GeV)",100,2.5,3.5);
//     m_vee    = BASF_Histogram->histogram("Mass(Vee - GeV)",100,0.,1.);
//     m_b      = BASF_Histogram->histogram("Mass(B - GeV)",50,5.0,5.5);
//     m_b_beam = BASF_Histogram->histogram("Mass(B.Beam - GeV)",50,5.0,5.5);
}

void 
LeptonAnalysis::term(void) {
    const char * c = "";
    disp_stat(c);
}

void 
LeptonAnalysis::event(BelleEvent *evptr, int *status) {
    ++_nEvents;

    ParticleManager gm(Gen_hepevt_Manager::get_manager());
    ParticleManager cm(Mdst_charged_Manager::get_manager());
    gm.link(cm);

    //...Select muons...
    ParticleManager muons("muons");
    for (unsigned i = 0; i < cm.size(); i++) {
	Particle & p = * cm[i];
//  	if (p.relation().mdstCharged().klm().muon() >= 2)
//  	    muons.push_back(& p);
    }

    //...For debug...
    if (muons.size()) {
	gm.dump();
	cm.dump();
	muons.dump();
    }
}

void
LeptonAnalysis::disp_stat(const char *) {
    dout(Debugout::RESULT,"lana") << ">>> LeptonAnalysis Summary <<<" << std::endl;
    dout(Debugout::RESULT,"lana") << "    # of events processed : " << _nEvents << std::endl;
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
