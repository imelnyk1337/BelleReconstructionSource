//
// ParticleManager.cc
//
// A class to handle particle objects.
//
// Filename : ParticleManger.cc
// Author : Yoshihito Iwasaki
// E-mail : yoshihito.iwasaki@kek.jp
//
// $Id: ParticleManager.cc 10613 2008-09-03 11:36:24Z katayama $
//
// $Log$
// Revision 1.13  2004/04/23 05:56:42  katayama
// Use count
//
// Revision 1.12  2002/02/25 01:44:09  katayama
// For -W warning suppression.
//
// Revision 1.11  2002/02/22 06:37:39  katayama
// Use __sparc
//
// Revision 1.10  2002/02/21 23:46:03  katayama
// For -ansi and other warning flags
//
// Revision 1.9  2001/12/23 09:58:25  katayama
// removed Strings.h
//
// Revision 1.8  2001/12/19 02:40:02  katayama
// Use itostring()
//
// Revision 1.7  2001/12/12 07:10:51  jtanaka
// (1) compatibility for obsoleted Mdst_vee/mdst_sim_xref table
// (2) compatibility for gcc3
//
// Revision 1.6  2000/06/13 05:20:47  yiwasaki
// updates
//
// Revision 1.5  2000/05/18 05:41:20  katayama
// Added panther/panther.h
//
// Revision 1.4  2000/05/03 09:25:29  yiwasaki
// minor changes
//
// Revision 1.3  2000/04/13 12:41:55  katayama
// Added std:: to cout,cerr,endl etc.
//
// Revision 1.2  2000/03/07 11:13:56  katayama
// compatibility with CC5.0
//
// Revision 1.1  2000/01/05 06:13:43  jtanaka
// major updates:please see BELLE whiteboard.
//
// Revision 1.7  1999/01/16 10:31:39  katayama
// clean up includes
//
// Revision 1.6  1998/12/14 08:10:04  katayama
// Because I removed static_Particle etc, I had to change these so that the
// library is consistent. There was another way to fix but I decided not to
// expose static_Particle etc. Authors, please check the changes I made.
//
// Revision 1.5  1998/10/26 08:41:18  jtanaka
// grandchildren -> finalStateParticles
//
// Revision 1.4  1998/10/12 19:24:09  jtanaka
// updated constructors Momentum.cc, Particle.cc, and Relation.cc. removed and added some functions in PID.cc. added (int*) to ParticleManager.cc(warning).
//
// Revision 1.3  1998/09/22 15:34:28  yiwasaki
// ParticleManager extended, Relation::mc modified
//
// Revision 1.2  1998/09/22 01:24:26  yiwasaki
// ParticleManager minor change
//
// Revision 1.1  1998/09/08 16:49:57  yiwasaki
// ParticleManager added
//
//
/* for isnan */
#if defined(__sparc)
#  if defined(__EXTENSIONS__)
#    include <cmath>
#  else
#    define __EXTENSIONS__
#    include <cmath>
#    undef __EXTENSIONS__
#  endif
#elif defined(__GNUC__)
#  if defined(_XOPEN_SOURCE)
#    include <cmath>
#  else
#    define _XOPEN_SOURCE
#    include <cmath>
#    undef _XOPEN_SOURCE
#  endif
#endif

#include "particle/ParticleManager.h"
#include "belle.h"
#include HEPEVT_H
#include MDST_H
#include "belleCLHEP/String/Strings.h"
#include "mdst/mdst.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


ParticleManager::ParticleManager()
: m_delete(false) {
}

ParticleManager::ParticleManager(const std::string & name)
: m_name(name), m_delete(false) {
}

ParticleManager::ParticleManager(Gen_hepevt_Manager & m) : m_delete(true) {
    m_name = "gen_hepevt manager";
//  unsigned n = m.size();
//  unsigned n = BsCouTab(GEN_HEPEVT);    
    unsigned n = m.count();
    for (unsigned i = 0; i < n; i++) {
	Gen_hepevt & g = m[i];
	Particle * p = new Particle(g);

	//	std::string nameTmp((int) i);
	//      std::string name = "mc";
	//	name += (const char *) nameTmp;
	//	p->name(name);
	std::string name = "mc";
	name += itostring(i);
	p->name(name);

	if (g.mother_ID()) {
	    Particle * m = (* this)[g.mother_ID() - 1];
	    p->relation().mother(* m);
	    m->relation().append(* p);
	}

	push_back(p);
    }
}

ParticleManager::ParticleManager(Mdst_charged_Manager & m) : m_delete(true) {
    m_name = "mdst_charged manager";
//  unsigned n = m.size();
//  unsigned n = BsCouTab(MDST_CHARGED);
    unsigned n = m.count();
    Ptype pip("PI+");
    Ptype pim("PI-");

// dout(Debugout::DDEBUG,"ParticleManager") << "ParticleManager ... n=" << n << std::endl;
    
    for (unsigned i = 0; i < n; i++) {
	Mdst_charged & g = m[i];

	Particle * p;
	if (g.charge() > 0.) p = new Particle(g, pip);
	else                 p = new Particle(g, pim);

	//	std::string nameTmp((int) i);
	//	std::string name = "c";
	//	name += (const char *) nameTmp;
	std::string name = "c";
	name += itostring(i);
	p->name(name);

	if (isnan(p->px()) || (isnan(p->pz())) || (isnan(p->py())))  {
	    Momentum & m = p->momentum();
	    m.momentum(HepLorentzVector(0., 0., 0., 0.));
	    m.position(HepPoint3D(0., 0., 0.));
	    p->usable(false);
	}
//	if (isnan(p->px())) {
//  	    delete p;
//  	    continue;
//  	}
//  	if (isnan(p->py())) {
//  	    delete p;
//  	    continue;
//  	}
//  	if (isnan(p->pz())) {
//  	    delete p;
//  	    continue;
//  	}

	push_back(p);
    }
}

ParticleManager::~ParticleManager() {
    if (! m_delete) return;

    unsigned n = size();
    for (unsigned i = 0; i < n; i++) {
	Particle * g = (*this)[i];
	delete g;
    }
}

void
ParticleManager::dump(const std::string &,
		      const std::string & prefix) {
    unsigned n = size();

    dout(Debugout::DUMP,"ParticleManager") << prefix;
    dout(Debugout::DUMP,"ParticleManager") << "list of " << m_name << " : #particles=" << n << std::endl;
    if (! n) return;

    dout(Debugout::DUMP,"ParticleManager") << "name type   mass     px     py     pz     vx     vy     vz";
    dout(Debugout::DUMP,"ParticleManager") << "   mo    mc" << std::endl;

//    bool hasNan = false;

    for (unsigned i = 0; i < n; i++) {
	Particle & g = * (*this)[i];

	dout(Debugout::DUMP,"ParticleManager") << prefix;
	dout(Debugout::DUMP,"ParticleManager").width(4);
	dout(Debugout::DUMP,"ParticleManager") << g.name().c_str() << " ";
	if (! g.pType()) {
	    dout(Debugout::DUMP,"ParticleManager") << "unknown" << std::endl;
	    continue;
	}
	dout(Debugout::DUMP,"ParticleManager").width(4);
	dout(Debugout::DUMP,"ParticleManager") << g.pType().name().c_str() << " ";
	dout(Debugout::DUMP,"ParticleManager").setf(std::ios::fixed,std::ios::floatfield);
	dout(Debugout::DUMP,"ParticleManager").precision(3);
	dout(Debugout::DUMP,"ParticleManager").width(6);
	dout(Debugout::DUMP,"ParticleManager") << g.mass() << " ";
	dout(Debugout::DUMP,"ParticleManager").width(6);
	dout(Debugout::DUMP,"ParticleManager") << g.px() << " ";
	dout(Debugout::DUMP,"ParticleManager").width(6);
	dout(Debugout::DUMP,"ParticleManager") << g.py() << " ";
	dout(Debugout::DUMP,"ParticleManager").width(6);
	dout(Debugout::DUMP,"ParticleManager") << g.pz() << " ";
	dout(Debugout::DUMP,"ParticleManager").precision(3);
	dout(Debugout::DUMP,"ParticleManager").width(6);
	dout(Debugout::DUMP,"ParticleManager") << g.x().x() << " ";
	dout(Debugout::DUMP,"ParticleManager").width(6);
	dout(Debugout::DUMP,"ParticleManager") << g.x().y() << " ";
	dout(Debugout::DUMP,"ParticleManager").width(6);
	dout(Debugout::DUMP,"ParticleManager") << g.x().z() << " ";

	//	if ((& g.mc()) != (& static_Particle)) {
	if (g.mother()) {
	    dout(Debugout::DUMP,"ParticleManager").width(4);
	    dout(Debugout::DUMP,"ParticleManager") << g.mother().name().c_str();
	}
	else {
	    dout(Debugout::DUMP,"ParticleManager") << "     ";
	}
	if (g.mc()) {
	    dout(Debugout::DUMP,"ParticleManager").width(4);
	    dout(Debugout::DUMP,"ParticleManager") << g.mc().name().c_str();
	}

	dout(Debugout::DUMP,"ParticleManager") << std::endl;

//	if (isnan(g.px())) hasNan = true;
    }
//    if (hasNan) BsShwDat(MDST_CHARGED);
}

void
ParticleManager::link(const ParticleManager & m) const {

    //...Check an argument...
    if (name() == m.name()) {
	dout(Debugout::ERR,"ParticleManager") << "ParticleManager::link !!! failed to link";
	dout(Debugout::ERR,"ParticleManager") << " because both managers are same." << std::endl;
	return;
    }

    //...Search gen_hepevt manager...
    const ParticleManager * gm = 0;
    if (name() == "gen_hepevt manager") gm = this;
    if (m.name() == "gen_hepevt manager") gm = & m;
    if (! gm) {
	dout(Debugout::ERR,"ParticleManager") << "ParticleManager::link !!! failed to link";
	dout(Debugout::ERR,"ParticleManager") << " : gen_hepevt manager is not found." << std::endl;
	return;
    }
    if (! gm->size()) return;

    //...Search mdst_charged manager...
    const ParticleManager * cm = 0;
    if (name() == "mdst_charged manager") cm = this;
    if (m.name() == "mdst_charged manager") cm = & m;
    if (! cm) {
	dout(Debugout::ERR,"ParticleManager") << "ParticleManager::link !!! failed to link";
	dout(Debugout::ERR,"ParticleManager") << " : mdst_charged manager is not found." << std::endl;
    }

    //...Link by xref...
    Gen_hepevt_Manager & gmm = Gen_hepevt_Manager::get_manager();
    Mdst_charged_Manager & cmm = Mdst_charged_Manager::get_manager();

    // Following code assumes that gen_hepevt manager and mdst_charged
    // manager are not modified by users.

    unsigned n = cmm.size();

//      dout(Debugout::INFO,"ParticleManager") << "cm,gm,xm:n = " << cm->size() << "," << gm->size() << ",";
//      dout(Debugout::INFO,"ParticleManager") << xm.size() << std::endl;

    for (unsigned i = 0; i < n; i++) {
      const Gen_hepevt& gh = get_hepevt(cmm[i]);
	int cid = cmm[i].get_ID();
	int gid = gh ? (int)(gh.get_ID()) : 0;
	if (cid && gid) {
	    --cid;
	    --gid;

//  	    dout(Debugout::INFO,"ParticleManager") << "cid,gid=" << cid << "," << gid << std::endl;

	    Particle & c = * (* cm)[cid];
	    Particle & g = * (* gm)[gid];

//  	    c.dump();
//  	    g.dump();

	    c.relation().mc(g);
	    c.relation().genHepevt(gmm[gid]);
	    g.relation().mc(c);
	    g.relation().mdstCharged(cmm[cid]);
	}
    }
}

void
ParticleManager::setPType(const Ptype & a) {
    unsigned n = size();
    if (! n) return;

    double m2 = a.mass() * a.mass();
    static HepLorentzVector mom;
    for (unsigned i = 0; i < n; i++) {
	Particle & p = * (* this)[i];
	p.pType() = a;
	mom = p.p();
	mom.setE(sqrt(mom.vect().mag2() + m2));
	p.momentum().momentum(mom);
    }
}

void
ParticleManager::setPType(const Ptype & a, Particle & p) {
    double m2 = a.mass() * a.mass();
    static HepLorentzVector mom;
    p.pType() = a;
    mom = p.p();
    mom.setE(sqrt(mom.vect().mag2() + m2));
    p.momentum().momentum(mom);
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
