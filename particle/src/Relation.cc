// Relation.cc
// 
// A class implimentation for <Relation> class 
// for the BELLE standard (hopefully!) Particle object liblary
// 
// <Relation> class supplies you interfaces to other objects,
// such as mother particle, MC particle, and MDST banks.
//
// Filename : Particle.cc
// Author : YOKOyama Masashi
// e-mail:yokoyama@belaxp1.phys.s.u-tokyo.ac.jp
// --> e-mail:jtanaka@belaxp1.phys.s.u-tokyo.ac.jp
//
// $Id: Relation.cc 9948 2006-12-05 06:11:00Z katayama $ 
//
// Revision history
// 
// $Log$
// Revision 1.34  2003/09/15 13:24:34  katayama
// reduce warnings
//
// Revision 1.33  2001/12/13 15:31:55  katayama
// MDST_OBS
//
// Revision 1.32  2001/10/22 02:21:08  jtanaka
// fix : memory leak when pi0,vee,vee2 particles are removed from "vector/list etc".
//
// Revision 1.31  2001/07/20 07:39:07  jtanaka
// bug fix of the "remove" function.(from houjo-san)
//
// Revision 1.30  2000/05/29 11:31:53  jtanaka
// add const and mutable to some member functions and parameters.
//
// Revision 1.29  2000/05/16 14:47:01  jtanaka
// added constructor of Mdst_ecl
//
// Revision 1.28  2000/04/27 12:57:38  jtanaka
// 1. use Mdst_vee_daughters information to construct Mdst_vee2.
// 2. add new constructor Particle(Mdst_charged, string <-- "new", pivot).
//
// Revision 1.27  2000/04/14 12:40:42  jtanaka
// updated for new table "mdst_vee2".
//
// Revision 1.26  2000/04/13 12:41:58  katayama
// Added std:: to cout,cerr,endl etc.
//
// Revision 1.25  2000/03/07 11:14:04  katayama
// compatibility with CC5.0
//
// Revision 1.24  2000/01/11 12:28:23  jtanaka
// Modify constructor of Mdst_vee to make LAM, ALAM, GAMM.
//
// Revision 1.23  2000/01/05 06:09:39  jtanaka
// major updates: please see BELLE whiteboard.
//
// Revision 1.22  1999/04/13 20:57:39  jtanaka
// bug fix(isIdenticalWith)
//
// Revision 1.21  1999/04/09 15:02:10  jtanaka
// Added some member functions,"const", and "&" and removed some members
//
// Revision 1.20  1999/01/16 10:31:42  katayama
// clean up includes
//
// Revision 1.19  1999/01/09 12:51:37  jtanaka
// Bugs fix, added some members in Ptype Class, removed some members from Particle.cc because of functions for setting member of other classes.
//
// Revision 1.18  1998/12/14 08:10:05  katayama
// Because I removed static_Particle etc, I had to change these so that the
// library is consistent. There was another way to fix but I decided not to
// expose static_Particle etc. Authors, please check the changes I made.
//
// Revision 1.17  1998/12/12 23:27:04  katayama
// kludge for returning non-const particle
//
// Revision 1.16  1998/12/12 22:54:06  katayama
// prevent basf seg. fault???
//
// Revision 1.15  1998/10/26 08:41:18  jtanaka
// grandchildren -> finalStateParticles
//
// Revision 1.14  1998/10/15 09:52:31  jtanaka
// bug fix, and added deep_copy and deep_delete in Particle.h
//
// Revision 1.13  1998/10/12 19:24:09  jtanaka
// updated constructors Momentum.cc, Particle.cc, and Relation.cc. removed and added some functions in PID.cc. added (int*) to ParticleManager.cc(warning).
//
// Revision 1.12  1998/10/05 10:55:57  jtanaka
// bug fix in qqread.c and modified static objects in Relation.cc
//
// Revision 1.11  1998/09/27 18:38:29  jtanaka
// Relation::mc and mother modified(remove const -> no waring)
//
// Revision 1.10  1998/09/22 15:34:30  yiwasaki
// ParticleManager extended, Relation::mc modified
//
// Revision 1.9  1998/09/08 11:58:02  jtanaka
// added static objects in order to return the reference.
//
// Revision 1.8  1998/09/08 09:54:01  jtanaka
// Modify Mdst_charged constructor in Particle and some functions in Momentum.
//
// Revision 1.7  1998/09/07 13:26:12  jtanaka
// implement some fuctions in PID Class, and modify some function in Relation Class
//
// Revision 1.6  1998/07/22 13:04:16  jtanaka
// add some functions and modify const functions a little.
//
// Revision 1.5  1998/07/13 10:57:05  jtanaka
// modify a little about const functions and add vertex to Momentum
//
// Revision 1.4  1998/07/02 09:28:23  higuchit
// null flag -> `usable' flag
//
// Revision 1.3  1998/07/01 11:54:18  jtanaka
// add m_null.
//
// Revision 1.2  1998/06/19 09:09:05  jtanaka
// *** empty log message ***
//
// Revision 1.1  1998/06/15 07:33:43  yokoyamm
// First version.
//

#include "belle.h"
#include "particle/Relation.h"

#include "particle/Particle.h"

#include "panther/panther.h"

#include MDST_H
#include HEPEVT_H
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


static Mdst_charged *static_Mdst_charged(NULL);
static Gen_hepevt   *static_Gen_hepevt(NULL);
static Mdst_gamma   *static_Mdst_gamma(NULL);
static Mdst_trk     *static_Mdst_trk(NULL);
static Mdst_vee2    *static_Mdst_vee2(NULL);
static Mdst_klong   *static_Mdst_klong(NULL);
static Mdst_pi0     *static_Mdst_pi0(NULL);
static Mdst_ecl     *static_Mdst_ecl(NULL);
static Particle     *static_Particle(NULL);

/* static int m_n_static_new_pi0 = 1;
static int m_n_static_new_gamma = 1;
static int m_n_static_del_pi0 = 1;
static int m_n_static_del_gamma = 1;
static int m_n_static_del_gamma2 = 1; */

//Default constructor
Relation::Relation() 
  : m_flagChildModification(0),
    m_vee2ChildCounter(0), m_pi0ChildCounter(0)
{
  //...null
  m_mother  = NULL;
  m_mc      = NULL;
  m_charged = NULL;
  m_gamma   = NULL;
  m_trk     = NULL;
  m_vee2    = NULL;
  m_klong   = NULL;
  m_pi0     = NULL;
  m_ecl     = NULL;
  m_hep     = NULL;
}

//copy constructor
Relation::Relation(const Relation &a, Particle *s)
  : m_vee2ChildCounter(0), m_pi0ChildCounter(0)
{
  m_self     = s;
  //m_mother   = &((Particle&)a.mother());
  m_mother   = a.m_mother;
  //m_mc       = &((Particle&)a.mc());
  m_mc       = a.m_mc;
  m_children = a.children();
  //...m_finalStateParticles
  if(m_children.size() == 0){
    m_flagChildModification = 0;
  }else{
    m_flagChildModification = 1;  
  }
  m_charged  = a.m_charged;
  m_gamma    = a.m_gamma;
  m_trk      = a.m_trk;
  m_vee2     = a.m_vee2;
  m_klong    = a.m_klong;
  m_pi0      = a.m_pi0;
  m_ecl      = a.m_ecl;
  m_hep      = a.m_hep;

#if 1
  // Oct22,2001 by jtanaka
  if(m_vee2){
    for(unsigned int i=0;i<nChildren();++i){
      if(child(i).relation().m_vee2ChildCounter >= 1)
        ++(child(i).relation().m_vee2ChildCounter);
    }
  }
  if(m_pi0){
    //dout(Debugout::INFO,"Relation") << "call copy constructor of pi0 in Relation." << std::endl;
    for(unsigned int i=0;i<nChildren();++i){
      if(child(i).relation().m_pi0ChildCounter >= 1)
        ++(child(i).relation().m_pi0ChildCounter);
    }
  }
#endif
}

//Constructor with Mdst\_charged
Relation::Relation(const Mdst_charged &a, Particle *s)
  : m_vee2ChildCounter(0), m_pi0ChildCounter(0)
{
  m_self    = s;
  m_charged = &a;
  m_flagChildModification = 0;

  //...null
  m_mother  = NULL;
  m_mc      = NULL;
  m_gamma   = NULL;
  m_trk     = NULL;
  m_vee2    = NULL;
  m_klong   = NULL;
  m_pi0     = NULL;
  m_ecl     = NULL;
  m_hep     = NULL;
}

//Constructor with Mdst\_gamma
Relation::Relation(const Mdst_gamma &a, Particle *s)
  : m_vee2ChildCounter(0), m_pi0ChildCounter(0)
{
  m_self  = s;
  m_gamma = &a;
  m_flagChildModification = 0;

  //...null
  m_mother  = NULL;
  m_mc      = NULL;
  m_charged = NULL;
  m_trk     = NULL;
  m_vee2    = NULL;
  m_klong   = NULL;
  m_pi0     = NULL;
  m_ecl     = NULL;
  m_hep     = NULL;
#if 0
  // debug
  dout(Debugout::INFO,"Relation") << "create Mdst_gamma : " << m_n_static_new_gamma << std::endl;
  m_n_static_new_gamma++;
#endif
}

//recalculation of Vee daughter using more detail data.
static void
reCalVeeDaughter(Particle *dau, const Mdst_vee2 &vee)
{
  const HepPoint3D pivot(vee.vx(),vee.vy(),vee.vz());
  HepVector  a(5);
  HepSymMatrix Ea(5,0);
  if(dau->pType().charge() > 0.){
    a[0] = vee.daut().helix_p(0); a[1] = vee.daut().helix_p(1);
    a[2] = vee.daut().helix_p(2); a[3] = vee.daut().helix_p(3);
    a[4] = vee.daut().helix_p(4);
    Ea[0][0] = vee.daut().error_p(0);  Ea[1][0] = vee.daut().error_p(1);
    Ea[1][1] = vee.daut().error_p(2);  Ea[2][0] = vee.daut().error_p(3);
    Ea[2][1] = vee.daut().error_p(4);  Ea[2][2] = vee.daut().error_p(5);
    Ea[3][0] = vee.daut().error_p(6);  Ea[3][1] = vee.daut().error_p(7);
    Ea[3][2] = vee.daut().error_p(8);  Ea[3][3] = vee.daut().error_p(9);
    Ea[4][0] = vee.daut().error_p(10); Ea[4][1] = vee.daut().error_p(11);
    Ea[4][2] = vee.daut().error_p(12); Ea[4][3] = vee.daut().error_p(13);
    Ea[4][4] = vee.daut().error_p(14);
  }else{
    a[0] = vee.daut().helix_m(0); a[1] = vee.daut().helix_m(1);
    a[2] = vee.daut().helix_m(2); a[3] = vee.daut().helix_m(3);
    a[4] = vee.daut().helix_m(4);
    Ea[0][0] = vee.daut().error_m(0);  Ea[1][0] = vee.daut().error_m(1);
    Ea[1][1] = vee.daut().error_m(2);  Ea[2][0] = vee.daut().error_m(3);
    Ea[2][1] = vee.daut().error_m(4);  Ea[2][2] = vee.daut().error_m(5);
    Ea[3][0] = vee.daut().error_m(6);  Ea[3][1] = vee.daut().error_m(7);
    Ea[3][2] = vee.daut().error_m(8);  Ea[3][3] = vee.daut().error_m(9);
    Ea[4][0] = vee.daut().error_m(10); Ea[4][1] = vee.daut().error_m(11);
    Ea[4][2] = vee.daut().error_m(12); Ea[4][3] = vee.daut().error_m(13);
    Ea[4][4] = vee.daut().error_m(14);
  }
  Helix helix(pivot, a, Ea);
  HepPoint3D x;
  HepSymMatrix dpx(7,0);
  HepLorentzVector p = helix.momentum(0.,dau->pType().mass(),x,dpx);
  dau->momentum().momentumPosition(p,x,dpx);
}

//Constructor with Mdst\_vee2
Relation::Relation(const Mdst_vee2 &a, const bool makeRelation, Particle *s)
  : m_vee2ChildCounter(0), m_pi0ChildCounter(0)
{
  m_self = s;
  m_vee2 = &a;
  m_flagChildModification = 0;
  
  //...null
  m_mother  = NULL;
  m_mc      = NULL;
  m_charged = NULL;
  m_gamma   = NULL;
  m_trk     = NULL;
  m_klong   = NULL;
  m_pi0     = NULL;
  m_ecl     = NULL;
  m_hep     = NULL;

#if 1
  if(makeRelation){
    std::string pname("E+");
    std::string mname("E-");
    if(a.kind() == 1){
      pname = "PI+";
      mname = "PI-";
    }else if(a.kind() == 2){
      pname = "P+";
      mname = "PI-";
    }else if(a.kind() == 3){
      pname = "PI+";
      mname = "AP+";
    }
    //Ptype ptype_pion_plus("PI+");
    //Ptype ptype_pion_minus("PI-");
    Ptype ptype_pion_plus(pname.c_str());
    Ptype ptype_pion_minus(mname.c_str());
    int pionID[2] = {0, 0};
    unsigned nPion = 0;
    Mdst_charged_Manager &charged_mag = Mdst_charged_Manager::get_manager();
    for(std::vector<Mdst_charged>::iterator i = charged_mag.begin();
	i != charged_mag.end(); ++i){
      if(a.chgd(0).get_ID() >= 1 && pionID[0] == 0 &&
	 a.chgd(0).get_ID() == i->get_ID()){
	pionID[0] = (int)(i->get_ID()); //pi+
	++nPion;
    }
      if(a.chgd(1).get_ID() >= 1 && pionID[1] == 0 &&
	 a.chgd(1).get_ID() == i->get_ID()){
	pionID[1] = (int)(i->get_ID()); //pi-
	++nPion;
      }
      if(nPion == 2)break;
    }
    HepPoint3D dauPivot(a.vx(),a.vy(),a.vz());
    if(pionID[0] >= 1){
      Particle *tmp;
      if(a.daut()){
	tmp = new Particle(charged_mag[pionID[0]-1],ptype_pion_plus);
	reCalVeeDaughter(tmp,a);
      }else tmp = new Particle(charged_mag[pionID[0]-1],ptype_pion_plus,dauPivot);
      tmp->relation().m_vee2ChildCounter = 1;
      append(*tmp);
    }
    if(pionID[1] >= 1){
      Particle *tmp;
      if(a.daut()){
	tmp = new Particle(charged_mag[pionID[1]-1],ptype_pion_minus);
	reCalVeeDaughter(tmp,a);
      }else tmp = new Particle(charged_mag[pionID[1]-1],ptype_pion_minus,dauPivot);
      tmp->relation().m_vee2ChildCounter = 1;
      append(*tmp);
    }
  }
#endif
}

//Constructor with Mdst\_klong
Relation::Relation(const Mdst_klong &a, Particle *s)
  : m_vee2ChildCounter(0), m_pi0ChildCounter(0)
{
  m_self  = s;
  m_klong = &a;
  m_flagChildModification = 0;

  //...null
  m_mother  = NULL;
  m_mc      = NULL;
  m_charged = NULL;
  m_gamma   = NULL;
  m_trk     = NULL;
  m_vee2    = NULL;
  m_pi0     = NULL;
  m_ecl     = NULL;
  m_hep     = NULL;
}

//Constructor with Mdst\_pi0
Relation::Relation(const Mdst_pi0 &a, const bool makeRelation, Particle *s)
  : m_vee2ChildCounter(0), m_pi0ChildCounter(0)
{
  m_self = s;
  m_pi0  = &a;
  m_flagChildModification = 0;

  //...null
  m_mother  = NULL;
  m_mc      = NULL;
  m_charged = NULL;
  m_gamma   = NULL;
  m_trk     = NULL;
  m_vee2    = NULL;
  m_klong   = NULL;
  m_ecl     = NULL;
  m_hep     = NULL;

#if 1
  if(makeRelation){
    unsigned nGamma = 0;
    Mdst_gamma_Manager &gamma_mag = Mdst_gamma_Manager::get_manager();
    for(std::vector<Mdst_gamma>::iterator i = gamma_mag.begin();
	i != gamma_mag.end(); ++i){
      if(a.gamma(0).get_ID() >= 1 &&
	 a.gamma(0).get_ID() == i->get_ID()){
	Particle *tmp = new Particle(gamma_mag[(int)(a.gamma(0).get_ID())-1]);
	tmp->relation().m_pi0ChildCounter = 1;
	append(*tmp);
	++nGamma;
    }
      if(a.gamma(1).get_ID() >= 1 &&
	 a.gamma(1).get_ID() == i->get_ID()){
      Particle *tmp = new Particle(gamma_mag[(int)(a.gamma(1).get_ID())-1]);
      tmp->relation().m_pi0ChildCounter = 1;
      append(*tmp);
      ++nGamma;
    }
      if(nGamma == 2)break;
    }
  }
#endif
#if 0
  // debug
  dout(Debugout::INFO,"Relation") << "create Mdst_pi0 : " << m_n_static_new_pi0 << std::endl;
  m_n_static_new_pi0++;
#endif
}
//Constructor with Mdst\_ecl
Relation::Relation(const Mdst_ecl &a, Particle *s)
  : m_vee2ChildCounter(0), m_pi0ChildCounter(0)
{
  m_self = s;
  m_ecl  = &a;
  m_flagChildModification = 0;

  //...null
  m_mother  = NULL;
  m_mc      = NULL;
  m_charged = NULL;
  m_gamma   = NULL;
  m_trk     = NULL;
  m_vee2    = NULL;
  m_klong   = NULL;
  m_pi0     = NULL;
  m_hep     = NULL;
}

//Constructor with Gen\_hepevt
Relation::Relation(const Gen_hepevt &a, Particle *s)
  : m_vee2ChildCounter(0), m_pi0ChildCounter(0)
{
  m_self = s;
  m_hep  = &a;
  m_flagChildModification = 0;
  
  //...null
  m_mother  = NULL;
  m_mc      = NULL;
  m_charged = NULL;
  m_gamma   = NULL;
  m_trk     = NULL;
  m_vee2    = NULL;
  m_klong   = NULL;
  m_pi0     = NULL;
  m_ecl     = NULL;
}

// Destructor
Relation::~Relation()
{ 
#if 1
  // Oct22,2001 by jtanaka
  if(m_vee2){
    for(unsigned int i=0;i<nChildren();++i){
      if(child(i).relation().m_vee2ChildCounter >= 2){
        --(child(i).relation().m_vee2ChildCounter);
      }else if(child(i).relation().m_vee2ChildCounter == 1){
        delete &(child(i));
      }
    }
  }

  if(m_pi0){
    //dout(Debugout::INFO,"Relation") << "call destrector of pi0 in Relation." << std::endl;
    for(unsigned int i=0;i<nChildren();++i){
      if(child(i).relation().m_pi0ChildCounter >= 2){
	--(child(i).relation().m_pi0ChildCounter);
	//dout(Debugout::INFO,"Relation") << "on mem gamma : " << m_n_static_del_gamma2 << std::endl;
	//m_n_static_del_gamma2++;
      }else if(child(i).relation().m_pi0ChildCounter == 1){
	delete &(child(i));
	//dout(Debugout::INFO,"Relation") << "delete gamma : " << m_n_static_del_gamma << std::endl;
	//m_n_static_del_gamma++;
      }
    }
  }
#endif
}

// Interfaces for particles.
// returns a reference to mother.
const Particle & 
Relation::mother(void) const
{
  if(m_mother){
    return *m_mother;
  }else{
    if(static_Particle) return *static_Particle;
    else return *(static_Particle=new Particle);
  }
}

Particle & 
Relation::mother(void)
{
  if(m_mother){
    return *m_mother;
  }else{
    //return *(Particle *)NULL;
    if(static_Particle) return *static_Particle;
    else return *(static_Particle=new Particle);
  }
}
  
// sets a reference to mother.
const Particle & 
Relation::mother(Particle &a)
{
  return *(m_mother = &a);
}
  
// returns a number of children.
unsigned 
Relation::nChildren(void) const 
{
  return m_children.size();
}
  
// returns a const reference to i'th child.
const Particle & 
Relation::child(unsigned i) const 
{
  return *(m_children[i]);
}

// returns a reference to i'th child.
Particle & 
Relation::child(unsigned i)
{
  return *(m_children[i]);
}
  
// returns a list of children.
const std::vector<Particle *> &
Relation::children(void) const 
{
  return m_children;
}

// appends a child.
void 
Relation::append(Particle &a)
{
  Particle *tmp = &a;
  m_children.push_back(tmp);
  m_flagChildModification = 1;
}
  
// removes a child.
void 
Relation::remove(Particle &a)
{
  //Particle *tmp = &a;
  //m_children.erase(&tmp);
  for (std::vector<Particle*>::iterator i = m_children.begin();
       i != m_children.end(); ++i){
    if (&a == *i) {
      m_children.erase(i);
    }
  }
  m_flagChildModification = 1;
}

// removes all children and finalStateParticles.
void 
Relation::removeAll(void)
{
  m_children.erase(m_children.begin(),m_children.end());
  m_finalStateParticles.erase(m_finalStateParticles.begin(),m_finalStateParticles.end());
  m_flagChildModification = 0;
}

// returns a number of finalStateParticles.
unsigned 
Relation::nFinalStateParticles(void) const
{
  if(m_flagChildModification != 0){
    fillFinalStateParticles();
    m_flagChildModification = 0;
  }else{
    fillFinalStateParticles2();
  }
  return m_finalStateParticles.size();
}
  
// returns a reference to i'th finalStateParticle.
const Particle & 
Relation::finalStateParticle(unsigned i) const
{
  if(m_flagChildModification != 0){
    fillFinalStateParticles();
    m_flagChildModification = 0;
  }else{
    fillFinalStateParticles2();
  }
  return *(m_finalStateParticles[i]);
}

// returns a no const reference to i'th finalStateParticle.
Particle & 
Relation::noConstFinalStateParticle(unsigned i) const
{
  if(m_flagChildModification != 0){
    fillFinalStateParticles();
    m_flagChildModification = 0;
  }else{
    fillFinalStateParticles2();
  }
  return *(m_finalStateParticles[i]);
}
 
// returns a list of finalStateParticles.
const std::vector<Particle *> & 
Relation::finalStateParticles(void) const
{
  if(m_flagChildModification != 0){
    fillFinalStateParticles();
    m_flagChildModification = 0;
  }else{
    fillFinalStateParticles2();
  }
  return m_finalStateParticles;
}

// returns a reference to MC particle.
const Particle & 
Relation::mc(void) const
{
  if(m_mc){
    return *m_mc;
  }else{
    if(static_Particle) return *static_Particle;
    else return *(static_Particle=new Particle);
  }
}

Particle & 
Relation::mc(void)
{
  if(m_mc){
    return *m_mc;
  }else{
    //return *(Particle *)NULL;
    if(static_Particle) return *static_Particle;
    else return *(static_Particle=new Particle);
  }
}
  
// sets a reference to MC particle.
const Particle & 
Relation::mc(Particle &a)
{
  return *(m_mc = &a);
}
  
// Interfaces for MDST banks.
// returns a reference to Mdst\_charged.
const Mdst_charged & 
Relation::mdstCharged(void) const
{
  if(m_charged){
    return *m_charged;
  }else{
    if(static_Mdst_charged) return *static_Mdst_charged;
    else return *(static_Mdst_charged=new Mdst_charged);
  }
}
  
// sets a reference to Mdst\_charged.
const Mdst_charged & 
Relation::mdstCharged(const Mdst_charged &a)
{
  return *(m_charged = &a);
}
  
// returns a reference to Mdst\_gamma.
const Mdst_gamma & 
Relation::mdstGamma(void) const
{
  if(m_gamma){
    return *m_gamma;
  }else{
    if(static_Mdst_gamma) return *static_Mdst_gamma;
    else return *(static_Mdst_gamma=new Mdst_gamma);
  }
}

// sets a reference to Mdst\_gamma.
const Mdst_gamma & 
Relation::mdstGamma(const Mdst_gamma &a)
{
  return *(m_gamma = &a);
}

// returns a reference to Mdst\_trk.
const Mdst_trk & 
Relation::mdstTrk(void) const
{
  if(m_trk){
    return *m_trk;
  }else{
    if(static_Mdst_trk) return *static_Mdst_trk;
    else return *(static_Mdst_trk=new Mdst_trk);
  }
}
  
// sets a reference to Mdst\_trk.
const Mdst_trk & 
Relation::mdstTrk(const Mdst_trk &a)
{
  return *(m_trk = &a);
}
  
// returns a reference to Mdst\_vee2.
const Mdst_vee2 & 
Relation::mdstVee2(void) const
{
  if(m_vee2){
    return *m_vee2;
  }else{
    if(static_Mdst_vee2) return *static_Mdst_vee2;
    else return *(static_Mdst_vee2=new Mdst_vee2);
  }
}
  
// sets a reference to Mdst\_vee2.
const Mdst_vee2 & 
Relation::mdstVee2(const Mdst_vee2 &a)
{
 return *(m_vee2 = &a);
}

// returns a reference to Mdst\_klong.
const Mdst_klong & 
Relation::mdstKlong(void) const
{
  if(m_klong){
    return *m_klong;
  }else{
    if(static_Mdst_klong) return *static_Mdst_klong;
    else return *(static_Mdst_klong=new Mdst_klong);
  }
}

// sets a reference to Mdst\_klong.
const Mdst_klong & 
Relation::mdstKlong(const Mdst_klong &a)
{
  return *(m_klong = &a);
}

// returns a reference to Mdst\_pi0.
const Mdst_pi0 & 
Relation::mdstPi0(void) const
{
  if(m_pi0){
    return *m_pi0;
  }else{
    if(static_Mdst_pi0) return *static_Mdst_pi0;
    else return *(static_Mdst_pi0=new Mdst_pi0);
  }
}

// sets a reference to Mdst\_pi0.
const Mdst_pi0 & 
Relation::mdstPi0(const Mdst_pi0 &a)
{
  return *(m_pi0 = &a);
}


// returns a reference to Mdst\_ecl.
const Mdst_ecl & 
Relation::mdstEcl(void) const
{
  if(m_ecl){
    return *m_ecl;
  }else{
    if(static_Mdst_ecl) return *static_Mdst_ecl;
    else return *(static_Mdst_ecl=new Mdst_ecl);
  }
}

// sets a reference to Mdst\_ecl.
const Mdst_ecl & 
Relation::mdstEcl(const Mdst_ecl &a)
{
  return *(m_ecl = &a);
}

// returns a reference to Gen\_hepevt.
const Gen_hepevt & 
Relation::genHepevt(void) const
{
  if(m_hep){
    return *m_hep;
  }else{
    if(static_Gen_hepevt) return *static_Gen_hepevt;
    else return *(static_Gen_hepevt=new Gen_hepevt);
  }
}

void
Relation::resetGenHepevt(void)
{
  m_hep = NULL;
}

//void
//Relation::reset_genHepevt(void)
//{
//  //m_hep = NULL;
//  resetGenHepevt();
//}
  
// sets a reference to Gen\_hepevt.
const Gen_hepevt & 
Relation::genHepevt(const Gen_hepevt &a)
{
  return *(m_hep = &a);
}

bool 
Relation::isIdenticalWith(const Relation &x, const unsigned &type) const
{
  switch(type){
  case PC_ALL:
    if(m_charged && x.m_charged) 
      return (*m_charged == *(x.m_charged));

    if(m_gamma   && x.m_gamma) 
      return (*m_gamma   == *(x.m_gamma));
   
    if(m_vee2    && x.m_vee2) 
      return (*m_vee2    == *(x.m_vee2));

    if(m_pi0     && x.m_pi0) 
      return (*m_pi0     == *(x.m_pi0));

    if(m_klong   && x.m_klong) 
      return (*m_klong   == *(x.m_klong));

    if(m_trk     && x.m_trk) 
      return (*m_trk     == *(x.m_trk));

    if(m_ecl     && x.m_ecl) 
      return (*m_ecl     == *(x.m_ecl));

    if(!isReconstructedParticle() && !x.isReconstructedParticle() && m_hep && x.m_hep)     
      return (*m_hep     == *(x.m_hep)); 
    
    return false;
  case PC_CHARGED:
    if(m_charged && x.m_charged) return (*m_charged == *(x.m_charged));
    return false;
  case PC_GAMMA:
    if(m_gamma   && x.m_gamma)   return (*m_gamma   == *(x.m_gamma));
    return false;
  case PC_VEE2:
    if(m_vee2    && x.m_vee2)    return (*m_vee2    == *(x.m_vee2));
    return false;
  case PC_PI0:
    if(m_pi0     && x.m_pi0)     return (*m_pi0     == *(x.m_pi0));
    return false;
  case PC_ECL:
    if(m_ecl     && x.m_ecl)     return (*m_ecl     == *(x.m_ecl));
    return false;
  case PC_KLONG:
    if(m_klong   && x.m_klong)   return (*m_klong   == *(x.m_klong));
    return false;
  case PC_HEPEVT:
    if(m_hep     && x.m_hep)     return (*m_hep     == *(x.m_hep));
    return false;
  case PC_TRK:
    if(m_trk     && x.m_trk)     return (*m_trk     == *(x.m_trk));
    return false;
  default:
    return false;
  }
}

//Operators
// copy operator
Relation & 
Relation::operator = (const Relation &a)
{
  if(this == &a) return *this;

  m_self     = a.m_self;
  //m_mother   = &((Particle&)a.mother());
  m_mother   = a.m_mother;
  //m_mc       = &((Particle&)a.mc());
  m_mc       = a.m_mc;
  m_children = a.children();
  //...m_finalStateParticles
  if(m_children.size() == 0){
    m_flagChildModification = 0;
  }else{
    m_flagChildModification = 1;  
  }
  m_charged  = a.m_charged;
  m_gamma    = a.m_gamma;
  m_trk      = a.m_trk;
  m_vee2     = a.m_vee2;
  m_klong    = a.m_klong;
  m_pi0      = a.m_pi0;
  m_ecl      = a.m_ecl;
  m_hep      = a.m_hep;
  return *this;
} 

void 
Relation::fillFinalStateParticles(void) const
{
  unsigned size_children = m_children.size();
  if(size_children == 0)return;

  //...clears
  m_finalStateParticles.erase(m_finalStateParticles.begin(),m_finalStateParticles.end());
  //...refills
  for(unsigned i=0;i<size_children;i++){
    unsigned int size = m_children[i]->relation().nFinalStateParticles();
    for(unsigned j=0;j<size;j++){
	Particle *tmp = &(m_children[i]->relation().noConstFinalStateParticle(j));
	m_finalStateParticles.push_back(tmp);
    }
  }
}

void 
Relation::fillFinalStateParticles2(void) const
{
  if(m_children.size() != 0)return;
  if(m_charged == NULL && m_gamma == NULL && 
     m_klong   == NULL && m_ecl   == NULL)return;
  if(m_finalStateParticles.size() == 1)return;

  //...clears
  m_finalStateParticles.erase(m_finalStateParticles.begin(),m_finalStateParticles.end());
  //...refills
  m_finalStateParticles.push_back(m_self);
}

bool
Relation::isReconstructedParticle(void) const
{
  if(m_charged || m_gamma || m_trk || m_vee2 || m_klong || m_pi0 || m_ecl)
    return true;

  return false;
}
void 
Relation::dump(const std::string & keyword,
	       const std::string & prefix) const
{
  dout(Debugout::DUMP,"Relation") << prefix;
  m_self->pType().dump("nameonly");
  if(nChildren() == 0)dout(Debugout::DUMP,"Relation") << " --> not-decay" << std::endl;
  else{
    dout(Debugout::DUMP,"Relation") << " --> ";
    for(unsigned int i=0;i<nChildren();++i){
      child(i).pType().dump("nameonly");
      dout(Debugout::DUMP,"Relation") << "  ";
    }
    dout(Debugout::DUMP,"Relation") << std::endl;
  }

  if(keyword.find("recursive") != std::string::npos) {
    if(nChildren() != 0){
      for(unsigned int i=0;i<nChildren();++i)
	child(i).relation().dump(keyword, prefix+"       ");
    }
  }
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
