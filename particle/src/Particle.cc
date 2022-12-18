//
// Particle.cc
// 
// A class implimentation for <Particle> class 
// for the BELLE standard (hopefully!) Particle object liblary
// 
// <Particle> class supplies you interfaces to various particle information,
// such as momentum, pid, etc. through its private member objects,
// Momentum, PID, Ptype, Relation.
//
// Filename : Particle.cc
// Author : YOKOyama Masashi
// e-mail:yokoyama@belaxp1.phys.s.u-tokyo.ac.jp
// --> e-mail:jtanaka@belaxp1.phys.s.u-tokyo.ac.jp
// 
// $Id: Particle.cc 10079 2007-04-13 08:12:30Z katayama $ 
//
// Revision history
// 
// $Log$
// Revision 1.41  2003/09/15 13:24:34  katayama
// reduce warnings
//
// Revision 1.40  2002/04/29 06:11:38  jtanaka
// added removeAll() in Particle::deepDelete(void). (from Kakuno-san's info.)
//
// Revision 1.39  2001/12/14 02:16:32  katayama
// MDST_OBS removed
//
// Revision 1.38  2001/12/13 15:31:55  katayama
// MDST_OBS
//
// Revision 1.37  2001/10/22 02:21:08  jtanaka
// fix : memory leak when pi0,vee,vee2 particles are removed from "vector/list etc".
//
// Revision 1.36  2000/05/16 14:47:01  jtanaka
// added constructor of Mdst_ecl
//
// Revision 1.35  2000/04/27 12:57:38  jtanaka
// 1. use Mdst_vee_daughters information to construct Mdst_vee2.
// 2. add new constructor Particle(Mdst_charged, string <-- "new", pivot).
//
// Revision 1.34  2000/04/14 12:40:41  jtanaka
// updated for new table "mdst_vee2".
//
// Revision 1.33  2000/04/13 12:41:57  katayama
// Added std:: to cout,cerr,endl etc.
//
// Revision 1.32  2000/03/08 07:56:47  jtanaka
// modify a constructor of the Mdst_charged.
//
// Revision 1.31  2000/03/07 11:14:04  katayama
// compatibility with CC5.0
//
// Revision 1.30  2000/01/11 12:28:22  jtanaka
// Modify constructor of Mdst_vee to make LAM, ALAM, GAMM.
//
// Revision 1.29  2000/01/05 06:09:38  jtanaka
// major updates: please see BELLE whiteboard.
//
// Revision 1.28  1999/12/17 04:10:51  yiwasaki
// Particle::dump implemented
//
// Revision 1.27  1999/11/10 14:03:02  jtanaka
// bug fix: operator and copy constructor in Particle.cc, constructor in Momentum.cc, checkSame in combination.cc, (int) in utility.cc
//
// Revision 1.26  1999/11/08 15:24:45  higuchit
// Modified copy constructor (request from Nakadaira-san)
//
// Revision 1.25  1999/04/27 12:16:09  jtanaka
// added child() and mc().
//
// Revision 1.24  1999/04/09 15:02:10  jtanaka
// Added some member functions,"const", and "&" and removed some members
//
// Revision 1.23  1999/01/31 03:36:59  jtanaka
// for new xref in fundamentalFunctions.cc and bug fix in Particle.cc
//
// Revision 1.22  1999/01/16 10:31:38  katayama
// clean up includes
//
// Revision 1.21  1999/01/09 12:51:36  jtanaka
// Bugs fix, added some members in Ptype Class, removed some members from Particle.cc because of functions for setting member of other classes.
//
// Revision 1.20  1999/01/08 06:08:11  katayama
// Do not create PID objects as default
//
// Revision 1.19  1998/11/17 22:27:54  jtanaka
// id -> lund in Particle, and bug fix in PID.cc
//
// Revision 1.18  1998/11/10 06:50:58  jtanaka
// For new MDST table.
//
// Revision 1.17  1998/11/07 08:02:34  katayama
// Offset starts from zero, unlike mdst.tdf comment
//
// Revision 1.16  1998/11/07 07:54:39  katayama
// Modified to compile with the new mdst.
//
// Revision 1.15  1998/10/20 14:45:18  jtanaka
// added FittedMomentum.
//
// Revision 1.14  1998/10/15 09:52:31  jtanaka
// bug fix, and added deep_copy and deep_delete in Particle.h
//
// Revision 1.13  1998/10/12 19:24:08  jtanaka
// updated constructors Momentum.cc, Particle.cc, and Relation.cc. removed and added some functions in PID.cc. added (int*) to ParticleManager.cc(warning).
//
// Revision 1.12  1998/09/09 08:29:56  jtanaka
// added new constructor to Particle and modified some parts of Momentum.
//
// Revision 1.11  1998/09/08 11:58:02  jtanaka
// added static objects in order to return the reference.
//
// Revision 1.10  1998/09/08 09:54:01  jtanaka
// Modify Mdst_charged constructor in Particle and some functions in Momentum.
//
// Revision 1.9  1998/09/03 23:55:14  jtanaka
// We need to include "Particle.h" only. Comment out some members of Ptype.
//
// Revision 1.8  1998/07/22 13:04:15  jtanaka
// add some functions and modify const functions a little.
//
// Revision 1.7  1998/07/13 11:13:40  jtanaka
// modify a little about const functions and constructors
//
// Revision 1.6  1998/07/10 00:50:56  katayama
// Use .vect() for new CLHEP LorentzVec
//
// Revision 1.5  1998/07/03 14:57:15  jtanaka
// add some constructors for Ptype.
//
// Revision 1.4  1998/07/02 09:28:23  higuchit
// null flag -> `usable' flag
//
// Revision 1.3  1998/07/01 11:53:54  jtanaka
// add m_null.
//
// Revision 1.2  1998/06/19 09:08:56  jtanaka
// *** empty log message ***
//
// Revision 1.1  1998/06/15 07:33:17  yokoyamm
// First version.
//

#include "belle.h"
#include "particle/Particle.h"
#include "particle/PID.h"
#include "particle/ParticleUserInfo.h"

#include "panther/panther.h"
#include MDST_H
#include HEPEVT_H
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


//#include <typeinfo> // for DEBUG

static PID              *static_PID(NULL);
static ParticleUserInfo *static_ParticleUserInfo(NULL);

// constructors
Particle::Particle() 
  : /* m_usable(USABLE) */ m_usable(UNUSABLE), m_pId(NULL), m_userInfo(NULL)
{
  // m_usable must be UNUSABLE because of bool operator.
  m_momentum = new Momentum();
  // due to it's size, do not create now (nk)
  //  m_pId      = new PID();
  m_relation = new Relation();
  m_relation->m_self = this;
  m_pType    = new Ptype();
  //m_fittedMomentum = new FittedMomentum();
}

Particle::Particle(const Particle & a)
  : m_pId(NULL), m_userInfo(NULL)
{
  m_momentum = new Momentum(a.momentum());
  // due to it's size, do not create now (nk)
  if (a.m_pId) m_pId = new PID(a.pId());
  m_relation = new Relation(a.relation(), this);
  m_pType    = new Ptype(a.pType());
  //m_fittedMomentum = new FittedMomentum(a.fittedMomentum());
  m_name   = a.m_name;
  m_usable = a.m_usable;

  if (a.m_userInfo) {
    // dout(Debugout::DDEBUG,"Particle") << "Particle1: " << typeid(*(a.m_userInfo)).name() << std::endl;
    m_userInfo = a.userInfo().clone();
    // dout(Debugout::DDEBUG,"Particle") << "Particle2: " << typeid(*m_userInfo).name() << std::endl;
    // dout(Debugout::DDEBUG,"Particle") << "Particle3: " << typeid(a.userInfo()).name() << std::endl;
  }
  
}

Particle::Particle(const HepLorentzVector & a, const Ptype &ptype)
  : m_usable(USABLE), m_pId(NULL), m_userInfo(NULL)
{
  m_momentum = new Momentum(a);
  // due to it's size, do not create now (nk)
  //  m_pId      = new PID();
  m_relation = new Relation();
  m_relation->m_self = this;
  m_pType    = new Ptype(ptype);
  //m_fittedMomentum = new FittedMomentum();
}

Particle::Particle(const Momentum & a, const Ptype &ptype)
  : m_usable(USABLE), m_pId(NULL), m_userInfo(NULL)
{
  m_momentum = new Momentum(a);
  // due to it's size, do not create now (nk)
  //  m_pId      = new PID();
  m_relation = new Relation();
  m_relation->m_self = this;
  m_pType    = new Ptype(ptype);
  //m_fittedMomentum = new FittedMomentum();
}

Particle::Particle(const Mdst_charged & a, const Ptype &ptype,
		   const HepPoint3D &pivot)
  : m_usable(USABLE), m_pId(NULL), m_userInfo(NULL)
{
  //m_momentum = new Momentum(a, ptype.mass());
  m_momentum = new Momentum(a, ptype, pivot);
  // due to it's size, do not create now (nk)
  //  m_pId      = new PID(a);
  m_relation = new Relation(a, this);
  m_pType    = new Ptype(ptype);
  //m_fittedMomentum = new FittedMomentum();
}

Particle::Particle(const Mdst_charged & a, const std::string &pTypeName,
		   const HepPoint3D &pivot)
  : m_usable(USABLE), m_pId(NULL), m_userInfo(NULL)
{
  m_pType    = new Ptype(pTypeName.c_str());
  m_momentum = new Momentum(a, *m_pType, pivot);
  m_relation = new Relation(a, this);
  //m_fittedMomentum = new FittedMomentum();
}

Particle::Particle(const Mdst_gamma & a) 
  : m_usable(USABLE), m_pId(NULL), m_userInfo(NULL)
{
  m_momentum = new Momentum(a);
  // due to it's size, do not create now (nk)
  //  m_pId      = new PID(a);
  m_relation = new Relation(a, this);
  m_pType    = new Ptype("GAMM");
  //m_fittedMomentum = new FittedMomentum();
}


Particle::Particle(const Mdst_vee2 & a, const bool makeRelation) 
  : m_usable(USABLE), m_pId(NULL), m_userInfo(NULL)
{
  m_momentum = new Momentum(a);
  // due to it's size, do not create now (nk)
  //  m_pId      = new PID(a);
  m_relation = new Relation(a, makeRelation, this);
  std::string pname("GAMM");
  if(a.kind() == 1)pname = "K0S";
  else if(a.kind() == 2)pname = "LAM";
  else if(a.kind() == 3)pname = "ALAM";
  m_pType = new Ptype(pname.c_str());
  //m_fittedMomentum = new FittedMomentum();
}

Particle::Particle(const Mdst_klong & a)
 : m_usable(USABLE), m_pId(NULL), m_userInfo(NULL)
{
  m_momentum = new Momentum(a);
  // due to it's size, do not create now (nk)
  //  m_pId      = new PID(a);
  m_relation = new Relation(a, this);
  m_pType    = new Ptype("K0L");
  //m_fittedMomentum = new FittedMomentum();
}

Particle::Particle(const Mdst_pi0 & a, const bool makeRelation)
 : m_usable(USABLE), m_pId(NULL), m_userInfo(NULL)
{
  m_momentum = new Momentum(a);
  // due to it's size, do not create now (nk)
  //  m_pId      = new PID(a);
  m_relation = new Relation(a, makeRelation, this);
  m_pType    = new Ptype("PI0");
  //m_fittedMomentum = new FittedMomentum();
}
Particle::Particle(const Mdst_ecl & a) 
  : m_usable(USABLE), m_pId(NULL), m_userInfo(NULL)
{
  m_momentum = new Momentum(a);
  m_relation = new Relation(a, this);
  m_pType    = new Ptype("GAMM");
}

Particle::Particle(const Gen_hepevt & a, const Ptype &ptype)
  : m_usable(USABLE), m_pId(NULL), m_userInfo(NULL)
{
  m_momentum = new Momentum(a);
  // due to it's size, do not create now (nk)
  //  m_pId      = new PID(a);
  m_relation = new Relation(a, this);
  m_pType    = new Ptype(ptype);
  //m_fittedMomentum = new FittedMomentum();
}

Particle::Particle(const Gen_hepevt & a)
  : m_usable(USABLE), m_pId(NULL), m_userInfo(NULL)
{
  m_momentum = new Momentum(a);
  // due to it's size, do not create now (nk)
  //  m_pId      = new PID(a);
  m_relation = new Relation(a, this);
  m_pType    = new Ptype(a.idhep());
  //m_fittedMomentum = new FittedMomentum();
}

// destructor
Particle::~Particle() 
{
  delete m_momentum;
  if (m_pId) delete m_pId;
  delete m_relation;
  delete m_pType;
  //delete m_fittedMomentum;
  if (m_userInfo) delete m_userInfo;
}

// create PID object
void 
Particle::createPID(void)
{
  if(isUsable()) {
    //if (m_pId) return;
    if(mdstCharged()) {
      //m_pId = new PID(mdstCharged());
      m_pId = new PID(*this);
      return;
    }
    return;
  }
}

// append daughter
Particle &Particle::append_daughter(Particle &d) {
  relation().append(d);
  if (d.m_userInfo) {
    if (! m_userInfo) {
      userInfo(d.userInfo());
    } else {
      userInfo().append(d.userInfo());
    }
  }
  return *this;
}


// returns object name
const std::string &
Particle::name(void) const 
{
  return m_name;
}

// sets object name
const std::string &
Particle::name(const std::string & a) 
{
  return m_name = a;
}


// Momentum interfaces

// returns a reference to Momentum
Momentum &
Particle::momentum(void)
{
  return *m_momentum;
}

const Momentum &
Particle::momentum(void) const 
{
  return *m_momentum;
}

// sets a reference to Momentum
const Momentum &
Particle::momentum(const Momentum & a)
{
  return *m_momentum = a;
}

// returns 4-momentum
const HepLorentzVector &
Particle::p(void) const 
{
  return m_momentum->p();
}

const Hep3Vector
Particle::p3(void) const
{
  return m_momentum->p().vect();
}

// returns position vector
const Hep3Vector &
Particle::x(void) const
{
  return m_momentum->x();
}

// returns magnitude of momentum
double 
Particle::ptot(void) const 
{
  return (m_momentum->p()).vect().mag();
}

// returns x component of momentum
double 
Particle::px(void) const 
{
  return (m_momentum->p()).x();
}

// returns y component of momentum
double 
Particle::py(void) const 
{
  return (m_momentum->p()).y();
}

// returns z component of momentum
double 
Particle::pz(void) const 
{
  return (m_momentum->p()).z();
}

// returns energy
double 
Particle::e(void) const 
{
  return (m_momentum->p()).t();
}

// returns invariant mass
double 
Particle::mass(void) const 
{
  return m_momentum->mass();
}

// PID interfaces

// returns a reference to PID
const PID &
Particle::pId(void) const 
{
  if(m_pId)return *m_pId;
  else{
    if(static_PID)return *static_PID;
    else return *(static_PID = new PID);
  }
}

PID &
Particle::pId(void) 
{
  if(m_pId)return *m_pId;
  else{
    if(static_PID)return *static_PID;
    else return *(static_PID = new PID);
  }
}

// sets a reference to PID
const PID &
Particle::pId(const PID &a) 
{
  if(m_pId){
    return *m_pId = a;
  }else{
    m_pId = new PID(a);
    return *m_pId;
  }
}

// Ptype interfaces

/// returns a reference to Ptype.
const Ptype &
Particle::pType(void) const 
{
  return *m_pType;
}

Ptype &
Particle::pType(void) 
{
  return *m_pType;
}

/// sets a reference to Ptype.
const Ptype & 
Particle::pType(const Ptype &a) 
{
  return *m_pType = a;
}

/// returns charge. (in units of e)
double 
Particle::charge(void) const 
{
  return m_pType->charge();
}

/// returns LUND7 particle code.
int 
Particle::lund(void) const 
{
  return m_pType->lund();
}
  
// Relation interfaces

/// returns a reference to Relation.
const Relation & 
Particle::relation(void) const 
{
  return *m_relation;
}

Relation & 
Particle::relation(void) 
{
  return *m_relation;
}
  
/// sets a reference to Relation.
const Relation & 
Particle::relation(const Relation &a) 
{
  return *m_relation = a;
}

/// returns a reference to mother.
const Particle & 
Particle::mother(void) const 
{
  return m_relation->mother();
}
  
/// returns a number of children.
unsigned 
Particle::nChildren(void) const 
{
  return m_relation->nChildren();
}

/// returns a reference to i'th child.
const Particle & 
Particle::child(unsigned i) const
{
    return m_relation->child(i);
}

/// returns a reference to i'th child.
Particle & 
Particle::child(unsigned i)
{
    return m_relation->child(i);
}
  
/// returns a reference to MC particle.
const Particle &  
Particle::mc(void) const 
{
  return m_relation->mc();
}

/// returns a reference to MC particle.
Particle &  
Particle::mc(void)
{
  return m_relation->mc();
}

/// returns a reference to Mdst\_charged.
const Mdst_charged & 
Particle::mdstCharged(void) const 
{
  return m_relation->mdstCharged();
}
  
/// returns a reference to Mdst\_gamma.
const Mdst_gamma & 
Particle::mdstGamma(void) const 
{
  return m_relation->mdstGamma();
}
  
/// returns a reference to Mdst\_trk.
const Mdst_trk & 
Particle::mdstTrk(void) const 
{
  return m_relation->mdstTrk();
}


/// returns a reference to Mdst\_vee2.
const Mdst_vee2 & 
Particle::mdstVee2(void) const 
{
  return m_relation->mdstVee2();
}

/// returns a reference to Mdst\_klong.
const Mdst_klong & 
Particle::mdstKlong(void) const 
{
  return m_relation->mdstKlong();
}

/// returns a reference to Mdst\_pi0.
const Mdst_pi0 & 
Particle::mdstPi0(void) const 
{
  return m_relation->mdstPi0();
}

/// returns a reference to Mdst\_ecl.
const Mdst_ecl & 
Particle::mdstEcl(void) const 
{
  return m_relation->mdstEcl();
}
 
/// returns a reference to Gen\_hepevt.
const Gen_hepevt & 
Particle::genHepevt(void) const 
{
  return m_relation->genHepevt();
}

// Operators
/// Copy operator
Particle & 
Particle::operator=(const Particle &a)
{
  if(this == &a) return *this;
  
  delete m_momentum;
  m_momentum = new Momentum(a.momentum());

  if(m_pId) delete m_pId;
  if(a.m_pId){
    m_pId = new PID(a.pId());
  }else{
    m_pId = NULL;
  }
  
  delete m_relation;
  m_relation = new Relation(a.relation(), this);

  delete m_pType;
  m_pType    = new Ptype(a.pType());

  //delete m_fittedMomentum;
  //m_fittedMomentum = new FittedMomentum(a.fittedMomentum());

  m_name   = a.m_name;
  m_usable = a.m_usable;

  if(m_userInfo) delete m_userInfo;
  if (a.m_userInfo) m_userInfo = a.userInfo().clone();
  else m_userInfo = NULL;
  return *this;
}

bool 
Particle::usable(const bool &n)
{
  switch (n){
  case USABLE:
    m_usable = USABLE;
    break;
  case UNUSABLE:
    m_usable = UNUSABLE;
    break;
  }
  return m_usable;
}

//bool 
//Particle::unusable(bool n)
//{
//  switch (n){
//  case USABLE:
//    m_usable = USABLE;
//    break;
//  case UNUSABLE:
//    m_usable = UNUSABLE;
//    break;
//  }
//  return m_usable;
//}  

Particle
Particle::deepCopy(void)
{
  Particle copied(*this);

  //...m_children
  std::vector<Particle *> children;
  unsigned size  = copied.relation().nChildren();
  //std::vector<unsigned> fitted;
  for(unsigned i=0;i<size;++i){
    children.push_back(new Particle(copied.relation().child(i).deepCopy()));
    //...for fitted momentum class
    //fitted.push_back(copied.relation().child(i).fittedMomentum().nCoMatrix());
    //for(unsigned j=0;j<size;++j){
    //if(copied.relation().child(i).fittedMomentum().nCoMatrix() == 0)break;
    //unsigned error;
    //HepMatrix tmp = 
    //  copied.relation().child(i).fittedMomentum().coMatrix(copied.relation().child(j), &error);
    //if(&error != NULL){
    //  copied.relation().child(i).fittedMomentum().coMatrix(children[i], tmp);
    //}
    //}
  }
  //...for fitted momentum class
  //for(unsigned i=0;i<size;++i){
  //if(fitted[i] == 0)break;
  //for(unsigned j=0;j<fitted[i];++j){
  //copied.relation().child(i).fittedMomentum().removeCoMatrix(this->relation().child(i).fittedMomentum().coParticle(j));
  //}
  //}
  // It is necessary to decrease "index" of original "Particle".
  if(mdstVee2()){
    for(unsigned int i=0;i<m_relation->nChildren();++i){
      --(m_relation->child(i).relation().m_vee2ChildCounter);
    }
  }
  if(mdstPi0()){
    for(unsigned int i=0;i<m_relation->nChildren();++i){
      --(m_relation->child(i).relation().m_pi0ChildCounter);
    }
  }

  copied.relation().removeAll();
  for(unsigned i=0;i<size;++i){
    copied.relation().append(*children[i]);
  }

  return copied;
}

void 
Particle::deepDelete(void)
{
  if(this->relation().nChildren() == 0)return;
  for(unsigned i=0;i<this->relation().nChildren();++i){
    this->relation().child(i).deepDelete();
    delete &this->relation().child(i);
  }
  this->relation().removeAll();
}

const ParticleUserInfo &
Particle::userInfo(void) const
{
  if(m_userInfo)return *m_userInfo;
  else{
//      if(static_ParticleUserInfo)return *static_ParticleUserInfo;
//      else return *(static_ParticleUserInfo = new ParticleUserInfo);
    return *static_ParticleUserInfo;
  }
}

ParticleUserInfo &
Particle::userInfo(void)
{
  if(m_userInfo)return *m_userInfo;
  else{
//      if(static_ParticleUserInfo)return *static_ParticleUserInfo;
//      else return *(static_ParticleUserInfo = new ParticleUserInfo);
    return *static_ParticleUserInfo;
  }
}

const ParticleUserInfo &
Particle::userInfo(const ParticleUserInfo &info)
{
  if(m_userInfo){
    *m_userInfo = info;
  }else{
    // dout(Debugout::DDEBUG,"Particle") << "ParticleU1: " << typeid(info).name() << std::endl;
    m_userInfo = info.clone();
    // dout(Debugout::DDEBUG,"Particle") << "ParticleU2: " << typeid(*m_userInfo).name() << std::endl;
  }
  return *m_userInfo;
}

void
Particle::dump(const std::string & keyword, const std::string & prefix) const {
    bool full = false;
    if (keyword.find("full") != std::string::npos) full = true;

    dout(Debugout::DUMP,"Particle") << prefix;
    if(m_name == std::string("") && pType())pType().dump("nameonly");
    else dout(Debugout::DUMP,"Particle") << m_name;
    if (full || keyword.find("mass") != std::string::npos)     dout(Debugout::DUMP,"Particle") << " m=" << mass();
    if (full || keyword.find("momentum") != std::string::npos) dout(Debugout::DUMP,"Particle") << " p=" << p();
    if (full || keyword.find("position") != std::string::npos) dout(Debugout::DUMP,"Particle") << " x=" << x();
    dout(Debugout::DUMP,"Particle") << std::endl;

    if (full || keyword.find("recursive") != std::string::npos) {
        unsigned n = nChildren();
        if (n) {
            for (unsigned i = 0; i < n; i++)
                child(i).dump(keyword, prefix + "    ");
        }
    }
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
