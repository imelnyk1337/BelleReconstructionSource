// Relation.h
// 
// A class definition for <Relation> class 
// for the BELLE standard (hopefully!) Particle object liblary
// 
// <Relation> class supplies you interfaces to other objects,
// such as mother particle, MC particle, and MDST banks.
//
// Filename : Particle.h
// Author : YOKOyama Masashi
// e-mail:yokoyama@belaxp1.phys.s.u-tokyo.ac.jp
// --> e-mail:jtanaka@belaxp1.phys.s.u-tokyo.ac.jp
// 
// $Id: Relation.h 9932 2006-11-12 14:26:53Z katayama $ 
//
// Revision history
// 
// $Log$
// Revision 1.28  2001/12/13 15:31:54  katayama
// MDST_OBS
//
// Revision 1.27  2000/05/29 11:32:19  jtanaka
// add const and mutable to some member functions and parameters.
//
// Revision 1.26  2000/05/16 14:46:28  jtanaka
// added constructor of Mdst_ecl.
//
// Revision 1.25  2000/04/14 12:40:57  jtanaka
// updated for new table "mdst_vee2".
//
// Revision 1.24  2000/03/07 11:14:01  katayama
// compatibility with CC5.0
//
// Revision 1.23  2000/01/23 13:01:54  katayama
// fixed typo found with gcc 2.95.2
//
// Revision 1.22  2000/01/05 06:10:39  jtanaka
// major updates:please see BELLE whiteboard.
//
// Revision 1.21  1999/04/09 15:01:09  jtanaka
// Added new member functions, "const", "&"...etc. Removed some member functions.
//
// Revision 1.20  1999/01/16 10:31:32  katayama
// clean up includes
//
// Revision 1.19  1999/01/09 12:51:00  jtanaka
// Bugs fix, added some members in Ptype Class, removed some member from Particle.h because of functions for setting member of other classes.
//
// Revision 1.18  1998/12/14 08:09:55  katayama
// Because I removed static_Particle etc, I had to change these so that the
// library is consistent. There was another way to fix but I decided not to
// expose static_Particle etc. Authors, please check the changes I made.
//
// Revision 1.17  1998/10/26 08:41:28  jtanaka
// grandchildren -> finalStateParticles
//
// Revision 1.16  1998/10/15 09:52:14  jtanaka
// bug fix, and added deep_copy and deep_delete in Particle.h
//
// Revision 1.15  1998/10/12 19:25:45  jtanaka
// updated constructors Momentum.h, Particle.h, Relation.h. removed and added some functions in PID.h. added some static_PIDs to static_particle.h.
//
// Revision 1.14  1998/10/05 10:59:30  jtanaka
// added some comments for doc++ in some headers files and modified static objects in Relation.h, ParticleManager.h.
//
// Revision 1.13  1998/09/27 18:38:49  jtanaka
// Relation::mc and mother modified(remove const -> no warning)
//
// Revision 1.12  1998/09/08 11:58:13  jtanaka
// added static objects in order to return the reference.
//
// Revision 1.11  1998/09/07 13:25:49  jtanaka
// implement some fuctions in PID Class, and modify some function in Relation Class
//
// Revision 1.10  1998/09/03 23:55:36  jtanaka
// We need to include "Particle.h" only. Comment out some members of Ptype.
//
// Revision 1.9  1998/07/23 12:51:43  katayama
// conform ANSI
//
// Revision 1.8  1998/07/22 13:04:39  jtanaka
// add some functions and modify const functions a little.
//
// Revision 1.7  1998/07/13 10:55:16  jtanaka
// modify a little about const functions and add vertex to Momentum
//
// Revision 1.6  1998/07/03 12:10:14  jtanaka
// modified some parts
//
// Revision 1.5  1998/07/02 09:29:10  higuchit
// null flag -> `usable' flag
//
// Revision 1.4  1998/07/01 11:55:16  jtanaka
// add m_null.
//
// Revision 1.3  1998/06/19 09:09:18  jtanaka
// *** empty log message ***
//
// Revision 1.2  1998/06/15 07:29:45  yokoyamm
// First version..
//

#ifndef PARTICLE_CLASS_RELATION_H
#define PARTICLE_CLASS_RELATION_H

#include "belle.h"
#include <string>
#include <vector>

#include "particle/Particle.h"
#include "particle/constant.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//class Particle;
class Mdst_charged;
class Mdst_gamma;
class Mdst_trk;
class Mdst_vee2;
class Mdst_klong;
class Mdst_pi0;
class Mdst_ecl;
class Gen_hepevt;

/// Relation class supplies you interfaces to other objects, such as mother particle, MC particle, and MDST banks.
class Relation {

public:
  /// Default constructor
  Relation();

  /// Copy constructor
  Relation(const Relation &, Particle * = NULL);
  
  /// Constructor with Mdst\_charged
  Relation(const Mdst_charged &, Particle * = NULL);
  /// Constructor with Mdst\_gamma
  Relation(const Mdst_gamma &, Particle * = NULL);
  /// Constructor with Mdst\_vee2
  Relation(const Mdst_vee2 &, const bool makeRelation = true, Particle * = NULL);
  /// Constructor with Mdst\_klong
  Relation(const Mdst_klong &, Particle * = NULL);
  /// Constructor with Mdst\_pi0
  Relation(const Mdst_pi0 &, const bool makeRelation = true, Particle * = NULL);
  /// Constructor with Mdst\_ecl
  Relation(const Mdst_ecl &, Particle * = NULL);
  /// Constructor with Gen\_hepevt
  Relation(const Gen_hepevt &, Particle * = NULL);

  /// Destructor
  virtual ~Relation();

public:// General interfaces
  /// returns class name.
  virtual std::string className(void){ return std::string("Relation"); }
  
  /// dumps debug information. (not implement.)
  virtual void dump(const std::string & keyword = std::string(""),
		    const std::string & prefix  = std::string("")) const;
  
public:// Interfaces for particles.
  /// returns a const reference to mother.
  virtual const Particle & mother(void) const;

  /// returns a reference to mother.
  virtual Particle & mother(void);
  
  /// sets a reference to mother and returns it.
  virtual const Particle & mother(Particle &);
  
  /// returns a number of children.
  virtual unsigned int nChildren(void) const;
  
  /// returns a const reference to i'th child.
  virtual const Particle & child(unsigned int) const;

  /// returns a reference to i'th child.
  virtual Particle & child(unsigned int);
  
  /// returns a list of children.
  virtual const std::vector<Particle *> & children(void) const;
  
  /// append a child.
  virtual void append(Particle &);
  
  /// remove a child.
  virtual void remove(Particle &);

  /// remove all children and finalStateParticles
  virtual void removeAll(void);
  
  /// returns a number of finalStateParticles.
  virtual unsigned int nFinalStateParticles(void) const;
  
  /// returns a const reference to i'th finalStateParticle.
  virtual const Particle & finalStateParticle(unsigned int) const;

  /// returns a list of finalStateParticles.
  virtual const std::vector<Particle *> & finalStateParticles(void) const;

  /// returns a const reference to MC particle.
  virtual const Particle & mc(void) const;

  /// returns a reference to MC particle.
  virtual Particle & mc(void);
  
  /// sets a reference to MC particle and returns it.
  virtual const Particle & mc(Particle &);
  
public:// Interfaces for MDST banks.
  /// returns a reference to Mdst\_charged.
  virtual const Mdst_charged & mdstCharged(void) const;
  
  /// sets a reference to Mdst\_charged and returns it.
  virtual const Mdst_charged & mdstCharged(const Mdst_charged &);
  
  /// returns a reference to Mdst\_gamma.
  virtual const Mdst_gamma & mdstGamma(void) const;
  
  /// sets a reference to Mdst\_gamma and returns it.
  virtual const Mdst_gamma & mdstGamma(const Mdst_gamma &);
  
  /// returns a reference to Mdst\_trk.
  virtual const Mdst_trk & mdstTrk(void) const;
  
  /// sets a reference to Mdst\_trk and returns it.
  virtual const Mdst_trk & mdstTrk(const Mdst_trk &);
  
  /// returns a reference to Mdst\_vee2.
  virtual const Mdst_vee2 & mdstVee2(void) const;
  
  /// sets a reference to Mdst\_vee2 and returns it.
  virtual const Mdst_vee2 & mdstVee2(const Mdst_vee2 &);

  /// returns a reference to Mdst\_klong.
  virtual const Mdst_klong & mdstKlong(void) const;

  /// sets a reference to Mdst\_klong and returns it.
  virtual const Mdst_klong & mdstKlong(const Mdst_klong &);

  /// returns a reference to Mdst\_pi0.
  virtual const Mdst_pi0 & mdstPi0(void) const;

  /// sets a reference to Mdst\_pi0 and returns it.
  virtual const Mdst_pi0 & mdstPi0(const Mdst_pi0 &);

  /// returns a reference to Mdst\_ecl.
  virtual const Mdst_ecl & mdstEcl(void) const;

  /// sets a reference to Mdst\_ecl and returns it.
  virtual const Mdst_ecl & mdstEcl(const Mdst_ecl &);
  
  /// returns a reference to Gen\_hepevt.
  virtual const Gen_hepevt & genHepevt(void) const;
  
  /// sets a reference to Gen\_hepevt and returns it.
  virtual const Gen_hepevt & genHepevt(const Gen_hepevt &);

  /// resets a reference to Gen\_hepevt.
  virtual void resetGenHepevt(void);
  //virtual void reset_genHepevt(void);

  /// identifies to particles.
  virtual bool isIdenticalWith(const Relation &, const unsigned &type = PC_ALL) const;

  /// true if Particle is constructed from mdst_charged/gamma/trk_vee2/klong/pi0/ecl
  bool isReconstructedParticle(void) const;

public://operators
  /// copy operator.
  Relation & operator = (const Relation &);

private://Calculation
  /// calculates m_finalStateParticles.
  virtual void fillFinalStateParticles(void) const;

  /// calculates m_finalStateParticles(self).
  virtual void fillFinalStateParticles2(void) const;

  /// returns a no const reference to i'th finalStateParticle.
  virtual Particle & noConstFinalStateParticle(unsigned int) const;

private:// Protected members
  Particle * m_self;
  Particle * m_mother;
  std::vector<Particle *> m_children;
  mutable std::vector<Particle *> m_finalStateParticles;
  Particle * m_mc;
  const Mdst_charged * m_charged;
  const Mdst_gamma * m_gamma;
  const Mdst_trk * m_trk;
  const Mdst_vee2 * m_vee2;
  const Mdst_klong * m_klong;
  const Mdst_pi0 * m_pi0;
  const Mdst_ecl * m_ecl;
  const Gen_hepevt * m_hep;

  mutable unsigned int m_flagChildModification;

  int m_vee2ChildCounter;
  int m_pi0ChildCounter;
  
  friend Particle::Particle(const Particle & a);
  friend Particle::Particle();
  friend Particle::Particle(const Momentum & a, const Ptype &ptype);
  friend Particle::Particle(const HepLorentzVector & a, const Ptype &ptype);
  friend Particle & Particle::operator = (const Particle &);
  friend Particle::~Particle();
  friend Particle Particle::deepCopy(void);
};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* PARTICLE_CLASS_RELATION_H */
