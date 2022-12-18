// -*- C++ -*-
//
// Package:     <package>
// Module:      Ptype
// 
// Description: <one line class summary>
//
// Usage:
//    <usage>
//
// Author:      KATAYAMA Nobuhiko
// Created:     Mon Jan 26 20:09:42 JST 1998
// $Id: Ptype.h 10564 2008-07-30 21:21:59Z katayama $
//
// Revision history
//
// $Log$
// Revision 1.19  2002/01/03 11:04:39  katayama
// Point3D and other header files are cleaned
//
// Revision 1.18  2000/03/07 11:14:01  katayama
// compatibility with CC5.0
//
// Revision 1.17  2000/01/05 06:10:39  jtanaka
// major updates:please see BELLE whiteboard.
//
// Revision 1.16  1999/11/10 14:04:11  jtanaka
// bug fix: explicit in Ptype.h, checkSame in combination.h, & in utility.h
//
// Revision 1.15  1999/04/09 15:01:09  jtanaka
// Added new member functions, "const", "&"...etc. Removed some member functions.
//
// Revision 1.14  1999/01/16 10:31:31  katayama
// clean up includes
//
// Revision 1.13  1999/01/09 12:51:00  jtanaka
// Bugs fix, added some members in Ptype Class, removed some member from Particle.h because of functions for setting member of other classes.
//
// Revision 1.12  1998/11/17 22:28:15  jtanaka
// id -> lund
//
// Revision 1.11  1998/10/14 12:50:40  jtanaka
// uses qq2panther to get Ptype information.
//
// Revision 1.10  1998/10/05 10:59:30  jtanaka
// added some comments for doc++ in some headers files and modified static objects in Relation.h, ParticleManager.h.
//
// Revision 1.9  1998/09/09 08:07:11  jtanaka
// added Ptype(idhep).
//
// Revision 1.8  1998/09/03 23:55:36  jtanaka
// We need to include "Particle.h" only. Comment out some members of Ptype.
//
// Revision 1.7  1998/07/23 12:51:42  katayama
// conform ANSI
//
// Revision 1.6  1998/07/16 07:05:23  higuchit
// gave up of using PDT/CLHEP
//
// Revision 1.5  1998/07/03 12:10:14  jtanaka
// modified some parts
//
// Revision 1.4  1998/07/02 09:29:09  higuchit
// null flag -> `usable' flag
//
// Revision 1.3  1998/07/01 11:55:04  jtanaka
// add m_null.
//
// Revision 1.2  1998/06/16 05:13:50  katayama
// added friend class Particle so that Particle.cc compiles
//
// Revision 1.1  1998/06/14 14:11:56  higuchit
// *** empty log message ***
//
#ifndef PARTICLE_CLASS_PTYPE_H
#define PARTICLE_CLASS_PTYPE_H

#include <string>

#include "particle/constant.h"

#include "belle.h"
#include QQ_H
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

/// Ptype class supplies you interfaces for various properties of particle, such as its mass etc..
class Ptype
{
public:
  /// Default constructor
  Ptype(void) : m_usable(false), m_qq(NULL) {}
  /// Constructor with LUND7 particle code
  explicit Ptype(const int);
  /// Constructor with particle name of the decay.dec or user.dec
  explicit Ptype(const char *);
  
  /// Destructor
  virtual ~Ptype(void){}

  /// dumps debug information. (not implement.)
  virtual void dump(const std::string & keyword = std::string("name mass return"),
		    const std::string & prefix  = std::string("")) const;

  /// returns LUND7 particle code.
  inline const long lund(void) const;
  /// returns particle name of your decay.dec and user.dec.
  inline const std::string & name(void) const;

  /// returns mass(nominal mass).
  inline const double mass(void) const;

  /// returns c*tau.
  inline const double cTau(void) const;

  /// returns charge. (in units of e)
  inline const double charge(void) const;
  /// returns charge. (in units of e/3)
  inline const double iCharge(void) const;

  /// returns the spin. In the set method the spin should be given as
  /// 2J+1 in units of hbar/2. In the spin() method the spin is
  /// returned in units of hbar.
  inline const double spin(void) const;

  /// returs "Ptype" to the corresponding anti partner.
  Ptype chargeConjugation(void) const;

  /// returns whether particle is to be considered stable. If the decay
  /// table is empty, this method always returns true, even if the
  /// member variable is false.
  inline const bool stable(void) const;

  /// returns number of decay modes
  const unsigned int nDecay(void) const;

  /// == operator (using particle internal pointer to PANTHER)
  bool operator == ( const Ptype &p ) const { return m_qq == p.m_qq;}

  bool operator == ( const char *p ) const { return *this == Ptype(p);}

  /// != operator (using particle internal pointer to PANTHER)
  bool operator != ( const Ptype &p ) const { return !(m_qq == p.m_qq);}
  bool operator != ( const char *p ) const { return !(*this == Ptype(p));}

  /// Bool operator
  operator bool(void) const { return  m_usable; }

  /// ! operator
  bool operator!(void) const { return !m_usable; }
	
private:
  bool m_usable;
  struct qq_particle *m_qq;
};

inline const long Ptype::lund(void) const { return m_qq->m_pdgid; }

inline const std::string &Ptype::name(void) const 
{
  static std::string name;
  return name = m_qq->m_name;
}

inline const double Ptype::mass(void) const { return m_qq->m_mass; }

inline const double Ptype::cTau(void) const  { return m_qq->m_ctau; }

inline const double Ptype::charge(void) const { return m_qq->m_charge; }
inline const double Ptype::iCharge(void) const { return m_qq->m_charge*3.0; }

inline const double Ptype::spin(void) const { return m_qq->m_spin; }

inline const bool Ptype::stable(void) const { return !(m_qq->m_decay[0]); }

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* PARTICLE_CLASS_PTYPE_H */
