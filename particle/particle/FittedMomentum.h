//
// Fitted Momentum Class of Particle Class
//
// --- FittedMomentum.h ---
//
// $Id: FittedMomentum.h 9932 2006-11-12 14:26:53Z katayama $
//
// $Log$
// Revision 1.5  2000/03/07 11:13:59  katayama
// compatibility with CC5.0
//
// Revision 1.4  1999/04/15 20:04:02  jtanaka
// add statical info.
//
// Revision 1.3  1999/04/09 15:01:07  jtanaka
// Added new member functions, "const", "&"...etc. Removed some member functions.
//
// Revision 1.2  1998/10/22 11:38:16  jtanaka
// changed names of the member.  e.g) del_info --> delInfo
//
// Revision 1.1  1998/10/20 14:45:49  jtanaka
// added FittedMomentum.
//
//
#ifndef PARTICLE_CLASS_FITTED_MOMENTUM_H
#define PARTICLE_CLASS_FITTED_MOMENTUM_H

#include "belle.h"
#include "particle/Momentum.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
class Particle;

/* detele types */
//#define REMOVE_COMATRIX 1
//#define REMOVE_COVERTEX 2
//#define REMOVE_CODECAYVERTEX 4
//#define REMOVE_ALL (REMOVE_COMATRIX+REMOVE_COVERTEX+REMOVE_CODECAYVERTEX)
const unsigned REMOVE_COMATRIX = 1;
const unsigned REMOVE_COVERTEX = 2;
const unsigned REMOVE_CODECAYVERTEX = 4;
const unsigned REMOVE_ALL = REMOVE_COMATRIX+REMOVE_COVERTEX+REMOVE_CODECAYVERTEX;

/// FittedMometum class supplies you interfaces for various values of fitted particle, such as its momentum, position, and vertex info. etc..
class FittedMomentum : public Momentum
{

public:// Constructor
  /// Default constructor
  FittedMomentum();
  
  /// Copy constructor
  FittedMomentum(const FittedMomentum &);

  /// Destructor
  virtual ~FittedMomentum();

public:// General interfaces
  /// returns class name.
  std::string className(void){ return std::string("FittedMomentum"); }
  
public:// Selectors
  /// returns correlation matrix with other particles.
  /// When we use "error", "*error = 0" means that correlation matrix can not found.
  virtual const HepMatrix & coMatrix(Particle &, unsigned *error = NULL) const;

  /// returns correlation matrix with other particles.
  /// When we use "error", "*error = 0" means that correlation matrix can not found.
  virtual const HepMatrix & coMatrix(Particle *, unsigned *error = NULL) const;
  
  /// returns correlation matrix with other particles.
  /// When we use "error", "*error = 0" means that correlation matrix can not found.
  virtual const HepMatrix & coMatrix(unsigned, unsigned *error = NULL) const;

  /// returns correlated particle.
  /// When we use "error", "*error = 0" means that correlated particle can not found.
  virtual const Particle * coParticle(unsigned, unsigned *error = NULL) const;

  /// returns correlation matrix with product vertex position.
  virtual const HepMatrix & coVertex(void) const;

  /// returns correlation matrix with decay vertex position.
  virtual const HepMatrix & coDecayVertex(void) const;

  /// returns number of correlated particles
  virtual unsigned nCoMatrix(void) const;

  /// returns CL.
  virtual double cl(void) const;

  /// returns chisq.
  virtual double chisq(void) const;

  /// returns degree of freedom.
  virtual unsigned dof(void) const;

public:// Modifiers
  /// sets correlation matrix with other particles and returns.
  virtual HepMatrix & coMatrix(Particle &, HepMatrix &);

  /// sets correlation matrix with other particles and returns.
  virtual HepMatrix & coMatrix(Particle *, HepMatrix &);

  /// sets correlation matrix with product vertex position and returns.
  virtual HepMatrix & coVertex(HepMatrix &);

  /// sets correlation matrix with decay vertex position and returns.
  virtual HepMatrix & coDecayVertex(HepMatrix &);

  /// removes info.
  virtual void remove(unsigned = REMOVE_ALL);

  /// removes correlated particles.
  virtual void removeCoMatrix(const Particle &);

  /// removes correlated particles.
  virtual void removeCoMatrix(const Particle *);

  /// sets CL.
  virtual double cl(const double &);

  /// sets chisq.
  virtual double chisq(const double &);

  /// sets degree of freedom.
  virtual unsigned dof(const unsigned &);

public:// Operators
  /// Copy operator
  FittedMomentum & operator = (const FittedMomentum &);

  /// Copy operator
  FittedMomentum & operator = (const Momentum &);

private:// Private members
  std::vector<Particle *> m_correlationParticles;
  std::vector<HepMatrix> m_correlationMatrix; // 7 x 7
  
  HepMatrix m_vertexAndSelfError; // 3 x 7
  HepMatrix m_decayVertexAndSelfError; // 3 x 7

  double m_cl;
  double m_chisq;
  unsigned m_dof;
};
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* PARTICLE_CLASS_FITTED_MOMENTUM_H */
