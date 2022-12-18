//
// Momentum Class of Particle Class
//
// --- Momentum.h ---
//
// $Id: Momentum.h 10002 2007-02-26 06:56:17Z katayama $
//
// $Log$
// Revision 1.27  2002/01/03 11:04:38  katayama
// Point3D and other header files are cleaned
//
// Revision 1.26  2001/12/14 02:16:24  katayama
// MDST_OBS removed
//
// Revision 1.25  2000/05/16 14:46:26  jtanaka
// added constructor of Mdst_ecl.
//
// Revision 1.24  2000/04/14 12:40:56  jtanaka
// updated for new table "mdst_vee2".
//
// Revision 1.23  2000/03/08 07:56:18  jtanaka
// modify a constructor of the Mdst_charged and
// add new functions(please use the lastest kfitter:May 8,2000)
//
// Revision 1.22  2000/03/07 11:14:00  katayama
// compatibility with CC5.0
//
// Revision 1.21  2000/01/05 06:10:37  jtanaka
// major updates:please see BELLE whiteboard.
//
// Revision 1.20  1999/04/09 15:01:08  jtanaka
// Added new member functions, "const", "&"...etc. Removed some member functions.
//
// Revision 1.19  1999/01/16 10:31:27  katayama
// clean up includes
//
// Revision 1.18  1998/11/10 06:50:39  jtanaka
// For new MDST table.
//
// Revision 1.17  1998/10/22 11:38:16  jtanaka
// changed names of the member.  e.g) del_info --> delInfo
//
// Revision 1.16  1998/10/20 14:45:49  jtanaka
// added FittedMomentum.
//
// Revision 1.15  1998/10/12 19:25:44  jtanaka
// updated constructors Momentum.h, Particle.h, Relation.h. removed and added some functions in PID.h. added some static_PIDs to static_particle.h.
//
// Revision 1.14  1998/10/05 10:59:29  jtanaka
// added some comments for doc++ in some headers files and modified static objects in Relation.h, ParticleManager.h.
//
// Revision 1.13  1998/09/08 09:54:35  jtanaka
// Modify some functions in Momentum.
//
// Revision 1.12  1998/07/23 12:51:39  katayama
// conform ANSI
//
// Revision 1.11  1998/07/22 13:04:39  jtanaka
// add some functions and modify const functions a little.
//
// Revision 1.10  1998/07/13 10:55:15  jtanaka
// modify a little about const functions and add vertex to Momentum
//
// Revision 1.9  1998/07/03 12:10:14  jtanaka
// modified some parts
//
// Revision 1.8  1998/07/02 09:29:08  higuchit
// null flag -> `usable' flag
//
// Revision 1.7  1998/07/01 11:54:34  jtanaka
// add m_null.
//
// Revision 1.6  1998/06/23 05:35:51  katayama
// Needs panther/panther.h
//
// Revision 1.5  1998/06/22 05:43:45  katayama
// Use XXXX_H instead of explicit file names
//
// Revision 1.4  1998/06/19 09:09:36  jtanaka
// *** empty log message ***
//
// 
#ifndef PARTICLE_CLASS_MOMENTUM_H
#define PARTICLE_CLASS_MOMENTUM_H

#include "belle.h"
#include <string>

#include "belleCLHEP/Vector/LorentzVector.h"
#include "belleCLHEP/Matrix/Vector.h"
#include "belleCLHEP/Matrix/SymMatrix.h"
#include "belleCLHEP/Matrix/Matrix.h"

#include "helix/Helix.h"
#include "particle/constant.h"
#include "particle/Ptype.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class Mdst_charged;
class Mdst_gamma;
class Mdst_trk;
class Mdst_vee2;
class Mdst_klong;
class Mdst_pi0;
class Mdst_ecl;
class Gen_hepevt;

/// Mometum class supplies you interfaces for various values of particle, such as its momentum, position, and vertex info. etc..
class Momentum 
{

public:// Constructor
  /// Default constructor
  Momentum();
  
  /// Copy constructor
  Momentum(const Momentum &);
  
  /// Constructor with momentum
  Momentum(const HepLorentzVector &, 
	   const HepSymMatrix &error = HepSymMatrix(4,0));

  /// Constructor with Mdst\_charged
  Momentum(const Mdst_charged &, const double mass,
	   const HepPoint3D & = HepPoint3D(0.,0.,0.));
  /// Constructor with Mdst\_charged
  Momentum(const Mdst_charged &, const Ptype &,
	   const HepPoint3D & = HepPoint3D(0.,0.,0.));
  /// Constructor with Mdst\_gamma
  Momentum(const Mdst_gamma &);
  /// Constructor with Mdst\_vee2
  Momentum(const Mdst_vee2 &);
  /// Constructor with Mdst\_klong 
  Momentum(const Mdst_klong &, double momentum_size = 1.);
  /// Constructor with Mdst\_pi0
  Momentum(const Mdst_pi0 &);
  /// Constructor with Mdst\_elid
  Momentum(const Mdst_ecl &,
	   const HepPoint3D & = HepPoint3D(0.,0.,0.));
  /// Constructor with Gen\_hepevt
  Momentum(const Gen_hepevt &);

  /// Destructor
  virtual ~Momentum();

public:// General interfaces
  /// returns class name.
  virtual std::string className(void){ return std::string("Momentum"); }

  /// dumps debug information. (not implement.)
  virtual void dump(const std::string &keyword = std::string("mass momentum return"),
		    const std::string &prefix  = std::string("")) const;
  
public:// Selectors
  /// returns momentum vector.
  virtual const HepLorentzVector & p(void) const { return m_momentum; }
  
  /// returns error matrix(4x4) of momentum vector. (not reference)
  virtual const HepSymMatrix dp(void) const { return m_error.sub(1,4); }
  
  /// returns position vector.
  virtual const HepPoint3D & x(void) const { return m_position; }
  
  /// returns error matrix(3x3) of position vector. (not reference)
  virtual const HepSymMatrix dx(void) const { return m_error.sub(5,7); }

  /// returns error matrix of momentum and position vector.
  virtual const HepSymMatrix & dpx(void) const { return m_error; }

  /// returns mass.
  virtual double mass(void) const { return m_momentum.mag(); }
  
  /// returns error of mass.
  virtual double dMass(void) const;
  
  /// retruns production vertex.
  virtual const HepPoint3D & vertex(void) const { return m_vertex; }

  /// retruns error matrix(3x3) of production vertex.
  virtual const HepSymMatrix & dVertex(void) const { return m_vertexError; }

  /// retruns decay vertex.
  virtual const HepPoint3D & decayVertex(void) const { return m_decayVertex; }

  /// retruns error matrix(3x3) of decay vertex.
  virtual const HepSymMatrix & dDecayVertex(void) const { return m_decayVertexError; }

public:// Modifiers
  /// sets momentum vector and its error matrix(4x4).
  virtual void momentum(const HepLorentzVector &,
			const HepSymMatrix     &error = HepSymMatrix(4,0));
  
  /// sets position vector and its error matrix(3x3).
  virtual void position(const HepPoint3D   &, 
			const HepSymMatrix &error = HepSymMatrix(3,0));

  /// sets momentum and position vector and its error matrix(7x7).
  virtual void momentumPosition(const HepLorentzVector &, const HepPoint3D &, 
				const HepSymMatrix &error = HepSymMatrix(7,0));

  /// sets production vertex and its error matrix(3x3).
  virtual HepPoint3D & vertex(const HepPoint3D &, 
			      const HepSymMatrix &error = HepSymMatrix(3,0)); 

  /// sets decay vertex and its error matrix(3x3).
  virtual HepPoint3D & decayVertex(const HepPoint3D &, 
				   const HepSymMatrix &error = HepSymMatrix(3,0)); 

public:// Operators
  /// copy operator.
  Momentum & operator = (const Momentum &);

  //protected:// Protected members
private:// Private members
  HepLorentzVector m_momentum;       // 4
  HepPoint3D       m_position;       // 3
  HepSymMatrix     m_error;          // 7 x 7

  HepPoint3D       m_vertex;         // 3      ... production vertex
  HepSymMatrix     m_vertexError;    // 3 x 3  ... production vertex
  HepPoint3D       m_decayVertex;         // 3
  HepSymMatrix     m_decayVertexError;    // 3 x 3
};
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* PARTICLE_CLASS_MOMENTUM_H */
