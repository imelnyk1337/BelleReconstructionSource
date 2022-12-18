//
//  PID.h
//
//  A class definition for the particle ID class <PID> for the Belle
//  particle class <Particle> structure.
//
//  An object of <PID> class is a private member of a <Particle>.
//
//  When a PID object is constructed, it retrieves relevant information
//  from Mdst_charged and other MDST tables.
//
//  First version written by M. Nakao (KEK)
//
// $Log$
// Revision 1.17  2001/12/13 15:31:54  katayama
// MDST_OBS
//
// Revision 1.16  2000/04/14 12:40:56  jtanaka
// updated for new table "mdst_vee2".
//
// Revision 1.15  2000/04/07 22:01:21  katayama
// muid_mdst is now in mdst
//
// Revision 1.14  2000/04/07 10:55:45  jtanaka
// add new functions to use new muon ID class
//
// Revision 1.13  2000/01/05 06:10:38  jtanaka
// major updates:please see BELLE whiteboard.
//
// Revision 1.12  1999/01/16 10:31:29  katayama
// clean up includes
//
// Revision 1.11  1998/11/10 06:50:40  jtanaka
// For new MDST table.
//
// Revision 1.10  1998/11/07 07:54:34  katayama
// Modified to compile with the new mdst.
//
// Revision 1.9  1998/10/22 11:38:17  jtanaka
// changed names of the member.  e.g) del_info --> delInfo
//
// Revision 1.8  1998/10/12 19:25:45  jtanaka
// updated constructors Momentum.h, Particle.h, Relation.h. removed and added some functions in PID.h. added some static_PIDs to static_particle.h.
//
// Revision 1.7  1998/10/05 10:59:29  jtanaka
// added some comments for doc++ in some headers files and modified static objects in Relation.h, ParticleManager.h.
//
// Revision 1.6  1998/09/07 13:25:49  jtanaka
// implement some fuctions in PID Class, and modify some function in Relation Class
//
// Revision 1.5  1998/07/23 12:51:40  katayama
// conform ANSI
//
// Revision 1.4  1998/07/03 12:10:13  jtanaka
// modified some parts
//
// Revision 1.3  1998/07/02 09:29:09  higuchit
// null flag -> `usable' flag
//
// Revision 1.2  1998/07/01 11:54:42  jtanaka
// add m_null.
//
// Revision 1.1  1998/06/05 05:35:15  nakao
// PID class for Particle class, first revision.
//
#ifndef PARTICLE_CLASS_PID_H
#define PARTICLE_CLASS_PID_H

#include "belle.h"
#include "particle/constant.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class Particle;
class Mdst_charged;
class atc_pid;
class eid;
class Muid_mdst;

class PID {
public:
  /// Default constructor
  PID() : m_self(NULL), m_kId(NULL), m_eId(NULL) {};

  /// Copy constructor
  PID(const PID &pid);
  
  /// Constructor with Particle
  PID(const Particle &p,
      const int &accq0 = 0,
      const int &tofq0 = 1,
      const int &cdcq0 = 0,
      const int &ids0  = 3,
      const int &idb0  = 2);

  /// Destructor
  ~PID(){};

  /// returns atc\_pid object.
  virtual const atc_pid & kId(void) const;

  /// returns probability of kaon.
  virtual double kIdProb(void) const;

  /// returns probability of kaon using given atc\_pid object.
  virtual double kIdProb(const atc_pid&) const;

  /// returns eid object.
  virtual const eid & eId(void) const;

  /// returns muon id.
  virtual unsigned muonId(void) const;

  /// returns muid object.
  virtual const Muid_mdst & muId(void) const;

  /// returns "usable" info.
  virtual bool usable(void) const;

  /// copy operator
  PID & operator = (const PID &);

protected:
  /// pointer to Particle Class(self)
  const Particle * m_self;
  /// Kaon ID
  atc_pid *m_kId;
  /// Electron ID
  eid *m_eId;
  /// Muon ID;
  unsigned m_muonId;
  Muid_mdst *m_muId;
};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* PARTICLE_CLASS_PID_H */
