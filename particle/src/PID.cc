//
//  PID.cc
//
//  A class implimentation for the particle ID class <PID> for the Belle
//  particle class <Particle> structure.
//
//  An object of <PID> class is a private member of a <Particle>.
//
//  When a PID object is constructed, it retrieves relevant information
//  from Mdst_charged and other MDST tables.
//
//  First version is written by M. Nakao (KEK)
//
// $Log$
// Revision 1.19  2003/09/15 13:24:34  katayama
// reduce warnings
//
// Revision 1.18  2001/12/13 15:31:55  katayama
// MDST_OBS
//
// Revision 1.17  2000/04/07 10:56:03  jtanaka
// add new functions to use new muon ID class
//
// Revision 1.16  2000/01/05 06:09:38  jtanaka
// major updates: please see BELLE whiteboard.
//
// Revision 1.15  1999/01/16 10:31:37  katayama
// clean up includes
//
// Revision 1.14  1998/11/17 22:27:54  jtanaka
// id -> lund in Particle, and bug fix in PID.cc
//
// Revision 1.13  1998/11/10 06:50:57  jtanaka
// For new MDST table.
//
// Revision 1.12  1998/11/07 07:54:38  katayama
// Modified to compile with the new mdst.
//
// Revision 1.11  1998/10/22 11:38:01  jtanaka
// changed names of the member.  e.g) del_info --> delInfo
//
// Revision 1.10  1998/10/12 19:24:08  jtanaka
// updated constructors Momentum.cc, Particle.cc, and Relation.cc. removed and added some functions in PID.cc. added (int*) to ParticleManager.cc(warning).
//
// Revision 1.9  1998/09/29 10:10:31  katayama
// sigma_tof is not an array anymore
//
// Revision 1.8  1998/09/28 06:31:52  jtanaka
// modified PID.cc for the new mdst_klm.
//
// Revision 1.7  1998/09/07 13:26:12  jtanaka
// implement some fuctions in PID Class, and modify some function in Relation Class
//
// Revision 1.6  1998/07/31 01:56:15  nakao
// acc, trk also need to be checked before using quality()
//
// Revision 1.5  1998/07/30 06:41:59  nakao
// need to check if tof is available in PIDtof::init()
//
// Revision 1.4  1998/07/18 14:00:25  katayama
// Tof has been updated
//
// Revision 1.3  1998/07/02 09:28:22  higuchit
// null flag -> `usable' flag
//
// Revision 1.2  1998/07/01 11:53:43  jtanaka
// add m_null.
//
// Revision 1.1  1998/06/05 05:35:19  nakao
// PID class for Particle class, first revision.
//

#include "particle/Particle.h"
#include "particle/PID.h"
#include "belle.h"
#include MDST_H

#include "kid/atc_pid.h"
#include "eid/eid.h"
#include "mdst/Muid_mdst.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


PID::PID(const PID &pid)
{
  m_self = pid.m_self;
  if(m_self){
    m_kId = new atc_pid(*(pid.m_kId));
    m_eId = new eid(m_self->mdstCharged());
    m_muId = new Muid_mdst(m_self->mdstCharged());
  }else{
    m_kId = NULL;
    m_eId = NULL;
    m_muId = NULL;
  }
  m_muonId = pid.m_muonId;
}

PID::PID(const Particle &p,
	 const int &accq0,
	 const int &tofq0,
	 const int &cdcq0,
	 const int &ids0,
	 const int &idb0)
{
  if(p.mdstCharged()){
    m_self = &p;
    m_kId = new atc_pid(accq0,tofq0,cdcq0,ids0,idb0);
    m_eId = new eid(p.mdstCharged());
    m_muId = new Muid_mdst(p.mdstCharged());
    if(p.mdstCharged().muid())m_muonId = p.mdstCharged().muid().muon();
    else m_muonId = 0;
  }else{
    m_self = NULL;
  }
}

const atc_pid & 
PID::kId(void) const
{
  static atc_pid * null(NULL);
  if(m_self)return *m_kId;
  else return *null;
}

double
PID::kIdProb(void) const
{
  if(m_self)return (double)(m_kId->prob(&(m_self->mdstCharged())));
  else return 0.;
}

double
PID::kIdProb(const atc_pid &kid) const
{
  if(m_self)return (double)(kid.prob(&(m_self->mdstCharged())));
  else return 0.;
}

const eid &
PID::eId(void) const
{
  static eid * null(NULL);
  if(m_self)return *m_eId;
  return *null;
}

unsigned
PID::muonId(void) const
{
  if(m_self)return m_muonId;
  else return 0;
}

const Muid_mdst & 
PID::muId(void) const
{
  static Muid_mdst * null(NULL);
  if(m_self)return *m_muId;
  return *null;
}

bool
PID::usable(void) const
{
  if(m_self)return true;
  else return false;
}

PID &
PID::operator=(const PID &pid)
{
  if(this == &pid)return *this;

  m_self = pid.m_self;
  if(m_self){
    delete m_kId;
    delete m_eId;
    delete m_muId;
    m_kId = new atc_pid(*(pid.m_kId));
    m_eId = new eid(m_self->mdstCharged());
    m_muId = new Muid_mdst(m_self->mdstCharged());
  }else{
    m_kId = NULL;
    m_eId = NULL;
    m_muId = NULL;
  }
  m_muonId = pid.m_muonId;
  return *this;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
