//
// ParticleManager.h
//
// A class to handle particle objects.
//
// Filename : ParticleManger.h
// Author : Yoshihito Iwasaki
// E-mail : yoshihito.iwasaki@kek.jp
//
// $Id: ParticleManager.h 9932 2006-11-12 14:26:53Z katayama $
//
// $Log$
// Revision 1.7  2000/06/13 05:20:48  yiwasaki
// updates
//
// Revision 1.6  2000/05/03 09:25:30  yiwasaki
// minor changes
//
// Revision 1.5  2000/03/07 11:14:00  katayama
// compatibility with CC5.0
//
// Revision 1.4  1998/10/05 10:59:30  jtanaka
// added some comments for doc++ in some headers files and modified static objects in Relation.h, ParticleManager.h.
//
// Revision 1.3  1998/09/22 15:34:23  yiwasaki
// ParticleManager extended, Relation::mc modified
//
// Revision 1.2  1998/09/22 01:24:22  yiwasaki
// ParticleManager minor change
//
// Revision 1.1  1998/09/08 16:49:54  yiwasaki
// ParticleManager added
//
//

#ifndef PARTICLE_CLASS_PARTICLEMANAGER_H
#define PARTICLE_CLASS_PARTICLEMANAGER_H

#include "belle.h"
#include "particle/Particle.h"
#include "particle/static_particle.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class Gen_hepevt_Manager;
class Mdst_charged_Manager;

/// Particle manager class handles particle objects.
class ParticleManager : public std::vector<Particle *> {

  public:
    /// Default constructor for general use.
    ParticleManager();

    /// Constructor for general use.
    ParticleManager(const std::string & name);

    /// Constructor with Gen_hepevt. Particles will be created and stored.
    ParticleManager(Gen_hepevt_Manager &);

    /// Constructor with Mdst_charged. Particles will be created and stored.
    ParticleManager(Mdst_charged_Manager &);

    /// Destructor
    virtual ~ParticleManager();

  public:
    /// returns/sets manager name.
    const std::string & name(void) const;
    const std::string & name(const std::string & name);

    /// dumps debug information.
    void dump(const std::string & keyword = std::string(""),
	      const std::string & prefix = std::string(""));

    /// links particle objects and mc particle objects.
    void link(const ParticleManager &) const;

    /// sets switch to delete objects when a manager destructed.
    bool deleteParticles(bool);

    /// sets all members' Ptype.
    void setPType(const Ptype &);

    /// sets Ptype.
    static void setPType(const Ptype &, Particle &);

    /// returns true if manager has a member p.
    bool hasMember(const Particle & p) const;
    
  private:
    std::string m_name;
    bool m_delete;
};

inline
const std::string &
ParticleManager::name(void) const {
    return m_name;
}

inline
const std::string &
ParticleManager::name(const std::string & name) {
    return m_name = name;
}

inline
bool
ParticleManager::deleteParticles(bool a) {
    return m_delete = a;
}

inline
bool
ParticleManager::hasMember(const Particle & p) const {
    for (unsigned i = 0; i < size(); i++)
	if ((* this)[i] == & p) return true;
    return false;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* PARTICLE_CLASS_PARTICLEMANAGER_H */
