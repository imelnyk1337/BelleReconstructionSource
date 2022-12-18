#ifndef PARTICLE_CLASS_PARTICLE_USER_INFO_H
#define PARTICLE_CLASS_PARTICLE_USER_INFO_H
#include "belle.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


class ParticleUserInfo {
  
public:
  /// Default constructor
  ParticleUserInfo(){};

  /// Copy constructor
  ParticleUserInfo(const ParticleUserInfo&){};

  /// Destructor
  virtual ~ParticleUserInfo(){};

  /// constructs self object.
  virtual ParticleUserInfo * clone(void) const = 0;

  /// Copy operator
  virtual ParticleUserInfo & operator = (const ParticleUserInfo &x){ return *this; }

  /// append daughter
  virtual void append(ParticleUserInfo &){}
};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* PARTICLE_CLASS_PARTICLE_USER_INFO_H */
