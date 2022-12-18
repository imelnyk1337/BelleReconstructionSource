//================================================================
//  File Name : beame.cc
//----------------------------------------------------------------
//  Creation ; 2004.01.10
//  Author   ; T.Shibata
//----------------------------------------------------------------
#include "belle.h"
#include "event/BelleEvent.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "benergy/Beame.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//----------------------------------------------------------------
extern "C" Module_descr *mdcl_beame() {

  beame *module   = new beame;
  Module_descr *dscr = new Module_descr("beame", module);

  BeamEnergy::define_global(dscr);

  return dscr;

}
//----------------------------------------------------------------
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
