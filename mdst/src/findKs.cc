//
//  findKs.cc ver.0
//  F. Fang (Univ. of Hawaii)
//  select good ks from MDST_vee
//
//  Revision 1.0 2000/04/11 F. Fang
//  modified according to the changes in mdst tables

#include <string>
#include <iostream>
#include <cmath>

#include "belle.h"
#include "mdst/findKs.h"

#include "belleCLHEP/Vector/ThreeVector.h"
#include "belleCLHEP/Matrix/Vector.h"
#include "helix/Helix.h"


#include "panther/panther.h"
#include MDST_H
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


double
FindKs::_impactParameter( const Mdst_vee_daughters& dau, HepPoint3D runIP, 
			  HepPoint3D vtx, int charge ) {
  double  dr = 1000.0;
#if defined(BELLE_DEBUG)
    try {
#endif
  if(dau){
    if (charge==1) {   //+
      HepVector a(5,0);
      a[0] = dau.helix_p(0);
      a[1] = dau.helix_p(1);
      a[2] = dau.helix_p(2);
      a[3] = dau.helix_p(3);
      a[4] = dau.helix_p(4);
      Helix ltrk(vtx, a);       
      ltrk.pivot(runIP);              
      dr=ltrk.dr();  
    }
    else if (charge==-1) {
      HepVector a(5,0);
      a[0] = dau.helix_m(0);
      a[1] = dau.helix_m(1);
      a[2] = dau.helix_m(2);
      a[3] = dau.helix_m(3);
      a[4] = dau.helix_m(4);
      Helix ltrk(vtx, a);       
      ltrk.pivot(runIP);              
      dr=ltrk.dr();  
    }
  }   
#if defined(BELLE_DEBUG)
    }
    catch(const std::string &e) {
      dout(Debugout::ERR,"findKs") << "findKs:_InpactParameter" << e << std::endl;
	Mdst_charged_Manager::get_manager().dump();
	Mdst_trk_Manager::get_manager().dump();
	Mdst_vee2_Manager::get_manager().dump();
      return 1000.0;
    }
#endif    
  return dr;
}

double 
FindKs::_impactParameter( const Mdst_charged& Charged, HepPoint3D runIP ){

  Mdst_trk& trk = Charged.trk();
  Mdst_trk_fit& trkFit = trk.mhyp(2);  /* assume pion mass */

  double dr = 1000.0;
#if defined(BELLE_DEBUG)
    try {
#endif
  if(trk){
    if(trkFit){
      HepVector a(5,0);
    
      a[0] = trkFit.helix(0);
      a[1] = trkFit.helix(1);
      a[2] = trkFit.helix(2);
      a[3] = trkFit.helix(3);
      a[4] = trkFit.helix(4);

      HepPoint3D pivot(trkFit.pivot(0), trkFit.pivot(1), trkFit.pivot(2));
      Helix ltrk(pivot, a);         
      ltrk.pivot(runIP);          
      dr=ltrk.dr();      
    }
  }
#if defined(BELLE_DEBUG)
    }
    catch(const std::string &e) {
      dout(Debugout::ERR,"findKs") << "findKs:_InpactParameter" << e << std::endl;
	Mdst_charged_Manager::get_manager().dump();
	Mdst_trk_Manager::get_manager().dump();
	Mdst_vee2_Manager::get_manager().dump();
      return 1000.0;
    }
#endif    
  return dr;
}

void
FindKs::candidates ( const Mdst_vee2& ksCand, HepPoint3D runIP )
{
  _IP = runIP;
  _kind = ksCand.kind();
  _chisq = ksCand.chisq();

  Hep3Vector ksP( ksCand.px(),ksCand.py(),ksCand.pz() );
  _pmag = ksP.mag();
  HepPoint3D vtx ( ksCand.vx(),ksCand.vy(),ksCand.vz() ) ;
  HepPoint3D v2IP = vtx - runIP; 
  //fl
  _fl = sqrt(v2IP.x()*v2IP.x()+v2IP.y()*v2IP.y());
  //dphi
  double dot =  v2IP.x()*ksP.x() + v2IP.y()*ksP.y();
  double atot = _fl*sqrt(ksP.x()*ksP.x() + ksP.y()*ksP.y());
  if(abs(dot)>atot) _dphi=0.0;
  else _dphi = atot==0.0 ? 0.0:acos(dot/atot);
  
  //z_dist
  _zdist=ksCand.z_dist();
  
  //dr1 and dr2
  Mdst_charged& pi1= ksCand.chgd(0);
  Mdst_charged& pi2= ksCand.chgd(1);
  double dr1 = _impactParameter(pi1, _IP);
  double dr2 = _impactParameter(pi2, _IP);
  _dr=std::min(abs(dr1),abs(dr2));
}

int
FindKs::goodKs()
{    //square cuts
  if ( _kind == 1 ) {
    bool low=_pmag<0.5 && _zdist<0.8 && _dr>0.05 && _dphi<0.3;
    bool mid=_pmag<1.5 && _pmag>0.5 && _zdist<1.8 && _dr>0.03 && _dphi<0.1 && _fl>.08;
    bool high=_pmag>1.5 && _zdist<2.4 && _dr>0.02 && _dphi<0.03 && _fl>.22;
    if ( low || mid || high ) {
      return 1;
    }
    else return 0;
  }
  else return 0;
}






#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
