//
// Momentum Class of Particle Class
//
// --- Momentum.cc ---
//
// $Id: Momentum.cc 9944 2006-11-29 07:36:07Z katayama $
//
// $Log$
// Revision 1.31  2002/02/25 02:13:39  katayama
// Use npos
//
// Revision 1.30  2001/12/14 02:16:32  katayama
// MDST_OBS removed
//
// Revision 1.29  2001/12/12 07:10:53  jtanaka
// (1) compatibility for obsoleted Mdst_vee/mdst_sim_xref table
// (2) compatibility for gcc3
//
// Revision 1.28  2000/06/30 10:37:16  jtanaka
// bug fix : Constructor of Mdst_ecl(from Yamaga-san).
//
// Revision 1.27  2000/05/16 14:47:00  jtanaka
// added constructor of Mdst_ecl
//
// Revision 1.26  2000/04/27 12:57:37  jtanaka
// 1. use Mdst_vee_daughters information to construct Mdst_vee2.
// 2. add new constructor Particle(Mdst_charged, string <-- "new", pivot).
//
// Revision 1.25  2000/04/14 12:40:41  jtanaka
// updated for new table "mdst_vee2".
//
// Revision 1.24  2000/04/13 12:41:57  katayama
// Added std:: to cout,cerr,endl etc.
//
// Revision 1.23  2000/03/08 07:56:46  jtanaka
// modify a constructor of the Mdst_charged.
//
// Revision 1.22  2000/03/07 11:14:03  katayama
// compatibility with CC5.0
//
// Revision 1.21  2000/01/05 06:09:37  jtanaka
// major updates: please see BELLE whiteboard.
//
// Revision 1.20  1999/11/10 14:03:01  jtanaka
// bug fix: operator and copy constructor in Particle.cc, constructor in Momentum.cc, checkSame in combination.cc, (int) in utility.cc
//
// Revision 1.19  1999/04/09 15:02:10  jtanaka
// Added some member functions,"const", and "&" and removed some members
//
// Revision 1.18  1999/01/16 10:31:35  katayama
// clean up includes
//
// Revision 1.17  1998/11/10 06:50:57  jtanaka
// For new MDST table.
//
// Revision 1.16  1998/11/07 08:02:33  katayama
// Offset starts from zero, unlike mdst.tdf comment
//
// Revision 1.15  1998/11/07 07:54:36  katayama
// Modified to compile with the new mdst.
//
// Revision 1.14  1998/10/22 11:38:00  jtanaka
// changed names of the member.  e.g) del_info --> delInfo
//
// Revision 1.13  1998/10/20 14:45:17  jtanaka
// added FittedMomentum.
//
// Revision 1.12  1998/10/15 09:52:31  jtanaka
// bug fix, and added deep_copy and deep_delete in Particle.h
//
// Revision 1.11  1998/10/12 19:24:07  jtanaka
// updated constructors Momentum.cc, Particle.cc, and Relation.cc. removed and added some functions in PID.cc. added (int*) to ParticleManager.cc(warning).
//
// Revision 1.10  1998/09/09 08:29:56  jtanaka
// added new constructor to Particle and modified some parts of Momentum.
//
// Revision 1.9  1998/09/08 09:54:00  jtanaka
// Modify Mdst_charged constructor in Particle and some functions in Momentum.
//
// Revision 1.8  1998/07/13 10:57:04  jtanaka
// modify a little about const functions and add vertex to Momentum
//
// Revision 1.7  1998/07/03 12:08:29  jtanaka
// modified some parts
//
// Revision 1.6  1998/07/02 09:28:22  higuchit
// null flag -> `usable' flag
//
// Revision 1.5  1998/07/01 11:53:28  jtanaka
// add m_null.
//
// Revision 1.4  1998/06/19 09:08:44  jtanaka
// *** empty log message ***
//
// 
#include "belle.h"
#include "particle/Momentum.h"

#include "panther/panther.h"
#include MDST_H
#include HEPEVT_H

// Calculation of Gamma(ECL)
#include "particle/gammac.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


Momentum::Momentum()
  : m_momentum(),
    m_position(),
    m_error(7,0),
    m_vertex(),
    m_vertexError(3,0),
    m_decayVertex(),
    m_decayVertexError(3,0)
{
}
  
Momentum::Momentum(const Momentum &p)
  : m_momentum(p.m_momentum),
    m_position(p.m_position),
    m_error(p.m_error),
    m_vertex(p.m_vertex),
    m_vertexError(p.m_vertexError),
    m_decayVertex(p.m_decayVertex),
    m_decayVertexError(p.m_decayVertexError)
{
}

Momentum::Momentum(const HepLorentzVector &p, 
		   const HepSymMatrix &error)
  : m_momentum(p),
    m_position(),
    m_error(7,0),
    m_vertex(),
    m_vertexError(3,0),
    m_decayVertex(),
    m_decayVertexError(3,0)
{
  if(error.num_row() == 4){
    m_error.sub(1,error);
    //sets 0.0 in momentum position correlation elements.  
    for(unsigned i=0;i<4;i++){
      for(unsigned j=4;j<7;j++){
	m_error[i][j] = 0.0;
      }
    }
    for(unsigned i=4;i<7;i++){
      for(unsigned j=i;j<7;j++){
	m_error[i][j] = 0.0;
      }
    }
  }else{
    dout(Debugout::WARN,"Momentum") << "Error in momentum(const HepLorentzVector&," << std::endl;
    dout(Debugout::WARN,"Momentum") << "                  const HepSymMatrix&) of Mometum Class" << std::endl;
    dout(Debugout::WARN,"Momentum") << "                        ~~~~~~~~~~~~~--- 4 x 4 error matrix!!! " << std::endl;
    //sets 0.0 in all elements.  
    for(unsigned i=0;i<7;i++){
      for(unsigned j=i;j<7;j++){
	m_error[i][j] = 0.0;
      }
    }
  }  
}

Momentum::Momentum(const Mdst_charged &p, const double thisMass,
		   const HepPoint3D &newPivot)
  : m_momentum(p.px(), p.py(), p.pz(),
	       sqrt(p.px()*p.px()+p.py()*p.py()+p.pz()*p.pz()+thisMass*thisMass)),
    m_error(7,0),
    m_vertex(),
    m_vertexError(3,0),
    m_decayVertex(),
    m_decayVertexError(3,0)
{
  int hyp = 4;
  if(thisMass < 0.005){ // e = 0.000511
    hyp = 0;
  }else if(thisMass < 0.110){ // mu = 0.1056
    hyp = 1;
  }else if(thisMass < 0.200){ // pi = 0.13956
    hyp = 2;
  }else if(thisMass < 0.5){ // K = 0.4936
    hyp = 3;
  }
  const HepPoint3D pivot(p.trk().mhyp(hyp).pivot_x(),
			 p.trk().mhyp(hyp).pivot_y(),
			 p.trk().mhyp(hyp).pivot_z());
  //HepMatrix  tmp_a(5,1);
  //tmp_a[0][0] = p.trk().mhyp(hyp).helix(0);
  //tmp_a[1][0] = p.trk().mhyp(hyp).helix(1);
  //tmp_a[2][0] = p.trk().mhyp(hyp).helix(2);
  //tmp_a[3][0] = p.trk().mhyp(hyp).helix(3);
  //tmp_a[4][0] = p.trk().mhyp(hyp).helix(4);
  //HepVector  a(tmp_a);
  HepVector  a(5);
  a[0] = p.trk().mhyp(hyp).helix(0);
  a[1] = p.trk().mhyp(hyp).helix(1);
  a[2] = p.trk().mhyp(hyp).helix(2);
  a[3] = p.trk().mhyp(hyp).helix(3);
  a[4] = p.trk().mhyp(hyp).helix(4);
  HepSymMatrix Ea(5,0);
  Ea[0][0] = p.trk().mhyp(hyp).error(0);
  Ea[1][0] = p.trk().mhyp(hyp).error(1);
  Ea[1][1] = p.trk().mhyp(hyp).error(2);
  Ea[2][0] = p.trk().mhyp(hyp).error(3);
  Ea[2][1] = p.trk().mhyp(hyp).error(4);
  Ea[2][2] = p.trk().mhyp(hyp).error(5);
  Ea[3][0] = p.trk().mhyp(hyp).error(6);
  Ea[3][1] = p.trk().mhyp(hyp).error(7);
  Ea[3][2] = p.trk().mhyp(hyp).error(8);
  Ea[3][3] = p.trk().mhyp(hyp).error(9);
  Ea[4][0] = p.trk().mhyp(hyp).error(10);
  Ea[4][1] = p.trk().mhyp(hyp).error(11);
  Ea[4][2] = p.trk().mhyp(hyp).error(12);
  Ea[4][3] = p.trk().mhyp(hyp).error(13);
  Ea[4][4] = p.trk().mhyp(hyp).error(14);
  Helix helix(pivot, a, Ea);

  if(newPivot.x() != 0. ||
     newPivot.y() != 0. ||
     newPivot.z() != 0.){
    helix.pivot(newPivot);
    m_momentum = helix.momentum(0.,thisMass,m_position,m_error);
    return;
  }else{
    //...finds ref. point.
    if(pivot.x() != 0. ||
       pivot.y() != 0. ||
       pivot.z() != 0.){
      //Point3D tmp(0.,0.,0.);
      //helix.pivot(tmp);
      helix.pivot(HepPoint3D(0.,0.,0.));
      m_momentum = helix.momentum(0.,thisMass,m_position,m_error);
      return;
    }else{
      m_momentum = helix.momentum(0.,thisMass,m_position,m_error);
      return;
    }
  }
}

Momentum::Momentum(const Mdst_charged &p, const Ptype &ptype,
		   const HepPoint3D &newPivot)
  : m_error(7,0),
    m_vertex(),
    m_vertexError(3,0),
    m_decayVertex(),
    m_decayVertexError(3,0)
{
  int hyp = 2;
  if(ptype.name() == "K+" || ptype.name() == "K-"){
    hyp = 3;
  }else if(ptype.name() == "E+" || ptype.name() == "E-"){
    hyp = 0;
  }else if(ptype.name() == "MU+" || ptype.name() == "MU-"){
    hyp = 1;
  }else if(ptype.name() == "P+" || ptype.name() == "AP+"){
    hyp = 4;
  }
  const HepPoint3D pivot(p.trk().mhyp(hyp).pivot_x(),
			 p.trk().mhyp(hyp).pivot_y(),
			 p.trk().mhyp(hyp).pivot_z());
  //HepMatrix  tmp_a(5,1);
  //tmp_a[0][0] = p.trk().mhyp(hyp).helix(0);
  //tmp_a[1][0] = p.trk().mhyp(hyp).helix(1);
  //tmp_a[2][0] = p.trk().mhyp(hyp).helix(2);
  //tmp_a[3][0] = p.trk().mhyp(hyp).helix(3);
  //tmp_a[4][0] = p.trk().mhyp(hyp).helix(4);
  //HepVector  a(tmp_a);
  HepVector  a(5);
  a[0] = p.trk().mhyp(hyp).helix(0);
  a[1] = p.trk().mhyp(hyp).helix(1);
  a[2] = p.trk().mhyp(hyp).helix(2);
  a[3] = p.trk().mhyp(hyp).helix(3);
  a[4] = p.trk().mhyp(hyp).helix(4);
  HepSymMatrix Ea(5,0);
  Ea[0][0] = p.trk().mhyp(hyp).error(0);
  Ea[1][0] = p.trk().mhyp(hyp).error(1);
  Ea[1][1] = p.trk().mhyp(hyp).error(2);
  Ea[2][0] = p.trk().mhyp(hyp).error(3);
  Ea[2][1] = p.trk().mhyp(hyp).error(4);
  Ea[2][2] = p.trk().mhyp(hyp).error(5);
  Ea[3][0] = p.trk().mhyp(hyp).error(6);
  Ea[3][1] = p.trk().mhyp(hyp).error(7);
  Ea[3][2] = p.trk().mhyp(hyp).error(8);
  Ea[3][3] = p.trk().mhyp(hyp).error(9);
  Ea[4][0] = p.trk().mhyp(hyp).error(10);
  Ea[4][1] = p.trk().mhyp(hyp).error(11);
  Ea[4][2] = p.trk().mhyp(hyp).error(12);
  Ea[4][3] = p.trk().mhyp(hyp).error(13);
  Ea[4][4] = p.trk().mhyp(hyp).error(14);
  Helix helix(pivot, a, Ea);

  if(newPivot.x() != 0. ||
     newPivot.y() != 0. ||
     newPivot.z() != 0.){
    helix.pivot(newPivot);
    m_momentum = helix.momentum(0.,ptype.mass(),m_position,m_error);
    return;
  }else{
    if(pivot.x() != 0. ||
       pivot.y() != 0. ||
       pivot.z() != 0.){
      //Point3D tmp(0.,0.,0.);
      //helix.pivot(tmp);
      helix.pivot(HepPoint3D(0.,0.,0.));
      m_momentum = helix.momentum(0.,ptype.mass(),m_position,m_error);
      return;
    }else{
      m_momentum = helix.momentum(0.,ptype.mass(),m_position,m_error);
      return;
    }
  }
}

Momentum::Momentum(const Mdst_gamma &p)
  : m_momentum(p.px(), p.py(), p.pz(), 
	       sqrt(p.px()*p.px()+p.py()*p.py()+p.pz()*p.pz())),
    m_position(),
    m_error(7,0),
    m_vertex(),
    m_vertexError(3,0),
    m_decayVertex(),
    m_decayVertexError(3,0)
{
}

Momentum::Momentum(const Mdst_vee2 &p)
  : m_momentum(p.px(), p.py(), p.pz(), p.energy()),
    m_position(p.vx(), p.vy(), p.vz()),
    m_error(7,0),
    m_vertex(),
    m_vertexError(3,0),
    m_decayVertex(),
    m_decayVertexError(3,0)
{
}

Momentum::Momentum(const Mdst_klong &p, double m)
  : m_position(),
    m_error(7,0),
    m_vertex(),
    m_vertexError(3,0),
    m_decayVertex(),
    m_decayVertexError(3,0)
{
  Ptype kl("K0L");
  //HepLorentzVector tmp(m*p.cos_x(),m*p.cos_y(),m*p.cos_z(),sqrt(m*m+kl.mass()*kl.mass()));
  //m_momentum = tmp;
  m_momentum = HepLorentzVector(m*p.cos_x(),m*p.cos_y(),m*p.cos_z(),sqrt(m*m+kl.mass()*kl.mass()));
}

Momentum::Momentum(const Mdst_pi0 &p)
  : m_momentum(p.px(),p.py(),p.pz(),p.energy()),
    m_position(),
    m_error(7,0),
    m_vertex(),
    m_vertexError(3,0),
    m_decayVertex(),
    m_decayVertexError(3,0)
{
}

Momentum::Momentum(const Mdst_ecl &p,
		   const HepPoint3D &org)
  : m_position(org), // Production Point
    m_error(7,0),
    m_vertex(org), // Production Point
    m_vertexError(3,0),
    m_decayVertex(),
    m_decayVertexError(3,0)
{
  HepSymMatrix err(3,0);
  err[0][0] = p.error(0);
  err[1][0] = p.error(1);
  err[1][1] = p.error(2);
  err[2][0] = p.error(3);
  err[2][1] = p.error(4);
  err[2][2] = p.error(5);
  GammaParticle g(p.energy(), p.phi(), p.theta(),
		  p.r(), err);
  if(org.x() != 0. ||
     org.y() != 0. ||
     org.z() != 0.)g.vertex(org);
  m_momentum = g.momentumEnergy();
  //m_error    = g.errorMomentumEnergy();
  m_error.sub(1,g.errorMomentumEnergy());
}

Momentum::Momentum(const Gen_hepevt &p)
  : m_momentum(p.PX(), p.PY(), p.PZ(), p.E()),
    m_position(p.VX()*0.1, p.VY()*0.1, p.VZ()*0.1),
    m_error(7,0),
    m_vertex(p.VX()*0.1,p.VY()*0.1,p.VZ()*0.1),
    m_vertexError(3,0),
    m_decayVertex(),
    m_decayVertexError(3,0)
{  
}
  
Momentum::~Momentum()
{
}

//void 
//Momentum::dump(const string &keyword, const string &prefix) const {
//} 


double 
Momentum::dMass(void) const 
{
  // mass  = sqrt(E^2 - p^2)
  // dmass = (EdE - pxdpx - pydpy -pzdpz)/mass
  // dMass = dmass^2
  double E  = m_momentum.t();
  double px = m_momentum.x();
  double py = m_momentum.y();
  double pz = m_momentum.z();
  double main = E*E*m_error[3][3] + px*px*m_error[0][0] + 
              py*py*m_error[1][1] + pz*pz*m_error[2][2];
  double sub  = E  * (px*m_error[0][3]+py*m_error[1][3]+pz*m_error[2][3]) -
                px * (py*m_error[0][1]+pz*m_error[0][2]) -
                py * pz *m_error[1][2];
  double mass = this->mass();
  if(mass != 0.)
    return (main-2.*sub)/mass/mass;
  else
    return (main-2.*sub)*1.0e10; //temporary
}
   
void 
Momentum::momentum(const HepLorentzVector &p, const HepSymMatrix &dp)
{
  m_momentum = p;
  if(dp.num_row() == 4){
    m_error.sub(1,dp);
    //sets 0.0 in momentum position correlation elements.  
    for(unsigned i=0;i<4;i++){
      for(unsigned j=4;j<7;j++){
	m_error[i][j] = 0.0;
      }
    }
  }else{
    dout(Debugout::WARN,"Momentum") << "Error in momentum(const HepLorentzVector&," << std::endl;
    dout(Debugout::WARN,"Momentum") << "                  const HepSymMatrix&) of Mometum Class" << std::endl;
    dout(Debugout::WARN,"Momentum") << "                        ~~~~~~~~~~~~~--- 4 x 4 error matrix!!! " << std::endl;
    //sets 0.0 in all elements.  
    for(unsigned i=0;i<7;i++){
      for(unsigned j=i;j<7;j++){
	m_error[i][j] = 0.0;
      }
    }
  }
}
  
void 
Momentum::position(const HepPoint3D &x, const HepSymMatrix &dx)
{
  m_position = x;
  if(dx.num_row() == 3){
    m_error.sub(5,dx);
    //sets 0.0 in momentum position correlation elements.  
    for(unsigned i=0;i<4;i++){
      for(unsigned j=4;j<7;j++){
	m_error[i][j] = 0.0;
      }
    }
  }else{
    dout(Debugout::WARN,"Momentum") << "Error in position(const HepPoint3D&," << std::endl;
    dout(Debugout::WARN,"Momentum") << "                  const HepSymMatrix&) of Mometum Class" << std::endl;
    dout(Debugout::WARN,"Momentum") << "                        ~~~~~~~~~~~~~--- 3 x 3 error matrix!!! " << std::endl;
    //sets 0.0 in all elements.  
    for(unsigned i=0;i<7;i++){
      for(unsigned j=i;j<7;j++){
	m_error[i][j] = 0.0;
      }
    }
  }
}

void 
Momentum::momentumPosition(const HepLorentzVector &p, const HepPoint3D &x, 
			   const HepSymMatrix &dpx)
{
  m_momentum = p;
  m_position = x;
  if(dpx.num_row() == 7){
    m_error = dpx;
  }else{
    dout(Debugout::WARN,"Momentum") << "Error in momentumPosition(const HepLorentzVector&," << std::endl;
    dout(Debugout::WARN,"Momentum") << "                          const HepPoint3D&,      " << std::endl;
    dout(Debugout::WARN,"Momentum") << "                          const HepSymMatrix&) of Mometum Class" << std::endl;
    dout(Debugout::WARN,"Momentum") << "                                ~~~~~~~~~~~~~--- 7 x 7 error matrix!!! " << std::endl;
    //sets 0.0 in all elements.  
    for(unsigned i=0;i<7;i++){
      for(unsigned j=i;j<7;j++){
	m_error[i][j] = 0.0;
      }
    }
  }
}

HepPoint3D &
Momentum::vertex(const HepPoint3D   &vertex, 
		 const HepSymMatrix &error)
{
  m_vertex = vertex;
  if(error.num_row() == 3){
    m_vertexError = error;
  }else{
    dout(Debugout::WARN,"Momentum") << "Error in vertex(const HepPoint3D&," << std::endl;
    dout(Debugout::WARN,"Momentum") << "                const HepSymMatrix&) of Mometum Class" << std::endl;
    dout(Debugout::WARN,"Momentum") << "                      ~~~~~~~~~~~~~--- 3 x 3 error matrix!!! " << std::endl;
    //sets 0.0 in all elements.  
    for(unsigned i=0;i<3;i++){
      for(unsigned j=i;j<3;j++){
	m_vertexError[i][j] = 0.0;
      }
    }
  }
  return m_vertex;
}


HepPoint3D &
Momentum::decayVertex(const HepPoint3D   &vertex, 
		      const HepSymMatrix &error)
{
  m_decayVertex = vertex;
  if(error.num_row() == 3){
    m_decayVertexError = error;
  }else{
    dout(Debugout::WARN,"Momentum") << "Error in decayVertex(const HepPoint3D&," << std::endl;
    dout(Debugout::WARN,"Momentum") << "                     const HepSymMatrix&) of Mometum Class" << std::endl;
    dout(Debugout::WARN,"Momentum") << "                           ~~~~~~~~~~~~~--- 3 x 3 error matrix!!! " << std::endl;
    //sets 0.0 in all elements.  
    for(unsigned i=0;i<3;i++){
      for(unsigned j=i;j<3;j++){
	m_decayVertexError[i][j] = 0.0;
      }
    }
  }
  return m_decayVertex;
}


Momentum & 
Momentum::operator = (const Momentum &p)
{
  if(this == &p)return *this;

  m_momentum = p.m_momentum;
  m_position = p.m_position;
  m_error    = p.m_error;
  m_vertex   = p.m_vertex;
  m_vertexError      = p.m_vertexError;
  m_decayVertex      = p.m_decayVertex;
  m_decayVertexError = p.m_decayVertexError;
  return *this;
}

void
Momentum::dump(const std::string & keyword, const std::string & prefix) const {
    bool full = false;
    if (keyword.find("full") != std::string::npos) full = true;

    dout(Debugout::DUMP,"Momentum") << prefix;
    dout(Debugout::DUMP,"Momentum") << "Momentum:";
    if (full || keyword.find("mass")       != std::string::npos)dout(Debugout::DUMP,"Momentum") << " m=" << mass();
    if (full || keyword.find("momentum")   != std::string::npos)dout(Debugout::DUMP,"Momentum") << " p=" << p();
    if (full || keyword.find("position")   != std::string::npos)dout(Debugout::DUMP,"Momentum") << " x=" << x();
    if (full || keyword.find("production") != std::string::npos)dout(Debugout::DUMP,"Momentum") << " ProductionV=" << vertex();
    if (full || keyword.find("decay")      != std::string::npos)dout(Debugout::DUMP,"Momentum") << " decayV=" << decayVertex();
    if (full || keyword.find("return")     != std::string::npos)dout(Debugout::DUMP,"Momentum") << std::endl;
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
