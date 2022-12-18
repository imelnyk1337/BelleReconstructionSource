//
// $Id: khelix2xyz.cc 9944 2006-11-29 07:36:07Z katayama $
//
// $Log$
// Revision 1.7  2000/04/13 12:35:57  katayama
// Added std:: to cout,cerr,endl etc.
//
// Revision 1.6  1999/03/29 05:39:14  jtanaka
// new class structure and new+old interfaces
//
// Revision 1.5  1998/10/18 12:55:55  jtanaka
// updates
//
// Revision 1.4  1998/08/24 14:02:58  jtanaka
// returns proper momentum etc...
//
// Revision 1.3  1998/06/30 21:06:51  jtanaka
// modification : this class uses helix class.
//
// Revision 1.2  1998/01/22 03:20:12  jtanaka
// Updated from Tanaka san. New Interface etc.
//
// Revision 1.1  1997/10/04 05:30:05  katayama
// New from Tanaka san
//
//
#include "belle.h"
#include "kfitter/khelix2xyz.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

helix2Xyz::helix2Xyz(void)
  : m_errorFlag(0)
{
}
//
// set helix
//
void 
helix2Xyz::helix(const double &h1,
		 const double &h2,
		 const double &h3,
		 const double &h4,
		 const double &h5)
{
  m_helix[0] = h1;
  m_helix[1] = h2;
  m_helix[2] = h3;
  m_helix[3] = h4;
  m_helix[4] = h5;
}
//
// set helix
//
void 
helix2Xyz::helix(const unsigned &i, const double &h)
{
  if(i <= 4)m_helix[i] = h;
  else dout(Debugout::ERR,"khelix2xyz") << "(khelix2xyz): Out of Range!" << std::endl;
}
//
// set helix
//
void 
helix2Xyz::helix(const double *h)
{
  for(unsigned int i=0;i<5;++i)
    m_helix[i] = h[i];
}
//
// set error of helix
//
void 
helix2Xyz::errHelix(const unsigned &i, const double &h)
{
  if(i <= 14){
    if(i == 0)m_errHelix[0][0] = h;
    else if(i ==  1)m_errHelix[1][0] = h;
    else if(i ==  2)m_errHelix[1][1] = h;
    else if(i ==  3)m_errHelix[2][0] = h;
    else if(i ==  4)m_errHelix[2][1] = h;
    else if(i ==  5)m_errHelix[2][2] = h;
    else if(i ==  6)m_errHelix[3][0] = h;
    else if(i ==  7)m_errHelix[3][1] = h;
    else if(i ==  8)m_errHelix[3][2] = h;
    else if(i ==  9)m_errHelix[3][3] = h;
    else if(i == 10)m_errHelix[4][0] = h;
    else if(i == 11)m_errHelix[4][1] = h;
    else if(i == 12)m_errHelix[4][2] = h;
    else if(i == 13)m_errHelix[4][3] = h;
    else m_errHelix[4][4] = h;
  }else{
    dout(Debugout::ERR,"khelix2xyz") << "(khelix2xyz): Out of Range!" <<std::endl;
  }
}
//
// set error of helix
//
void 
helix2Xyz::errHelix(const double *h)
{
  m_errHelix[0][0] = h[0];
  for(unsigned i=0;i<2;++i)
    m_errHelix[1][i] = h[1+i];
  for(unsigned i=0;i<3;++i)
    m_errHelix[2][i] = h[3+i];
  for(unsigned i=0;i<4;++i)
    m_errHelix[3][i] = h[6+i];
  for(unsigned i=0;i<5;++i)
    m_errHelix[4][i] = h[10+i];
}
//
// set mass
//
void 
helix2Xyz::mass(const double &m)
{
  m_mass = m;
}
//
// set pivot
//
void 
helix2Xyz::pivot(const HepPoint3D &p)
{
  m_pivot = p;
}
//
// get momentum
//
HepLorentzVector 
helix2Xyz::momentum(void)
{
  return m_momentum;
}
//
// get position
//
HepPoint3D
helix2Xyz::position(void)
{
  return m_position;
}
//
// get error of momentum and position
//
HepSymMatrix 
helix2Xyz::errMomentumPosition(void)
{
  return m_errMomentumPosition;
}
//
// get momentum and position
//
double 
helix2Xyz::xyz(const unsigned &i)
{
  if(i <= 6){
    if(i == 0)return m_momentum.x();
    else if(i == 1)return m_momentum.y();
    else if(i == 2)return m_momentum.z();
    else if(i == 3)return m_momentum.t();
    else if(i == 4)return m_position.x();
    else if(i == 5)return m_position.y();
    else return m_position.z();
  }else{
    dout(Debugout::ERR,"khelix2xyz") << "(khelix2xyz): Out of Range!" << std::endl;
    return 0.;
  }
}
//
// get error of momentum and postion
//
double 
helix2Xyz::errXyz(const unsigned &i,
		  const unsigned &j)
{
  if(i <= 6 && j <= 6){
    return m_errMomentumPosition[i][j];
  }else{
    dout(Debugout::ERR,"khelix2xyz") << "(khelix2xyz): Out of Range!" << std::endl;
    return 0.;
  }
}
//
// update
//
unsigned
helix2Xyz::update(void)
{
  //...returns momentum etc.. which are near the origin(0, 0, 0)
  Helix helixclass(m_pivot, m_helix, m_errHelix);
  
  HepPoint3D ori(0.,0.,0.);
  if(m_pivot.x() != 0. ||
     m_pivot.y() != 0. ||
     m_pivot.z() != 0.)helixclass.pivot(ori);
  
  m_momentum = 
    helixclass.momentum(0., m_mass, m_position, m_errMomentumPosition);

  return m_errorFlag;
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
