//
// $Id: kfitterparticle.cc 9944 2006-11-29 07:36:07Z katayama $
//
// $Log$
// Revision 1.9  2002/06/11 04:46:46  jtanaka
// Bug fix : kfitterparticle::momentum() (mass was not updated.) from Kai-Feng Chen-san.
//
// Revision 1.8  2002/03/27 23:41:57  jtanaka
// Add new fitter to the mass-constraint fit.
//
// Revision 1.7  2002/02/23 17:56:53  katayama
// Added (int) for int/unsigned comparison
//
// Revision 1.6  2000/04/13 12:35:56  katayama
// Added std:: to cout,cerr,endl etc.
//
// Revision 1.5  2000/03/07 17:09:00  jtanaka
// bug fixed, use "similarity" etc.
//
// Revision 1.4  1999/03/29 05:39:14  jtanaka
// new class structure and new+old interfaces
//
// Revision 1.3  1998/09/08 11:51:48  jtanaka
// non-const --> const
//
// Revision 1.2  1998/07/23 12:51:05  katayama
// Removed control-M (probably added from a PC)
//
// Revision 1.1  1998/01/22 03:20:12  jtanaka
// Updated from Tanaka san. New Interface etc.
//
//
#include "belle.h"
#include "kfitter/kfitterparticle.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


kfitterparticle::kfitterparticle(void)
  : m_momentum_b(),
    m_position_b(),
    m_error_b(KF_NUM7,0),
    m_momentum_a(),
    m_position_a(),
    m_error_a(KF_NUM7,0),
    m_charge(0.),
    m_mass(0.),
    m_vertex(),
    m_errVertex(3,0)
{
}

kfitterparticle::kfitterparticle(const HepLorentzVector &p,
				 const HepPoint3D &x,
				 const HepSymMatrix &e,
				 const double q,
				 const unsigned f)
  : m_charge(q),
    m_mass(p.mag()),
    m_vertex(),
    m_errVertex(3,0)
{
  if(f == KF_BEFORE_FIT){
    m_momentum_b = p;
    m_position_b = x;
    if(e.num_row() == (int)KF_NUM7){
      m_error_b = e;
    }else dout(Debugout::ERR,"kfitterparticle") << "(kfitterparticle): Out of Range!" << std::endl;
    m_momentum_a = HepLorentzVector();
    m_position_a = HepPoint3D();
    m_error_a = HepSymMatrix(KF_NUM7,0);
  }else if(f == KF_AFTER_FIT){
    m_momentum_b = HepLorentzVector();
    m_position_b = HepPoint3D();
    m_error_b = HepSymMatrix(KF_NUM7,0);
    m_momentum_a = p;
    m_position_a = x;
    if(e.num_row() == (int)KF_NUM7){
      m_error_a = e;
    }else dout(Debugout::ERR,"kfitterparticle") << "(kfitterparticle): Out of Range!" << std::endl;
  }
}

#if KF_WITH_OLD_INTERFACE
kfitterparticle::kfitterparticle(const HepLorentzVector &p,
				 const HepPoint3D &x,
				 const HepSymMatrix &e,
				 const double q,
				 const double m,
				 const unsigned f)
  : m_charge(q),
    m_mass(p.mag()),
    m_vertex(),
    m_errVertex(3,0)
{
  if(f == KF_BEFORE_FIT){
    m_momentum_b = p;
    m_position_b = x;
    if(e.num_row() == (int)KF_SUBPARAMETER_NUMBER){
      m_error_b = e;
    }else dout(Debugout::ERR,"kfitterparticle") << "(kfitterparticle): Out of Range!" << std::endl;
    m_momentum_a = HepLorentzVector();
    m_position_a = HepPoint3D();
    m_error_a = HepSymMatrix(KF_SUBPARAMETER_NUMBER,0);
  }else if(f == KF_AFTER_FIT){
    m_momentum_b = HepLorentzVector();
    m_position_b = HepPoint3D();
    m_error_b = HepSymMatrix(KF_SUBPARAMETER_NUMBER,0);
    m_momentum_a = p;
    m_position_a = x;
    if(e.num_row() == (int)KF_SUBPARAMETER_NUMBER){
      m_error_a = e;
    }else dout(Debugout::ERR,"kfitterparticle") << "(kfitterparticle): Out of Range!" << std::endl;
  }
}
#endif

kfitterparticle::kfitterparticle(const kfitterparticle &p)
  : m_momentum_b(p.m_momentum_b),
    m_position_b(p.m_position_b),
    m_error_b(p.m_error_b),
    m_momentum_a(p.m_momentum_a),
    m_position_a(p.m_position_a),
    m_error_a(p.m_error_a),
    m_charge(p.m_charge),
    m_mass(p.m_mass),
    m_vertex(p.m_vertex),
    m_errVertex(p.m_errVertex)
{
}

void
kfitterparticle::momentum(const HepLorentzVector &p,
			  const unsigned n)
{
  m_mass = p.mag();
  if(n == KF_BEFORE_FIT)m_momentum_b = p;
  else if(n == KF_AFTER_FIT)m_momentum_a = p;
}

void
kfitterparticle::position(const HepPoint3D &x,
			  const unsigned n)
{
  if(n == KF_BEFORE_FIT)m_position_b = x;
  else if(n == KF_AFTER_FIT)m_position_a = x;
}

void
kfitterparticle::error(const HepSymMatrix &e,
		       const unsigned n)
{
  if(e.num_row() == (int)KF_NUM7){
    if(n == KF_BEFORE_FIT)m_error_b = e;
    else if(n == KF_AFTER_FIT)m_error_a = e;
  }else dout(Debugout::ERR,"kfitterparticle") << "(kfitterparticle): Out of Range!" << std::endl;
}

void
kfitterparticle::charge(const double q)
{
  m_charge = q;
}

//void
//kfitterparticle::mass(const double m)
//{
//  m_mass = m;
//}

void
kfitterparticle::vertex(const HepPoint3D &v)
{
  m_vertex = v;
}

void
kfitterparticle::errVertex(const HepSymMatrix &e)
{
  if(e.num_row() == 3)m_errVertex = e;
  else dout(Debugout::ERR,"kfitterparticle") << "(kfitterparticle): Out of Range!" << std::endl;
}

HepLorentzVector 
kfitterparticle::momentum(const unsigned n)
{
  if(n == KF_AFTER_FIT)return m_momentum_a;
  else if(n == KF_BEFORE_FIT)return m_momentum_b;
  dout(Debugout::ERR,"kfitterparticle") << "(kfitterparticle): Out of Range!" << std::endl;
  return HepLorentzVector();
}

HepLorentzVector 
kfitterparticle::momentum(const unsigned n) const
{
  if(n == KF_AFTER_FIT)return m_momentum_a;
  else if(n == KF_BEFORE_FIT)return m_momentum_b;
  dout(Debugout::ERR,"kfitterparticle") << "(kfitterparticle): Out of Range!" << std::endl;
  return HepLorentzVector();
}

HepPoint3D
kfitterparticle::position(const unsigned n)
{
  if(n == KF_AFTER_FIT)return m_position_a;
  else if(n == KF_BEFORE_FIT)return m_position_b;
  dout(Debugout::ERR,"kfitterparticle") << "(kfitterparticle): Out of Range!" << std::endl;
  return HepPoint3D();
}

HepPoint3D
kfitterparticle::position(const unsigned n) const
{
  if(n == KF_AFTER_FIT)return m_position_a;
  else if(n == KF_BEFORE_FIT)return m_position_b;
  dout(Debugout::ERR,"kfitterparticle") << "(kfitterparticle): Out of Range!" << std::endl;
  return HepPoint3D();
}

HepSymMatrix
kfitterparticle::error(const unsigned n)
{
  if(n == KF_AFTER_FIT)return m_error_a;
  else if(n == KF_BEFORE_FIT)return m_error_b;
  dout(Debugout::ERR,"kfitterparticle") << "(kfitterparticle): Out of Range!" << std::endl;
  return HepSymMatrix(KF_NUM7,0);
}

HepSymMatrix
kfitterparticle::error(const unsigned n) const
{
  if(n == KF_AFTER_FIT)return m_error_a;
  else if(n == KF_BEFORE_FIT)return m_error_b;
  dout(Debugout::ERR,"kfitterparticle") << "(kfitterparticle): Out of Range!" << std::endl;
  return HepSymMatrix(KF_NUM7,0);
}

double
kfitterparticle::charge(void)
{
  return m_charge;
}

double
kfitterparticle::charge(void) const
{
  return m_charge;
}

double
kfitterparticle::mass(void)
{
  return m_mass;
}

double
kfitterparticle::mass(void) const
{
  return m_mass;
}

HepPoint3D
kfitterparticle::vertex(void)
{
  return m_vertex;
}

HepPoint3D
kfitterparticle::vertex(void) const
{
  return m_vertex;
}

HepSymMatrix
kfitterparticle::errVertex(void)
{
  return m_errVertex;
}

HepSymMatrix
kfitterparticle::errVertex(void) const
{
  return m_errVertex;
}

double
kfitterparticle::getFitParameter(const unsigned n,
				 const unsigned f)
{
  if(f == KF_BEFORE_FIT){
    if(n == 0){
      return m_momentum_b.x();
    }else if(n == 1){
      return m_momentum_b.y();
    }else if(n == 2){
      return m_momentum_b.z();
    }else if(n == 3){
      return m_position_b.x();
    }else if(n == 4){
      return m_position_b.y();
    }else if(n == 5){
      return m_position_b.z();
    }else{
      dout(Debugout::ERR,"kfitterparticle") << "(kfitterparticle): Out of Range!" << std::endl;
      return 0;
    }
  }else if(f == KF_AFTER_FIT){
    if(n == 0){
      return m_momentum_a.x();
    }else if(n == 1){
      return m_momentum_a.y();
    }else if(n == 2){
      return m_momentum_a.z();
    }else if(n == 3){
      return m_position_a.x();
    }else if(n == 4){
      return m_position_a.y();
    }else if(n == 5){
      return m_position_a.z();
    }else{
      dout(Debugout::ERR,"kfitterparticle") << "(kfitterparticle): Out of Range!" << std::endl;
      return 0;
    }
  }
  dout(Debugout::ERR,"kfitterparticle") << "(kfitterparticle): Out of Range!" << std::endl;
  return 0;
}

double
kfitterparticle::getFitParameter(const unsigned n,
				 const unsigned f) const
{
  if(f == KF_BEFORE_FIT){
    if(n == 0){
      return m_momentum_b.x();
    }else if(n == 1){
      return m_momentum_b.y();
    }else if(n == 2){
      return m_momentum_b.z();
    }else if(n == 3){
      return m_position_b.x();
    }else if(n == 4){
      return m_position_b.y();
    }else if(n == 5){
      return m_position_b.z();
    }else{
      dout(Debugout::ERR,"kfitterparticle") << "(kfitterparticle): Out of Range!" << std::endl;
      return 0;
    }
  }else if(f == KF_AFTER_FIT){
    if(n == 0){
      return m_momentum_a.x();
    }else if(n == 1){
      return m_momentum_a.y();
    }else if(n == 2){
      return m_momentum_a.z();
    }else if(n == 3){
      return m_position_a.x();
    }else if(n == 4){
      return m_position_a.y();
    }else if(n == 5){
      return m_position_a.z();
    }else{
      dout(Debugout::ERR,"kfitterparticle") << "(kfitterparticle): Out of Range!" << std::endl;
      return 0;
    }
  }
  dout(Debugout::ERR,"kfitterparticle") << "(kfitterparticle): Out of Range!" << std::endl;
  return 0;
}

HepMatrix
kfitterparticle::getFitParameter(const unsigned f)
{
  HepMatrix a(KF_NUM6,1,0);
  if(f == KF_BEFORE_FIT){
    a[0][0] = m_momentum_b.x();
    a[1][0] = m_momentum_b.y();
    a[2][0] = m_momentum_b.z();
    a[3][0] = m_position_b.x();
    a[4][0] = m_position_b.y();
    a[5][0] = m_position_b.z();
  }else if(f == KF_AFTER_FIT){
    a[0][0] = m_momentum_a.x();
    a[1][0] = m_momentum_a.y();
    a[2][0] = m_momentum_a.z();
    a[3][0] = m_position_a.x();
    a[4][0] = m_position_a.y();
    a[5][0] = m_position_a.z();
  }
  return a;
}

HepMatrix
kfitterparticle::getFitParameter(const unsigned f) const
{
  HepMatrix a(KF_NUM6,1,0);
  if(f == KF_BEFORE_FIT){
    a[0][0] = m_momentum_b.x();
    a[1][0] = m_momentum_b.y();
    a[2][0] = m_momentum_b.z();
    a[3][0] = m_position_b.x();
    a[4][0] = m_position_b.y();
    a[5][0] = m_position_b.z();
  }else if(f == KF_AFTER_FIT){
    a[0][0] = m_momentum_a.x();
    a[1][0] = m_momentum_a.y();
    a[2][0] = m_momentum_a.z();
    a[3][0] = m_position_a.x();
    a[4][0] = m_position_a.y();
    a[5][0] = m_position_a.z();
  }
  return a;
}

HepSymMatrix 
kfitterparticle::getFitError(const unsigned f)
{
  HepSymMatrix err(KF_NUM6,0);
  if(f == KF_BEFORE_FIT){
    for(unsigned i=0;i<3;++i){
      for(unsigned j=i;j<3;++j){
	err[i][j]     = m_error_b[i][j];
	err[3+i][3+j] = m_error_b[4+i][4+j];
      }
    }
    for(unsigned i=0;i<3;++i){
      for(unsigned j=0;j<3;++j){
	err[i][3+j]   = m_error_b[i][4+j];
      }
    }
  }else if(f == KF_AFTER_FIT){
    for(unsigned i=0;i<3;++i){
      for(unsigned j=i;j<3;++j){
	err[i][j]     = m_error_a[i][j];
	err[3+i][3+j] = m_error_a[4+i][4+j];
      }
    }
    for(unsigned i=0;i<3;++i){
      for(unsigned j=0;j<3;++j){
	err[i][3+j]   = m_error_a[i][4+j];
      }
    }
  }
  return err;
}

HepSymMatrix 
kfitterparticle::getFitError(const unsigned f) const
{
  HepSymMatrix err(KF_NUM6,0);
  if(f == KF_BEFORE_FIT){
    for(unsigned i=0;i<3;++i){
      for(unsigned j=i;j<3;++j){
	err[i][j]     = m_error_b[i][j];
	err[3+i][3+j] = m_error_b[4+i][4+j];
      }
    }
    for(unsigned i=0;i<3;++i){
      for(unsigned j=0;j<3;++j){
	err[i][3+j]   = m_error_b[i][4+j];
      }
    }
  }else if(f == KF_AFTER_FIT){
    for(unsigned i=0;i<3;++i){
      for(unsigned j=i;j<3;++j){
	err[i][j]     = m_error_a[i][j];
	err[3+i][3+j] = m_error_a[4+i][4+j];
      }
    }
    for(unsigned i=0;i<3;++i){
      for(unsigned j=0;j<3;++j){
	err[i][3+j]   = m_error_a[i][4+j];
      }
    }
  }
  return err;
}

kfitterparticle & 
kfitterparticle::operator=(const kfitterparticle &a)
{
  if(this == &a)return *this;
  
  m_momentum_b = a.m_momentum_b;
  m_position_b = a.m_position_b;
  m_error_b    = a.m_error_b;
  
  m_momentum_a = a.m_momentum_a;
  m_position_a = a.m_position_a;
  m_error_a    = a.m_error_a;
  
  m_charge = a.m_charge;
  m_mass   = a.m_mass;

  m_vertex = a.m_vertex;
  m_errVertex = a.m_errVertex;
  return *this;
}

HepMatrix
kfitterparticle::mompos(const unsigned f)
{
  HepMatrix a(KF_NUM7,1,0);
  if(f == KF_BEFORE_FIT){
    a[0][0] = m_momentum_b.x();
    a[1][0] = m_momentum_b.y();
    a[2][0] = m_momentum_b.z();
    a[3][0] = m_momentum_b.t();
    a[4][0] = m_position_b.x();
    a[5][0] = m_position_b.y();
    a[6][0] = m_position_b.z();
  }else if(f == KF_AFTER_FIT){
    a[0][0] = m_momentum_a.x();
    a[1][0] = m_momentum_a.y();
    a[2][0] = m_momentum_a.z();
    a[3][0] = m_momentum_a.t();
    a[4][0] = m_position_a.x();
    a[5][0] = m_position_a.y();
    a[6][0] = m_position_a.z();
  }
  return a;
}

HepMatrix
kfitterparticle::mompos(const unsigned f) const
{
  HepMatrix a(KF_NUM7,1,0);
  if(f == KF_BEFORE_FIT){
    a[0][0] = m_momentum_b.x();
    a[1][0] = m_momentum_b.y();
    a[2][0] = m_momentum_b.z();
    a[3][0] = m_momentum_b.t();
    a[4][0] = m_position_b.x();
    a[5][0] = m_position_b.y();
    a[6][0] = m_position_b.z();
  }else if(f == KF_AFTER_FIT){
    a[0][0] = m_momentum_a.x();
    a[1][0] = m_momentum_a.y();
    a[2][0] = m_momentum_a.z();
    a[3][0] = m_momentum_a.t();
    a[4][0] = m_position_a.x();
    a[5][0] = m_position_a.y();
    a[6][0] = m_position_a.z();
  }
  return a;
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
