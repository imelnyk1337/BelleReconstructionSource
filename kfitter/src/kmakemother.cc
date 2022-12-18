//
// $Id: kmakemother.cc 9932 2006-11-12 14:26:53Z katayama $
//
// $Log$
// Revision 1.11  2002/03/27 23:41:57  jtanaka
// Add new fitter to the mass-constraint fit.
//
// Revision 1.10  2002/02/23 17:56:53  katayama
// Added (int) for int/unsigned comparison
//
// Revision 1.9  2000/03/07 17:09:01  jtanaka
// bug fixed, use "similarity" etc.
//
// Revision 1.8  2000/03/07 10:52:39  katayama
// compatibility with CC5.0
//
// Revision 1.7  1999/03/29 05:39:15  jtanaka
// new class structure and new+old interfaces
//
// Revision 1.6  1998/09/28 13:25:18  jtanaka
// bug fix
//
// Revision 1.5  1998/09/08 11:51:49  jtanaka
// non-const --> const
//
// Revision 1.4  1998/07/29 00:00:01  katayama
// Less warnings with g++ -Wall
//
// Revision 1.3  1998/01/22 03:20:13  jtanaka
// Updated from Tanaka san. New Interface etc.
//
// Revision 1.2  1997/10/24 07:05:15  katayama
// Updated from Tanaka san. Bug fixes etc.
//
// Revision 1.1  1997/10/04 05:30:06  katayama
// New from Tanaka san
//
//
#include "belle.h"
#include "kfitter/kmakemother.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


kmakemother::kmakemother(void)
  : m_errorFlag(0),
    m_atDecayPoint(1),
    m_trackNum(0),
    m_errVertexFlag(0),
    m_correlationFlag(0),
    m_errVertexTrackFlag(0),
    m_beforeAfterFlag(KF_AFTER_FIT),
    m_magneticField(KF_MAGNETIC_FIELD),
    m_vertex(),
    m_errVertex(3,0),
    m_charge(0.)
{
}

kmakemother::~kmakemother(void)
{
  m_plist.erase(m_plist.begin(),m_plist.end());
  m_correlation.erase(m_correlation.begin(),m_correlation.end());
  m_errVertexTrack.erase(m_errVertexTrack.begin(),m_errVertexTrack.end());
}

void 
kmakemother::addTrack(const HepLorentzVector &p, 
		      const HepPoint3D &x,
		      const HepSymMatrix &e,
		      const double q,
		      const unsigned f)
{
  if(e.num_row() == (int)KF_NUM7){
    // (not need mass information)
    kfitterparticle pdata(p,x,e,q,f);
    m_plist.push_back(pdata);
    m_charge += q;
    m_trackNum = m_plist.size();
  }else m_errorFlag = 1;
}

void 
kmakemother::addTrack(const kfitterparticle &pdata)
{
  m_plist.push_back(pdata);
  m_charge += pdata.charge();
  m_trackNum = m_plist.size();
}
 
void 
kmakemother::magneticField(const double mf)
{
  m_magneticField = mf;
}

void 
kmakemother::vertex(const HepPoint3D &v)
{
  m_vertex = v;
}

void 
kmakemother::errVertex(const HepSymMatrix &e)
{
  if(e.num_row() == 3){
    m_errVertex = e;
    m_errVertexFlag = 1;
  }else m_errorFlag = 1;
}

void 
kmakemother::correlation(const HepMatrix &e)
{
  if(e.num_row() == (int)KF_NUM7 && 
     e.num_col() == (int)KF_NUM7){
    m_correlation.push_back(e);
    m_correlationFlag = 1;
  }else m_errorFlag = 1;
}
 
void 
kmakemother::correlation(void)
{
  m_correlation.push_back(HepMatrix(KF_NUM7,
				    KF_NUM7,0));
  m_correlationFlag = 1;
}

void 
kmakemother::errVertexTrack(const HepMatrix &e)
{
  if(e.num_row() == 3 && e.num_col() == (int)KF_NUM7){
    m_errVertexTrack.push_back(e);
    m_errVertexTrackFlag = 1;
  }else m_errorFlag = 1;
}

void 
kmakemother::errVertexTrack(void)
{
  m_errVertexTrack.push_back(HepMatrix(3,KF_NUM7,0));
  m_errVertexTrackFlag = 1;
}

void 
kmakemother::atDecayPoint(void)
{
  m_atDecayPoint = 1;
}

void 
kmakemother::notDecayPoint(void)
{
  m_atDecayPoint = 0;
}

void 
kmakemother::beforeAfter(const unsigned f)
{
  m_beforeAfterFlag = f;
}

kfitterparticle
kmakemother::track(void)
{
  return m_mother;
}

HepLorentzVector
kmakemother::momentum(void)
{
  return m_mother.momentum(KF_BEFORE_FIT);
}

HepPoint3D
kmakemother::position(void)
{
  return m_mother.position(KF_BEFORE_FIT);
}

HepSymMatrix
kmakemother::error(void)
{
  return m_mother.error(KF_BEFORE_FIT);
}

void
kmakemother::m_error(HepSymMatrix &Ec)
{
  //...error matrix of tracks
  unsigned i(0);
  for(std::vector<kfitterparticle>::iterator it = m_plist.begin(), endIt = m_plist.end();
      it != endIt; ++it){
    Ec.sub(i*KF_NUM7+1, it->error(m_beforeAfterFlag));
    ++i;
  }

  //...error matrix between tracks
  if(m_correlationFlag == 1){
    if( m_correlation.size() != static_cast<unsigned>((m_trackNum*m_trackNum-m_trackNum)/2) ){
      m_errorFlag = 1;
      return;
    }
    unsigned i(0);
    unsigned j(1);
    for(std::vector<HepMatrix>::iterator it = m_correlation.begin(), endIt = m_correlation.end();
	it != endIt; ++it){
      HepMatrix &hm = *it;
      for(unsigned k=0;k<KF_NUM7;++k){
	  for(unsigned l=0;l<KF_NUM7;++l){
	    Ec[k+i*KF_NUM7][l+j*KF_NUM7] = hm[k][l];
	  }
      }
      if(j != m_trackNum - 1){
	++j;
      }else if(i != m_trackNum - 2){
	++i;
	j = i+1;
      }else{
	break;
      }
    }
  }

  //...error of vertex
  if(m_errVertexFlag == 1){
    Ec.sub(m_trackNum*KF_NUM7+1, m_errVertex);
  }

  //...error matrix between vertex and tracks
  if(m_errVertexTrackFlag == 1){
    if( m_errVertexTrack.size() != m_trackNum ){
      m_errorFlag = 1;
      return;
    }
    unsigned i(0);
    for(std::vector<HepMatrix>::iterator it = m_errVertexTrack.begin(), endIt = m_errVertexTrack.end();
	it != endIt; ++it){
      HepMatrix &hm = *it;
      for(unsigned j=0;j<3;++j){
	for(unsigned k=0;k<KF_NUM7;++k){
	  Ec[j+m_trackNum*KF_NUM7][k+i*KF_NUM7] = hm[j][k];
	}
      }
      ++i;
    }
  }  
  return;
}

void 
kmakemother::m_delMdelC(HepMatrix &dMdC)
{
  //...local parameters
  double sum_a(0.);

  for(unsigned i=0;i<m_trackNum;++i){
    double a(-KF_PHOTON_VELOCITY*m_magneticField*m_plist[i].charge());
    sum_a += a;
    
    //...sets "a" in dMdC.
    dMdC[0][5+i*KF_NUM7] =  a;
    dMdC[1][4+i*KF_NUM7] = -a;

    //...sets "1" in dMdC.
    dMdC[0][0+i*KF_NUM7] = 1.;
    dMdC[1][1+i*KF_NUM7] = 1.;
    dMdC[2][2+i*KF_NUM7] = 1.;
    dMdC[3][3+i*KF_NUM7] = 1.;
  }
  
  //...sets "1" in dMdC.
  dMdC[4][0+m_trackNum*KF_NUM7] = 1.;
  dMdC[5][1+m_trackNum*KF_NUM7] = 1.;
  dMdC[6][2+m_trackNum*KF_NUM7] = 1.;

  //...sets "sum_a" in dMdC.
  dMdC[0][1+m_trackNum*KF_NUM7] = - sum_a;
  dMdC[1][0+m_trackNum*KF_NUM7] =   sum_a;
}

unsigned
kmakemother::make(void)
{
  //...makes matrix.
  HepMatrix dMdC(KF_NUM7,KF_NUM7*m_trackNum+3,0);
  HepSymMatrix Ec(KF_NUM7*m_trackNum+3,0);
  //...sets error matrix
  m_error(Ec);

  //...checks error about input data.
  if(m_errorFlag != 0)return m_errorFlag;

  //...makes delMdelC to calculate error matrix of mother particle.
  m_delMdelC(dMdC);

  //...calculates error matrix of mother particle.
  HepSymMatrix Em(Ec.similarity(dMdC));
  //...makes mother particle
  double px(0.), py(0.), pz(0.), E(0.);
  for(unsigned i=0;i<m_trackNum;++i){
    double a(0.), dx(0.), dy(0.);
    if(m_atDecayPoint == 1){
      a = -KF_PHOTON_VELOCITY*m_magneticField*m_plist[i].charge();
      dx = m_vertex.x() - m_plist[i].position(m_beforeAfterFlag).x();
      dy = m_vertex.y() - m_plist[i].position(m_beforeAfterFlag).y();
    }
    px += m_plist[i].momentum(m_beforeAfterFlag).x() - a*dy;
    py += m_plist[i].momentum(m_beforeAfterFlag).y() + a*dx;
    pz += m_plist[i].momentum(m_beforeAfterFlag).z();
    E  += m_plist[i].momentum(m_beforeAfterFlag).t();
  }
  //...momentum
  m_mother.momentum(HepLorentzVector(px,py,pz,E));
  //...position
  m_mother.position(m_vertex);
  //...error
  m_mother.error(Em);
  m_mother.charge(m_charge);
  //m_mother.mass(0.);

  return m_errorFlag;
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
