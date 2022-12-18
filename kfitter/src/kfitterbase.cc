#include "kfitter/kfitterbase.h"

#if KF_WITH_BELLE_INTERFACE
#include "belle.h"
#include MDST_H
#include TRK_H
#include "helix/Helix.h"
#endif
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

kfitterbase::kfitterbase(void)
: m_magneticField(KF_MAGNETIC_FIELD),
  m_dgf(0),
  m_chisq(-1.),
  m_cl(0.),
  m_trackNum(0),
  m_errorFlag(KF_NO_ERROR),
  m_correlationFlag(0),
  m_overIterationFlag(KF_NO_OVER_ITERATION),
  m_errorMsgFlag(1)
{
}

kfitterbase::~kfitterbase(void){
  m_plist.erase(m_plist.begin(),m_plist.end());
  m_correlation_b.erase(m_correlation_b.begin(),m_correlation_b.end());
}

void 
kfitterbase::addTrack(const HepLorentzVector &p, 
		      const HepPoint3D &x,
		      const HepSymMatrix &e, 
		      const double q){
  m_plist.push_back(kfitterparticle(p, x, e, q, p.mag()));
  m_trackNum = m_plist.size();
  if(e.num_row() != (int)KF_NUM7)m_errorFlag = KF_INPUT_MATRIX_SIZE;
}

#if KF_WITH_OLD_INTERFACE
void 
kfitterbase::addTrack(const HepLorentzVector &p, 
		      const HepPoint3D &x,
		      const HepSymMatrix &e, 
		      const double q, 
		      const double m){
  addTrack(p,x,e,q);
  //m_plist.push_back(kfitterparticle(p, x, e, q, m));
  //m_trackNum = m_plist.size();
  //if(e.num_row() != KF_NUM7)m_errorFlag = KF_INPUT_MATRIX_SIZE;
}
#endif

void 
kfitterbase::addTrack(const kfitterparticle &p){
  m_plist.push_back(p);
  m_trackNum = m_plist.size();
}

#if KF_WITH_BELLE_INTERFACE
void
kfitterbase::addTrack(const Mdst_charged &p,
		      const double mass,
		      const unsigned id){
  const HepPoint3D pivot(p.trk().mhyp(id).pivot_x(),p.trk().mhyp(id).pivot_y(),p.trk().mhyp(id).pivot_z());
  HepMatrix  tmp_a(5,1);
  tmp_a[0][0] = p.trk().mhyp(id).helix(0);
  tmp_a[1][0] = p.trk().mhyp(id).helix(1);
  tmp_a[2][0] = p.trk().mhyp(id).helix(2);
  tmp_a[3][0] = p.trk().mhyp(id).helix(3);
  tmp_a[4][0] = p.trk().mhyp(id).helix(4);
  HepVector  a(tmp_a);
  HepSymMatrix Ea(5,0);
  Ea[0][0] = p.trk().mhyp(id).error(0);
  Ea[1][0] = p.trk().mhyp(id).error(1);
  Ea[1][1] = p.trk().mhyp(id).error(2);
  Ea[2][0] = p.trk().mhyp(id).error(3);
  Ea[2][1] = p.trk().mhyp(id).error(4);
  Ea[2][2] = p.trk().mhyp(id).error(5);
  Ea[3][0] = p.trk().mhyp(id).error(6);
  Ea[3][1] = p.trk().mhyp(id).error(7);
  Ea[3][2] = p.trk().mhyp(id).error(8);
  Ea[3][3] = p.trk().mhyp(id).error(9);
  Ea[4][0] = p.trk().mhyp(id).error(10);
  Ea[4][1] = p.trk().mhyp(id).error(11);
  Ea[4][2] = p.trk().mhyp(id).error(12);
  Ea[4][3] = p.trk().mhyp(id).error(13);
  Ea[4][4] = p.trk().mhyp(id).error(14);
  Helix helix(pivot, a, Ea);
  HepLorentzVector momentum;
  HepPoint3D position;
  HepSymMatrix error(7,0);
  //...finds ref. point.
  HepPoint3D tmp;
  if(pivot.x() == 0. &&
     pivot.y() == 0. &&
     pivot.z() == 0.)goto easy;
  helix.pivot(tmp);
  momentum = helix.momentum(0.,mass,position,error);
  m_plist.push_back(kfitterparticle(momentum, position, error, p.charge(), mass));
  m_trackNum = m_plist.size();
  return;

 easy:
  momentum = helix.momentum(0.,mass,position,error);
  
  m_plist.push_back(kfitterparticle(momentum, position, error, p.charge(), mass));
  m_trackNum = m_plist.size();
}

void
kfitterbase::addTrack(const Rectrk &p,
		      const double mass,
		      const unsigned id){

  const HepPoint3D pivot(p.zero(id).pivot(0),p.zero(id).pivot(1),p.zero(id).pivot(2));
  HepMatrix  tmp_a(5,1);
  tmp_a[0][0] = p.zero(id).helix(0);
  tmp_a[1][0] = p.zero(id).helix(1);
  tmp_a[2][0] = p.zero(id).helix(2);
  tmp_a[3][0] = p.zero(id).helix(3);
  tmp_a[4][0] = p.zero(id).helix(4);
  HepVector  a(tmp_a);
  HepSymMatrix Ea(5,0);
  Ea[0][0] = p.zero(id).error(0);
  Ea[1][0] = p.zero(id).error(1);
  Ea[1][1] = p.zero(id).error(2);
  Ea[2][0] = p.zero(id).error(3);
  Ea[2][1] = p.zero(id).error(4);
  Ea[2][2] = p.zero(id).error(5);
  Ea[3][0] = p.zero(id).error(6);
  Ea[3][1] = p.zero(id).error(7);
  Ea[3][2] = p.zero(id).error(8);
  Ea[3][3] = p.zero(id).error(9);
  Ea[4][0] = p.zero(id).error(10);
  Ea[4][1] = p.zero(id).error(11);
  Ea[4][2] = p.zero(id).error(12);
  Ea[4][3] = p.zero(id).error(13);
  Ea[4][4] = p.zero(id).error(14);

  Helix helix(pivot, a, Ea);
  HepLorentzVector momentum;
  HepPoint3D position;
  HepSymMatrix error(7,0);
  //...finds ref. point.
  HepPoint3D tmp;
  if(pivot.x() == 0. &&
     pivot.y() == 0. &&
     pivot.z() == 0.)goto easy;
  helix.pivot(tmp);
  momentum = helix.momentum(0.,mass,position,error);
  m_plist.push_back(kfitterparticle(momentum, position, error, a[2]>0 ? 1 : -1, mass));
  m_trackNum = m_plist.size();
  return;

 easy:
  momentum = helix.momentum(0.,mass,position,error);
  
  m_plist.push_back(kfitterparticle(momentum, position, error, a[2]>0 ? 1 : -1, mass));
  m_trackNum = m_plist.size();
}
#endif

void 
kfitterbase::magneticField(const double mf){
  m_magneticField = mf;
}

void 
kfitterbase::correlation(const HepMatrix &e){
  m_correlation_b.push_back(e);
  m_correlationFlag = 1;
  if(e.num_row() != (int)KF_NUM7)
    m_errorFlag = KF_INPUT_MATRIX_SIZE;
}

void 
kfitterbase::correlation(void){
  m_correlation_b.push_back(HepMatrix(KF_NUM7,
				      KF_NUM7,0));
  m_correlationFlag = 1;
}

unsigned
kfitterbase::error(void){
  return m_errorFlag;
}

unsigned
kfitterbase::tracks(void){
  return m_trackNum;
}

unsigned
kfitterbase::iteration(void){
  return m_overIterationFlag;
}

unsigned
kfitterbase::dgf(void){
  return m_dgf;
}

double
kfitterbase::chisq(void){
  return m_chisq;
}
 
double
kfitterbase::cl(void){
  return m_cl;
}

double
kfitterbase::chisq(const unsigned n){
  if(n < m_trackNum){
    HepMatrix da(m_plist[n].getFitParameter(KF_BEFORE_FIT)-
		 m_plist[n].getFitParameter(KF_AFTER_FIT));
    int errInverse(0);
    double chiSq((da.T()*(m_plist[n].getFitError(KF_BEFORE_FIT).inverse(errInverse))*da)[0][0]);
    if(errInverse){
      m_errorFlag = KF_OUTPUT_INVERSE;
      return -1.;
    }else return chiSq;
  }else{
    m_errorFlag = KF_OUTPUT_TRACK_NUM;
    if(m_errorMsgFlag)dout(Debugout::ERR,"kfitterbase") << "(kfitterbase::chisq): Out of Range!" << std::endl;
    return -1.;
  }
}

double
kfitterbase::magneticField(void){
  return m_magneticField;
}

void 
kfitterbase::errorMsg(unsigned int f){
  m_errorMsgFlag = f;
}

HepLorentzVector
kfitterbase::momentum(const unsigned n){
  if(n < m_trackNum){
    return m_plist[n].momentum();
  }else{
    m_errorFlag = KF_OUTPUT_TRACK_NUM;
    if(m_errorMsgFlag)dout(Debugout::ERR,"kfitterbase") << "(kfitterbase::momentum): Out of Range!" << std::endl;
    return HepLorentzVector();
  }
}

HepPoint3D
kfitterbase::position(const unsigned n){
  if(n < m_trackNum){
    return m_plist[n].position();
  }else{
    m_errorFlag = KF_OUTPUT_TRACK_NUM;
    if(m_errorMsgFlag)dout(Debugout::ERR,"kfitterbase") << "(kfitterbase::position): Out of Range!" << std::endl;
    return HepPoint3D();
  }
}

HepSymMatrix
kfitterbase::error(const unsigned n){
  if(n < m_trackNum){
    return m_plist[n].error();
  }else{
    m_errorFlag = KF_OUTPUT_TRACK_NUM;
    if(m_errorMsgFlag)dout(Debugout::ERR,"kfitterbase") << "(kfitterbase::error): Out of Range!" << std::endl;
    return HepSymMatrix(KF_NUM7,0);
  }
}

kfitterparticle
kfitterbase::track(const unsigned n){
  if(n < m_trackNum){
    return m_plist[n];
  }else{
    m_errorFlag = KF_OUTPUT_TRACK_NUM;
    if(m_errorMsgFlag)dout(Debugout::ERR,"kfitterbase") << "(kfitterbase::track): Out of Range!" << std::endl;
    return kfitterparticle();
  }
}

HepMatrix
kfitterbase::correlation(const unsigned n, const unsigned m, const unsigned flag){
  if(flag == KF_AFTER_FIT){
    if(n < m_trackNum && m < m_trackNum){
      return m_makeError(momentum(n),
			 momentum(m),
			 m_V_al_1.sub(KF_NUM6*n+1,KF_NUM6*(n+1),
				      KF_NUM6*m+1,KF_NUM6*(m+1)));
    }else{
      m_errorFlag = KF_OUTPUT_TRACK_NUM;
      if(m_errorMsgFlag)dout(Debugout::ERR,"kfitterbase") << "(kfitterbase::correlation,1): Out of Range!" << std::endl;
      return HepMatrix(KF_NUM7,KF_NUM7,0);
    }
  }else{
    if(n < m_trackNum && m < m_trackNum){
      if(n == m){
	return static_cast<HepMatrix>(m_plist[n].error(KF_BEFORE_FIT));
      }else{
	// 2000/03/07
	unsigned n2(n), m2(m);
	if(n2 > m2){
	  unsigned tmp(n2);
	  n2 = m2;
	  m2 = tmp;
	}
	unsigned index(0);
	for(unsigned i=0;i<n2;++i)index += m_trackNum-1-i;
	index -= n2+1;
	index += m2;
	if(n == n2)return m_correlation_b[index+m2];
	else return m_correlation_b[index+m2].T();
      }
    }else{
      m_errorFlag = KF_OUTPUT_TRACK_NUM;
      if(m_errorMsgFlag)dout(Debugout::ERR,"kfitterbase") << "(kfitterbase::correlation,2): Out of Range!" << std::endl;
      return HepMatrix(KF_NUM7,KF_NUM7,0);
    }
  }
}

HepSymMatrix 
kfitterbase::m_makeError(const HepLorentzVector &p,
			 const HepMatrix &e){
  // self track
  // Error(6x6,e) ==> Error(7x7,output(hsm)) using Momentum(p).
  HepSymMatrix hsm(KF_NUM7,0);
  if(p.t() == 0.){
    m_errorFlag = KF_OUTPUT_DIV_ZERO;
    if(m_errorMsgFlag)dout(Debugout::ERR,"kfitterbase") << "(kfitterbase::m_makeError,1): Divide by ZERO!" << std::endl;
    return hsm;
  }
  for(unsigned i=0;i<3;++i){
    for(unsigned j=i;j<3;++j){
      hsm[i][j]     = e[i][j];
      hsm[4+i][4+j] = e[3+i][3+j];
    }
  }
  for(unsigned i=0;i<3;++i){
    for(unsigned j=0;j<3;++j){
      hsm[i][4+j] = e[i][3+j];
    }
  }
  double invE = 1./p.t();
  hsm[0][3] = (p.x()*hsm[0][0]+p.y()*hsm[0][1]+p.z()*hsm[0][2])*invE;
  hsm[1][3] = (p.x()*hsm[0][1]+p.y()*hsm[1][1]+p.z()*hsm[1][2])*invE;
  hsm[2][3] = (p.x()*hsm[0][2]+p.y()*hsm[1][2]+p.z()*hsm[2][2])*invE;
  hsm[3][3] = (p.x()*p.x()*hsm[0][0]+p.y()*p.y()*hsm[1][1]+p.z()*p.z()*hsm[2][2]
               +2.0*p.x()*p.y()*hsm[0][1]
               +2.0*p.x()*p.z()*hsm[0][2]
               +2.0*p.y()*p.z()*hsm[1][2])*invE*invE;
  hsm[3][4] = (p.x()*hsm[0][4]+p.y()*hsm[1][4]+p.z()*hsm[2][4])*invE;
  hsm[3][5] = (p.x()*hsm[0][5]+p.y()*hsm[1][5]+p.z()*hsm[2][5])*invE;
  hsm[3][6] = (p.x()*hsm[0][6]+p.y()*hsm[1][6]+p.z()*hsm[2][6])*invE;
  return hsm;
}

HepMatrix 
kfitterbase::m_makeError(const HepLorentzVector &p1,
			 const HepLorentzVector &p2,
			 const HepMatrix &e){
  // track and track
  // Error(6x6,e) ==> Error(7x7,output(hm)) using Momentum(p1&p2).
  HepMatrix hm(KF_NUM7,KF_NUM7,0);
  if(p1.t() == 0. || p2.t() == 0.){
    m_errorFlag = KF_OUTPUT_DIV_ZERO;
    if(m_errorMsgFlag)dout(Debugout::ERR,"kfitterbase") << "(kfitterbase::m_makeError,2): Divide by ZERO!" << std::endl;
    return hm;
  }
  for(unsigned i=0;i<3;++i){
    for(unsigned j=0;j<3;++j){
      hm[i][j]     = e[i][j];
      hm[4+i][4+j] = e[3+i][3+j];
      hm[4+i][j]   = e[3+i][j];
      hm[i][4+j]   = e[i][3+j];
    }
  }
  double invE1 = 1./p1.t();
  double invE2 = 1./p2.t();
  hm[0][3] = (p2.x()*hm[0][0]+p2.y()*hm[0][1]+p2.z()*hm[0][2])*invE2;
  hm[1][3] = (p2.x()*hm[1][0]+p2.y()*hm[1][1]+p2.z()*hm[1][2])*invE2;
  hm[2][3] = (p2.x()*hm[2][0]+p2.y()*hm[2][1]+p2.z()*hm[2][2])*invE2;
  hm[4][3] = (p2.x()*hm[4][0]+p2.y()*hm[4][1]+p2.z()*hm[4][2])*invE2;
  hm[5][3] = (p2.x()*hm[5][0]+p2.y()*hm[5][1]+p2.z()*hm[5][2])*invE2;
  hm[6][3] = (p2.x()*hm[6][0]+p2.y()*hm[6][1]+p2.z()*hm[6][2])*invE2;
  hm[3][3] = (p1.x()*p2.x()*hm[0][0]+p1.y()*p2.y()*hm[1][1]+p1.z()*p2.z()*hm[2][2]+
              p1.x()*p2.y()*hm[0][1]+p2.x()*p1.y()*hm[1][0]+
              p1.x()*p2.z()*hm[0][2]+p2.x()*p1.z()*hm[2][0]+
              p1.y()*p2.z()*hm[1][2]+p2.y()*p1.z()*hm[2][1])*invE1*invE2;
  hm[3][0] = (p1.x()*hm[0][0]+p1.y()*hm[1][0]+p1.z()*hm[2][0])*invE1;
  hm[3][1] = (p1.x()*hm[0][1]+p1.y()*hm[1][1]+p1.z()*hm[2][1])*invE1;
  hm[3][2] = (p1.x()*hm[0][2]+p1.y()*hm[1][2]+p1.z()*hm[2][2])*invE1;
  hm[3][4] = (p1.x()*hm[0][4]+p1.y()*hm[1][4]+p1.z()*hm[2][4])*invE1;
  hm[3][5] = (p1.x()*hm[0][5]+p1.y()*hm[1][5]+p1.z()*hm[2][5])*invE1;
  hm[3][6] = (p1.x()*hm[0][6]+p1.y()*hm[1][6]+p1.z()*hm[2][6])*invE1;

  return hm;
}

HepMatrix
kfitterbase::m_makeError2(const HepLorentzVector &p,
			  const HepMatrix &e){
  // vertex and track
  // Error(3x6,e) ==> Error(3x7,output(hm)) using Momentum(p).
  HepMatrix hm(3,KF_NUM7,0);
  if(p.t() == 0.){
    m_errorFlag = KF_OUTPUT_DIV_ZERO;
    if(m_errorMsgFlag)dout(Debugout::ERR,"kfitterbase") << "(kfitterbase::m_makeError2): Divide by ZERO!" << std::endl;
    return hm;
  }
  for(unsigned i=0;i<3;++i){
    for(unsigned j=0;j<3;++j){
      hm[i][j]   = e[i][j];
      hm[i][4+j] = e[i][3+j];
    }
  }
  double invE = 1./p.t();
  hm[0][3] = (p.x()*hm[0][0]+p.y()*hm[0][1]+p.z()*hm[0][2])*invE;
  hm[1][3] = (p.x()*hm[1][0]+p.y()*hm[1][1]+p.z()*hm[1][2])*invE;
  hm[2][3] = (p.x()*hm[2][0]+p.y()*hm[2][1]+p.z()*hm[2][2])*invE;
  return hm;
}

HepSymMatrix 
kfitterbase::m_makeError3(const HepLorentzVector &p,
			  const HepMatrix &e,
			  const int isFixMass){
  // self track
  // Error(7x7,e) ==> Error(7x7,output(hsm)) using Momentum(p).
  // isFixMass = 1 : Energy term is recalculated from the other parameters.
  // isFixMass = 0 : hsm = e.
  HepSymMatrix hsm(KF_NUM7,0);
  if(isFixMass == KF_UNFIX_MASS){
    for(unsigned i=0;i<7;++i){
      for(unsigned j=i;j<7;++j){
	hsm[i][j] = e[i][j];
      }
    }
    return hsm;
  }
  if(p.t() == 0.){
    m_errorFlag = KF_OUTPUT_DIV_ZERO;
    if(m_errorMsgFlag)dout(Debugout::ERR,"kfitterbase") << "(kfitterbase::m_makeError3,1): Divide by ZERO!" << std::endl;
    return hsm;
  }
  for(unsigned i=0;i<7;++i){
    if(i != 3){
      for(unsigned j=i;j<7;++j){
	hsm[i][j] = e[i][j];
      }
    }
  }
  double invE = 1./p.t();
  hsm[0][3] = (p.x()*hsm[0][0]+p.y()*hsm[0][1]+p.z()*hsm[0][2])*invE;
  hsm[1][3] = (p.x()*hsm[0][1]+p.y()*hsm[1][1]+p.z()*hsm[1][2])*invE;
  hsm[2][3] = (p.x()*hsm[0][2]+p.y()*hsm[1][2]+p.z()*hsm[2][2])*invE;
  hsm[3][3] = (p.x()*p.x()*hsm[0][0]+p.y()*p.y()*hsm[1][1]+p.z()*p.z()*hsm[2][2]
               +2.0*p.x()*p.y()*hsm[0][1]
               +2.0*p.x()*p.z()*hsm[0][2]
               +2.0*p.y()*p.z()*hsm[1][2])*invE*invE;
  hsm[3][4] = (p.x()*hsm[0][4]+p.y()*hsm[1][4]+p.z()*hsm[2][4])*invE;
  hsm[3][5] = (p.x()*hsm[0][5]+p.y()*hsm[1][5]+p.z()*hsm[2][5])*invE;
  hsm[3][6] = (p.x()*hsm[0][6]+p.y()*hsm[1][6]+p.z()*hsm[2][6])*invE;
  return hsm;
}

HepMatrix 
kfitterbase::m_makeError3(const HepLorentzVector &p1,
			  const HepLorentzVector &p2,
			  const HepMatrix &e,
			  const int isFixMass1,
			  const int isFixMass2){
  // track and track
  // Error(7x7,e) ==> Error(7x7,output(hm)) using Momentum(p1&p2).
  // isFixMass = 1 : Energy term is recalculated from the other parameters.
  // isFixMass = 0 : not.
  if(isFixMass1 == KF_UNFIX_MASS && isFixMass2 == KF_UNFIX_MASS){
    return e;
  }else if(isFixMass1 == KF_FIX_MASS && isFixMass2 == KF_FIX_MASS){
    HepMatrix hm(e);
    if(p1.t() == 0. || p2.t() == 0.){
      m_errorFlag = KF_OUTPUT_DIV_ZERO;
      if(m_errorMsgFlag)dout(Debugout::ERR,"kfitterbase") << "(kfitterbase::m_makeError3,2): Divide by ZERO!" << std::endl;
      return hm;
    }
    double invE1 = 1./p1.t();
    double invE2 = 1./p2.t();
    hm[0][3] = (p2.x()*hm[0][0]+p2.y()*hm[0][1]+p2.z()*hm[0][2])*invE2;
    hm[1][3] = (p2.x()*hm[1][0]+p2.y()*hm[1][1]+p2.z()*hm[1][2])*invE2;
    hm[2][3] = (p2.x()*hm[2][0]+p2.y()*hm[2][1]+p2.z()*hm[2][2])*invE2;
    hm[4][3] = (p2.x()*hm[4][0]+p2.y()*hm[4][1]+p2.z()*hm[4][2])*invE2;
    hm[5][3] = (p2.x()*hm[5][0]+p2.y()*hm[5][1]+p2.z()*hm[5][2])*invE2;
    hm[6][3] = (p2.x()*hm[6][0]+p2.y()*hm[6][1]+p2.z()*hm[6][2])*invE2;
    hm[3][0] = (p1.x()*hm[0][0]+p1.y()*hm[1][0]+p1.z()*hm[2][0])*invE1;
    hm[3][1] = (p1.x()*hm[0][1]+p1.y()*hm[1][1]+p1.z()*hm[2][1])*invE1;
    hm[3][2] = (p1.x()*hm[0][2]+p1.y()*hm[1][2]+p1.z()*hm[2][2])*invE1;
    hm[3][3] = (p1.x()*p2.x()*hm[0][0]+p1.y()*p2.y()*hm[1][1]+p1.z()*p2.z()*hm[2][2]+
		p1.x()*p2.y()*hm[0][1]+p2.x()*p1.y()*hm[1][0]+
		p1.x()*p2.z()*hm[0][2]+p2.x()*p1.z()*hm[2][0]+
		p1.y()*p2.z()*hm[1][2]+p2.y()*p1.z()*hm[2][1])*invE1*invE2;
    hm[3][4] = (p1.x()*hm[0][4]+p1.y()*hm[1][4]+p1.z()*hm[2][4])*invE1;
    hm[3][5] = (p1.x()*hm[0][5]+p1.y()*hm[1][5]+p1.z()*hm[2][5])*invE1;
    hm[3][6] = (p1.x()*hm[0][6]+p1.y()*hm[1][6]+p1.z()*hm[2][6])*invE1;

    return hm;
  }else if(isFixMass1 == KF_FIX_MASS && isFixMass2 == KF_UNFIX_MASS){
    HepMatrix hm(e);
    if(p1.t() == 0.){
      m_errorFlag = KF_OUTPUT_DIV_ZERO;
      if(m_errorMsgFlag)dout(Debugout::ERR,"kfitterbase") << "(kfitterbase::m_makeError3,3): Divide by ZERO!" << std::endl;
      return hm;
    }
    double invE1 = 1./p1.t();
    hm[3][0] = (p1.x()*hm[0][0]+p1.y()*hm[1][0]+p1.z()*hm[2][0])*invE1;
    hm[3][1] = (p1.x()*hm[0][1]+p1.y()*hm[1][1]+p1.z()*hm[2][1])*invE1;
    hm[3][2] = (p1.x()*hm[0][2]+p1.y()*hm[1][2]+p1.z()*hm[2][2])*invE1;
    hm[3][3] = (p1.x()*hm[0][3]+p1.y()*hm[1][3]+p1.z()*hm[2][3])*invE1;
    hm[3][4] = (p1.x()*hm[0][4]+p1.y()*hm[1][4]+p1.z()*hm[2][4])*invE1;
    hm[3][5] = (p1.x()*hm[0][5]+p1.y()*hm[1][5]+p1.z()*hm[2][5])*invE1;
    hm[3][6] = (p1.x()*hm[0][6]+p1.y()*hm[1][6]+p1.z()*hm[2][6])*invE1;

    return hm;    
  }else{
    HepMatrix hm(e);
    if(p2.t() == 0.){
      m_errorFlag = KF_OUTPUT_DIV_ZERO;
      if(m_errorMsgFlag)dout(Debugout::ERR,"kfitterbase") << "(kfitterbase::m_makeError3,4): Divide by ZERO!" << std::endl;
      return hm;
    }
    double invE2 = 1./p2.t();
    hm[0][3] = (p2.x()*hm[0][0]+p2.y()*hm[0][1]+p2.z()*hm[0][2])*invE2;
    hm[1][3] = (p2.x()*hm[1][0]+p2.y()*hm[1][1]+p2.z()*hm[1][2])*invE2;
    hm[2][3] = (p2.x()*hm[2][0]+p2.y()*hm[2][1]+p2.z()*hm[2][2])*invE2;
    hm[3][3] = (p2.x()*hm[3][0]+p2.y()*hm[3][1]+p2.z()*hm[3][2])*invE2;
    hm[4][3] = (p2.x()*hm[4][0]+p2.y()*hm[4][1]+p2.z()*hm[4][2])*invE2;
    hm[5][3] = (p2.x()*hm[5][0]+p2.y()*hm[5][1]+p2.z()*hm[5][2])*invE2;
    hm[6][3] = (p2.x()*hm[6][0]+p2.y()*hm[6][1]+p2.z()*hm[6][2])*invE2;

    return hm;    
  }
}

HepMatrix
kfitterbase::m_makeError4(const HepLorentzVector &p,
			  const HepMatrix &e){
  // vertex and track
  // Error(3x7,e) ==> Error(3x7,output(hm)) using Momentum(p).
  // Energy term is recalculated from the other parameters.
  HepMatrix hm(e);
  if(p.t() == 0.){
    m_errorFlag = KF_OUTPUT_DIV_ZERO;
    if(m_errorMsgFlag)dout(Debugout::ERR,"kfitterbase") << "(kfitterbase::m_makeError4): Divide by ZERO!" << std::endl;
    return hm;
  }
  double invE = 1./p.t();
  hm[0][3] = (p.x()*hm[0][0]+p.y()*hm[0][1]+p.z()*hm[0][2])*invE;
  hm[1][3] = (p.x()*hm[1][0]+p.y()*hm[1][1]+p.z()*hm[1][2])*invE;
  hm[2][3] = (p.x()*hm[2][0]+p.y()*hm[2][1]+p.z()*hm[2][2])*invE;
  return hm;
}

unsigned
kfitterbase::m_setCorrelation(void){
  if(m_correlation_b.size() != static_cast<unsigned>((m_trackNum*(m_trackNum-1))/2)){
    m_errorFlag = KF_INPUT_CORRE_SIZE;
    return m_errorFlag;
  }

  HepMatrix tmp_hm(KF_NUM6,KF_NUM6,0);
  unsigned row(0), col(0);

  for(std::vector<HepMatrix>::iterator it = m_correlation_b.begin(),
      endIt = m_correlation_b.end();
      it != endIt; ++it){
    HepMatrix &hm = *it;
    // 2000/03/07
    //counter
    ++row;
    if(row == m_trackNum){
      ++col;
      row = col + 1;
    }

    //7x7 --> 6x6
    for(unsigned i=0;i<3;++i){
      for(unsigned j=0;j<3;++j){
        tmp_hm[i][j]     = hm[i][j];
        tmp_hm[3+i][3+j] = hm[4+i][4+j];
        tmp_hm[3+i][j]   = hm[4+i][j];
        tmp_hm[i][3+j]   = hm[i][4+j];
      }
    }

    unsigned ii(0), jj(0);
    for(unsigned i=KF_NUM6*row; i<KF_NUM6*(row+1); ++i){
      for(unsigned j=KF_NUM6*col; j<KF_NUM6*(col+1); ++j){
	m_V_al_0[i][j] = tmp_hm[ii][jj];
	++jj;
      }
      jj = 0;
      ++ii;
    }
  }
  return m_errorFlag;
}

unsigned 
kfitterbase::fit(void){
  if(m_errorFlag != KF_NO_ERROR)return m_errorFlag;
  if(m_trackNum < m_necessaryTrackNum){
    m_errorFlag = KF_TRACK_SIZE;
    return m_errorFlag;
  }
  if(m_setInputMatrix() != KF_NO_ERROR)return m_errorFlag;
  if(m_calDgf() != KF_NO_ERROR)return m_errorFlag;

  double chiSq(0.);
  int errInverse(0);
  double tmp_chiSq(KF_INIT_CHI2);
  HepMatrix tmp_al_1(m_al_1);
  HepMatrix tmp_V_al_1(m_V_al_1);

  m_al_a = m_al_0;
  HepMatrix tmp_al_a(m_al_a);

  for(unsigned i=0;i<KF_MAX_ITERATION_NUMBER;++i){
    if(m_makeCoreMatrix() != KF_NO_ERROR)return m_errorFlag;
    //m_V_D = (m_D*m_V_al_0*(m_D.T())).inverse(errInverse);
    m_V_D = (m_V_al_0.similarity(m_D)).inverse(errInverse); // 2000/03/06
    if(errInverse != 0){
      m_errorFlag = KF_INVERSE;
      return m_errorFlag;
    }
    m_lam    = m_V_D*(m_D*(m_al_0-m_al_1)+m_d);
    chiSq    = ((m_lam.T())*(m_D*(m_al_0-m_al_1)+m_d))(1,1);
    //m_al_1   = m_al_0 - m_V_al_0*(m_D.T())*(m_lam.T());
    m_al_1   = m_al_0 - m_V_al_0*(m_D.T())*m_lam;
    m_V_al_1 = m_V_al_0 - m_V_al_0*(m_D.T())*m_V_D*m_D*m_V_al_0;
    if(tmp_chiSq > chiSq){
      tmp_chiSq  = chiSq;
      tmp_al_a   = tmp_al_1;
      tmp_al_1   = m_al_1;
      tmp_V_al_1 = m_V_al_1;
      if(i == KF_MAX_ITERATION_NUMBER-1){
	m_al_a = tmp_al_1;
	m_overIterationFlag = KF_OVER_ITERATION;
      }else continue;
    }else if(i != 0){
      chiSq    = tmp_chiSq;
      m_al_1   = tmp_al_1;
      m_al_a   = tmp_al_a;
      m_V_al_1 = tmp_V_al_1;
      break;
    }else if(i == 0){
      m_errorFlag = KF_INIT_CHISQ;
      return m_errorFlag;
    }
  }
  if(m_errorFlag != KF_NO_ERROR)return m_errorFlag;
  if(m_setOutputMatrix() != KF_NO_ERROR)return m_errorFlag;

  m_chisq = chiSq;
  m_cl    = chisq2Confi(m_dgf, m_chisq);
  return m_errorFlag;
}

//double 
//kfitterbase::estimation(void){
//rough estimation...
//This # is almost 0, I think..... 
//return 2.*((m_lam.T()*(m_D*(m_al_1-m_al_a)+m_d))[0][0]);
//}

void 
kfitterbase::dump(const unsigned flag){
  if(flag == KF_DUMP_MEASUREMENT){
    dout(Debugout::DUMP,"kfitterbase") << "Track#=" << m_trackNum << ", ErrorFlag=" << m_errorFlag << std::endl;
    dout(Debugout::DUMP,"kfitterbase") << "m_al_0=" << m_al_0 << std::endl;
    dout(Debugout::DUMP,"kfitterbase") << "m_V_al_0=" << m_V_al_0 << std::endl;
  }else if(flag == KF_DUMP_CORE_MATRIX){
    dout(Debugout::DUMP,"kfitterbase") << "Track#=" << m_trackNum << ", ErrorFlag=" << m_errorFlag << std::endl;
    dout(Debugout::DUMP,"kfitterbase") << "m_D=" << m_D << std::endl;
    dout(Debugout::DUMP,"kfitterbase") << "m_d=" << m_d << std::endl;
  }else if(flag == KF_DUMP_FITTED){
    dout(Debugout::DUMP,"kfitterbase") << "Track#=" << m_trackNum << ", ErrorFlag=" << m_errorFlag << std::endl;
    dout(Debugout::DUMP,"kfitterbase") << "Chi2=" << m_chisq << ", Dgf=" << m_dgf << std::endl;
    dout(Debugout::DUMP,"kfitterbase") << "m_al_1=" << m_al_1 << std::endl;
    dout(Debugout::DUMP,"kfitterbase") << "m_V_al_1=" << m_V_al_1 << std::endl;
  }
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
