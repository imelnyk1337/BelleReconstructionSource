//
// $Id: kvertexfitter.cc 11161 2010-06-06 18:20:44Z hitoshi $
//
// $Log$
// Revision 1.25  2003/04/08 03:35:24  sumisawa
// modify to fix OS dependence in the calculation of xi^2
//
// Revision 1.24  2002/06/20 16:17:53  sumisawa
// Nakadaira modified to avoid "Segmentaion violation"
//
// Revision 1.23  2002/06/04 12:56:06  sumisawa
// add functions which return modified chi2
//
// Revision 1.22  2002/03/27 23:41:57  jtanaka
// Add new fitter to the mass-constraint fit.
//
// Revision 1.21  2002/03/12 14:55:09  jtanaka
// bug fix, which doesn't affect the result of the physics.
//
// Revision 1.20  2002/02/23 17:56:53  katayama
// Added (int) for int/unsigned comparison
//
// Revision 1.19  2000/06/07 12:32:05  jtanaka
// add m_errorMsgFlag to control "dump error message".
// Default is ON. errorMsg(0) --> OFF.
//
// Revision 1.18  2000/05/30 08:42:29  jtanaka
// avoid overestimation of vertex-z error at fit3.
//
// Revision 1.17  2000/04/24 09:29:36  jtanaka
// bug fix : calculation of each tracks in fit3() of kvertexfitter.
// add vertexChisq() function.
//
// Revision 1.16  2000/04/13 12:35:57  katayama
// Added std:: to cout,cerr,endl etc.
//
// Revision 1.15  2000/03/26 06:27:57  jtanaka
// bug fix : eachChisq
//
// Revision 1.14  2000/03/07 17:09:01  jtanaka
// bug fixed, use "similarity" etc.
//
// Revision 1.13  2000/03/07 10:52:40  katayama
// compatibility with CC5.0
//
// Revision 1.12  1999/12/10 15:50:26  jtanaka
// add refitter for slowpion to kvertexfitter
//
// Revision 1.11  1999/07/19 18:06:49  jtanaka
// bug fix : chisq
//
// Revision 1.10  1999/07/19 15:58:23  jtanaka
// added Mdst_charged
//
// Revision 1.9  1999/06/20 06:30:23  jtanaka
// added beam position fit
//
// Revision 1.8  1999/05/26 12:02:17  jtanaka
// speed up etc.
//
// Revision 1.7  1999/04/08 12:48:07  jtanaka
// remove "exit(1);"
//
// Revision 1.6  1999/03/29 05:39:19  jtanaka
// new class structure and new+old interfaces
//
// Revision 1.5  1998/07/29 00:00:11  katayama
// Less warnings with g++ -Wall
//
// Revision 1.4  1998/01/22 03:20:16  jtanaka
// Updated from Tanaka san. New Interface etc.
//
// Revision 1.3  1997/11/12 08:03:50  katayama
// Minor mods from Tanaka san
//
// Revision 1.2  1997/10/24 07:05:17  katayama
// Updated from Tanaka san. Bug fixes etc.
//
// Revision 1.1  1997/10/04 05:30:08  katayama
// New from Tanaka san
//
//

#include <stdlib.h>
#include "belle.h"
#include "kfitter/kvertexfitter.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

kvertexfitter::kvertexfitter(void)
: m_mode(KF_MODE_WO_CORRELATION),
  m_vertexChisq(0.),
  m_vertex_b(0.,0.,0.),
  m_errVertex_a(3,0),
  m_errBeam(3,0),
  m_beam(0),
  m_knownVertex(0),
  m_tube(0),      // added by T.H 2006/05/11
  m_itrk_tube(-1) // added by T.H 2006/05/11
	{
  m_eachChisq = new double[KF_MAX_TRACK_NUMBER2];
  m_necessaryTrackNum = 2;
  m_V_E = HepMatrix(3,3,0);
  m_v   = HepMatrix(3,1,0);
  m_v_a = HepMatrix(3,1,0);
}

kvertexfitter::~kvertexfitter(void){
  delete [] m_eachChisq;
  m_errVertexTrack_a.erase(m_errVertexTrack_a.begin(),m_errVertexTrack_a.end());
}

void
kvertexfitter::initialVertex(const HepPoint3D &v){
  m_vertex_b = v;
}

unsigned
kvertexfitter::mode(const unsigned m){
  m_mode = m;
  return m_mode;
}

void
kvertexfitter::beamProfile(const HepSymMatrix &e){
  vertexProfile(e);
}

void
kvertexfitter::vertexProfile(const HepSymMatrix &e){
	if( m_tube ){
	  dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::vertexProfile] already constrained to IPtube" << std::endl;
		exit(1);
	}

  m_beam = 1;
  m_errBeam = e;
}

void
kvertexfitter::tubeProfile(const kfitterparticle &p){
	if( m_beam ){
	  dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::tubeProfile] already constrained to IP"<< std::endl;
		exit(1);
	}

  m_tube = 1;
  m_tube_track = p;
}

// added by T.H 2006/05/11
void
kvertexfitter::addTube(void)
{
	if( !m_tube ) return;

	if( m_itrk_tube!=-1 ){
	  dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::addTube] internal error, duplicated addTube() call?" << std::endl;
		exit(1);
	}

	m_plist.push_back(m_tube_track);
	m_trackNum  = m_plist.size();
	m_itrk_tube = m_trackNum;
}

// added by T.H 2006/05/11
void
kvertexfitter::delTube(void)
{
	if( !m_tube ) return;

	if( m_itrk_tube==-1 ){
	  dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::addTube] internal error, duplicated delTube() call?" << std::endl;
		exit(1);
	}

	m_plist.pop_back();
	m_trackNum  = m_plist.size();
	m_itrk_tube = -1;
}

unsigned
kvertexfitter::mode(void){
  return m_mode;
}

unsigned
kvertexfitter::m_setInputMatrix(void){
  if(m_mode == KF_MODE_WO_CORRELATION || m_beam == 1 || m_knownVertex == 1){
    if(m_trackNum > KF_MAX_TRACK_NUMBER2){
      m_errorFlag = KF_INPUT_TRACK_SIZE;
      return m_errorFlag;
    }
  }else{
    if(m_trackNum > KF_MAX_TRACK_NUMBER){
      m_errorFlag = KF_INPUT_TRACK_SIZE;
      return m_errorFlag;
    }
  }

  unsigned index(0);
  HepMatrix    tmp_al_0(KF_NUM6*m_trackNum,1,0);
  HepSymMatrix tmp_V_al_0(KF_NUM6*m_trackNum,0);
  HepMatrix    tmp_property(m_trackNum,3,0);

  for(std::vector<kfitterparticle>::iterator it = m_plist.begin(), endIt = m_plist.end();
      it != endIt; ++it){
    //momentum x,y,z and position x,y,z
    for(unsigned j=0;j<KF_NUM6;++j)
      tmp_al_0[index*KF_NUM6+j][0] = it->getFitParameter(j, KF_BEFORE_FIT);
    //these error
    tmp_V_al_0.sub(index*KF_NUM6+1, it->getFitError(KF_BEFORE_FIT));
    //charge , mass , a
    tmp_property[index][0] =  it->charge();
    tmp_property[index][1] =  it->mass();
    tmp_property[index][2] = -KF_PHOTON_VELOCITY*m_magneticField*it->charge();
    ++index;
  }

  //error between tarck and track
  m_V_al_0 = tmp_V_al_0;
  if(m_knownVertex == 0 && m_beam == 0 && m_mode != KF_MODE_WO_CORRELATION){
    if(m_correlationFlag == 1){
      m_errorFlag = m_setCorrelation();
      if(m_errorFlag != KF_NO_ERROR)return m_errorFlag;
    }
  }

  //vertex
  m_v_a[0][0] = m_vertex_b.x();
  m_v_a[1][0] = m_vertex_b.y();
  m_v_a[2][0] = m_vertex_b.z();

  //set member matrix
  m_al_0     = tmp_al_0;
  m_al_1     = m_al_0;
  m_property = tmp_property;

  //define size of matrix
  m_V_al_1 = HepMatrix(KF_NUM6*m_trackNum,KF_NUM6*m_trackNum,0);
  m_D    = m_V_al_1.sub(1,m_trackNum*2,1,KF_NUM6*m_trackNum);
  m_E    = m_V_al_1.sub(1,m_trackNum*2,1,3);
  m_d    = m_V_al_1.sub(1,m_trackNum*2,1,1);
  m_V_D  = m_V_al_1.sub(1,m_trackNum*2,1,m_trackNum*2);
  m_lam  = m_V_al_1.sub(1,m_trackNum*2,1,1);
  m_lam0 = m_V_al_1.sub(1,m_trackNum*2,1,1);
  m_V_Dt = m_V_al_1.sub(1,m_trackNum*2,1,m_trackNum*2);
  m_Cov_v_al_1 = m_V_al_1.sub(1,3,1,KF_NUM6*m_trackNum);

  return m_errorFlag;
}

unsigned
kvertexfitter::m_setInputSubMatrix(void){
  //vertex
  for(unsigned i=0;i<3;++i){
    m_v[i][0] = m_v_a[i][0];
  }
  return m_errorFlag;
}

unsigned
kvertexfitter::m_makeCoreMatrix(void){
#if 0
  double px,py,pz,x,y,z,a;
  double pt,pt2,dlx,dly,dlz,a1,a2,r2d2,B;
  double sininv,sqrtag;
  double tmp0,tmp1;

  for(unsigned i=0;i<m_trackNum;++i){
    px = m_al_1[i*KF_NUM6+0][0];
    py = m_al_1[i*KF_NUM6+1][0];
    pz = m_al_1[i*KF_NUM6+2][0];
    x  = m_al_1[i*KF_NUM6+3][0];
    y  = m_al_1[i*KF_NUM6+4][0];
    z  = m_al_1[i*KF_NUM6+5][0];
    a  = m_property[i][2];

    pt   = sqrt(px*px+py*py);
    pt2  =  pt*pt;
    dlx  =  m_v_a[0][0]-x;
    dly  =  m_v_a[1][0]-y;
    dlz  =  m_v_a[2][0]-z;
    a1   = -dlx*py+dly*px;
    a2   =  dlx*px+dly*py;
    r2d2 =  dlx*dlx+dly*dly;
    if(a != 0.){
      if(pt2 == 0.){
	m_errorFlag = KF_DIV_ZERO;
	return m_errorFlag;
      }
      B = a*a2/pt2;
      if(fabs(B) > 1.){
	m_errorFlag = KF_ARCSIN;
	return m_errorFlag;
      }
      //sin^(-1)(B)
      sininv = asin(B);
      tmp0 = 1.0-B*B;
      if(tmp0 == 0.){
	m_errorFlag = KF_DIV_ZERO;
	return m_errorFlag;
      }
      //1/sqrt(1-B^2)
      sqrtag = 1.0/sqrt(tmp0);
    }else{
      sininv = 0.0;
      sqrtag = 1.0;
    }
    if(a != 0.){
      if(pt2 == 0.){
        m_errorFlag = KF_DIV_ZERO;
        return m_errorFlag;
      }
      //charged track
      tmp1 = dlz-pz*sininv/a;
      //d
      m_d[i*2+0][0] = a1-0.5*a*r2d2;
      m_d[i*2+1][0] = tmp1*pt2;

      //D
      m_D[i*2+0][i*KF_NUM6+0] =  dly;
      m_D[i*2+0][i*KF_NUM6+1] = -dlx;
      m_D[i*2+0][i*KF_NUM6+2] =  0.0;
      m_D[i*2+0][i*KF_NUM6+3] =  py+a*dlx;
      m_D[i*2+0][i*KF_NUM6+4] = -px+a*dly;
      m_D[i*2+0][i*KF_NUM6+5] =  0.0;
      m_D[i*2+1][i*KF_NUM6+0] = -pz*sqrtag*(dlx-2.0*px*a2/pt2)+2.0*tmp1*px;
      m_D[i*2+1][i*KF_NUM6+1] = -pz*sqrtag*(dly-2.0*py*a2/pt2)+2.0*tmp1*py;
      m_D[i*2+1][i*KF_NUM6+2] = -sininv/a*pt2;
      m_D[i*2+1][i*KF_NUM6+3] =  px*pz*sqrtag;
      m_D[i*2+1][i*KF_NUM6+4] =  py*pz*sqrtag;
      m_D[i*2+1][i*KF_NUM6+5] = -pt2;

      //E
      m_E[i*2+0][0] = -py-a*dlx;
      m_E[i*2+0][1] =  px-a*dly;
      m_E[i*2+0][2] =  0.0;
      m_E[i*2+1][0] = -px*pz*sqrtag;
      m_E[i*2+1][1] = -py*pz*sqrtag;
      m_E[i*2+1][2] =  pt2;
    }else{
      //neutral track
      //d
      m_d[i*2+0][0] = a1;
      m_d[i*2+1][0] = dlz*pt2-pz*a2;

      //D
      m_D[i*2+0][i*KF_NUM6+0] =  dly;
      m_D[i*2+0][i*KF_NUM6+1] = -dlx;
      m_D[i*2+0][i*KF_NUM6+2] =  0.0;
      m_D[i*2+0][i*KF_NUM6+3] =  py;
      m_D[i*2+0][i*KF_NUM6+4] = -px;
      m_D[i*2+0][i*KF_NUM6+5] =  0.0;
      m_D[i*2+1][i*KF_NUM6+0] =  2.0*px*dlz-pz*dlx;
      m_D[i*2+1][i*KF_NUM6+1] =  2.0*py*dlz-pz*dly;
      m_D[i*2+1][i*KF_NUM6+2] = -a2;
      m_D[i*2+1][i*KF_NUM6+3] =  px*pz;
      m_D[i*2+1][i*KF_NUM6+4] =  py*pz;
      m_D[i*2+1][i*KF_NUM6+5] = -pt2;

      //E
      m_E[i*2+0][0] = -py;
      m_E[i*2+0][1] =  px;
      m_E[i*2+0][2] =  0.0;
      m_E[i*2+1][0] = -px*pz;
      m_E[i*2+1][1] = -py*pz;
      m_E[i*2+1][2] =  pt2;
    }
  }
#else
  // 2000/03/07
  //vertex fit
  double px,py,pz,x,y,z,a;
  double pt,invPt,invPt2,dlx,dly,dlz,a1,a2,r2d2,B,Rx,Ry,S,U;
  double sininv,sqrtag;

  for(unsigned i=0;i<m_trackNum;++i){
    px = m_al_1[i*KF_NUM6+0][0];
    py = m_al_1[i*KF_NUM6+1][0];
    pz = m_al_1[i*KF_NUM6+2][0];
    x  = m_al_1[i*KF_NUM6+3][0];
    y  = m_al_1[i*KF_NUM6+4][0];
    z  = m_al_1[i*KF_NUM6+5][0];
    a  = m_property[i][2];

    pt   =  sqrt(px*px+py*py);

    if(pt == 0.){
      m_errorFlag = KF_DIV_ZERO;
      return m_errorFlag;
    }

    invPt   = 1./pt;
    invPt2  = invPt*invPt;
    dlx  =  m_v_a[0][0]-x;
    dly  =  m_v_a[1][0]-y;
    dlz  =  m_v_a[2][0]-z;
    a1   = -dlx*py+dly*px;
    a2   =  dlx*px+dly*py;
    r2d2 =  dlx*dlx+dly*dly;
    Rx   =  dlx-2.*px*a2*invPt2;
    Ry   =  dly-2.*py*a2*invPt2;
    if(a != 0.){
      // charged
      B = a*a2*invPt2;
      if(fabs(B) > 1.){
	m_errorFlag = KF_ARCSIN;
	return m_errorFlag;
      }
      //sin^(-1)(B)
      sininv = asin(B);
      double tmp0 = 1.0-B*B;
      if(tmp0 == 0.){
	m_errorFlag = KF_DIV_ZERO;
	return m_errorFlag;
      }
      //1/sqrt(1-B^2)
      sqrtag = 1.0/sqrt(tmp0);
      S = sqrtag*invPt2;
      U = dlz-pz*sininv/a;
    }else{
      // neutral
      B      = 0.0;
      sininv = 0.0;
      sqrtag = 1.0;
      S = invPt2;
      U = dlz-pz*a2*invPt2;
    }

    //d
    m_d[i*2+0][0] = a1-0.5*a*r2d2;
    m_d[i*2+1][0] = U*pt;

    //D
    m_D[i*2+0][i*KF_NUM6+0] =  dly;
    m_D[i*2+0][i*KF_NUM6+1] = -dlx;
    m_D[i*2+0][i*KF_NUM6+2] =  0.0;
    m_D[i*2+0][i*KF_NUM6+3] =  py+a*dlx;
    m_D[i*2+0][i*KF_NUM6+4] = -px+a*dly;
    m_D[i*2+0][i*KF_NUM6+5] =  0.0;
    m_D[i*2+1][i*KF_NUM6+0] = -pz*pt*S*Rx+U*px*invPt;
    m_D[i*2+1][i*KF_NUM6+1] = -pz*pt*S*Ry+U*py*invPt;
    if(a != 0.){
      m_D[i*2+1][i*KF_NUM6+2] = -sininv*pt/a;
    }else{
      m_D[i*2+1][i*KF_NUM6+2] = -a2*invPt;
    }
    m_D[i*2+1][i*KF_NUM6+3] =  px*pz*pt*S;
    m_D[i*2+1][i*KF_NUM6+4] =  py*pz*pt*S;
    m_D[i*2+1][i*KF_NUM6+5] = -pt;

    //E
    m_E[i*2+0][0] = -py-a*dlx;
    m_E[i*2+0][1] =  px-a*dly;
    m_E[i*2+0][2] =  0.0;
    m_E[i*2+1][0] = -px*pz*pt*S;
    m_E[i*2+1][1] = -py*pz*pt*S;
    m_E[i*2+1][2] =  pt;
  }
#endif
  return m_errorFlag;
}

HepPoint3D
kvertexfitter::vertex(const unsigned f){
  if(f == KF_AFTER_FIT){
    return m_vertex_a;
  }else if(f == KF_BEFORE_FIT){
    return m_vertex_b;
  }
  m_errorFlag = KF_OUTPUT_OUT_RANGE;
  if(m_errorMsgFlag)dout(Debugout::ERR,"kvertexfitter") << "(kvertexfitter): Out of Range!" << std::endl;
  return HepPoint3D();
}

HepSymMatrix
kvertexfitter::errVertex(void){
  return m_errVertex_a;
}

HepMatrix
kvertexfitter::errVertexTrack(const unsigned n){
  if(n < m_trackNum){
    return m_errVertexTrack_a[n];
  }else{
    m_errorFlag = KF_OUTPUT_TRACK_NUM;
    if(m_errorMsgFlag)dout(Debugout::ERR,"kvertexfitter") << "(kvertexfitter): Out of Range!" << std::endl;
    return HepMatrix(3,KF_NUM7,0);
  }
}

unsigned
kvertexfitter::m_setOutputMatrix(void){
  Hep3Vector h3v;
  unsigned index(0);
  for(std::vector<kfitterparticle>::iterator it = m_plist.begin(),
      endIt = m_plist.end();
      it != endIt; ++it){
    kfitterparticle &pdata = *it;
    //tracks
    //momentum
    h3v.setX(m_al_1[index*KF_NUM6+0][0]);
    h3v.setY(m_al_1[index*KF_NUM6+1][0]);
    h3v.setZ(m_al_1[index*KF_NUM6+2][0]);
    pdata.momentum(HepLorentzVector(h3v, sqrt(h3v.mag2()+pdata.mass()*pdata.mass())),KF_AFTER_FIT);
    //position
    pdata.position(HepPoint3D(m_al_1[index*KF_NUM6+3][0],
                              m_al_1[index*KF_NUM6+4][0],
                              m_al_1[index*KF_NUM6+5][0]),KF_AFTER_FIT);
    //error of the tracks
    pdata.error(m_makeError(pdata.momentum(),
			    m_V_al_1.sub(index*KF_NUM6+1,
					 (index+1)*KF_NUM6,
					 index*KF_NUM6+1,
					 (index+1)*KF_NUM6)),KF_AFTER_FIT);
    if(m_errorFlag != KF_NO_ERROR)break;
    ++index;
  }

  //vertex
  m_vertex_a.setX(m_v_a[0][0]);
  m_vertex_a.setY(m_v_a[1][0]);
  m_vertex_a.setZ(m_v_a[2][0]);

  //error of the vertex
  for(unsigned i=0;i<3;++i){
    for(unsigned j=i;j<3;++j){
      m_errVertex_a[i][j] = m_V_E[i][j];
    }
  }

  //error between vertex and tracks
  for(unsigned i=0;i<m_trackNum;++i){
    HepMatrix hm(3,KF_NUM6,0);
    for(unsigned j=0;j<3;++j){
      for(unsigned k=0;k<KF_NUM6;++k){
	hm[j][k] = m_Cov_v_al_1[j][KF_NUM6*i+k];
      }
    }
    m_errVertexTrack_a.push_back(m_makeError2(m_plist[i].momentum(),hm));
  }
  return m_errorFlag;
}

unsigned
kvertexfitter::m_calDgf(void){
  m_dgf = 2*m_trackNum-3;
  if(m_beam == 1)m_dgf = 2*m_trackNum;

	// T.H
	// m_calDgf is called with in _fit().
	// In case of IP tube, m_trackNum==ntrk+1 within _fit()
	// due to the additional IP-tube-track.
	// real ndf = 2ntrk-1, but we have to compute ndf by 2*(m_trackNum-1)-1
  if(m_tube     )m_dgf = 2*(m_trackNum-1)-1;

  if(m_knownVertex == 1)m_dgf = 2*m_trackNum;
  return m_errorFlag;
}

unsigned
kvertexfitter::fit(void){
	unsigned int ret;
	if( m_tube ) this->addTube();
	ret = this->_fit();
	if( m_tube ) this->delTube();
	return ret;
}

unsigned
kvertexfitter::_fit(void){
  if(m_beam == 1)return _fit3();
  if(m_knownVertex == 1)return _fit4();
  if(m_mode != KF_MODE_WO_CORRELATION){
    return kfitterbase2::fit();
  }else{
    return _fit2();
  }
}

double
kvertexfitter::chisq(void){
  return kfitterbase::chisq();
}

double
kvertexfitter::chisq(const unsigned n){
  if(m_knownVertex == 0 && m_beam == 0 &&
     m_mode != KF_MODE_WO_CORRELATION){
    return kfitterbase::chisq(n);
  }else{
    if(n < m_trackNum){
      return m_eachChisq[n];
    }else{
      m_errorFlag = KF_OUTPUT_TRACK_NUM;
      if(m_errorMsgFlag)dout(Debugout::ERR,"kvertexfitter") << "(kvertexfitter): Out of Range!" << std::endl;
      return -1.;
    }
  }
}

double
kvertexfitter::vertexChisq(void){
  // only m_beam = 1
  return m_vertexChisq;
}

unsigned
kvertexfitter::fit2(void){
	unsigned int ret;
	if( m_tube ) this->addTube();
	ret = this->_fit2();
	if( m_tube ) this->delTube();
	return ret;
}

unsigned
kvertexfitter::_fit2(void){
  // use small Matrix --> No Correlation
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
  double tmp2_chiSq(KF_INIT_CHI2);

  m_al_a = m_al_0;
  HepMatrix tmp_al_a(m_al_a);

  HepMatrix tmp_D(m_D),tmp_E(m_E);
  HepMatrix tmp_V_D(m_V_D),tmp_V_E(m_V_E);
  HepMatrix tmp_lam0(m_lam0),tmp_v_a(m_v_a);

  HepMatrix tmp2_D(m_D),tmp2_E(m_E);
  HepMatrix tmp2_V_D(m_V_D),tmp2_V_E(m_V_E);
  HepMatrix tmp2_lam0(m_lam0),tmp2_v_a(m_v_a);

  double * tmp_eachChisq  = new double[m_trackNum];
  double * tmp2_eachChisq = new double[m_trackNum];

  for(unsigned j=0;j<KF_MAX_ITERATION_NUMBER;++j){
    tmp_chiSq = KF_INIT_CHI2;
    for(unsigned i=0;i<KF_MAX_ITERATION_NUMBER;++i){
      if(m_setInputSubMatrix() != KF_NO_ERROR){
	delete [] tmp_eachChisq;
	delete [] tmp2_eachChisq;
	return m_errorFlag;
      }
      if(m_makeCoreMatrix() != KF_NO_ERROR){
	delete [] tmp_eachChisq;
	delete [] tmp2_eachChisq;
	return m_errorFlag;
      }
      HepMatrix tV_Ein(3,3,0);
      chiSq = 0.;
      for(unsigned k=0;k<m_trackNum;++k){
	HepMatrix tD = m_D.sub(2*k+1,2*(k+1),KF_NUM6*k+1,KF_NUM6*(k+1));//2x6
#if 0
	HepMatrix tV_D = (tD*m_V_al_0.sub(KF_NUM6*k+1,(int)(KF_NUM6*(k+1)))*
			  (tD.T())).inverse(errInverse);//2x2
#else
	// 2000/03/06
	HepMatrix tV_D = ((m_V_al_0.sub(KF_NUM6*k+1,(int)(KF_NUM6*(k+1)))).similarity(tD)).inverse(errInverse);//2x2
#endif
	if(errInverse != 0){
	  delete [] tmp_eachChisq;
	  delete [] tmp2_eachChisq;
	  m_errorFlag = KF_INVERSE;
	  return m_errorFlag;
	}
	m_V_D.sub(2*k+1,2*k+1,tV_D);
	HepMatrix tE = m_E.sub(2*k+1,2*(k+1),1,3);//2x3
	tV_Ein += (tE.T())*tV_D*tE;//3x3
	HepMatrix tDeltaAl = (m_al_0-m_al_1).sub(KF_NUM6*k+1,KF_NUM6*(k+1),1,1);//6x1
	HepMatrix td = m_d.sub(2*k+1,2*(k+1),1,1);//2x1
	HepMatrix tlam0 = tV_D*(tD*tDeltaAl+td);//2x2x(2x6x6x1+2x1) = 2x1
	m_lam0.sub(2*k+1,1,tlam0);
	m_eachChisq[k] = ((tlam0.T())*(tD*tDeltaAl+tE*(m_v-m_v_a)+td))(1,1);//1x2x(2x6x6x1+2x3x3x1+2x1)
	chiSq  += m_eachChisq[k];
      }
      m_V_E = tV_Ein.inverse(errInverse);
      if(errInverse != 0){
	delete [] tmp_eachChisq;
	delete [] tmp2_eachChisq;
	m_errorFlag = KF_INVERSE;
	return m_errorFlag;
      }
      m_v_a = m_v_a - m_V_E*(m_E.T())*m_lam0;
      if(tmp_chiSq > chiSq){
	for(unsigned k=0;k<m_trackNum;++k)tmp_eachChisq[k] = m_eachChisq[k];
	tmp_chiSq = chiSq;
	tmp_v_a   = m_v_a;
	tmp_V_E   = m_V_E;
	tmp_V_D   = m_V_D;
	tmp_lam0  = m_lam0;
	tmp_E     = m_E;
	tmp_D     = m_D;
	if(i == KF_MAX_ITERATION_NUMBER-1)
	  m_overIterationFlag = KF_OVER_ITERATION;
	else continue;
      }else if(i != 0){
	for(unsigned k=0;k<m_trackNum;++k)m_eachChisq[k] = tmp_eachChisq[k];
	chiSq   = tmp_chiSq;
	m_v_a   = tmp_v_a;
	m_V_E   = tmp_V_E;
	m_V_D   = tmp_V_D;
	m_lam0  = tmp_lam0;
	m_E     = tmp_E;
	m_D     = tmp_D;
	break;
      }else if(i == 0){
	m_errorFlag = KF_INIT_CHISQ;
	break;
      }
    }
    m_al_a = m_al_1;
    m_lam  = m_lam0 - m_V_D*m_E*m_V_E*(m_E.T())*m_lam0;
    m_al_1 = m_al_0 - m_V_al_0*(m_D.T())*m_lam;
    if(j == 0){
      for(unsigned k=0;k<m_trackNum;++k)tmp2_eachChisq[k] = m_eachChisq[k];
      tmp2_chiSq = chiSq;
      tmp2_v_a   = m_v_a;
      tmp2_V_E   = m_V_E;
      tmp2_V_D   = m_V_D;
      tmp2_lam0  = m_lam0;
      tmp2_E     = m_E;
      tmp2_D     = m_D;
      tmp_al_a   = m_al_a;
      continue;
    }else{
      if(tmp2_chiSq > chiSq){
	for(unsigned k=0;k<m_trackNum;++k)tmp2_eachChisq[k] = m_eachChisq[k];
	tmp2_chiSq = chiSq;
	tmp2_v_a   = m_v_a;
	tmp2_V_E   = m_V_E;
	tmp2_V_D   = m_V_D;
	tmp2_lam0  = m_lam0;
	tmp2_E     = m_E;
	tmp2_D     = m_D;
	tmp_al_a   = m_al_a;
	if(j == KF_MAX_ITERATION_NUMBER-1)
	  m_overIterationFlag = KF_OVER_ITERATION;
	else continue;
      }else{
	for(unsigned k=0;k<m_trackNum;++k)m_eachChisq[k] = tmp2_eachChisq[k];
	chiSq   = tmp2_chiSq;
	m_v_a   = tmp2_v_a;
	m_V_E   = tmp2_V_E;
	m_V_D   = tmp2_V_D;
	m_lam0  = tmp2_lam0;
	m_E     = tmp2_E;
	m_D     = tmp2_D;
	m_al_a  = tmp_al_a;
	break;
      }
    }
  }
  if(m_errorFlag != KF_NO_ERROR){
    delete [] tmp_eachChisq;
    delete [] tmp2_eachChisq;
    return m_errorFlag;
  }
  m_lam    = m_lam0 - m_V_D*m_E*m_V_E*(m_E.T())*m_lam0;
  m_al_1   = m_al_0 - m_V_al_0*(m_D.T())*m_lam;
  m_V_Dt   = m_V_D  - m_V_D*m_E*m_V_E*(m_E.T())*m_V_D;
  m_V_al_1 = m_V_al_0 - m_V_al_0*(m_D.T())*m_V_Dt*m_D*m_V_al_0;
  m_Cov_v_al_1 = -m_V_E*(m_E.T())*m_V_D*m_D*m_V_al_0;
  if(m_setOutputMatrix() != KF_NO_ERROR){
    delete [] tmp_eachChisq;
    delete [] tmp2_eachChisq;
    return m_errorFlag;
  }

  m_chisq = chiSq;
  m_cl    = chisq2Confi(m_dgf, m_chisq);

  delete [] tmp_eachChisq;
  delete [] tmp2_eachChisq;
  return m_errorFlag;
}

//double
//kvertexfitter::estimation(void){
//rough estimation...
//This # is almost 0, I think.....
//if(m_mode != KF_MODE_WO_CORRELATION){
//  return kfitterbase2::estimation();
//}else{
//  return 0.;
//}
//}

unsigned
kvertexfitter::fit3(void){
	unsigned int ret;
	if( m_tube ) this->addTube();
	ret = this->_fit3();
	if( m_tube ) this->delTube();
	return ret;
}

unsigned
kvertexfitter::_fit3(void){
  // included beam position constraint (only no correlation)
  if(m_errorFlag != KF_NO_ERROR)return m_errorFlag;
  if(m_trackNum < m_necessaryTrackNum-1){
    m_errorFlag = KF_TRACK_SIZE;
    return m_errorFlag;
  }
  if(m_setInputMatrix() != KF_NO_ERROR)return m_errorFlag;
  if(m_calDgf() != KF_NO_ERROR)return m_errorFlag;

  double chiSq(0.);
  int errInverse(0);
  double tmp_chiSq(KF_INIT_CHI2);

  m_al_a = m_al_0;
  HepMatrix tmp_al_a(m_al_a);

  HepMatrix tmp_D(m_D),tmp_E(m_E);
  HepMatrix tmp_lam(m_lam);

  //vertex
  m_v[0][0] = m_vertex_b.x();
  m_v[1][0] = m_vertex_b.y();
  m_v[2][0] = m_vertex_b.z();

  double * tmp_eachChisq = new double[m_trackNum];
  double tmp_vertexChisq = 1.e+30; // An init-value is not needed but the C++ complier requires the init-value.

  // to avoid overestimation of vertex-z error.
  int itFlag = 0;

  for(unsigned i=0;i<KF_MAX_ITERATION_NUMBER;++i){
    if(m_makeCoreMatrix() != KF_NO_ERROR){
      delete [] tmp_eachChisq;
      return m_errorFlag;
    }
    chiSq = 0.;
    //HepMatrix tV_Dtin = m_D*m_V_al_0*(m_D.T())+m_E*m_errBeam*(m_E.T());
    HepMatrix tV_Dtin = m_V_al_0.similarity(m_D)+m_errBeam.similarity(m_E); // 2000/03/06
    HepMatrix tV_Dt = tV_Dtin.inverse(errInverse);
    if(errInverse != 0){
      delete [] tmp_eachChisq;
      m_errorFlag = KF_INVERSE;
      return m_errorFlag;
    }
    m_lam = tV_Dt*(m_D*(m_al_0-m_al_1)+m_E*(m_v-m_v_a)+m_d);//(2*nTrk)x1
    for(unsigned k=0;k<m_trackNum;++k){
      HepMatrix tD = m_D.sub(2*k+1,2*(k+1),
			     KF_NUM6*k+1,
			     KF_NUM6*(k+1));//2x6
      HepMatrix tDeltaAl = (m_al_0-m_al_1).sub(KF_NUM6*k+1,
					       KF_NUM6*(k+1),1,1);//6x1
      HepMatrix td = m_d.sub(2*k+1,2*(k+1),1,1);//2x1
      HepMatrix tE = m_E.sub(2*k+1,2*(k+1),1,3);//2x3
      /* m_eachChisq[k] = ((m_lam.sub(2*k+1,2*(k+1),1,1).T())*
	 (tD*tDeltaAl+tE*(m_v-m_v_a)+td))(1,1);//1x2x(2x6x6x1+2x3x3x1+2x1) */
      chiSq += ((m_lam.sub(2*k+1,2*(k+1),1,1).T())*
		(tD*tDeltaAl+tE*(m_v-m_v_a)+td))(1,1);//1x2x(2x6x6x1+2x3x3x1+2x1)
      m_eachChisq[k] = (m_lam.sub(2*k+1,2*(k+1),1,1).T()*tD*m_V_al_0.sub(KF_NUM6*k+1,
									 KF_NUM6*(k+1))*
			(tD.T())*m_lam.sub(2*k+1,2*(k+1),1,1))(1,1);
    }
    // dout(Debugout::DDEBUG,"kvertexfitter") << chiSq << ", <--> " << (m_lam.T())*(m_D*(m_al_0-m_al_1)+m_E*(m_v-m_v_a)+m_d) << std::endl;
    m_vertexChisq = (m_lam.T()*m_E*m_errBeam*(m_E.T())*m_lam)(1,1);
    m_al_a = m_al_1;
    m_v_a  = m_v - m_errBeam*(m_E.T())*m_lam;
    m_al_1 = m_al_0 - m_V_al_0*(m_D.T())*m_lam;
    // dout(Debugout::DDEBUG,"kvertexfitter") << tmp_chiSq << " --> " << chiSq << ", it = " << itFlag << std::endl;
    if(tmp_chiSq > chiSq || (tmp_chiSq <= chiSq && itFlag == 0)){
      if(tmp_chiSq <= chiSq)itFlag = 1;
      for(unsigned k=0;k<m_trackNum;++k)tmp_eachChisq[k] = m_eachChisq[k];
      tmp_vertexChisq = m_vertexChisq;
      tmp_chiSq = chiSq;
      tmp_lam   = m_lam;
      tmp_E     = m_E;
      tmp_D     = m_D;
      tmp_al_a  = m_al_a;
      if(i == KF_MAX_ITERATION_NUMBER-1)
	m_overIterationFlag = KF_OVER_ITERATION;
      else continue;
    }else if(i != 0){
      for(unsigned k=0;k<m_trackNum;++k)m_eachChisq[k] = tmp_eachChisq[k];
      m_vertexChisq = tmp_vertexChisq;
      chiSq   = tmp_chiSq;
      m_lam   = tmp_lam;
      m_E     = tmp_E;
      m_D     = tmp_D;
      m_al_a  = tmp_al_a;
      break;
    }else if(i == 0){
      m_errorFlag = KF_INIT_CHISQ;
      break;
    }
  }
  if(m_errorFlag != KF_NO_ERROR){
    delete [] tmp_eachChisq;
    return m_errorFlag;
  }
  m_al_1 = m_al_0 - m_V_al_0*(m_D.T())*m_lam;
  m_v_a  = m_v - m_errBeam*(m_E.T())*m_lam;
  //HepMatrix tV_Dtin = m_D*m_V_al_0*(m_D.T())+m_E*m_errBeam*(m_E.T());
  HepMatrix tV_Dtin = m_V_al_0.similarity(m_D)+m_errBeam.similarity(m_E); // 2000/03/06
  m_V_Dt = tV_Dtin.inverse(errInverse);
  if(errInverse != 0){
    delete [] tmp_eachChisq;
    m_errorFlag = KF_INVERSE;
    return m_errorFlag;
  }
  m_V_al_1 = m_V_al_0 - m_V_al_0*(m_D.T())*m_V_Dt*m_D*m_V_al_0;
  m_Cov_v_al_1 = -m_errBeam*(m_E.T())*m_V_Dt*m_D*m_V_al_0;
  // m_V_v is m_V_E
  // --> need to replace m_V_E for my implimentaion.
  m_V_E = m_errBeam-m_errBeam*(m_E.T())*m_V_Dt*m_E*m_errBeam;
  if(m_setOutputMatrix() != KF_NO_ERROR){
    delete [] tmp_eachChisq;
    return m_errorFlag;
  }

  m_chisq = chiSq;
  m_cl    = chisq2Confi(m_dgf, m_chisq);

  delete [] tmp_eachChisq;
  return m_errorFlag;
}

void
kvertexfitter::knownVertex(void)
{
  m_knownVertex = 1;
}

unsigned
kvertexfitter::fit4(void){
	unsigned int ret;
	if( m_tube ) this->addTube();
	ret = this->_fit4();
	if( m_tube ) this->delTube();
	return ret;
}

unsigned
kvertexfitter::_fit4(void){
  // known vertex --> do not find vertex. (only no correlation)
  if(m_errorFlag != KF_NO_ERROR)return m_errorFlag;
  if(m_trackNum < m_necessaryTrackNum-1){
    m_errorFlag = KF_TRACK_SIZE;
    return m_errorFlag;
  }
  if(m_setInputMatrix() != KF_NO_ERROR)return m_errorFlag;
  if(m_calDgf() != KF_NO_ERROR)return m_errorFlag;

  double chiSq(0.);
  int errInverse(0);
  double tmp_chiSq(KF_INIT_CHI2);

  m_al_a = m_al_0;
  HepMatrix tmp_al_a(m_al_a);

  HepMatrix tmp_al_0(m_al_1);
  HepMatrix tmp_V_al_0(m_V_al_1);

  double * tmp_eachChisq  = new double[m_trackNum];

  for(unsigned i=0;i<KF_MAX_ITERATION_NUMBER;++i){
    if(m_makeCoreMatrix() != KF_NO_ERROR){
      delete [] tmp_eachChisq;
      return m_errorFlag;
    }
    chiSq = 0.;
    for(unsigned k=0;k<m_trackNum;++k){
      HepMatrix tD = m_D.sub(2*k+1,2*(k+1),KF_NUM6*k+1,KF_NUM6*(k+1));//2x6
#if 0
      HepMatrix tV_D = (tD*m_V_al_0.sub(KF_NUM6*k+1,(int)(KF_NUM6*(k+1)))*
			(tD.T())).inverse(errInverse);//2x2
#else
      // 2000/03/06
      HepMatrix tV_D = ((m_V_al_0.sub(KF_NUM6*k+1,(int)(KF_NUM6*(k+1)))).similarity(tD)).inverse(errInverse);//2x2
#endif
      if(errInverse != 0){
	delete [] tmp_eachChisq;
	m_errorFlag = KF_INVERSE;
	return m_errorFlag;
      }
      m_V_D.sub(2*k+1,2*k+1,tV_D);

      HepMatrix tDeltaAl = (m_al_0-m_al_1).sub(KF_NUM6*k+1,KF_NUM6*(k+1),1,1);//6x1
      HepMatrix td = m_d.sub(2*k+1,2*(k+1),1,1);//2x1
      HepMatrix tlam = tV_D*(tD*tDeltaAl+td);//2x2x(2x6x6x1+2x1) = 2x1
      m_lam.sub(2*k+1,1,tlam);
      m_eachChisq[k] = ((tlam.T())*(tD*tDeltaAl+td))(1,1);//1x2x(2x6x6x1+2x1)
      chiSq  += m_eachChisq[k];
    }
    m_al_a   = m_al_1;
    m_al_1   = m_al_0 - m_V_al_0*(m_D.T())*m_lam;
    m_V_al_1 = m_V_al_0 - m_V_al_0*(m_D.T())*m_V_D*m_D*m_V_al_0;
    if(tmp_chiSq > chiSq){
      for(unsigned k=0;k<m_trackNum;++k)tmp_eachChisq[k] = m_eachChisq[k];
      tmp_chiSq  = chiSq;
      tmp_al_0   = m_al_1;
      tmp_V_al_0 = m_V_al_1;
      tmp_al_a   = m_al_a;
      if(i == KF_MAX_ITERATION_NUMBER-1)
	m_overIterationFlag = KF_OVER_ITERATION;
      else continue;
    }else if(i != 0){
      for(unsigned k=0;k<m_trackNum;++k)m_eachChisq[k] = tmp_eachChisq[k];
      chiSq    = tmp_chiSq;
      m_al_1   = tmp_al_0;
      m_V_al_1 = tmp_V_al_0;
      m_al_a   = tmp_al_a;
      break;
    }else if(i == 0){
      m_errorFlag = KF_INIT_CHISQ;
      break;
    }
  }
  if(m_errorFlag != KF_NO_ERROR){
    delete [] tmp_eachChisq;
    return m_errorFlag;
  }
  if(m_setOutputMatrix() != KF_NO_ERROR){
    delete [] tmp_eachChisq;
    return m_errorFlag;
  }

  m_chisq = chiSq;
  m_cl    = chisq2Confi(m_dgf, m_chisq);

  delete [] tmp_eachChisq;
  return m_errorFlag;
}

// Add by T.N (2002/05/14)
double kvertexfitter::trackchisq(const unsigned int n,
                               const bool usex, const bool usey, const bool usez){
  if(n>=m_trackNum){
    dout(Debugout::ERR,"kvertexfitter") << "Track # is out of range.( "
              <<  n << "/ " << m_trackNum << ")" << std::endl;
    m_errorFlag = KF_OUTPUT_TRACK_NUM;
    return -1.0;
  }
  static const double hfact(1.e1);

  HepMatrix pbf = m_plist[n].getFitParameter(KF_BEFORE_FIT);
  HepMatrix paf = m_plist[n].getFitParameter(KF_AFTER_FIT);
  HepSymMatrix pbe = m_plist[n].getFitError(KF_BEFORE_FIT);
  HepSymMatrix inv_pbe;

  int inv_error = 0;
  inv_pbe = pbe.inverse(inv_error);
  while(pbe.fast(6,6) < 1.e10 && (inv_error ||inv_pbe.fast(3,3) < 1.)){
    pbf *= hfact;
    paf *= hfact;
    pbe *= hfact * hfact;
    inv_pbe = pbe.inverse(inv_error);
  }

  pbf[0][0] = 0.0;
  paf[0][0] = 0.0;
  pbf[1][0] = 0.0;
  paf[1][0] = 0.0;
  pbf[2][0] = 0.0;
  paf[2][0] = 0.0;

  if(!usex){
    pbf[3][0] = 0.0;
    paf[3][0] = 0.0;
  }
  if(!usey){
    pbf[4][0] = 0.0;
    paf[4][0] = 0.0;
  }
  if(!usez){
    pbf[5][0] = 0.0;
    paf[5][0] = 0.0;
  }
  HepMatrix dpa(pbf-paf);
  //double chisq = (dpa.T()*pbe.inverse(inv_error)*dpa)[0][0];
  double chisq = (dpa.T()*inv_pbe*dpa)[0][0];
  if(inv_error){
    m_errorFlag = KF_OUTPUT_INVERSE;
    dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::trackchisq] Can not calculate matrix inverse." << std::endl;
    return -1.0;
  }
  return chisq;
}


// following three functions are modified by T.H on 2006/05/11
double
kvertexfitter::chisq_wo_ip(const bool usex, const bool usey, const bool usez, const bool dont_use_chash)
{
  static double chisq = -1.0;
  static bool _usex=false, _usey=false, _usez=false;
  if((!dont_use_chash)&&chisq>=0.0&&_usex==usex&&_usey==usey&&_usez==usez){
    return chisq; /* use chash */
  }

	if( !(usex==false && usey==false && usez==true) ){
		dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::chisq_wo_ip] usex/usey/usez must be F/F/T" << std::endl;
		exit(1);
	}

  chisq = 0.0;
  if(m_trackNum==0){
    dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::chisq_wo_ip] Error : # of Tracks is 0." << std::endl;
    m_errorFlag = KF_TRACK_SIZE;
    return -1.0;
  }
  if(m_errorFlag!=KF_NO_ERROR){
    dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::chisq_wo_ip] Fit is failed.("
	      << m_errorFlag << ")"<< std::endl;
    return -1.0;
  }
  for(unsigned int itr=0;itr<m_trackNum;++itr){
    const double i_chisq = trackchisq(itr, usex, usey, usez);
    if(i_chisq<0.0){
      dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::chisq_wo_ip] Error : Can not get chisq" << std::endl;
      m_errorFlag = KF_OUTPUT_INVERSE;
      dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::chisq_wo_ip] Can not calculate matrix inverse." << std::endl;
      chisq = -1.0;
      return -1.0;
    }
    chisq += i_chisq;
  }
  return chisq;
}

unsigned
kvertexfitter::dgf_wo_ip(const bool usex, const bool usey, const bool usez)
{
	if( !(usex==false && usey==false && usez==true) ){
		dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::dgf_wo_ip] usex/usey/usez must be F/F/T" << std::endl;
		exit(1);
	}

  if(m_trackNum==0){
    dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::dgf_wo_ip] Error : # of Tracks is 0." << std::endl;
    m_errorFlag = KF_TRACK_SIZE;
    return 0;
  }
  if(m_errorFlag!=KF_NO_ERROR){
    dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::dgf_wo_ip] Fit is failed.("
	      << m_errorFlag << ")"<< std::endl;
    return 0;
  }

  return m_trackNum - 1;
}

double
kvertexfitter::cl_wo_ip(const bool usex, const bool usey, const bool usez)
{
	if( !(usex==false && usey==false && usez==true) ){
		dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::cl_wo_ip] usex/usey/usez must be F/F/T" << std::endl;
		exit(1);
	}
  if(m_trackNum==0){
    dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::cl_wo_ip] Error : # of Tracks is 0." << std::endl;
    m_errorFlag = KF_TRACK_SIZE;
    return -1.0;
  }
  if(m_errorFlag!=KF_NO_ERROR){
    dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::cl_wo_ip] Fit is failed.("
	      << m_errorFlag << ")"<< std::endl;
    return -1.0;
  }
  const double _chisq = chisq_wo_ip(usex, usey, usez, false);
  const int    _ndf   = dgf_wo_ip(usex, usey, usez);
  return chisq2Confi(_ndf, _chisq);
}

// following three functions are for backward compatibility up to 2005summer (T.H 2006/06/16)
double
kvertexfitter::chisq_wo_ip_2005summer(const bool usex, const bool usey, const bool usez, const bool dont_use_chash)
{
	return chisq_wo_ip(usex, usey, usez, dont_use_chash);
}

unsigned
kvertexfitter::dgf_wo_ip_2005summer(const bool usex, const bool usey, const bool usez)
{
	if( !(usex==false && usey==false && usez==true) ){
		dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::dgf_wo_ip_2005summer] usex/usey/usez must be F/F/T" << std::endl;
		exit(1);
	}

  if(m_trackNum==0){
    dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::dgf_wo_ip_2005summer] Error : # of Tracks is 0." << std::endl;
    m_errorFlag = KF_TRACK_SIZE;
    return 0;
  }
  if(m_errorFlag!=KF_NO_ERROR){
    dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::dgf_wo_ip_2005summer] Fit is failed.("
	      << m_errorFlag << ")"<< std::endl;
    return 0;
  }

  if( m_tube ){
		dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::dgf_wo_ip_2005summer] You request dgf w/o IP with 2005summer version" << std::endl;
		dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::dgf_wo_ip_2005summer] But no IP tube was prepared then" << std::endl;
		dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::dgf_wo_ip_2005summer] Are you sure being correct?" << std::endl;
	}
  return 2 * m_trackNum;
}

double
kvertexfitter::cl_wo_ip_2005summer(const bool usex, const bool usey, const bool usez)
{
	if( !(usex==false && usey==false && usez==true) ){
		dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::cl_wo_ip_2005summer] usex/usey/usez must be F/F/T" << std::endl;
		exit(1);
	}
  if(m_trackNum==0){
    dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::cl_wo_ip_2005summer] Error : # of Tracks is 0." << std::endl;
    m_errorFlag = KF_TRACK_SIZE;
    return -1.0;
  }
  if(m_errorFlag!=KF_NO_ERROR){
    dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::cl_wo_ip_2005summer] Fit is failed.("
	      << m_errorFlag << ")"<< std::endl;
    return -1.0;
  }
  const double _chisq = chisq_wo_ip(usex, usey, usez, false);
  const int    _ndf   = dgf_wo_ip_2005summer(usex, usey, usez);
  return chisq2Confi(_ndf, _chisq);
}

// following two functions are for 2010summer (T.H 2010/05/17)
double
kvertexfitter::chisq_tracks(void)
{
  double chisq = 0.0;

  if(m_trackNum==0){
    dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::chisq_wo_ip] Error : # of Tracks is 0." << std::endl;
    m_errorFlag = KF_TRACK_SIZE;
    return -1.0;
  }
  if(m_errorFlag!=KF_NO_ERROR){
    dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::chisq_wo_ip] Fit is failed.("
	      << m_errorFlag << ")"<< std::endl;
    return -1.0;
  }

  for(unsigned int itr=0;itr<m_trackNum;++itr){
    const double i_chisq = this->chisq(itr);
    chisq += i_chisq;
  }

  return chisq;
}

unsigned
kvertexfitter::dgf_tracks(void)
{
  if(m_trackNum==0){
    dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::dgf_wo_ip] Error : # of Tracks is 0." << std::endl;
    m_errorFlag = KF_TRACK_SIZE;
    return 0;
  }
  if(m_errorFlag!=KF_NO_ERROR){
    dout(Debugout::ERR,"kvertexfitter") << "[kvertexfitter::dgf_wo_ip] Fit is failed.("
	      << m_errorFlag << ")"<< std::endl;
    return 0;
  }

  return m_trackNum*2 - 2;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
