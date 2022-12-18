//
// $Id: kmassvertexfitter.cc 9997 2007-02-23 23:59:43Z katayama $
//
// $Log$
// Revision 1.15  2002/04/22 11:08:06  jtanaka
// bug fix for virtual functions (Kataoka-san's info.).
//
// Revision 1.14  2002/03/27 23:41:57  jtanaka
// Add new fitter to the mass-constraint fit.
//
// Revision 1.13  2000/06/07 12:32:05  jtanaka
// add m_errorMsgFlag to control "dump error message".
// Default is ON. errorMsg(0) --> OFF.
//
// Revision 1.12  2000/04/13 12:35:57  katayama
// Added std:: to cout,cerr,endl etc.
//
// Revision 1.11  2000/03/07 17:09:01  jtanaka
// bug fixed, use "similarity" etc.
//
// Revision 1.10  2000/03/07 10:52:40  katayama
// compatibility with CC5.0
//
// Revision 1.9  1999/04/08 12:48:07  jtanaka
// remove "exit(1);"
//
// Revision 1.8  1999/03/29 05:39:18  jtanaka
// new class structure and new+old interfaces
//
// Revision 1.7  1998/07/29 00:00:08  katayama
// Less warnings with g++ -Wall
//
// Revision 1.6  1998/07/23 12:31:36  jtanaka
// delete parameter --> delete [] parameter
//
// Revision 1.5  1998/07/23 11:56:31  jtanaka
// modify a bug(use new and delete).
//
// Revision 1.4  1998/01/22 03:20:15  jtanaka
// Updated from Tanaka san. New Interface etc.
//
// Revision 1.3  1997/11/12 08:03:49  katayama
// Minor mods from Tanaka san
//
// Revision 1.2  1997/10/24 07:05:16  katayama
// Updated from Tanaka san. Bug fixes etc.
//
// Revision 1.1  1997/10/04 05:30:07  katayama
// New from Tanaka san
//
//
#include "belle.h"
#include "kfitter/kmassvertexfitter.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

kmassvertexfitter::kmassvertexfitter(void)
: m_vertex_b(0.,0.,0.),
  m_errVertex_a(3,0)
{
  m_necessaryTrackNum = 2;
  m_V_E = HepMatrix(3,3,0);
  m_v   = HepMatrix(3,1,0);
  m_v_a = HepMatrix(3,1,0);
}

kmassvertexfitter::~kmassvertexfitter()
{ 
  m_errVertexTrack_a.erase(m_errVertexTrack_a.begin(),m_errVertexTrack_a.end());
  m_isFixMass.erase(m_isFixMass.begin(),m_isFixMass.end());
}

void
kmassvertexfitter::initialVertex(const HepPoint3D &v)
{
  m_vertex_b = v;
}

void 
kmassvertexfitter::invariantMass(const double m)
{
  m_invariantMass = m;
}

unsigned
kmassvertexfitter::m_setInputMatrix(void)
{
  if(m_trackNum > KF_MAX_TRACK_NUMBER){
    m_errorFlag = KF_INPUT_TRACK_SIZE;
    return m_errorFlag;
    // dout(Debugout::DDEBUG,"kmassvertexfitter") << "(kmassvertexfitter): Too Many Tracks!" << std::endl;
    // dout(Debugout::DDEBUG,"kmassvertexfitter") << "E-mail to jtanaka@hep.phys.s.u-tokyo.ac.jp." << std::endl;
    //exit(1);
  }

  if(m_isFixMass.size() == 0){
    // If no information whether mass is fixed or free,
    // we consider all tracks as those with the fixed mass.
    for(unsigned int i=0;i<m_trackNum;++i){
      fixMass();
    }
  }else if(m_trackNum != m_isFixMass.size()){
    m_errorFlag = KF_INPUT_TRACK_SIZE;
    return m_errorFlag;    
  }

  unsigned index(0);
  //HepMatrix    tmp_al_0(KF_NUM6*m_trackNum,1,0);
  //HepSymMatrix tmp_V_al_0(KF_NUM6*m_trackNum,0);
  //HepMatrix    tmp_property(m_trackNum,3,0);
  m_al_0     = HepMatrix(KF_NUM7*m_trackNum,1,0);
  m_property = HepMatrix(m_trackNum,3,0);
  m_V_al_0   = HepSymMatrix(KF_NUM7*m_trackNum,0);
  for(std::vector<kfitterparticle>::iterator it = m_plist.begin(), endIt = m_plist.end();
      it != endIt; ++it){
    //momentum x,y,z and position x,y,z
    //for(unsigned j=0;j<KF_NUM6;++j)
    //tmp_al_0[index*KF_NUM6+j][0] = it->getFitParameter(j,KF_BEFORE_FIT);
    m_al_0[index*KF_NUM7+0][0] = it->momentum(KF_BEFORE_FIT).x();
    m_al_0[index*KF_NUM7+1][0] = it->momentum(KF_BEFORE_FIT).y();
    m_al_0[index*KF_NUM7+2][0] = it->momentum(KF_BEFORE_FIT).z();
    m_al_0[index*KF_NUM7+3][0] = it->momentum(KF_BEFORE_FIT).t();
    m_al_0[index*KF_NUM7+4][0] = it->position(KF_BEFORE_FIT).x();
    m_al_0[index*KF_NUM7+5][0] = it->position(KF_BEFORE_FIT).y();
    m_al_0[index*KF_NUM7+6][0] = it->position(KF_BEFORE_FIT).z();
    //these error
    //tmp_V_al_0.sub(index*KF_NUM6+1, it->getFitError(KF_BEFORE_FIT));
    m_V_al_0.sub(index*KF_NUM7+1, it->error(KF_BEFORE_FIT));
    //charge , mass , a
    m_property[index][0] =  it->charge();
    m_property[index][1] =  it->mass();
    m_property[index][2] = -KF_PHOTON_VELOCITY*m_magneticField*it->charge();
    ++index;
  }

  //error between track and track
  //m_V_al_0 = tmp_V_al_0;
  if(m_correlationFlag == 1){
    m_errorFlag = m_setCorrelation();
    if(m_errorFlag != KF_NO_ERROR)return m_errorFlag;
  }
  
  //vertex
  m_v_a[0][0] = m_vertex_b.x();
  m_v_a[1][0] = m_vertex_b.y();
  m_v_a[2][0] = m_vertex_b.z();

  //set member matrix 
  //m_al_0     = tmp_al_0;
  m_al_1     = m_al_0;
  //m_property = tmp_property;

  //define size of matrix
  //HepMatrix tmp_Matrix(KF_NUM6*m_trackNum,KF_NUM6*m_trackNum,0);
  //m_V_al_1     = tmp_Matrix;
  //m_D          = tmp_Matrix.sub(1,m_trackNum*2+1,1,KF_NUM6*m_trackNum);
  //m_E          = tmp_Matrix.sub(1,m_trackNum*2+1,1,3);
  //m_d          = tmp_Matrix.sub(1,m_trackNum*2+1,1,1);
  //m_V_D        = tmp_Matrix.sub(1,m_trackNum*2+1,1,m_trackNum*2+1);
  //m_lam        = tmp_Matrix.sub(1,m_trackNum*2+1,1,1);
  //m_lam0       = tmp_Matrix.sub(1,m_trackNum*2+1,1,1);
  //m_V_Dt       = tmp_Matrix.sub(1,m_trackNum*2+1,1,m_trackNum*2+1);
  //m_Cov_v_al_1 = tmp_Matrix.sub(1,3,1,KF_NUM6*m_trackNum);
  
  m_V_al_1     = HepMatrix(KF_NUM7*m_trackNum,KF_NUM7*m_trackNum,0);
  m_D          = m_V_al_1.sub(1,m_trackNum*2+1,1,KF_NUM7*m_trackNum);
  m_E          = m_V_al_1.sub(1,m_trackNum*2+1,1,3);
  m_d          = m_V_al_1.sub(1,m_trackNum*2+1,1,1);
  m_V_D        = m_V_al_1.sub(1,m_trackNum*2+1,1,m_trackNum*2+1);
  m_lam        = m_V_al_1.sub(1,m_trackNum*2+1,1,1);
  m_lam0       = m_V_al_1.sub(1,m_trackNum*2+1,1,1);
  m_V_Dt       = m_V_al_1.sub(1,m_trackNum*2+1,1,m_trackNum*2+1);
  m_Cov_v_al_1 = m_V_al_1.sub(1,3,1,KF_NUM7*m_trackNum);

  return m_errorFlag;
}

unsigned
kmassvertexfitter::m_setInputSubMatrix(void)
{
  //vertex
  for(unsigned i=0;i<3;++i){
    m_v[i][0] = m_v_a[i][0];
  }
  return m_errorFlag;
}

unsigned
kmassvertexfitter::m_makeCoreMatrix(void)
{
  //Mass Constraint
  HepMatrix al_1_prime(m_al_1);
  HepMatrix Sum_al_1(4,1,0);
  double   *energy = new double [m_trackNum];
  double    a;

  for(unsigned i=0;i<m_trackNum;++i){
    a = m_property[i][2];
    
    al_1_prime[i*KF_NUM7+0][0] -= a*(m_v_a[1][0]-al_1_prime[i*KF_NUM7+5][0]);
    al_1_prime[i*KF_NUM7+1][0] += a*(m_v_a[0][0]-al_1_prime[i*KF_NUM7+4][0]);

    energy[i] = sqrt(al_1_prime[i*KF_NUM7+0][0]*al_1_prime[i*KF_NUM7+0][0]+
		     al_1_prime[i*KF_NUM7+1][0]*al_1_prime[i*KF_NUM7+1][0]+
		     al_1_prime[i*KF_NUM7+2][0]*al_1_prime[i*KF_NUM7+2][0]+
		     m_property[i][1]*m_property[i][1]); 
  }
  for(unsigned i=0;i<m_trackNum;++i){
    if(m_isFixMass[i] == KF_FIX_MASS){
      Sum_al_1[3][0] += energy[i];
    }else{
      Sum_al_1[3][0] += al_1_prime[i*KF_NUM7+3][0];
    }
    for(unsigned j=0;j<3;++j){
      Sum_al_1[j][0] += al_1_prime[i*KF_NUM7+j][0];
    }
  }

  m_d[2*m_trackNum][0] = Sum_al_1[3][0]*Sum_al_1[3][0]-Sum_al_1[0][0]*Sum_al_1[0][0]
    -                    Sum_al_1[1][0]*Sum_al_1[1][0]-Sum_al_1[2][0]*Sum_al_1[2][0]
    -                    m_invariantMass*m_invariantMass;

  double Sum_a(0.), Sum_tmpx(0.), Sum_tmpy(0.);
  for(unsigned i=0;i<m_trackNum;++i){
    if(energy[i] == 0.){
      m_errorFlag = KF_DIV_ZERO;
      delete [] energy;
      return m_errorFlag;
    }
    a = m_property[i][2];

    if(m_isFixMass[i] == KF_FIX_MASS){
      double invE = 1./energy[i];
      m_D[2*m_trackNum][i*KF_NUM7+0] = 2.*(Sum_al_1[3][0]*al_1_prime[i*KF_NUM7+0][0]*invE-Sum_al_1[0][0]);      
      m_D[2*m_trackNum][i*KF_NUM7+1] = 2.*(Sum_al_1[3][0]*al_1_prime[i*KF_NUM7+1][0]*invE-Sum_al_1[1][0]);      
      m_D[2*m_trackNum][i*KF_NUM7+2] = 2.*(Sum_al_1[3][0]*al_1_prime[i*KF_NUM7+2][0]*invE-Sum_al_1[2][0]);
      m_D[2*m_trackNum][i*KF_NUM7+3] = 0.;
      m_D[2*m_trackNum][i*KF_NUM7+4] =-2.*(Sum_al_1[3][0]*al_1_prime[i*KF_NUM7+1][0]*invE-Sum_al_1[1][0])*a;
      m_D[2*m_trackNum][i*KF_NUM7+5] = 2.*(Sum_al_1[3][0]*al_1_prime[i*KF_NUM7+0][0]*invE-Sum_al_1[0][0])*a;
      m_D[2*m_trackNum][i*KF_NUM7+6] = 0.;
      Sum_tmpx += al_1_prime[i*KF_NUM7+0][0]*invE*a;
      Sum_tmpy += al_1_prime[i*KF_NUM7+1][0]*invE*a;
    }else{
      m_D[2*m_trackNum][i*KF_NUM7+0] = -2.*Sum_al_1[0][0];
      m_D[2*m_trackNum][i*KF_NUM7+1] = -2.*Sum_al_1[1][0];
      m_D[2*m_trackNum][i*KF_NUM7+2] = -2.*Sum_al_1[2][0];
      m_D[2*m_trackNum][i*KF_NUM7+3] =  2.*Sum_al_1[3][0];
      m_D[2*m_trackNum][i*KF_NUM7+4] =  2.*Sum_al_1[1][0]*a;
      m_D[2*m_trackNum][i*KF_NUM7+5] = -2.*Sum_al_1[0][0]*a;
      m_D[2*m_trackNum][i*KF_NUM7+6] =  0.;
    }
    Sum_a    += a;
  }
  
  //m_E
  m_E[2*m_trackNum][0] = -2.*Sum_al_1[1][0]*Sum_a + 2.*Sum_al_1[3][0]*Sum_tmpy;
  m_E[2*m_trackNum][1] =  2.*Sum_al_1[0][0]*Sum_a - 2.*Sum_al_1[3][0]*Sum_tmpx;
  m_E[2*m_trackNum][2] =  0.;

  // 2000/03/07
  //vertex fit
  double px,py,pz,x,y,z;
  double pt,invPt,invPt2,dlx,dly,dlz,a1,a2,r2d2,B,Rx,Ry,S,U;
  double sininv,sqrtag;
  
  for(unsigned i=0;i<m_trackNum;++i){
    px = m_al_1[i*KF_NUM7+0][0];
    py = m_al_1[i*KF_NUM7+1][0];
    pz = m_al_1[i*KF_NUM7+2][0];
    x  = m_al_1[i*KF_NUM7+4][0];
    y  = m_al_1[i*KF_NUM7+5][0];
    z  = m_al_1[i*KF_NUM7+6][0];
    a  = m_property[i][2];

    pt =  sqrt(px*px+py*py);

    if(pt == 0.){
      m_errorFlag = KF_DIV_ZERO;
      delete [] energy;
      return m_errorFlag;
    }    

    invPt   =  1./pt;
    invPt2  =  invPt*invPt;
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
	delete [] energy;
	return m_errorFlag;
      }
      //sin^(-1)(B)
      sininv = asin(B);
      double tmp0 = 1.0-B*B;
      if(tmp0 == 0.){
	m_errorFlag = KF_DIV_ZERO;
	delete [] energy;
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
    m_D[i*2+0][i*KF_NUM7+0] =  dly;
    m_D[i*2+0][i*KF_NUM7+1] = -dlx;
    m_D[i*2+0][i*KF_NUM7+2] =  0.0;
    m_D[i*2+0][i*KF_NUM7+4] =  py+a*dlx;
    m_D[i*2+0][i*KF_NUM7+5] = -px+a*dly;
    m_D[i*2+0][i*KF_NUM7+6] =  0.0;
    m_D[i*2+1][i*KF_NUM7+0] = -pz*pt*S*Rx+U*px*invPt;
    m_D[i*2+1][i*KF_NUM7+1] = -pz*pt*S*Ry+U*py*invPt;
    if(a != 0.){
      m_D[i*2+1][i*KF_NUM7+2] = -sininv*pt/a;
    }else{
      m_D[i*2+1][i*KF_NUM7+2] = -a2*invPt;
    }
    m_D[i*2+1][i*KF_NUM7+4] =  px*pz*pt*S;
    m_D[i*2+1][i*KF_NUM7+5] =  py*pz*pt*S;
    m_D[i*2+1][i*KF_NUM7+6] = -pt;
      
    //E
    m_E[i*2+0][0] = -py-a*dlx;
    m_E[i*2+0][1] =  px-a*dly;
    m_E[i*2+0][2] =  0.0;
    m_E[i*2+1][0] = -px*pz*pt*S;
    m_E[i*2+1][1] = -py*pz*pt*S;
    m_E[i*2+1][2] =  pt;
  }

  delete [] energy;
  return m_errorFlag;
}

HepPoint3D
kmassvertexfitter::vertex(const unsigned f)
{
  if(f == KF_AFTER_FIT){
    return m_vertex_a;
  }else if(f == KF_BEFORE_FIT){
    return m_vertex_b;
  }
  m_errorFlag = KF_OUTPUT_OUT_RANGE;
  if(m_errorMsgFlag)dout(Debugout::ERR,"kmassvertexfitter") << "(kmassvertexfitter): Out of Range!" << std::endl;
  return HepPoint3D();
}

HepSymMatrix 
kmassvertexfitter::errVertex(void)
{
  return m_errVertex_a;
}

HepMatrix
kmassvertexfitter::errVertexTrack(const unsigned n)
{
  if(n < m_trackNum){
    return m_errVertexTrack_a[n];
  }else{
    m_errorFlag = KF_OUTPUT_TRACK_NUM;
    if(m_errorMsgFlag)dout(Debugout::ERR,"kmassvertexfitter") << "(kmassvertexfitter): Out of Range!" << std::endl;
    return HepMatrix(3,KF_NUM7,0);
  }
}

double
kmassvertexfitter::invariantMass(void)
{
  return m_invariantMass ;
}

unsigned
kmassvertexfitter::m_setOutputMatrix(void)
{
  Hep3Vector h3v;
  unsigned index(0);
  for(std::vector<kfitterparticle>::iterator it = m_plist.begin(),
      endIt = m_plist.end();
      it != endIt; ++it){
    kfitterparticle &pdata = *it;
    //tracks
    //momentum
    h3v.setX(m_al_1[index*KF_NUM7+0][0]);
    h3v.setY(m_al_1[index*KF_NUM7+1][0]);
    h3v.setZ(m_al_1[index*KF_NUM7+2][0]);
    if(m_isFixMass[index] == KF_FIX_MASS){
      pdata.momentum(HepLorentzVector(h3v, sqrt(h3v.mag2()+pdata.mass()*pdata.mass())),KF_AFTER_FIT);
    }else{
      pdata.momentum(HepLorentzVector(h3v, m_al_1[index*KF_NUM7+3][0]), KF_AFTER_FIT);
    }
    //position
    pdata.position(HepPoint3D(m_al_1[index*KF_NUM7+4][0],
                              m_al_1[index*KF_NUM7+5][0],
                              m_al_1[index*KF_NUM7+6][0]),KF_AFTER_FIT);
    //error of the tracks
    pdata.error(m_makeError3(pdata.momentum(),
			     m_V_al_1.sub(index    *KF_NUM7+1,
					  (index+1)*KF_NUM7,
					  index    *KF_NUM7+1,
					  (index+1)*KF_NUM7),
			     m_isFixMass[index]),KF_AFTER_FIT);

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
    HepMatrix hm(3,KF_NUM7,0);
    for(unsigned j=0;j<3;++j){
      for(unsigned k=0;k<KF_NUM7;++k){
	hm[j][k] = m_Cov_v_al_1[j][KF_NUM7*i+k];
      }
    }
    if(m_isFixMass[i] == KF_FIX_MASS){
      m_errVertexTrack_a.push_back(m_makeError4(m_plist[i].momentum(),hm));
    }else{
      m_errVertexTrack_a.push_back(hm);
    }
  }

  return m_errorFlag;
}

unsigned
kmassvertexfitter::m_calDgf(void)
{
  m_dgf = 2*m_trackNum-3+1;
  return m_errorFlag;
}

void
kmassvertexfitter::fixMass(void)
{
  m_isFixMass.push_back(KF_FIX_MASS);
}

void
kmassvertexfitter::unfixMass(void)
{
  m_isFixMass.push_back(KF_UNFIX_MASS);
}

unsigned
kmassvertexfitter::m_setCorrelation(void){
  if(m_correlation_b.size() != static_cast<unsigned>((m_trackNum*(m_trackNum-1))/2)){
    m_errorFlag = KF_INPUT_CORRE_SIZE;
    return m_errorFlag;
  }

  unsigned row(0), col(0);

  for(std::vector<HepMatrix>::iterator it = m_correlation_b.begin(),
	endIt = m_correlation_b.end();
      it != endIt; ++it){
    HepMatrix &hm = *it;

    //counter
    ++row;
    if(row == m_trackNum){
      ++col;
      row = col + 1;
    }

    unsigned ii(0), jj(0);    
    for(unsigned i=KF_NUM7*row; i<KF_NUM7*(row+1); ++i){
      for(unsigned j=KF_NUM7*col; j<KF_NUM7*(col+1); ++j){
        m_V_al_0[i][j] = hm[ii][jj];
        ++jj;
      }
      jj = 0;
      ++ii;
    }
  }

  return m_errorFlag;
}

double
kmassvertexfitter::chisq(void){
  return m_chisq;
}

double
kmassvertexfitter::chisq(const unsigned n){
  if(n < m_trackNum){
    if(m_isFixMass[n] == KF_FIX_MASS){
      HepMatrix da(m_plist[n].getFitParameter(KF_BEFORE_FIT)-
		   m_plist[n].getFitParameter(KF_AFTER_FIT));
      int errInverse(0);
      double chiSq((da.T()*(m_plist[n].getFitError(KF_BEFORE_FIT).inverse(errInverse))*da)[0][0]);
      if(errInverse){
	m_errorFlag = KF_OUTPUT_INVERSE;
	return -1.;
      }else return chiSq;
    }else{
      HepMatrix da(m_plist[n].mompos(KF_BEFORE_FIT)-
		   m_plist[n].mompos(KF_AFTER_FIT));
      int errInverse(0);
      double chiSq((da.T()*(m_plist[n].error(KF_BEFORE_FIT).inverse(errInverse))*da)[0][0]);
      if(errInverse){
	m_errorFlag = KF_OUTPUT_INVERSE;
	return -1.;
      }else return chiSq;
    }
  }else{
    m_errorFlag = KF_OUTPUT_TRACK_NUM;
    if(m_errorMsgFlag)dout(Debugout::ERR,"kmassvertexfitter") << "(kmassvertexfitter): Out of Range!" << std::endl;
    return -1.;
  }
}

//  void 
//  kmassvertexfitter::correlation(const HepMatrix &e){
//    m_correlation_b.push_back(e);
//    m_correlationFlag = 1;
//    if(e.num_row() != (int)KF_NUM7)
//      m_errorFlag = KF_INPUT_MATRIX_SIZE;
//  }

//  void 
//  kmassvertexfitter::correlation(void){
//    m_correlation_b.push_back(HepMatrix(KF_NUM7,
//  				      KF_NUM7,0));
//    m_correlationFlag = 1;
//  }

HepMatrix
kmassvertexfitter::correlation(const unsigned n, const unsigned m, const unsigned flag){
  if(flag == KF_AFTER_FIT){
    if(n < m_trackNum && m < m_trackNum){
      return m_makeError3(momentum(n),
			  momentum(m),
			  m_V_al_1.sub(KF_NUM7*n+1,KF_NUM7*(n+1),
				       KF_NUM7*m+1,KF_NUM7*(m+1)),
			  m_isFixMass[n],
			  m_isFixMass[m]);
    }else{
      m_errorFlag = KF_OUTPUT_TRACK_NUM;
      if(m_errorMsgFlag)dout(Debugout::ERR,"kmassvertexfitter") << "(kmassvertexfitter): Out of Range!" << std::endl;
      return HepMatrix(KF_NUM7,KF_NUM7,0);
    }
  }else{
    return kfitterbase::correlation(n,m,flag);
  }
}

void
kmassvertexfitter::correlation(const HepMatrix &m)
{
  kfitterbase::correlation(m);
}

void
kmassvertexfitter::correlation(void)
{
  kfitterbase::correlation();
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
