//
// $Id: kmassfitter.cc 9997 2007-02-23 23:59:43Z katayama $
//
// $Log$
// Revision 1.20  2002/04/22 11:08:05  jtanaka
// bug fix for virtual functions (Kataoka-san's info.).
//
// Revision 1.19  2002/03/27 23:41:57  jtanaka
// Add new fitter to the mass-constraint fit.
//
// Revision 1.18  2002/03/12 14:55:09  jtanaka
// bug fix, which doesn't affect the result of the physics.
//
// Revision 1.17  2002/02/23 17:56:53  katayama
// Added (int) for int/unsigned comparison
//
// Revision 1.16  2000/07/01 17:59:33  jtanaka
// bug : A "fitWithError" mode did not work well
//       because an error of vertex was not set.
//
// Revision 1.15  2000/06/07 12:32:05  jtanaka
// add m_errorMsgFlag to control "dump error message".
// Default is ON. errorMsg(0) --> OFF.
//
// Revision 1.14  2000/04/13 12:35:57  katayama
// Added std:: to cout,cerr,endl etc.
//
// Revision 1.13  2000/03/08 07:44:12  jtanaka
// add new member function.
//
// Revision 1.12  2000/03/07 17:09:01  jtanaka
// bug fixed, use "similarity" etc.
//
// Revision 1.11  2000/03/07 10:52:40  katayama
// compatibility with CC5.0
//
// Revision 1.10  1999/04/08 12:48:07  jtanaka
// remove "exit(1);"
//
// Revision 1.9  1999/03/29 05:39:16  jtanaka
// new class structure and new+old interfaces
//
// Revision 1.8  1998/07/29 00:00:03  katayama
// Less warnings with g++ -Wall
//
// Revision 1.7  1998/07/23 12:31:35  jtanaka
// delete parameter --> delete [] parameter
//
// Revision 1.6  1998/07/23 11:56:32  jtanaka
// modify a bug(use new and delete).
//
// Revision 1.5  1998/05/07 04:11:21  jtanaka
// version 2.1
//
// Revision 1.4  1998/01/22 03:20:14  jtanaka
// Updated from Tanaka san. New Interface etc.
//
// Revision 1.3  1997/11/12 08:03:48  katayama
// Minor mods from Tanaka san
//
// Revision 1.2  1997/10/24 07:05:15  katayama
// Updated from Tanaka san. Bug fixes etc.
//
// Revision 1.1  1997/10/04 05:30:07  katayama
// New from Tanaka san
//
//
#include "belle.h"
#include "kfitter/kmassfitter.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

kmassfitter::kmassfitter(void)
: m_errVertexTrackFlag(0),
  m_fitIncludingVertex(0),
  m_atDecayPoint(KF_AT_DECAY_POINT)
{
  m_necessaryTrackNum = 2;
  m_d   = HepMatrix(1,1,0);
  m_V_D = HepMatrix(1,1,0);
  m_lam = HepMatrix(1,1,0);
  m_errVertex_a = HepSymMatrix(3,0);
}

kmassfitter::~kmassfitter(void)
{
  m_errVertexTrack_b.erase(m_errVertexTrack_b.begin(),m_errVertexTrack_b.end());
  m_errVertexTrack_a.erase(m_errVertexTrack_a.begin(),m_errVertexTrack_a.end());
  m_isFixMass.erase(m_isFixMass.begin(),m_isFixMass.end());
}

void 
kmassfitter::vertex(const HepPoint3D &v)
{
  m_vertex_b = v;
}

void 
kmassfitter::errVertex(const HepSymMatrix &eV)
{
  if(eV.num_row() != 3){
    m_errorFlag = KF_INPUT_MATRIX_SIZE;
  }else{
    m_errVertex_b = eV;
    m_fitIncludingVertex = 1;
  }
}

void 
kmassfitter::errVertexTrack(const HepMatrix &eVT)
{
  m_errVertexTrack_b.push_back(eVT);
  m_errVertexTrackFlag = 1;
  m_fitIncludingVertex = 1;
  if(eVT.num_row() != 3 || eVT.num_col() != (int)KF_NUM7)
    m_errorFlag = KF_INPUT_MATRIX_SIZE;
}

void 
kmassfitter::errVertexTrack(void)
{
  m_errVertexTrack_b.push_back(HepMatrix(3,KF_NUM7,0));
  m_errVertexTrackFlag = 1;
  m_fitIncludingVertex = 1;
}

void 
kmassfitter::invariantMass(const double m)
{
  m_invariantMass = m;
}

void 
kmassfitter::atDecayPoint(void)
{
  m_atDecayPoint = KF_AT_DECAY_POINT;
}

void 
kmassfitter::notDecayPoint(void)
{
  m_atDecayPoint = KF_NOT_DECAY_POINT;
}

unsigned
kmassfitter::m_setInputMatrix(void)
{
  if(m_trackNum > KF_MAX_TRACK_NUMBER){
    m_errorFlag = KF_INPUT_TRACK_SIZE;
    return m_errorFlag;
    // dout(Debugout::DDEBUG,"kmassfitter") << "(kmassfitter): Too Many Tracks!" << std::endl;
    // dout(Debugout::DDEBUG,"kmassfitter") << "E-mail to jtanaka@hep.phys.s.u-tokyo.ac.jp." << std::endl;
    //exit(1);
  }

  if(m_isFixMass.size() == 0){
    // If no information whether mass is fixed or free,
    // we consider all tracks as those with the fixed mass.
    for(unsigned i=0;i<m_trackNum;++i){
      fixMass();
    }
  }else if(m_trackNum != m_isFixMass.size()){
    m_errorFlag = KF_INPUT_TRACK_SIZE;
    return m_errorFlag;    
  }

  if(m_fitIncludingVertex == 0){
    unsigned index(0);
    //HepMatrix    tmp_al_0(KF_NUM6*m_trackNum,1,0);
    //HepSymMatrix tmp_V_al_0(KF_NUM6*m_trackNum,0);
    m_al_0     = HepMatrix(KF_NUM7*m_trackNum,1,0);
    m_property = HepMatrix(m_trackNum,3,0);
    m_V_al_0   = HepSymMatrix(KF_NUM7*m_trackNum,0);
    for(std::vector<kfitterparticle>::iterator it = m_plist.begin(), endIt = m_plist.end();
	it != endIt; ++it){
      //momentum x,y,z and position x,y,z
      //for(unsigned j=0;j<KF_NUM6;++j)
      //m_al_0[index*KF_NUM6+j][0] = it->getFitParameter(j,KF_BEFORE_FIT);
      m_al_0[index*KF_NUM7+0][0] = it->momentum(KF_BEFORE_FIT).x();
      m_al_0[index*KF_NUM7+1][0] = it->momentum(KF_BEFORE_FIT).y();
      m_al_0[index*KF_NUM7+2][0] = it->momentum(KF_BEFORE_FIT).z();
      m_al_0[index*KF_NUM7+3][0] = it->momentum(KF_BEFORE_FIT).t();
      m_al_0[index*KF_NUM7+4][0] = it->position(KF_BEFORE_FIT).x();
      m_al_0[index*KF_NUM7+5][0] = it->position(KF_BEFORE_FIT).y();
      m_al_0[index*KF_NUM7+6][0] = it->position(KF_BEFORE_FIT).z();
      //these error
      //m_V_al_0.sub(index*KF_NUM6+1, it->getFitError(KF_BEFORE_FIT));
      m_V_al_0.sub(index*KF_NUM7+1, it->error(KF_BEFORE_FIT));
      //charge , mass , a
      m_property[index][0] =  it->charge();
      m_property[index][1] =  it->mass();
      m_property[index][2] = -KF_PHOTON_VELOCITY*m_magneticField*it->charge();
      ++index;
      /*
      dout(Debugout::INFO,"kmassfitter") << index << " particle " << std::endl;
      dout(Debugout::INFO,"kmassfitter") << m_al_0 << std::endl;
      dout(Debugout::INFO,"kmassfitter") << m_V_al_0 << std::endl;
      dout(Debugout::INFO,"kmassfitter") << std::endl;
      */
    }
    
    //error between track and track
    //m_V_al_0 = tmp_V_al_0;
    if(m_correlationFlag == 1){
      m_errorFlag = m_setCorrelation();
      if(m_errorFlag != KF_NO_ERROR)return m_errorFlag;
    }
    

    //set member matrix 
    //m_al_0      = tmp_al_0;
    m_al_1 = m_al_0;
    //m_property  = tmp_property;
    
    //define size of matrix
    //HepMatrix  tmp_Matrix(KF_NUM6*m_trackNum,KF_NUM6*m_trackNum,0);
    //m_V_al_1 = tmp_Matrix;
    //m_D      = tmp_Matrix.sub(1,1,1,KF_NUM6*m_trackNum);
    //HepMatrix  tmp_Matrix(KF_NUM6*m_trackNum,KF_NUM6*m_trackNum,0);
    m_V_al_1 = HepMatrix(KF_NUM7*m_trackNum,KF_NUM7*m_trackNum,0);
    m_D      = m_V_al_1.sub(1,1,1,KF_NUM7*m_trackNum);
  }else{
    //m_fitIncludingVertex == 1
    unsigned index(0);
    //HepMatrix    tmp_al_0(KF_NUM6*m_trackNum+3,1,0);
    //HepSymMatrix tmp_V_al_0(KF_NUM6*m_trackNum+3,0);
    //HepMatrix    tmp_property(m_trackNum,3,0);
    m_al_0     = HepMatrix(KF_NUM7*m_trackNum+3,1,0);
    m_property = HepMatrix(m_trackNum,3,0);
    m_V_al_0   = HepSymMatrix(KF_NUM7*m_trackNum+3,0);

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
      //tmp_property[index][0] =  it->charge();
      //tmp_property[index][1] =  it->mass();
      //tmp_property[index][2] = -KF_PHOTON_VELOCITY*m_magneticField*it->charge();
      ++index;
    }
    // vertex
    //tmp_al_0[KF_NUM6*m_trackNum+0][0] = m_vertex_b.x();
    //tmp_al_0[KF_NUM6*m_trackNum+1][0] = m_vertex_b.y();
    //tmp_al_0[KF_NUM6*m_trackNum+2][0] = m_vertex_b.z();
    //tmp_V_al_0.sub(KF_NUM6*m_trackNum+1,m_errVertex_b);
    m_al_0[KF_NUM7*m_trackNum+0][0] = m_vertex_b.x();
    m_al_0[KF_NUM7*m_trackNum+1][0] = m_vertex_b.y();
    m_al_0[KF_NUM7*m_trackNum+2][0] = m_vertex_b.z();
    m_V_al_0.sub(KF_NUM7*m_trackNum+1,m_errVertex_b);

    //error between track and track
    //m_V_al_0 = tmp_V_al_0;
    if(m_correlationFlag == 1){
      m_errorFlag = m_setCorrelation();
      if(m_errorFlag != KF_NO_ERROR)return m_errorFlag;
    }
    
    //set member matrix 
    //m_al_0      = tmp_al_0;
    m_al_1      = m_al_0;
    //m_property  = tmp_property;
    
    //define size of matrix
    //HepMatrix  tmp_Matrix(KF_NUM6*m_trackNum+3,KF_NUM6*m_trackNum+3,0);
    //m_V_al_1 = tmp_Matrix;
    //m_D      = tmp_Matrix.sub(1,1,1,KF_NUM6*m_trackNum+3);
    m_V_al_1 = HepMatrix(KF_NUM7*m_trackNum+3,KF_NUM7*m_trackNum+3,0);
    m_D      = m_V_al_1.sub(1,1,1,KF_NUM7*m_trackNum+3);
  }
  return m_errorFlag;
}

unsigned
kmassfitter::m_setCorrelation(void)
{
  // 2000/03/07
  //m_errorFlag = kfitterbase::m_setCorrelation();

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

  if(m_fitIncludingVertex == 1){
    //...error of vertex
    m_V_al_0.sub(KF_NUM7*m_trackNum+1, m_errVertex_b);
    
    //...error matrix between vertex and tracks
    if(m_errVertexTrackFlag == 1){
      if(m_errVertexTrack_b.size() != m_trackNum){
	m_errorFlag = KF_INPUT_CORRE_SIZE;
	return m_errorFlag;
      }
      unsigned i(0);
      for(std::vector<HepMatrix>::iterator it = m_errVertexTrack_b.begin(), endIt = m_errVertexTrack_b.end();
	  it != endIt; ++it){
	HepMatrix &hm = *it;
	for(unsigned j=0;j<3;++j){
	  for(unsigned k=0;k<KF_NUM7;++k){
	    m_V_al_0[j+KF_NUM7*m_trackNum][k+i*KF_NUM7] = hm[j][k];
	  }
	}
	++i;
      }
    }  
  }
  return m_errorFlag;
}

unsigned
kmassfitter::m_makeCoreMatrix(void)
{
  if(m_fitIncludingVertex == 0){
    HepMatrix al_1_prime(m_al_1);
    HepMatrix Sum_al_1(4,1,0);
    double   *energy = new double[m_trackNum];
    double    a;
    
    for(unsigned i=0;i<m_trackNum;++i){
      a = m_property[i][2];
      if(m_atDecayPoint == KF_NOT_DECAY_POINT)a = 0.;
      al_1_prime[i*KF_NUM7+0][0] -= a*(m_vertex_b.y()-al_1_prime[i*KF_NUM7+5][0]);
      al_1_prime[i*KF_NUM7+1][0] += a*(m_vertex_b.x()-al_1_prime[i*KF_NUM7+4][0]);
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
    
    m_d[0][0] = Sum_al_1[3][0]*Sum_al_1[3][0]-Sum_al_1[0][0]*Sum_al_1[0][0]
      -         Sum_al_1[1][0]*Sum_al_1[1][0]-Sum_al_1[2][0]*Sum_al_1[2][0]
      -         m_invariantMass*m_invariantMass;
    
    for(unsigned i=0;i<m_trackNum;++i){
      if(energy[i] == 0.){
	m_errorFlag = KF_DIV_ZERO;
	break;
      }
      
      a = m_property[i][2];
      if(m_atDecayPoint == KF_NOT_DECAY_POINT)a = 0.;
      
      if(m_isFixMass[i] == KF_FIX_MASS){
	double invE = 1./energy[i];
	m_D[0][i*KF_NUM7+0] = 2.*(Sum_al_1[3][0]*al_1_prime[i*KF_NUM7+0][0]*invE-Sum_al_1[0][0]);
	m_D[0][i*KF_NUM7+1] = 2.*(Sum_al_1[3][0]*al_1_prime[i*KF_NUM7+1][0]*invE-Sum_al_1[1][0]);
	m_D[0][i*KF_NUM7+2] = 2.*(Sum_al_1[3][0]*al_1_prime[i*KF_NUM7+2][0]*invE-Sum_al_1[2][0]);
	m_D[0][i*KF_NUM7+3] = 0.;
	m_D[0][i*KF_NUM7+4] =-2.*(Sum_al_1[3][0]*al_1_prime[i*KF_NUM7+1][0]*invE-Sum_al_1[1][0])*a;
	m_D[0][i*KF_NUM7+5] = 2.*(Sum_al_1[3][0]*al_1_prime[i*KF_NUM7+0][0]*invE-Sum_al_1[0][0])*a;
	m_D[0][i*KF_NUM7+6] = 0.;
      }else{
	m_D[0][i*KF_NUM7+0] = -2.*Sum_al_1[0][0];
	m_D[0][i*KF_NUM7+1] = -2.*Sum_al_1[1][0];
	m_D[0][i*KF_NUM7+2] = -2.*Sum_al_1[2][0];
	m_D[0][i*KF_NUM7+3] =  2.*Sum_al_1[3][0];
	m_D[0][i*KF_NUM7+4] =  2.*Sum_al_1[1][0]*a;
	m_D[0][i*KF_NUM7+5] = -2.*Sum_al_1[0][0]*a;
	m_D[0][i*KF_NUM7+6] =  0.;
      }
    }

    delete [] energy;
  }else{
    //m_fitIncludingVertex == 1
    HepMatrix al_1_prime(m_al_1);
    HepMatrix Sum_al_1(7,1,0);
    // Sum_al_1[0][0] = Sigma(px^prime)
    // Sum_al_1[1][0] = Sigma(py^prime)
    // Sum_al_1[2][0] = Sigma(pz^prime)
    // Sum_al_1[3][0] = Sigma(energy)
    // Sum_al_1[4][0] = Sigma(py^prime*a/energy)
    // Sum_al_1[5][0] = Sigma(px^prime*a/energy)
    // Sum_al_1[6][0] = Sigma(a)
    double   *energy = new double[m_trackNum];
    double    a;

    for(unsigned i=0;i<m_trackNum;++i){
      a = m_property[i][2];
      //if(m_atDecayPoint == KF_NOT_DECAY_POINT)a = 0.;
      al_1_prime[i*KF_NUM7+0][0] -= a*(al_1_prime[KF_NUM7*m_trackNum+1][0]-
				       al_1_prime[i*KF_NUM7+5][0]);
      al_1_prime[i*KF_NUM7+1][0] += a*(al_1_prime[KF_NUM7*m_trackNum+0][0]-
				       al_1_prime[i*KF_NUM7+4][0]);
      energy[i] = sqrt(al_1_prime[i*KF_NUM7+0][0]*al_1_prime[i*KF_NUM7+0][0]+
		       al_1_prime[i*KF_NUM7+1][0]*al_1_prime[i*KF_NUM7+1][0]+
		       al_1_prime[i*KF_NUM7+2][0]*al_1_prime[i*KF_NUM7+2][0]+
		       m_property[i][1]*m_property[i][1]);
      Sum_al_1[6][0] =+ a;
    }
    for(unsigned i=0;i<m_trackNum;++i){
      if(energy[i] == 0.){
	m_errorFlag = KF_DIV_ZERO;
	break;
      }
      if(m_isFixMass[i] == KF_FIX_MASS){
	double invE = 1./energy[i];
	Sum_al_1[3][0] += energy[i];
	Sum_al_1[4][0] += al_1_prime[i*KF_NUM7+1][0]*m_property[i][2]*invE;
	Sum_al_1[5][0] += al_1_prime[i*KF_NUM7+0][0]*m_property[i][2]*invE;
      }else{
        Sum_al_1[3][0] += al_1_prime[i*KF_NUM7+3][0];	
      }
      for(unsigned j=0;j<3;++j){
	Sum_al_1[j][0] += al_1_prime[i*KF_NUM7+j][0];
      }
    }
    
    m_d[0][0] = Sum_al_1[3][0]*Sum_al_1[3][0]-Sum_al_1[0][0]*Sum_al_1[0][0]
      -         Sum_al_1[1][0]*Sum_al_1[1][0]-Sum_al_1[2][0]*Sum_al_1[2][0]
      -         m_invariantMass*m_invariantMass;
    
    for(unsigned i=0;i<m_trackNum;++i){
      if(energy[i] == 0.){
	m_errorFlag = KF_DIV_ZERO;
	break;
      }
      
      a = m_property[i][2];
      //if(m_atDecayPoint == KF_NOT_DECAY_POINT)a = 0.;
      
      if(m_isFixMass[i] == KF_FIX_MASS){
	double invE = 1./energy[i];
	m_D[0][i*KF_NUM7+0] = 2.*(Sum_al_1[3][0]*al_1_prime[i*KF_NUM7+0][0]*invE-Sum_al_1[0][0]);
	m_D[0][i*KF_NUM7+1] = 2.*(Sum_al_1[3][0]*al_1_prime[i*KF_NUM7+1][0]*invE-Sum_al_1[1][0]);
	m_D[0][i*KF_NUM7+2] = 2.*(Sum_al_1[3][0]*al_1_prime[i*KF_NUM7+2][0]*invE-Sum_al_1[2][0]);
	m_D[0][i*KF_NUM7+3] = 0.;
	m_D[0][i*KF_NUM7+4] =-2.*(Sum_al_1[3][0]*al_1_prime[i*KF_NUM7+1][0]*invE-Sum_al_1[1][0])*a;
	m_D[0][i*KF_NUM7+5] = 2.*(Sum_al_1[3][0]*al_1_prime[i*KF_NUM7+0][0]*invE-Sum_al_1[0][0])*a;
	m_D[0][i*KF_NUM7+6] = 0.;
      }else{
	m_D[0][i*KF_NUM7+0] = -2.*Sum_al_1[0][0];
        m_D[0][i*KF_NUM7+1] = -2.*Sum_al_1[1][0];
        m_D[0][i*KF_NUM7+2] = -2.*Sum_al_1[2][0];
        m_D[0][i*KF_NUM7+3] =  2.*Sum_al_1[3][0];
        m_D[0][i*KF_NUM7+4] =  2.*Sum_al_1[1][0]*a;
        m_D[0][i*KF_NUM7+5] = -2.*Sum_al_1[0][0]*a;
        m_D[0][i*KF_NUM7+6] =  0.;
      }
    }    
    m_D[0][KF_NUM7*m_trackNum+0] = 2.*(Sum_al_1[3][0]*Sum_al_1[4][0]-
				       Sum_al_1[1][0]*Sum_al_1[6][0]);
    m_D[0][KF_NUM7*m_trackNum+1] =-2.*(Sum_al_1[3][0]*Sum_al_1[5][0]-
				       Sum_al_1[0][0]*Sum_al_1[6][0]);
    m_D[0][KF_NUM7*m_trackNum+2] = 0.;

    delete [] energy;
  }
  return m_errorFlag;
}

HepPoint3D
kmassfitter::vertex(const unsigned f)
{
  if(f == KF_AFTER_FIT/* && m_fitIncludingVertex == 1*/){
    return m_vertex_a;
  }else if(f == KF_BEFORE_FIT){
    return m_vertex_b;
  }
  m_errorFlag = KF_OUTPUT_OUT_RANGE;
  if(m_errorMsgFlag)dout(Debugout::ERR,"kmassfitter") << "(kmassfitter): Out of Range!" << std::endl;
  return HepPoint3D();
}

HepSymMatrix 
kmassfitter::errVertex(const unsigned f)
{
  if(f == KF_AFTER_FIT && m_fitIncludingVertex == 1){
    return m_errVertex_a;
  }else if(f == KF_BEFORE_FIT){
    return m_errVertex_b;
  }
  m_errorFlag = KF_OUTPUT_OUT_RANGE;
  if(m_errorMsgFlag)dout(Debugout::ERR,"kmassfitter") << "(kmassfitter): Out of Range!" << std::endl;
  return HepSymMatrix(3,0);
}

HepMatrix
kmassfitter::errVertexTrack(const unsigned n, const unsigned f)
{
  if(n < m_trackNum){
    if(f == KF_AFTER_FIT && m_fitIncludingVertex == 1){
      return m_errVertexTrack_a[n];
    }else if(f == KF_BEFORE_FIT){
      return m_errVertexTrack_b[n];
    }
    m_errorFlag = KF_OUTPUT_OUT_RANGE;
    if(m_errorMsgFlag)dout(Debugout::ERR,"kmassfitter") << "(kmassfitter): Out of Range!" << std::endl;
    return HepMatrix(3,KF_NUM7,0);
  }else{
    m_errorFlag = KF_OUTPUT_TRACK_NUM;
    if(m_errorMsgFlag)dout(Debugout::ERR,"kmassfitter") << "(kmassfitter): Out of Range!" << std::endl;
    return HepMatrix(3,KF_NUM7,0);
  }
}

double
kmassfitter::invariantMass(void)
{
  return m_invariantMass;
}

unsigned
kmassfitter::decayPoint(void)
{
  return m_atDecayPoint;
}

unsigned
kmassfitter::fitWithVertex(void)
{
  return m_fitIncludingVertex;
}

unsigned
kmassfitter::m_setOutputMatrix(void)
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
      pdata.momentum(HepLorentzVector(h3v, sqrt(h3v.mag2()+pdata.mass()*pdata.mass())), KF_AFTER_FIT);
    }else{
      pdata.momentum(HepLorentzVector(h3v, m_al_1[index*KF_NUM7+3][0]), KF_AFTER_FIT);
    }
    //position
    pdata.position(HepPoint3D(m_al_1[index*KF_NUM7+4][0],
			      m_al_1[index*KF_NUM7+5][0],
			      m_al_1[index*KF_NUM7+6][0]), KF_AFTER_FIT);
    //error of the tracks
    pdata.error(m_makeError3(pdata.momentum(),
			     m_V_al_1.sub(index    *KF_NUM7+1,
					  (index+1)*KF_NUM7,
					  index    *KF_NUM7+1,
					  (index+1)*KF_NUM7),
			     m_isFixMass[index]), KF_AFTER_FIT);

    if(m_errorFlag != KF_NO_ERROR)break;
    ++index;
  } 
  if(m_fitIncludingVertex == 1){
    //vertex
    m_vertex_a.setX(m_al_1[KF_NUM7*m_trackNum+0][0]);
    m_vertex_a.setY(m_al_1[KF_NUM7*m_trackNum+1][0]);
    m_vertex_a.setZ(m_al_1[KF_NUM7*m_trackNum+2][0]);
    //error of the vertex
    for(unsigned i=0;i<3;++i){
      for(unsigned j=i;j<3;++j){
	m_errVertex_a[i][j] = m_V_al_1[KF_NUM7*m_trackNum+i][KF_NUM7*m_trackNum+j];
      }
    }
    //error between vertex and tracks
    for(unsigned i=0;i<m_trackNum;++i){
      HepMatrix hm(3,KF_NUM7,0);
      for(unsigned j=0;j<3;++j){
	for(unsigned k=0;k<KF_NUM7;++k){
	  hm[j][k] = m_V_al_1[KF_NUM7*m_trackNum+j][KF_NUM7*i+k];
	}
      }
      if(m_isFixMass[i] == KF_FIX_MASS){
	m_errVertexTrack_a.push_back(m_makeError4(m_plist[i].momentum(),hm));
      }else{
	m_errVertexTrack_a.push_back(hm);
      }
    }
  }else{
    // not fit
    m_vertex_a = m_vertex_b;
  }
  return m_errorFlag;
}

unsigned
kmassfitter::m_calDgf(void)
{
  m_dgf = 1;
  return m_errorFlag;
}

void
kmassfitter::fixMass(void)
{
  m_isFixMass.push_back(KF_FIX_MASS);
}

void
kmassfitter::unfixMass(void)
{
  m_isFixMass.push_back(KF_UNFIX_MASS);
}

double
kmassfitter::chisq(void){
  return m_chisq;
}

double
kmassfitter::chisq(const unsigned n){
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
    if(m_errorMsgFlag)dout(Debugout::ERR,"kmassfitter") << "(kmassfitter): Out of Range!" << std::endl;
    return -1.;
  }
}

//  void 
//  kmassfitter::correlation(const HepMatrix &e){
//    m_correlation_b.push_back(e);
//    m_correlationFlag = 1;
//    if(e.num_row() != (int)KF_NUM7)
//      m_errorFlag = KF_INPUT_MATRIX_SIZE;
//  }

//  void 
//  kmassfitter::correlation(void){
//    m_correlation_b.push_back(HepMatrix(KF_NUM7,
//  				      KF_NUM7,0));
//    m_correlationFlag = 1;
//  }

HepMatrix
kmassfitter::correlation(const unsigned n, const unsigned m, const unsigned flag){
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
      if(m_errorMsgFlag)dout(Debugout::ERR,"kmassfitter") << "(kmassfitter): Out of Range!" << std::endl;
      return HepMatrix(KF_NUM7,KF_NUM7,0);
    }
  }else{
    return kfitterbase::correlation(n,m,flag);
  }
}

void
kmassfitter::correlation(const HepMatrix &m)
{
  kfitterbase::correlation(m);
}

void
kmassfitter::correlation(void)
{ 
  kfitterbase::correlation();
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
