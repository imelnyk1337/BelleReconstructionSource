#include "belle.h"
#include <cmath>
#include "particle/gammac.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


GammaParticle::GammaParticle(const double energy,
                             const double phi,
                             const double theta,
                             const double R,
                             const HepSymMatrix & Eg)
 : m_energy(energy), m_theta(theta), m_phi(phi),
   m_errorGamma(Eg), m_vertex(), m_errorVertex(3,0), m_errorGammaVertex(6,0),
   m_momentum(), m_errorMomentum(3,0),
   m_momentumEnergy(), m_errorMomentumEnergy(4,0), m_R(R),
   m_ipFlag(1)
{
  m_errorGammaVertex.sub(1,m_errorGamma);
  calcu();
}

GammaParticle::~GammaParticle(){
}

const HepPoint3D & 
GammaParticle::vertex(const HepPoint3D &vertex)
{
  m_ipFlag = 0;
  m_vertex = vertex;
  calcu();
  return m_vertex;
}

const HepPoint3D & 
GammaParticle::vertex(const HepPoint3D &vertex, const HepSymMatrix &error)
{
  m_ipFlag = 0;
  m_vertex = vertex;
  if(error.num_row() == 3){
    m_errorVertex = error;
    m_errorGammaVertex.sub(4,m_errorVertex);
  }else{
    dout(Debugout::WARN,"gammac") << "Error in vertex(const HepPoint3D&," << std::endl;
    dout(Debugout::WARN,"gammac") << "                const HepSymMatrix&) of Gamma Class" << std::endl;
    dout(Debugout::WARN,"gammac") << "                      ~~~~~~~~~~~~~--- 3 x 3 error matrix!!! " <<std::endl;
    //sets 0.0 in all elements.  
    for(unsigned i=0;i<3;i++){
      for(unsigned j=i;j<3;j++){
        m_errorVertex[i][j] = 0.0;
      }
    }
  }
  calcu();
  return m_vertex;
}

void
GammaParticle::calcu(void)
{
  if(m_ipFlag == 1){
    double sint = sin(m_theta);
    double cost = cos(m_theta);
    double sinp = sin(m_phi);
    double cosp = cos(m_phi);
    m_momentum.setX(m_energy*sint*cosp);
    m_momentum.setY(m_energy*sint*sinp);
    m_momentum.setZ(m_energy*cost);
    m_momentumEnergy = HepLorentzVector(m_momentum,m_energy);
    m_errorMomentumEnergy = m_errorGamma.similarity(delPEdelG());
    m_errorMomentum = m_errorMomentumEnergy.sub(1,3);
  }else{
    calcuVector();
    calcuError();
  }
}

HepMatrix 
GammaParticle::delPEdelG(void)
{
  HepMatrix dPEdG(4,3,0);
  
  double sint = sin(m_theta);
  double cost = cos(m_theta);
  double sinp = sin(m_phi);
  double cosp = cos(m_phi);
  
  //@* = @px, @py, @pz, @E

  //@*/@E
  dPEdG[0][0] = sint*cosp;
  dPEdG[1][0] = sint*sinp;
  dPEdG[2][0] = cost;
  dPEdG[3][0] = 1.;

  //@*/@phi
  dPEdG[0][1] = -m_energy*sint*sinp;
  dPEdG[1][1] =  m_energy*sint*cosp;

  //@*/@theta
  dPEdG[0][2] =  m_energy*cost*cosp;
  dPEdG[1][2] =  m_energy*cost*sinp;
  dPEdG[2][2] = -m_energy*sint;

  return dPEdG;
}

unsigned
GammaParticle::calcuVector(void){
  double sint = sin(m_theta);
  double cost = cos(m_theta);
  double sinp = sin(m_phi);
  double cosp = cos(m_phi);

  double vx = m_vertex.x();
  double vy = m_vertex.y();
  double vz = m_vertex.z();

  double lxy = sqrt((m_R*sint*cosp-vx)*(m_R*sint*cosp-vx)+
		    (m_R*sint*sinp-vy)*(m_R*sint*sinp-vy));
  double l = sqrt(lxy*lxy+(m_R*cost-vz)*(m_R*cost-vz));

  if(l == 0. || lxy == 0.)return 0;
  
  double cosT = (m_R*cost-vz)/l;
  double sinT = sqrt(1.-cosT*cosT);
  double cosP = (m_R*sint*cosp-vx)/lxy;
  double sinP = (m_R*sint*sinp-vy)/lxy;

  m_momentum.setX(m_energy*sinT*cosP);
  m_momentum.setY(m_energy*sinT*sinP);
  m_momentum.setZ(m_energy*cosT);
  m_momentumEnergy = HepLorentzVector(m_momentum,m_energy);
  
  return 1;
}

unsigned
GammaParticle::calcuError(void)
{
  if(m_ipFlag == 0){
    m_errorMomentumEnergy = m_errorGammaVertex.similarity(delPEdelGV());
    m_errorMomentum = m_errorMomentumEnergy.sub(1,3);
    return 1;
  }else{
    return 0;
  }
}

HepMatrix 
GammaParticle::delPEdelGV(void)
{
  double sint = sin(m_theta);
  double cost = cos(m_theta);
  double sinp = sin(m_phi);
  double cosp = cos(m_phi);
  
  double vx = m_vertex.x();
  double vy = m_vertex.y();
  double vz = m_vertex.z();

  double lxy = sqrt((m_R*sint*cosp-vx)*(m_R*sint*cosp-vx)+
		    (m_R*sint*sinp-vy)*(m_R*sint*sinp-vy));
  double l = sqrt(lxy*lxy+(m_R*cost-vz)*(m_R*cost-vz));

  HepMatrix dPEdGV(4,6,0);
  
  if(l == 0. || lxy == 0. || sint == 0.)return dPEdGV;

  double inLxy = 1./lxy;
  double inL   = 1./l;

  double cosT = (m_R*cost-vz)*inL;
  double sinT = sqrt(1.-cosT*cosT);
  double cosP = (m_R*sint*cosp-vx)*inLxy;
  double sinP = (m_R*sint*sinp-vy)*inLxy;

  if(sinT == 0.)return dPEdGV;
  double inTanT = cosT/sinT;

  double DlxyDphi   = inLxy*(-m_R*sint*sinp*(m_R*sint*cosp-vx)+
			     m_R*sint*cosp*(m_R*sint*sinp-vy));
  double DlxyDtheta = inLxy*(m_R*cost*cosp*(m_R*sint*cosp-vx)+
			     m_R*cost*sinp*(m_R*sint*sinp-vy));
  double DlxyDvx    = -inLxy*(m_R*sint*cosp-vx);
  double DlxyDvy    = -inLxy*(m_R*sint*sinp-vy);

  double DlDphi   = lxy*inL*DlxyDphi;
  double DlDtheta = lxy*inL*DlxyDtheta-(m_R*cost-vz)*inL*m_R*sint;
  double DlDvx    = lxy*inL*DlxyDvx;
  double DlDvy    = lxy*inL*DlxyDvy;
  double DlDvz    = -(m_R*cost-vz)*inL;

  double DcosTDphi   = -cosT*inL*DlDphi;
  double DcosTDtheta = -cosT*inL*DlDtheta-m_R*sint*inL;
  double DcosTDvx    = -cosT*inL*DlDvx;
  double DcosTDvy    = -cosT*inL*DlDvy;
  double DcosTDvz    = -cosT*inL*DlDvz-inL;

  double DsinTDphi   = -inTanT*DcosTDphi;
  double DsinTDtheta = -inTanT*DcosTDtheta;
  double DsinTDvx    = -inTanT*DcosTDvx;
  double DsinTDvy    = -inTanT*DcosTDvy;
  double DsinTDvz    = -inTanT*DcosTDvz;

  double DsinPDphi   =  inLxy*(m_R*sint*cosp-sinT*DlxyDphi);
  double DsinPDtheta =  inLxy*(m_R*cost*sinp-sinP*DlxyDtheta);
  double DsinPDvx    = -inLxy*sinT*DlxyDvx;
  double DsinPDvy    = -inLxy*(1.+sinT*DlxyDvy);

  double DcosPDphi   = -inLxy*(m_R*sint*sinp+cosT*DlxyDphi);
  double DcosPDtheta =  inLxy*(m_R*cost*cosp-cosP*DlxyDtheta);
  double DcosPDvx    = -inLxy*(1.+cosT*DlxyDvx);
  double DcosPDvy    = -inLxy*cosT*DlxyDvy;

  //@* = @px, @py, @pz, @E (new parameter)

  //@*/@E
  dPEdGV[0][0] = sinT*cosP;
  dPEdGV[1][0] = sinT*sinP;
  dPEdGV[2][0] = cosT;
  dPEdGV[3][0] = 1.;

  //@*/@phi
  dPEdGV[0][1] = m_energy*(DsinTDphi*cosP+sinT*DcosPDphi);
  dPEdGV[1][1] = m_energy*(DsinTDphi*sinP+sinT*DsinPDphi);
  dPEdGV[2][1] = m_energy*DcosTDphi;

  //@*/@theta
  dPEdGV[0][2] = m_energy*(DsinTDtheta*cosP+sinT*DcosPDtheta);
  dPEdGV[1][2] = m_energy*(DsinTDtheta*sinP+sinT*DsinPDtheta);
  dPEdGV[2][2] = m_energy*DcosTDtheta;

  //@*/@vx
  dPEdGV[0][3] = m_energy*(DsinTDvx*cosP+sinT*DcosPDvx);
  dPEdGV[1][3] = m_energy*(DsinTDvx*sinP+sinT*DsinPDvx);
  dPEdGV[2][3] = m_energy*DcosTDvx;

  //@*/@vy
  dPEdGV[0][4] = m_energy*(DsinTDvy*cosP+sinT*DcosPDvy);
  dPEdGV[1][4] = m_energy*(DsinTDvy*sinP+sinT*DsinPDvy);
  dPEdGV[2][4] = m_energy*DcosTDvy;

  //@*/@vz
  dPEdGV[0][5] = m_energy*DsinTDvz*cosP;
  dPEdGV[1][5] = m_energy*DsinTDvz*sinP;
  dPEdGV[2][5] = m_energy*DcosTDvz;

  return dPEdGV;
}

GammaParticle & 
GammaParticle::operator = (const GammaParticle &g){
  if(this == &g)return *this;

  m_energy              = g.m_energy;
  m_theta               = g.m_theta;
  m_phi                 = g.m_phi;
  m_errorGamma          = g.m_errorGamma;
  m_vertex              = g.m_vertex;
  m_errorVertex         = g.m_errorVertex;
  m_errorGammaVertex    = g.m_errorGammaVertex;
  m_momentum            = g.m_momentum;
  m_errorMomentum       = g.m_errorMomentum;
  m_momentumEnergy      = g.m_momentumEnergy;
  m_errorMomentumEnergy = g.m_errorMomentumEnergy;
  m_R                   = g.m_R;
  m_ipFlag              = g.m_ipFlag;
  
  return *this;
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
