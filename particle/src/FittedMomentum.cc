//
// Fitted Momentum Class of Particle Class
//
// --- FittedMomentum.cc ---
//
// $Id: FittedMomentum.cc 9932 2006-11-12 14:26:53Z katayama $
//
// $Log$
// Revision 1.6  2002/02/27 01:38:20  katayama
// int/unsigned comparison.
//
// Revision 1.5  2001/12/12 07:10:53  jtanaka
// (1) compatibility for obsoleted Mdst_vee/mdst_sim_xref table
// (2) compatibility for gcc3
//
// Revision 1.4  1999/04/15 20:03:47  jtanaka
// add statical info.
//
// Revision 1.3  1999/04/09 15:02:09  jtanaka
// Added some member functions,"const", and "&" and removed some members
//
// Revision 1.2  1998/10/22 11:38:00  jtanaka
// changed names of the member.  e.g) del_info --> delInfo
//
// Revision 1.1  1998/10/20 14:45:17  jtanaka
// added FittedMomentum.
//
// 
#include "belle.h"
#include "particle/FittedMomentum.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


static HepMatrix *static_zeroMatrix(NULL);

FittedMomentum::FittedMomentum()
  : m_vertexAndSelfError(3,7,0),
    m_decayVertexAndSelfError(3,7,0),
    m_cl(0.0),
    m_chisq(0.0),
    m_dof(0)
{
}

FittedMomentum::FittedMomentum(const FittedMomentum &p)
  : Momentum(p),
    m_correlationParticles(p.m_correlationParticles),
    m_correlationMatrix(p.m_correlationMatrix),
    m_vertexAndSelfError(p.m_vertexAndSelfError),
    m_decayVertexAndSelfError(p.m_decayVertexAndSelfError),
    m_cl(p.m_cl),
    m_chisq(p.m_chisq),
    m_dof(p.m_dof)
{
}

FittedMomentum::~FittedMomentum()
{
  m_correlationMatrix.erase(m_correlationMatrix.begin(),m_correlationMatrix.end());
  m_correlationParticles.erase(m_correlationParticles.begin(),m_correlationParticles.end());
}

FittedMomentum & 
FittedMomentum::operator = (const FittedMomentum &p)
{
  if(this == &p)return *this;
  Momentum::operator=(p);
  m_vertexAndSelfError       = p.m_vertexAndSelfError;
  m_decayVertexAndSelfError = p.m_decayVertexAndSelfError;
  m_correlationMatrix    = p.m_correlationMatrix;
  m_correlationParticles = p.m_correlationParticles;
  m_cl = p.m_cl;
  m_chisq = p.m_chisq;
  m_dof = p.m_dof;
  return *this;
}

FittedMomentum &
FittedMomentum::operator = (const Momentum &p)
{
  if(this == &p)return *this;
  Momentum::operator=(p);
  return *this;
}

// Selectors
const HepMatrix &
FittedMomentum::coMatrix(Particle &p, unsigned *error) const
{
  Particle *t = &p;
  return coMatrix(t, error);
}

const HepMatrix &
FittedMomentum::coMatrix(Particle *p, unsigned *error) const
{
  int size = m_correlationParticles.size();
  for(int i=0;i<size;i++){
    if(m_correlationParticles[i] == p){
	if(error != NULL)*error = 1;
	return m_correlationMatrix[i];
    }
  }
  if(error != NULL)*error = 0;
  if(static_zeroMatrix) return *static_zeroMatrix;
  else return *(static_zeroMatrix=new HepMatrix(7,7,0));
}

const HepMatrix &
FittedMomentum::coMatrix(unsigned index, unsigned *error) const
{
  int size = m_correlationParticles.size();
  //  if(index >= 0 && size > (int)index){
  if(size > (int)index){
    if(error != NULL)*error = 1;
    return m_correlationMatrix[index];
  }
  if(error != NULL)*error = 0;
  if(static_zeroMatrix) return *static_zeroMatrix;
  else return *(static_zeroMatrix=new HepMatrix(7,7,0));
}

const Particle *
FittedMomentum::coParticle(unsigned index, unsigned *error) const
{
  int size = m_correlationParticles.size();
  //  if(index >= 0 && size > (int)index){
  if(size > (int)index){
    if(error != NULL)*error = 1;
    return m_correlationParticles[index];
  }
  if(error != NULL)*error = 0;
  return NULL;
}

const HepMatrix &
FittedMomentum::coVertex(void) const
{
  return m_vertexAndSelfError;
}

const HepMatrix &
FittedMomentum::coDecayVertex(void) const
{
  return m_decayVertexAndSelfError;
}

double 
FittedMomentum::cl(void) const
{
  return m_cl;
}

double 
FittedMomentum::chisq(void) const
{
  return m_chisq;
}

unsigned 
FittedMomentum::dof(void) const
{
  return m_dof;
}

// Modifiers
HepMatrix & 
FittedMomentum::coMatrix(Particle &p, HepMatrix &h)
{
  Particle *t = &p;
  return coMatrix(t,h);
}

HepMatrix & 
FittedMomentum::coMatrix(Particle *p, HepMatrix &h)
{
  //...checks
  int size = m_correlationParticles.size();
  for(int i=0;i<size;i++){
    if(m_correlationParticles[i] == p){
	m_correlationMatrix[i] = h;
	return m_correlationMatrix[i];
    }
  }
  //...makes
  m_correlationParticles.push_back(p);
  m_correlationMatrix.push_back(h);
  return m_correlationMatrix[size+1];
}

HepMatrix & 
FittedMomentum::coVertex(HepMatrix &h)
{
  m_vertexAndSelfError = h;
  return m_vertexAndSelfError;
}

HepMatrix & 
FittedMomentum::coDecayVertex(HepMatrix &h)
{
  m_decayVertexAndSelfError = h;
  return m_decayVertexAndSelfError;
}

unsigned 
FittedMomentum::nCoMatrix(void) const
{
  return m_correlationMatrix.size();
}

void 
FittedMomentum::remove(unsigned type)
{
  if(type & REMOVE_COMATRIX){
    m_correlationMatrix.erase(m_correlationMatrix.begin(),m_correlationMatrix.end());
    m_correlationParticles.erase(m_correlationParticles.begin(),m_correlationParticles.end());
  }
  if(type & REMOVE_COVERTEX){
    m_vertexAndSelfError = HepMatrix(3,7,0);
  }
  if(type & REMOVE_CODECAYVERTEX){
    m_decayVertexAndSelfError = HepMatrix(3,7,0);
  }
}

void 
FittedMomentum::removeCoMatrix(const Particle &p)
{
  const Particle *t = &p;
  removeCoMatrix(t);
}

void 
FittedMomentum::removeCoMatrix(const Particle *p)
{
  int size = m_correlationParticles.size();
  for(int i=0;i<size;i++){
    if(m_correlationParticles[i] == p){
      m_correlationParticles.erase(m_correlationParticles.begin()+i);
      m_correlationMatrix.erase(m_correlationMatrix.begin()+i);
      break;
    }
  }
}

double 
FittedMomentum::cl(const double &cl)
{
  return m_cl = cl;
}

double 
FittedMomentum::chisq(const double &chisq)
{
  return m_chisq = chisq;
}

unsigned 
FittedMomentum::dof(const unsigned &dof)
{
  return m_dof = dof;
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
