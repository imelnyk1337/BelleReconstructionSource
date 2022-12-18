#include "belle.h"
#include "kfitter/kfitterbase2.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

unsigned 
kfitterbase2::fit(void)
{
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
  HepMatrix tmp2_lam0(m_lam0),tmp2_v_a(m_v_a),tmp2_v(m_v_a);

  for(unsigned j=0;j<KF_MAX_ITERATION_NUMBER;++j){
    tmp_chiSq = KF_INIT_CHI2;
    for(unsigned i=0;i<KF_MAX_ITERATION_NUMBER;++i){
      if(m_setInputSubMatrix() != KF_NO_ERROR)return m_errorFlag;
      if(m_makeCoreMatrix() != KF_NO_ERROR)return m_errorFlag;
      //m_V_D = (m_D*m_V_al_0*(m_D.T())).inverse(errInverse);
      m_V_D = (m_V_al_0.similarity(m_D)).inverse(errInverse); // 2000/03/06
      if(errInverse != 0){
	m_errorFlag = KF_INVERSE;
	return m_errorFlag;
      }
      m_V_E = ((m_E.T())*m_V_D*m_E).inverse(errInverse);
      if(errInverse != 0){
	m_errorFlag = KF_INVERSE;
	return m_errorFlag;
      }
      m_lam0 = m_V_D*(m_D*(m_al_0-m_al_1)+m_d);
      chiSq  = ((m_lam0.T())*(m_D*(m_al_0-m_al_1)+m_E*(m_v-m_v_a)+m_d))(1,1); 
      m_v_a  = m_v_a - m_V_E*(m_E.T())*m_lam0;
      
      if(tmp_chiSq > chiSq){
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
	return m_errorFlag;
      }
    }
    m_al_a   = m_al_1;
    m_lam    = m_lam0 - m_V_D*m_E*m_V_E*(m_E.T())*m_lam0;
    m_al_1   = m_al_0 - m_V_al_0*(m_D.T())*m_lam;
    if(j == 0){
      tmp2_chiSq = chiSq;
      tmp2_v_a   = m_v_a;
      tmp2_v     = m_v;
      tmp2_V_E   = m_V_E;
      tmp2_V_D   = m_V_D;
      tmp2_lam0  = m_lam0;
      tmp2_E     = m_E;
      tmp2_D     = m_D;
      tmp_al_a   = m_al_a;
      continue;
    }else{
      if(tmp2_chiSq > chiSq){
	tmp2_chiSq = chiSq;
	tmp2_v_a   = m_v_a;
	tmp2_v     = m_v;
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
	chiSq   = tmp2_chiSq;
	m_v_a   = tmp2_v_a;
	m_v     = tmp2_v;
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
  if(m_errorFlag != KF_NO_ERROR)return m_errorFlag;
  m_lam    = m_lam0 - m_V_D*m_E*m_V_E*(m_E.T())*m_lam0;
  m_al_1   = m_al_0 - m_V_al_0*(m_D.T())*m_lam;
  m_V_Dt   = m_V_D  - m_V_D*m_E*m_V_E*(m_E.T())*m_V_D;
  m_V_al_1 = m_V_al_0 - m_V_al_0*(m_D.T())*m_V_Dt*m_D*m_V_al_0;
  m_Cov_v_al_1 = -m_V_E*(m_E.T())*m_V_D*m_D*m_V_al_0;
  if(m_setOutputMatrix() != KF_NO_ERROR)return m_errorFlag;

  m_chisq = chiSq;
  m_cl    = chisq2Confi(m_dgf, m_chisq);

  return m_errorFlag;
}

//double 
//kfitterbase2::estimation(void){
//rough estimation...
//This # is almost 0, I think..... 
//  return 2.*((m_lam.T()*(m_D*(m_al_1-m_al_a)+m_E*(m_v_a-m_v)+m_d))[0][0]);
//}

void 
kfitterbase2::dump(const unsigned flag)
{
  if(flag == KF_DUMP_MEASUREMENT){
    dout(Debugout::DUMP,"kfitterbase2") << "Track#=" << m_trackNum << ", ErrorFlag=" << m_errorFlag << std::endl;
    dout(Debugout::DUMP,"kfitterbase2") << "m_al_0=" << m_al_0 << std::endl;
    dout(Debugout::DUMP,"kfitterbase2") << "m_v=" << m_v << std::endl;
    dout(Debugout::DUMP,"kfitterbase2") << "m_V_al_0=" << m_V_al_0 << std::endl;
  }else if(flag == KF_DUMP_CORE_MATRIX){
    dout(Debugout::DUMP,"kfitterbase2") << "Track#=" << m_trackNum << ", ErrorFlag=" << m_errorFlag << std::endl;
    dout(Debugout::DUMP,"kfitterbase2") << "m_D=" << m_D << std::endl;
    dout(Debugout::DUMP,"kfitterbase2") << "m_E=" << m_E << std::endl;
    dout(Debugout::DUMP,"kfitterbase2") << "m_d=" << m_d << std::endl;
  }else if(flag == KF_DUMP_FITTED){
    dout(Debugout::DUMP,"kfitterbase2") << "Track#=" << m_trackNum << ", ErrorFlag=" << m_errorFlag << std::endl;
    dout(Debugout::DUMP,"kfitterbase2") << "Chi2=" << m_chisq << ", Dgf=" << m_dgf << std::endl;
    dout(Debugout::DUMP,"kfitterbase2") << "m_al_1=" << m_al_1 << std::endl;
    dout(Debugout::DUMP,"kfitterbase2") << "m_v_a=" << m_v_a << std::endl;
    dout(Debugout::DUMP,"kfitterbase2") << "m_V_al_1=" << m_V_al_1 << std::endl;
  }
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
