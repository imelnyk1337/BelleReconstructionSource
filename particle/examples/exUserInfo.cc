#include "belle.h"
#include "exUserInfo.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//#include <iostream>
//#include <typeinfo>

UserInfo::UserInfo()
  : m_chisq(-7.), m_cl(-7.), m_test(0)
{
  //dout(Debugout::INFO,"exUserInfo") << "CONS1" << std::endl;
}

UserInfo::~UserInfo()
{
  //dout(Debugout::INFO,"exUserInfo") << "DEST" << std::endl;
  if(m_test)delete m_test;
}

UserInfo::UserInfo(const UserInfo &x)
  : m_chisq(x.m_chisq), m_cl(x.m_cl), m_test(0)
{
  //dout(Debugout::INFO,"exUserInfo") << "CONS2" << std::endl;

  if(x.m_test)m_test = new double(*(x.m_test));
}

UserInfo*
UserInfo::clone(void) const
{
  //dout(Debugout::INFO,"exUserInfo") << "CONS3" << std::endl;
  UserInfo *x = new UserInfo(*this);
  return x;
}

UserInfo &
UserInfo::operator = (const UserInfo &x)
{
  //dout(Debugout::INFO,"exUserInfo") << "CONS4" << std::endl;
  m_chisq = x.m_chisq;
  m_cl    = x.m_cl;
  
  if(m_test)delete m_test;
  if(x.m_test){
    m_test = new double(*(x.m_test));
  }else{
    m_test = 0;
  }
}

void 
UserInfo::test(const double &v)
{
  if(m_test){
    *m_test = v;
  }else{
    m_test = new double(v);
  }
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
