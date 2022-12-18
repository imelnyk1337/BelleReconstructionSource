#include "belle.h"
#include "particle/ParticleUserInfo.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


class UserInfo : public ParticleUserInfo
{
public:
  /// Default constructor
  UserInfo();

  /// Copy constructor
  UserInfo(const UserInfo &);

  /// Destructor
  virtual ~UserInfo();

  /// constructs self object.
  UserInfo * clone(void) const;

  /// Copy operator
  UserInfo & operator = (const UserInfo &);

public:
  void chisq(const double &v) { m_chisq = v; }
  void cl(const double &v)    { m_cl = v; }
  void test(const double &);

  const double & chisq(void) const { return m_chisq; }
  const double & cl(void)    const { return m_cl; }
  const double & test(void)  const { return *m_test; }

private:
  double m_chisq;
  double m_cl;

  double *m_test;
};
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
