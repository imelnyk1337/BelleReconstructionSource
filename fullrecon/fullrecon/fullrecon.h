#if ! defined(PANTHER_H)
#include "panther/panther.h"
#endif
#if defined(__cplusplus)
extern "C" {
#endif   /* __cplusplus */
#ifndef PANTHER_FULLRECON_
#define PANTHER_FULLRECON_
struct fullrecon {
  int m_panther_dummy_;
  int m_ID;
  int m_FullBrec;
  int m_kind;
  int m_b_mode;
  int m_sub0_mode;
  int m_sub1_mode;
  float m_deltE;
  float m_Mbc;
  float m_Aux0;
  float m_Aux1;
  float m_Aux2;
  float m_Aux3;
  float m_Aux4;
  float m_Aux5;
  float m_Aux6;
  float m_Aux7;
  float m_Aux8;
  float m_Aux9;
  int m_flag0;
  int m_flag1;
  int m_flag2;
  int m_flag3;
  int m_flag4;
};
extern PNT_Common *FULLRECON;
#endif  /* PANTHER_FULLRECON_ */

#if defined(__cplusplus)
}       /* extern "C" */
#if defined(BELLE_NAMESPACE)
  namespace Belle {
#endif
#ifndef PANTHER_FULLRECON_CLASS_
#define PANTHER_FULLRECON_CLASS_
class Brecon;
class Fullrecon;
class Fullrecon_Manager;
class Fullrecon {

  public:
     Fullrecon():m_p(0){}
     Fullrecon(const Fullrecon&e):m_p(e.m_p){}
  private:
  friend class Table_Manager<Fullrecon>;
  friend class Fullrecon_Manager;
     Fullrecon(const Panther_ID&, PNT_TableIndex*);
     void newent(void);
  public:
    static Fullrecon_Manager & get_manager();
    ~Fullrecon(){}

  public: // Extractors
    Brecon& FullBrec(void) const;
    int FullBrec_ID(void) const { return m_p->m_FullBrec; }
    int kind(void) const;
    int b_mode(void) const;
    int sub0_mode(void) const;
    int sub1_mode(void) const;
    float deltE(void) const;
    float Mbc(void) const;
    float Aux0(void) const;
    float Aux1(void) const;
    float Aux2(void) const;
    float Aux3(void) const;
    float Aux4(void) const;
    float Aux5(void) const;
    float Aux6(void) const;
    float Aux7(void) const;
    float Aux8(void) const;
    float Aux9(void) const;
    int flag0(void) const;
    int flag1(void) const;
    int flag2(void) const;
    int flag3(void) const;
    int flag4(void) const;
    Panther_ID get_ID(void) const { return Panther_ID(m_p?m_p->m_ID:0); }

  public: // Modifiers
    void FullBrec(const Brecon&);
    void reset_FullBrec(void);
    int kind(int);
    int b_mode(int);
    int sub0_mode(int);
    int sub1_mode(int);
    float deltE(float);
    float Mbc(float);
    float Aux0(float);
    float Aux1(float);
    float Aux2(float);
    float Aux3(float);
    float Aux4(float);
    float Aux5(float);
    float Aux6(float);
    float Aux7(float);
    float Aux8(float);
    float Aux9(float);
    int flag0(int);
    int flag1(int);
    int flag2(int);
    int flag3(int);
    int flag4(int);

  public:
    bool operator!(void) const { return m_p == 0; }
    operator  bool(void) const { return m_p != 0; }

    inline friend bool operator==(const Fullrecon&, const Fullrecon&);
    inline friend bool operator!=(const Fullrecon&, const Fullrecon&);
  private:
    struct fullrecon* m_p;
};
    inline bool operator==(const Fullrecon&a, const Fullrecon&b){ return a.m_p==b.m_p;} 
    inline bool operator!=(const Fullrecon&a, const Fullrecon&b){ return a.m_p!=b.m_p;}
//-----------------------------------------------------
class Fullrecon_Index;
class Fullrecon_Selector;
class Fullrecon_Manager : public Table_Manager<Fullrecon> {
  private:
    static Fullrecon_Manager *manager;
    static Fullrecon         *null;

  protected:
    Fullrecon_Manager();
   ~Fullrecon_Manager();

  public:
typedef std::vector<Fullrecon>::iterator iterator;
typedef std::vector<Fullrecon>::const_iterator const_iterator;

    static Fullrecon_Manager& get_manager(void);
    static  Fullrecon& get_NULL(void);
    void delete_NULL(void);
    Fullrecon_Index    index(const char*);
    Fullrecon_Selector& selector(void);

  private:
    Fullrecon_Selector *m_selector;
    static Fullrecon_Manager& get_manager_impl(void);
    Fullrecon& get_NULL(const int&){ return get_NULL(); }
    virtual void remove_wrapper(int id=0);

  friend class Fullrecon_Index;
  friend class Fullrecon_Selector;
  friend std::vector<Fullrecon> point_from(const Panther_ID&,Fullrecon_Index&);
};

class Fullrecon_Index : public Panther_Index<Fullrecon> {
  private:
    Fullrecon_Index();

    Fullrecon_Index(const char *name);

  public:
   ~Fullrecon_Index();

  public:
    void update(void);
    void update(Fullrecon_Index&);
    void update_reverse(Fullrecon_Index&);

    Fullrecon& operator()(const Panther_ID&) const;
    Fullrecon& operator[](const int&) const;

  friend class Fullrecon_Manager;
  friend std::vector<Fullrecon> point_from(const Panther_ID&,Fullrecon_Index&);
};

class Fullrecon_Selector {
  private:
    Fullrecon_Selector();
    PNT_Selector *selector_address;

  public:
   ~Fullrecon_Selector();

  public:
    void insert(Fullrecon&);
    int  count(void) const;
    void erase(int i=0);
    int  check(Panther_ID) const;
    void dump(void) const;

    Fullrecon& operator()(const Panther_ID&) const;

  friend class Fullrecon_Manager;
};
#endif /* PANTHER_FULLRECON_CLASS_ */
//-----------------------------------------------------
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif  /* __cplusplus */
#if defined(__cplusplus)
#if defined(BELLE_NAMESPACE)
  namespace Belle {
#endif
#ifndef PANTHER_FULLRECON_CLASS_DEFINITION_
#define PANTHER_FULLRECON_CLASS_DEFINITION_
// Extractors
inline int Fullrecon::kind(void) const { return m_p->m_kind; }
inline int Fullrecon::b_mode(void) const { return m_p->m_b_mode; }
inline int Fullrecon::sub0_mode(void) const { return m_p->m_sub0_mode; }
inline int Fullrecon::sub1_mode(void) const { return m_p->m_sub1_mode; }
inline float Fullrecon::deltE(void) const { return m_p->m_deltE; }
inline float Fullrecon::Mbc(void) const { return m_p->m_Mbc; }
inline float Fullrecon::Aux0(void) const { return m_p->m_Aux0; }
inline float Fullrecon::Aux1(void) const { return m_p->m_Aux1; }
inline float Fullrecon::Aux2(void) const { return m_p->m_Aux2; }
inline float Fullrecon::Aux3(void) const { return m_p->m_Aux3; }
inline float Fullrecon::Aux4(void) const { return m_p->m_Aux4; }
inline float Fullrecon::Aux5(void) const { return m_p->m_Aux5; }
inline float Fullrecon::Aux6(void) const { return m_p->m_Aux6; }
inline float Fullrecon::Aux7(void) const { return m_p->m_Aux7; }
inline float Fullrecon::Aux8(void) const { return m_p->m_Aux8; }
inline float Fullrecon::Aux9(void) const { return m_p->m_Aux9; }
inline int Fullrecon::flag0(void) const { return m_p->m_flag0; }
inline int Fullrecon::flag1(void) const { return m_p->m_flag1; }
inline int Fullrecon::flag2(void) const { return m_p->m_flag2; }
inline int Fullrecon::flag3(void) const { return m_p->m_flag3; }
inline int Fullrecon::flag4(void) const { return m_p->m_flag4; }

// Modifiers
inline int Fullrecon::kind(int i){ return m_p->m_kind = i; }
inline int Fullrecon::b_mode(int i){ return m_p->m_b_mode = i; }
inline int Fullrecon::sub0_mode(int i){ return m_p->m_sub0_mode = i; }
inline int Fullrecon::sub1_mode(int i){ return m_p->m_sub1_mode = i; }
inline float Fullrecon::deltE(float i){ return m_p->m_deltE = i; }
inline float Fullrecon::Mbc(float i){ return m_p->m_Mbc = i; }
inline float Fullrecon::Aux0(float i){ return m_p->m_Aux0 = i; }
inline float Fullrecon::Aux1(float i){ return m_p->m_Aux1 = i; }
inline float Fullrecon::Aux2(float i){ return m_p->m_Aux2 = i; }
inline float Fullrecon::Aux3(float i){ return m_p->m_Aux3 = i; }
inline float Fullrecon::Aux4(float i){ return m_p->m_Aux4 = i; }
inline float Fullrecon::Aux5(float i){ return m_p->m_Aux5 = i; }
inline float Fullrecon::Aux6(float i){ return m_p->m_Aux6 = i; }
inline float Fullrecon::Aux7(float i){ return m_p->m_Aux7 = i; }
inline float Fullrecon::Aux8(float i){ return m_p->m_Aux8 = i; }
inline float Fullrecon::Aux9(float i){ return m_p->m_Aux9 = i; }
inline int Fullrecon::flag0(int i){ return m_p->m_flag0 = i; }
inline int Fullrecon::flag1(int i){ return m_p->m_flag1 = i; }
inline int Fullrecon::flag2(int i){ return m_p->m_flag2 = i; }
inline int Fullrecon::flag3(int i){ return m_p->m_flag3 = i; }
inline int Fullrecon::flag4(int i){ return m_p->m_flag4 = i; }
//-----------------------------------------------------
inline Fullrecon_Manager& Fullrecon::get_manager(void){
  return Fullrecon_Manager::get_manager();
}

inline Fullrecon_Manager& Fullrecon_Manager::get_manager(void){
  if(manager!=0) return *manager;
  else return get_manager_impl();
}

#endif /* PANTHER_FULLRECON_CLASS_DEFINITION_ */
//-----------------------------------------------------
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* __cplusplus */

