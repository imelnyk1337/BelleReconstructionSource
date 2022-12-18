#if ! defined(PANTHER_H)
#include "panther/panther.h"
#endif
#if defined(__cplusplus)
extern "C" {
#endif   /* __cplusplus */
#ifndef PANTHER_FILL_IP_VERTEX_
#define PANTHER_FILL_IP_VERTEX_
struct fill_ip_vertex {
  int m_panther_dummy_;
  int m_ID;
  int m_flag;
  float m_v[3];
  float m_error[6];
  int m_ndf;
  float m_chisq;
  float m_cl;
};
extern PNT_Common *FILL_IP_VERTEX;
#endif  /* PANTHER_FILL_IP_VERTEX_ */

#if defined(__cplusplus)
}       /* extern "C" */
#if defined(BELLE_NAMESPACE)
  namespace Belle {
#endif
#ifndef PANTHER_FILL_IP_VERTEX_CLASS_
#define PANTHER_FILL_IP_VERTEX_CLASS_
class Fill_ip_vertex;
class Fill_ip_vertex_Manager;
class Fill_ip_vertex {

  public:
     Fill_ip_vertex():m_p(0){}
     Fill_ip_vertex(const Fill_ip_vertex&e):m_p(e.m_p){}
  private:
  friend class Table_Manager<Fill_ip_vertex>;
  friend class Fill_ip_vertex_Manager;
     Fill_ip_vertex(const Panther_ID&, PNT_TableIndex*);
     void newent(void);
  public:
    static Fill_ip_vertex_Manager & get_manager();
    ~Fill_ip_vertex(){}

  public: // Extractors
    int flag(void) const;
    float v(int) const;
    float vx(void) const;
    float vy(void) const;
    float vz(void) const;
    float error(int) const;
    int ndf(void) const;
    float chisq(void) const;
    float cl(void) const;
    Panther_ID get_ID(void) const { return Panther_ID(m_p?m_p->m_ID:0); }

  public: // Modifiers
    int flag(int);
    float v(int,float);
    float vx(float);
    float vy(float);
    float vz(float);
    float error(int,float);
    int ndf(int);
    float chisq(float);
    float cl(float);

  public:
    bool operator!(void) const { return m_p == 0; }
    operator  bool(void) const { return m_p != 0; }

    inline friend bool operator==(const Fill_ip_vertex&, const Fill_ip_vertex&);
    inline friend bool operator!=(const Fill_ip_vertex&, const Fill_ip_vertex&);
  private:
    struct fill_ip_vertex* m_p;
};
    inline bool operator==(const Fill_ip_vertex&a, const Fill_ip_vertex&b){ return a.m_p==b.m_p;} 
    inline bool operator!=(const Fill_ip_vertex&a, const Fill_ip_vertex&b){ return a.m_p!=b.m_p;}
//-----------------------------------------------------
class Fill_ip_vertex_Index;
class Fill_ip_vertex_Selector;
class Fill_ip_vertex_Manager : public Table_Manager<Fill_ip_vertex> {
  private:
    static Fill_ip_vertex_Manager *manager;
    static Fill_ip_vertex         *null;

  protected:
    Fill_ip_vertex_Manager();
   ~Fill_ip_vertex_Manager();

  public:
typedef std::vector<Fill_ip_vertex>::iterator iterator;
typedef std::vector<Fill_ip_vertex>::const_iterator const_iterator;

    static Fill_ip_vertex_Manager& get_manager(void);
    static  Fill_ip_vertex& get_NULL(void);
    void delete_NULL(void);
    Fill_ip_vertex_Index    index(const char*);
    Fill_ip_vertex_Selector& selector(void);

  private:
    Fill_ip_vertex_Selector *m_selector;
    static Fill_ip_vertex_Manager& get_manager_impl(void);
    Fill_ip_vertex& get_NULL(const int&){ return get_NULL(); }
    virtual void remove_wrapper(int id=0);

  friend class Fill_ip_vertex_Index;
  friend class Fill_ip_vertex_Selector;
  friend std::vector<Fill_ip_vertex> point_from(const Panther_ID&,Fill_ip_vertex_Index&);
};

class Fill_ip_vertex_Index : public Panther_Index<Fill_ip_vertex> {
  private:
    Fill_ip_vertex_Index();

    Fill_ip_vertex_Index(const char *name);

  public:
   ~Fill_ip_vertex_Index();

  public:
    void update(void);
    void update(Fill_ip_vertex_Index&);
    void update_reverse(Fill_ip_vertex_Index&);

    Fill_ip_vertex& operator()(const Panther_ID&) const;
    Fill_ip_vertex& operator[](const int&) const;

  friend class Fill_ip_vertex_Manager;
  friend std::vector<Fill_ip_vertex> point_from(const Panther_ID&,Fill_ip_vertex_Index&);
};

class Fill_ip_vertex_Selector {
  private:
    Fill_ip_vertex_Selector();
    PNT_Selector *selector_address;

  public:
   ~Fill_ip_vertex_Selector();

  public:
    void insert(Fill_ip_vertex&);
    int  count(void) const;
    void erase(int i=0);
    int  check(Panther_ID) const;
    void dump(void) const;

    Fill_ip_vertex& operator()(const Panther_ID&) const;

  friend class Fill_ip_vertex_Manager;
};
#endif /* PANTHER_FILL_IP_VERTEX_CLASS_ */
//-----------------------------------------------------
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif  /* __cplusplus */
#if defined(__cplusplus)
#if defined(BELLE_NAMESPACE)
  namespace Belle {
#endif
#ifndef PANTHER_FILL_IP_VERTEX_CLASS_DEFINITION_
#define PANTHER_FILL_IP_VERTEX_CLASS_DEFINITION_
// Extractors
inline int Fill_ip_vertex::flag(void) const { return m_p->m_flag; }
inline float Fill_ip_vertex::v(int i) const { return m_p->m_v[i]; }
inline float Fill_ip_vertex::vx(void) const { return m_p->m_v[1-1]; }
inline float Fill_ip_vertex::vy(void) const { return m_p->m_v[2-1]; }
inline float Fill_ip_vertex::vz(void) const { return m_p->m_v[3-1]; }
inline float Fill_ip_vertex::error(int i) const { return m_p->m_error[i]; }
inline int Fill_ip_vertex::ndf(void) const { return m_p->m_ndf; }
inline float Fill_ip_vertex::chisq(void) const { return m_p->m_chisq; }
inline float Fill_ip_vertex::cl(void) const { return m_p->m_cl; }

// Modifiers
inline int Fill_ip_vertex::flag(int i){ return m_p->m_flag = i; }
inline float Fill_ip_vertex::v(int i, float j){ return m_p->m_v[i] = j; }
inline float Fill_ip_vertex::vx(float i){ return m_p->m_v[1-1] = i; }
inline float Fill_ip_vertex::vy(float i){ return m_p->m_v[2-1] = i; }
inline float Fill_ip_vertex::vz(float i){ return m_p->m_v[3-1] = i; }
inline float Fill_ip_vertex::error(int i, float j){ return m_p->m_error[i] = j; }
inline int Fill_ip_vertex::ndf(int i){ return m_p->m_ndf = i; }
inline float Fill_ip_vertex::chisq(float i){ return m_p->m_chisq = i; }
inline float Fill_ip_vertex::cl(float i){ return m_p->m_cl = i; }
//-----------------------------------------------------
inline Fill_ip_vertex_Manager& Fill_ip_vertex::get_manager(void){
  return Fill_ip_vertex_Manager::get_manager();
}

inline Fill_ip_vertex_Manager& Fill_ip_vertex_Manager::get_manager(void){
  if(manager!=0) return *manager;
  else return get_manager_impl();
}

#endif /* PANTHER_FILL_IP_VERTEX_CLASS_DEFINITION_ */
//-----------------------------------------------------
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* __cplusplus */

