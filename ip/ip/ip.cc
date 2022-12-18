#include <iostream>
#include "tables/ip.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
  namespace Belle {
#endif

// #ifdef __DECCXX
// #pragma define_template Table_Manager<Fill_ip_vertex>
// #else
// template class Table_Manager<Fill_ip_vertex>;
// #endif /* __DECCXX */

//---------- Table FILL_IP_VERTEX ----------
Fill_ip_vertex_Manager* Fill_ip_vertex_Manager::manager = 0;
Fill_ip_vertex*         Fill_ip_vertex_Manager::null    = 0;
template<> PNT_Table *Table_Manager<Fill_ip_vertex>::m_table_address = 0;
void Fill_ip_vertex::newent(void){
  m_p = (fill_ip_vertex *)::BsNewEnt(FILL_IP_VERTEX);
}
Fill_ip_vertex::Fill_ip_vertex(const Panther_ID& ID, PNT_TableIndex* index){
  m_p = (fill_ip_vertex *)::BsGetEnt(FILL_IP_VERTEX,(int)ID, index);
}
  extern "C" void *(*Panther_Get_Manager_Address)(void*);
  extern "C" void  (*Panther_Remove_Wrapper)(void*,int);
Fill_ip_vertex_Manager::Fill_ip_vertex_Manager() : Table_Manager<Fill_ip_vertex>(), m_selector(NULL) { }
Fill_ip_vertex_Manager::~Fill_ip_vertex_Manager() {
  delete m_selector;
}

Fill_ip_vertex_Manager& Fill_ip_vertex_Manager::get_manager_impl(void){
  if( manager == 0 ){
    manager = new Fill_ip_vertex_Manager;
    null    = new Fill_ip_vertex;
    m_table_address = pnt_get_table(FILL_IP_VERTEX);
    void* addr = FILL_IP_VERTEX;
    pnt_register_manager((Table_Manager_Base*)manager,addr);
    Panther_Get_Manager_Address = ::pnt_get_manager_address;
    Panther_Remove_Wrapper      = ::remove_wrapper;
  }
  return *manager;
}

Fill_ip_vertex& Fill_ip_vertex_Manager::get_NULL(void){
  if( manager == 0 )
    Fill_ip_vertex_Manager::get_manager();
  return *null;
}

void Fill_ip_vertex_Manager::delete_NULL(void){
  if( null != 0 ) {
    delete null;
    null = NULL;
  }
}

void Fill_ip_vertex_Manager::remove_wrapper(int id){
  if( id == 0 ){
    erase( begin(), end() );
  }else{
    iterator k;
    for( k=begin(); k != end(); k++ ){
      if( (*k).get_ID() == id ) break;
    }
    if( k == end() ){
      dout(Debugout::ERR,"Fill_ip_vertex") << "Remove_wrapper warning: unknown object (ID=" << id << ")." << std::endl;
      return;
    }
    erase(k);
  }
}

Fill_ip_vertex_Index::Fill_ip_vertex_Index() { }

Fill_ip_vertex_Index::Fill_ip_vertex_Index(const char *name) : Panther_Index<Fill_ip_vertex>(name){ }

Fill_ip_vertex_Index::~Fill_ip_vertex_Index(){  
  ::BsDelInd(index_address); }

Fill_ip_vertex_Selector::Fill_ip_vertex_Selector() { }
Fill_ip_vertex_Selector::~Fill_ip_vertex_Selector() {  
::BsDelSel(selector_address); }

Fill_ip_vertex_Index Fill_ip_vertex_Manager::index(const char* attr) {
  return Fill_ip_vertex_Index(attr);
}

Fill_ip_vertex_Selector& Fill_ip_vertex_Manager::selector(void) {
  if(m_selector) return *m_selector;
  m_selector = new Fill_ip_vertex_Selector;
  m_selector->selector_address = ::pnt_maksel(m_table_address);
  return *m_selector;
}
//-------------------------------------------
void Fill_ip_vertex_Index::update(void){
  ::pnt_updind(Fill_ip_vertex_Manager::m_table_address, index_address);
}

void Fill_ip_vertex_Index::update(Fill_ip_vertex_Index& Index){
  ::pnt_upind2(Fill_ip_vertex_Manager::m_table_address,
        Index.index_address,index_address);
}

void Fill_ip_vertex_Index::update_reverse(Fill_ip_vertex_Index& Index){
  ::pnt_upindr(Fill_ip_vertex_Manager::m_table_address,
        Index.index_address,index_address);
}

Fill_ip_vertex& Fill_ip_vertex_Index::operator()(const Panther_ID& ID) const {
  int i((int)ID);
  int id = ::pnt_getindid(Fill_ip_vertex_Manager::m_table_address,i,index_address);
#ifdef HAVE_EXCEPTION
  if( id == 0 ) throw Panther_Invalid_ID(i);
#else
  if( id == 0 ) return Fill_ip_vertex_Manager::get_NULL();
#endif
  Fill_ip_vertex_Manager& mgr = Fill_ip_vertex_Manager::get_manager();
  return mgr(Panther_ID(id));
}

Fill_ip_vertex& Fill_ip_vertex_Index::operator[](const int& nth) const {
  int id;
  if( (id = ::pnt_nth_id(Fill_ip_vertex_Manager::m_table_address, nth+1,
            index_address)) == 0 )
#ifdef HAVE_EXCEPTION
    throw Panther_Invalid_ID(nth);
#else
    return Fill_ip_vertex_Manager::get_NULL();
#endif
  Panther_ID ID(id);
  return (*this)(ID);
}
//-------------------------------------------
void Fill_ip_vertex_Selector::insert(Fill_ip_vertex& r){
  Panther_ID id = r.get_ID();
  ::pnt_inssel(Fill_ip_vertex_Manager::m_table_address,id, selector_address);
}

int Fill_ip_vertex_Selector::count(void) const {
  return ::BsCouSel(selector_address);
}

void Fill_ip_vertex_Selector::erase(int id){
  if( id == 0 )
    ::BsClrSel(selector_address);
  else
    ::BsEraSel(selector_address, id);
}

void Fill_ip_vertex_Selector::dump(void) const {
  ::pnt_shwsel(Fill_ip_vertex_Manager::m_table_address,selector_address);
}

int Fill_ip_vertex_Selector::check(Panther_ID id) const {
  int ID((int)id);
  return ::BsChkSel(selector_address, ID);
}

Fill_ip_vertex& Fill_ip_vertex_Selector::operator()(const Panther_ID& ID) const {
  int Id((int)ID), id;
  if( (id = ::pnt_getselid(Fill_ip_vertex_Manager::m_table_address,
              Id, selector_address)) == 0 )
#ifdef HAVE_EXCEPTION
    throw Panther_Invalid_ID(Id);
#else
    return Fill_ip_vertex_Manager::get_NULL();
#endif
  Fill_ip_vertex_Manager& mgr = Fill_ip_vertex_Manager::get_manager();
  return mgr(Panther_ID(id));
}
//-------------------------------------------
 std::vector<Fill_ip_vertex> point_from(const Panther_ID& Id, Fill_ip_vertex_Index& Index){
  Fill_ip_vertex_Manager& mgr = Fill_ip_vertex_Manager::get_manager();
  int n;
  int id((int)Id);
  int *v = ::pnt_ptfrom_ids(id,mgr.m_table_address,Index.index_address,&n);
  if( v != 0 ){
    std::vector<Fill_ip_vertex> vec;
    for( int i=0; i<n; i++ ){
      Panther_ID ID(*(v+i+1));
      vec.push_back( mgr(ID) );
    }
    free(v);
    return vec;
  }
  std::vector<Fill_ip_vertex> vec;
  return vec;
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
