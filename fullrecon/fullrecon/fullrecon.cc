#include <iostream>
#include "tables/fullrecon.h"
#include "belleutil/debugout.h"

#ifndef PANTHER_INCLUDE_BRECON_H_
#define PANTHER_INCLUDE_BRECON_H_
#include BRECON_H
#endif
#if defined(BELLE_NAMESPACE)
  namespace Belle {
#endif

// #ifdef __DECCXX
// #pragma define_template Table_Manager<Fullrecon>
// #else
// template class Table_Manager<Fullrecon>;
// #endif /* __DECCXX */

//---------- Table FULLRECON ----------
Fullrecon_Manager* Fullrecon_Manager::manager = 0;
Fullrecon*         Fullrecon_Manager::null    = 0;
template<> PNT_Table *Table_Manager<Fullrecon>::m_table_address = 0;
void Fullrecon::newent(void){
  m_p = (fullrecon *)::BsNewEnt(FULLRECON);
}
Fullrecon::Fullrecon(const Panther_ID& ID, PNT_TableIndex* index){
  m_p = (fullrecon *)::BsGetEnt(FULLRECON,(int)ID, index);
}
void Fullrecon::FullBrec(const Brecon &i){ m_p->m_FullBrec = (int)i.get_ID(); }
void Fullrecon::reset_FullBrec(void){ m_p->m_FullBrec = 0; }
  extern "C" void *(*Panther_Get_Manager_Address)(void*);
  extern "C" void  (*Panther_Remove_Wrapper)(void*,int);
Fullrecon_Manager::Fullrecon_Manager() : Table_Manager<Fullrecon>(), m_selector(NULL) { }
Fullrecon_Manager::~Fullrecon_Manager() {
  delete m_selector;
}

Fullrecon_Manager& Fullrecon_Manager::get_manager_impl(void){
  if( manager == 0 ){
    manager = new Fullrecon_Manager;
    null    = new Fullrecon;
    m_table_address = pnt_get_table(FULLRECON);
    void* addr = FULLRECON;
    pnt_register_manager((Table_Manager_Base*)manager,addr);
    Panther_Get_Manager_Address = ::pnt_get_manager_address;
    Panther_Remove_Wrapper      = ::remove_wrapper;
  }
  return *manager;
}

Fullrecon& Fullrecon_Manager::get_NULL(void){
  if( manager == 0 )
    Fullrecon_Manager::get_manager();
  return *null;
}

void Fullrecon_Manager::delete_NULL(void){
  if( null != 0 ) {
    delete null;
    null = NULL;
  }
}

void Fullrecon_Manager::remove_wrapper(int id){
  if( id == 0 ){
    erase( begin(), end() );
  }else{
    iterator k;
    for( k=begin(); k != end(); k++ ){
      if( (*k).get_ID() == id ) break;
    }
    if( k == end() ){
      dout(Debugout::ERR,"Fullrecon") << "Remove_wrapper warning: unknown object (ID=" << id << ")." << std::endl;
      return;
    }
    erase(k);
  }
}

Fullrecon_Index::Fullrecon_Index() { }

Fullrecon_Index::Fullrecon_Index(const char *name) : Panther_Index<Fullrecon>(name){ }

Fullrecon_Index::~Fullrecon_Index(){  
  ::BsDelInd(index_address); }

Fullrecon_Selector::Fullrecon_Selector() { }
Fullrecon_Selector::~Fullrecon_Selector() {  
::BsDelSel(selector_address); }

Fullrecon_Index Fullrecon_Manager::index(const char* attr) {
  return Fullrecon_Index(attr);
}

Fullrecon_Selector& Fullrecon_Manager::selector(void) {
  if(m_selector) return *m_selector;
  m_selector = new Fullrecon_Selector;
  m_selector->selector_address = ::pnt_maksel(m_table_address);
  return *m_selector;
}
//-------------------------------------------
void Fullrecon_Index::update(void){
  ::pnt_updind(Fullrecon_Manager::m_table_address, index_address);
}

void Fullrecon_Index::update(Fullrecon_Index& Index){
  ::pnt_upind2(Fullrecon_Manager::m_table_address,
        Index.index_address,index_address);
}

void Fullrecon_Index::update_reverse(Fullrecon_Index& Index){
  ::pnt_upindr(Fullrecon_Manager::m_table_address,
        Index.index_address,index_address);
}

Fullrecon& Fullrecon_Index::operator()(const Panther_ID& ID) const {
  int i((int)ID);
  int id = ::pnt_getindid(Fullrecon_Manager::m_table_address,i,index_address);
#ifdef HAVE_EXCEPTION
  if( id == 0 ) throw Panther_Invalid_ID(i);
#else
  if( id == 0 ) return Fullrecon_Manager::get_NULL();
#endif
  Fullrecon_Manager& mgr = Fullrecon_Manager::get_manager();
  return mgr(Panther_ID(id));
}

Fullrecon& Fullrecon_Index::operator[](const int& nth) const {
  int id;
  if( (id = ::pnt_nth_id(Fullrecon_Manager::m_table_address, nth+1,
            index_address)) == 0 )
#ifdef HAVE_EXCEPTION
    throw Panther_Invalid_ID(nth);
#else
    return Fullrecon_Manager::get_NULL();
#endif
  Panther_ID ID(id);
  return (*this)(ID);
}
//-------------------------------------------
void Fullrecon_Selector::insert(Fullrecon& r){
  Panther_ID id = r.get_ID();
  ::pnt_inssel(Fullrecon_Manager::m_table_address,id, selector_address);
}

int Fullrecon_Selector::count(void) const {
  return ::BsCouSel(selector_address);
}

void Fullrecon_Selector::erase(int id){
  if( id == 0 )
    ::BsClrSel(selector_address);
  else
    ::BsEraSel(selector_address, id);
}

void Fullrecon_Selector::dump(void) const {
  ::pnt_shwsel(Fullrecon_Manager::m_table_address,selector_address);
}

int Fullrecon_Selector::check(Panther_ID id) const {
  int ID((int)id);
  return ::BsChkSel(selector_address, ID);
}

Fullrecon& Fullrecon_Selector::operator()(const Panther_ID& ID) const {
  int Id((int)ID), id;
  if( (id = ::pnt_getselid(Fullrecon_Manager::m_table_address,
              Id, selector_address)) == 0 )
#ifdef HAVE_EXCEPTION
    throw Panther_Invalid_ID(Id);
#else
    return Fullrecon_Manager::get_NULL();
#endif
  Fullrecon_Manager& mgr = Fullrecon_Manager::get_manager();
  return mgr(Panther_ID(id));
}
//-------------------------------------------
 std::vector<Fullrecon> point_from(const Panther_ID& Id, Fullrecon_Index& Index){
  Fullrecon_Manager& mgr = Fullrecon_Manager::get_manager();
  int n;
  int id((int)Id);
  int *v = ::pnt_ptfrom_ids(id,mgr.m_table_address,Index.index_address,&n);
  if( v != 0 ){
    std::vector<Fullrecon> vec;
    for( int i=0; i<n; i++ ){
      Panther_ID ID(*(v+i+1));
      vec.push_back( mgr(ID) );
    }
    free(v);
    return vec;
  }
  std::vector<Fullrecon> vec;
  return vec;
}
Brecon& Fullrecon::FullBrec(void) const {
  Brecon_Manager& m = Brecon_Manager::get_manager();
  if( m_p->m_FullBrec < 1  ||  m_p->m_FullBrec > m.count() ){
    return m.get_NULL();
  }
  Panther_ID ID(m_p->m_FullBrec);
  Brecon& e = m(ID);
  return e;
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
