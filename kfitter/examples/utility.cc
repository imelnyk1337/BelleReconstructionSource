#include "utility.h"
typedef Hep3Vector Vector3;
#include "toolbox/FoxWolfr.h"
#include "toolbox/FuncPtr.h"
#include "toolbox/Thrust.h"
#include "mdst/mdst.h"
#include "kfitter/kmakemother.h"
#include MDST_H
#include HEPEVT_H
#include <algorithm>

#include "belleutil/debugout.h"

namespace exkfit {

void 
makeKPi(vector<Particle> &k_p, 
	vector<Particle> &k_m, 
	vector<Particle> &pi_p, 
	vector<Particle> &pi_m)
{
  Mdst_charged_Manager &charged_mag = Mdst_charged_Manager::get_manager();
  
  //...Particle Type
  Ptype ptype_pion_plus("PI+");
  Ptype ptype_pion_minus("PI-");
  Ptype ptype_kaon_plus("K+");
  Ptype ptype_kaon_minus("K-");
 
  //...Fills pion and kaon lists with MDST_Charged Data Base
  for(vector<Mdst_charged>::iterator i = charged_mag.begin();
      i != charged_mag.end(); ++i){
    if(!good_charged(*i))continue;
    if((*i).charge() > 0.){
      Particle tmp1(*i, ptype_kaon_plus);
      k_p.push_back(tmp1);
      Particle tmp2(*i, ptype_pion_plus);
      pi_p.push_back(tmp2);
    }else{
      Particle tmp1(*i, ptype_kaon_minus);
      k_m.push_back(tmp1);
      Particle tmp2(*i, ptype_pion_minus);
      pi_m.push_back(tmp2);
    }
  }
}

void 
makeLepton(vector<Particle> &e_p,
	   vector<Particle> &e_m,
	   vector<Particle> &mu_p,
	   vector<Particle> &mu_m)
{
  Mdst_charged_Manager &charged_mag = Mdst_charged_Manager::get_manager();
  
  //...Particle Type
  Ptype ptype_elec_plus("E+");
  Ptype ptype_elec_minus("E-");
  Ptype ptype_muon_plus("MU+");
  Ptype ptype_muon_minus("MU-");
 
  //...Fills elec and muon lists with MDST_Charged Data Base
  for(vector<Mdst_charged>::iterator i = charged_mag.begin();
      i != charged_mag.end(); ++i){
    if(!good_charged(*i))continue;
    if((*i).charge() > 0.){
      Particle tmp1(*i, ptype_elec_plus);
      e_p.push_back(tmp1);
      Particle tmp2(*i, ptype_muon_plus);
      mu_p.push_back(tmp2);
    }else{
      Particle tmp1(*i, ptype_elec_minus);
      e_m.push_back(tmp1);
      Particle tmp2(*i, ptype_muon_minus);
      mu_m.push_back(tmp2);
    }
  }
}

void 
withPCut(vector<Particle> &list, 
	 const double &pL,
	 const double &pR){
  for(int i=0;i<list.size();++i){
    if(list[i].momentum().p().vect().mag() < pL ||
       list[i].momentum().p().vect().mag() > pR) {
      list.erase(&list[i]);
      --i;
    }
  }
}

void 
withPCut(vector<Particle> &list, 
	 const double &p){
  for(int i=0;i<list.size();++i){
    if(list[i].momentum().p().vect().mag() < p){
      list.erase(&list[i]);
      --i;
    }
  }
}

HepLorentzVector
pStar(Particle &p,
      double elec, 
      double posi){
  return pStar(p.momentum().p(),elec,posi);
}

HepLorentzVector
pStar(HepLorentzVector p,
      double elec, 
      double posi){
  HepLorentzVector boost_vector(0., 0., elec-posi, elec+posi);
  p.boost(-boost_vector.boostVector());
  return p;
}

void 
withPSCut(vector<Particle> &list, 
	       const double &pL,
	       const double &pR){
  for(int i=0;i<list.size();++i){
    HepLorentzVector P(pStar(list[i]));
    if(P.vect().mag() < pL ||
       P.vect().mag() > pR){
      list.erase(&list[i]);
      --i;
    }
  }
}

void 
withPSCut(vector<Particle> &list, 
	       const double &p){
  for(int i=0;i<list.size();++i){
    HepLorentzVector P(pStar(list[i]));
    if(P.vect().mag() < p){
      list.erase(&list[i]);
      --i;
    }
  }
}

void 
withSVD(vector<Particle> &list, 
	     const unsigned &n){
  for(int i=0;i<list.size();++i){
    if(list[i].mdstCharged()){
      int mhyp = 2;//pion
      if(abs(list[i].pType().lund()) == 321) mhyp = 3;
      if(abs(list[i].pType().lund()) == 11)  mhyp = 0;
      if(abs(list[i].pType().lund()) == 13)  mhyp = 1;
      if(abs(list[i].pType().lund()) == 2212)mhyp = 4;
      if(list[i].mdstCharged().trk().mhyp(mhyp).nhits(4) < n){
	list.erase(&list[i]);
	--i;
      }
    }
  }
}

void 
withSVD(vector<Particle> &list, 
	     const unsigned &n,
	     const unsigned &childID){
  for(int i=0;i<list.size();++i){
    if(list[i].relation().child(childID).mdstCharged()){
      int mhyp = 2;//pion
      if(abs(list[i].relation().child(childID).pType().lund()) == 321) mhyp = 3;
      if(abs(list[i].relation().child(childID).pType().lund()) == 11)  mhyp = 0;
      if(abs(list[i].relation().child(childID).pType().lund()) == 13)  mhyp = 1;
      if(abs(list[i].relation().child(childID).pType().lund()) == 2212)mhyp = 4;
      if(list[i].relation().child(childID).mdstCharged().trk().mhyp(mhyp).nhits(4) < n){
	list.erase(&list[i]);
	--i;
      }
    }
  }
}

void 
withMassCut(vector<Particle> &list, 
		 const double &pL,
		 const double &pR){
  for(int i=0;i<list.size();++i){
    if(list[i].momentum().mass() < pL ||
       list[i].momentum().mass() > pR){
      list.erase(&list[i]);
      --i;
    }
  }
}

void 
withMassCut(vector<Particle> &list, 
		 vector<Particle> &nlist,
		 const double &pL,
		 const double &pR){
  for(int i=0;i<list.size();++i){
    if(list[i].momentum().mass() < pL ||
       list[i].momentum().mass() > pR){
      nlist.push_back(list[i]);
      list.erase(&list[i]);
      --i;
    }
  }
}

void 
withMassDifCut(vector<Particle> &list, 
		    const double &pL,
		    const double &pR,
		    const unsigned &child){
  for(int i=0;i<list.size();++i){
    if(!withMassDifCut(list[i],pL,pR,child)){
      list.erase(&list[i]);
      --i;
    }
  }
}

unsigned
withMassDifCut(Particle &p, 
		    const double &pL,
		    const double &pR,
		    const unsigned &child){
  double massdif(p.momentum().mass()-
		 p.relation().child(child).momentum().mass());
  if(massdif < pL || massdif > pR)return 0;
  return 1;
}

void
withMuonId(vector<Particle> &list, const unsigned &th){
  for(int i=0;i<list.size();++i){
    if(list[i].mdstCharged() && 
       list[i].mdstCharged().klm() &&
       list[i].mdstCharged().klm().muon() >= th){
      ;
    }else{
      list.erase(&list[i]);
      --i;
    }
  }
}

void
withEId(vector<Particle> &list, const double &th){
  for(int i=0;i<list.size();++i){
    if(list[i].mdstCharged() && 
       list[i].mdstCharged().elid() &&
       list[i].mdstCharged().elid().le() >= th){
      ;
    }else{
      list.erase(&list[i]);
      --i;
    }
  }
}

void
makeGenHepKPi(vector<Particle> &k_p, 
		   vector<Particle> &k_m, 
		   vector<Particle> &pi_p, 
		   vector<Particle> &pi_m){
  Gen_hepevt_Manager &genMgr  = Gen_hepevt_Manager::get_manager();

  //...Particle Type
  Ptype ptype_pion_plus("PI+");
  Ptype ptype_pion_minus("PI-");
  Ptype ptype_kaon_plus("K+");
  Ptype ptype_kaon_minus("K-");

  for(vector<Gen_hepevt>::iterator i = genMgr.begin();
      i != genMgr.end(); ++i){
    if((*i).idhep() == ptype_kaon_plus.lund()){
      k_p.push_back(Particle(*i));
      continue;
    }
    if((*i).idhep() == ptype_kaon_minus.lund()){
      k_m.push_back(Particle(*i));
      continue;
    }
    if((*i).idhep() == ptype_pion_plus.lund()){
      pi_p.push_back(Particle(*i));
      continue;      
    }
    if((*i).idhep() == ptype_pion_minus.lund()){
      pi_m.push_back(Particle(*i));
      continue;
    }
  }
}

void
makeGenHepLepton(vector<Particle> &e_p,
		      vector<Particle> &e_m,
		      vector<Particle> &mu_p,
		      vector<Particle> &mu_m){
  Gen_hepevt_Manager &genMgr  = Gen_hepevt_Manager::get_manager();

  //...Particle Type
  Ptype ptype_elec_plus("E+");
  Ptype ptype_elec_minus("E-");
  Ptype ptype_muon_plus("MU+");
  Ptype ptype_muon_minus("MU-");

  for(vector<Gen_hepevt>::iterator i = genMgr.begin();
      i != genMgr.end(); ++i){
    if((*i).idhep() == ptype_elec_plus.lund()){
      e_p.push_back(Particle(*i));
      continue;
    }
    if((*i).idhep() == ptype_elec_minus.lund()){
      e_m.push_back(Particle(*i));
      continue;
    }
    if((*i).idhep() == ptype_muon_plus.lund()){
      mu_p.push_back(Particle(*i));
      continue;      
    }
    if((*i).idhep() == ptype_muon_minus.lund()){
      mu_m.push_back(Particle(*i));
      continue;
    }
  }
}

void 
setGenHepInfoF(vector<Particle> &list){
  //...final state particles
  if(list.size() == 0)return;

  for(vector<Particle>::iterator i = list.begin();
      i != list.end(); ++i){
    //...mdst_charged
    if(i->mdstCharged()){
      Gen_hepevt & hep(get_hepevt(i->mdstCharged()));
      if (hep && i->pType().lund() == hep.idhep()){
	i->relation().genHepevt(hep);
      }
    }
  }
}

void 
setUniqueGenHepInfoFBySvdAndDeltaP(vector<Particle> &list)
{
  int size = list.size();
  if(size == 0 || size == 1)return;
#if 0
  dout(Debugout::INFO,"utility") << "setUniqueGenHepInfoFBySvd--init--" << std::endl;
  for(int i=0;i<size;++i){
    if(list[i].genHepevt()){
      dout(Debugout::INFO,"utility") << i << ": hepID=" << list[i].genHepevt().get_ID() << std::endl;
    }
  }
#endif
  vector<int> list_checked;
  vector<int> list_best;
  for(int i=0;i<size;++i){
    if(!list[i].genHepevt())continue;    
    int best_index = i;
    Hep3Vector diff;
    diff.setX(list[i].p().x()-list[i].genHepevt().PX());
    diff.setY(list[i].p().y()-list[i].genHepevt().PY());
    diff.setZ(list[i].p().z()-list[i].genHepevt().PZ());
    double min_frac = diff.mag();
    int mhyp = 2;//pion
    if(abs(list[i].pType().lund()) == 321)mhyp = 3;
    if(abs(list[i].pType().lund()) == 11)mhyp = 0;
    if(abs(list[i].pType().lund()) == 13)mhyp = 1;
    if(abs(list[i].pType().lund()) == 2212)mhyp = 4;
    int svdHits = list[i].mdstCharged().trk().mhyp(mhyp).nhits(4);
    for(vector<int>::iterator it = list_checked.begin();
	it != list_checked.end();++it){
      if(list[i].genHepevt().get_ID() == *it)goto kokodayo;
    }
    for(int j=i+1;j<size;j++){
      if(!list[j].genHepevt())continue;
      if(list[i].genHepevt().get_ID() == list[j].genHepevt().get_ID()){
	Hep3Vector diff2;
	diff2.setX(list[j].p().x()-list[j].genHepevt().PX());
	diff2.setY(list[j].p().y()-list[j].genHepevt().PY());
	diff2.setZ(list[j].p().z()-list[j].genHepevt().PZ());
	int svdHits2 = list[j].mdstCharged().trk().mhyp(mhyp).nhits(4);
	if(svdHits2 > svdHits){
	  svdHits = svdHits2;
	  min_frac   = diff2.mag();
          best_index = j;
	}else if(svdHits2 == svdHits && diff2.mag() < min_frac){
	  min_frac   = diff2.mag();
	  best_index = j;
	}
      }
    }
    list_checked.push_back(list[i].genHepevt().get_ID());
    list_best.push_back(best_index);
  kokodayo:;
  }

  vector<int> remove_list;
  for(int i=0;i<size;i++){
    for(vector<int>::iterator it = list_best.begin();
	it != list_best.end();++it){
      if(i == *it)goto kokokana;
    }
    remove_list.push_back(i);
  kokokana:;
  }
  
  for(vector<int>::iterator it = remove_list.begin();
      it != remove_list.end();++it){
    list[*it].relation().resetGenHepevt();
  }
#if 0
  dout(Debugout::INFO,"utility") << "setUniqueGenHepInfoFBySvd--result--" << std::endl;
  for(int i=0;i<size;++i){
    if(list[i].genHepevt()){
      dout(Debugout::INFO,"utility") << i << ": hepID=" << list[i].genHepevt().get_ID() << std::endl;
    }
  }
#endif
}

void
addTrack2fit(kvertexfitter &kv, vector<Particle> &list){
  for(unsigned i=0;i<list.size();++i){
    addTrack2fit(kv,list[i]);
  }
}

void
addTrack2fit(kvertexfitter &kv, Particle &p)
{
  //No Correlation.
  kv.addTrack(p.momentum().p(),
	      p.momentum().x(),
	      p.momentum().dpx(),
	      p.pType().charge(),
	      p.pType().mass());
}

void
addTrack2fit(kmassvertexfitter &kmv, vector<Particle> &list){
  for(unsigned i=0;i<list.size();++i){
    addTrack2fit(kmv,list[i]);
  }
}

void
addTrack2fit(kmassvertexfitter &kmv, Particle &p)
{
  //No Correlation.
  kmv.addTrack(p.momentum().p(),
	       p.momentum().x(),
	       p.momentum().dpx(),
	       p.pType().charge(),
	       p.pType().mass());
}

unsigned
removeParticle(vector<Particle> &list, const Particle &p){
  unsigned count(0);
  for(unsigned i=0;i<list.size();++i){
#if 0
    if(list[i].relation().isIdenticalWith(p.relation())){
      list.erase(&list[i]);
      --i;
      ++count;
    }
#else
    if(list[i].relation().mdstCharged() && p.relation().mdstCharged()){
      if(list[i].relation().mdstCharged().get_ID() ==
	 p.relation().mdstCharged().get_ID()){
	list.erase(&list[i]);
	--i;
	++count;
      }
    }
#endif
  }
  return count;  
}

void 
calcuFoxWolfram(vector<Particle> &list, double *r,
		     double e, double p){
  vector<Hep3Vector> vec;
  HepLorentzVector boost_vector(0., 0., e-p, e+p);
  for(int i=0;i<list.size();i++){
    HepLorentzVector b0(list[i].p());
    b0.boost(-(boost_vector.boostVector()));
    Hep3Vector tmp(b0.vect());
    vec.push_back(tmp);
  }
  
  FoxWolfram t;

  t = foxwolfram(vec.begin(), vec.end(), SelfFunc(Hep3Vector()));

  r[0] = t.R(0);
  r[1] = t.R(1);
  r[2] = t.R(2);
}

void 
massCut(vector<Particle> &p1, 
	     vector<Particle> &p2, 
	     const double &width)
{
  for(vector<Particle>::iterator i = p1.begin();
      i != p1.end(); ++i){
    double mass = (*i).momentum().mass();
    double nominal_mass = (*i).pType().mass();
    if(nominal_mass - width <= mass && mass <= nominal_mass + width){
      p2.push_back(*i);
    }
  }
}

void 
massDifCut(vector<Particle> &p1, 
		vector<Particle> &p2, 
		const double &width, 
		const unsigned &child)
{
  for(vector<Particle>::iterator i = p1.begin();
      i != p1.end(); ++i){
    double massDif = (*i).momentum().mass()-(*i).relation().child(child).momentum().mass();
    double nominal_massDif = (*i).pType().mass()-(*i).relation().child(child).pType().mass();
    if(nominal_massDif - width <= massDif && massDif <= nominal_massDif + width){
      p2.push_back(*i);
    }
  }
}

void 
deepCopy(vector<Particle> &p1, 
	      vector<Particle> &p2)
{
  for(vector<Particle>::iterator i = p1.begin();
      i != p1.end(); ++i){
    p2.push_back((*i).deepCopy());
  }
}

void 
deleteDeepCopiedObjects(vector<Particle> &p)
{
  for(vector<Particle>::iterator i = p.begin();
      i != p.end(); ++i){
    (*i).deepDelete();
  }
}

unsigned setGenHepInfoR(Particle &);
unsigned setGenHepInfoR_sub2(Particle &);
unsigned setGenHepInfoR_sub3(Particle &);
unsigned setGenHepInfoR_sub4(Particle &);
unsigned setGenHepInfoR_sub5(Particle &);

void 
setGenHepInfoR(vector<Particle> &list)
{
  //...reconstructed particles
  if(list.size() == 0)return;
  for(vector<Particle>::iterator i = list.begin();
      i != list.end(); ++i){
    setGenHepInfoR(*i);
  }
}

unsigned
setGenHepInfoR(Particle &p)
{
  //...supports 2, 3, 4, and 5
  
  //...returns = 0 : can NOT set.
  //...returns = 1 : can set.
  
  switch(p.relation().nChildren()){
  case 2:
    return setGenHepInfoR_sub2(p);
  case 3:
    return setGenHepInfoR_sub3(p);
  case 4:
    return setGenHepInfoR_sub4(p);
  case 5:
    return setGenHepInfoR_sub5(p);
  default:
    return 0;
  }
}

unsigned
setGenHepInfoR_sub2(Particle &p)
{
  if(!p.relation().child(0).genHepevt() || !p.relation().child(1).genHepevt()){
    return 0;
  }

  if(p.relation().child(0).genHepevt().idhep() == 0 ||
     p.relation().child(1).genHepevt().idhep() == 0){
    return 0;
  }

  if(p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(1).genHepevt().get_ID()){
    return 0;
  }

  if(p.relation().child(0).genHepevt().mother().get_ID() ==
     p.relation().child(1).genHepevt().mother().get_ID()){
    if(p.relation().child(0).genHepevt().mother().daLast() -
       p.relation().child(0).genHepevt().mother().daFirst() + 1 == 2){
      if(p.pType().lund() == p.relation().child(0).genHepevt().mother().idhep()){
	p.relation().genHepevt(p.relation().child(0).genHepevt().mother());
	return 1;
      }
    }
  }
  
  return 0;
}

unsigned
setGenHepInfoR_sub3(Particle &p)
{
  if(!p.relation().child(0).genHepevt()||
     !p.relation().child(1).genHepevt()||
     !p.relation().child(2).genHepevt()){
    return 0;
  }

  if(p.relation().child(0).genHepevt().idhep() == 0 ||
     p.relation().child(1).genHepevt().idhep() == 0 ||
     p.relation().child(2).genHepevt().idhep() == 0){
    return 0;
  }
  
  if(p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(1).genHepevt().get_ID() ||
     p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(2).genHepevt().get_ID() ||
     p.relation().child(1).genHepevt().get_ID() == 
     p.relation().child(2).genHepevt().get_ID()){
    return 0;
  }

  if(p.relation().child(0).genHepevt().mother().get_ID() ==
     p.relation().child(1).genHepevt().mother().get_ID() &&
     p.relation().child(0).genHepevt().mother().get_ID() ==
     p.relation().child(2).genHepevt().mother().get_ID()){
    if(p.relation().child(0).genHepevt().mother().daLast() -
       p.relation().child(0).genHepevt().mother().daFirst() + 1 == 3){
      if(p.pType().lund() == p.relation().child(0).genHepevt().mother().idhep()){
	p.relation().genHepevt(p.relation().child(0).genHepevt().mother());
	return 1;
      }
    }
  }
 
  return 0;
}


unsigned
setGenHepInfoR_sub4(Particle &p)
{
  if(!p.relation().child(0).genHepevt() ||
     !p.relation().child(1).genHepevt() ||
     !p.relation().child(2).genHepevt() ||
     !p.relation().child(3).genHepevt()){
    return 0;
  }
  if(p.relation().child(0).genHepevt().idhep() == 0 ||
     p.relation().child(1).genHepevt().idhep() == 0 ||
     p.relation().child(2).genHepevt().idhep() == 0 ||
     p.relation().child(3).genHepevt().idhep() == 0){
    return 0;
  }
  if(p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(1).genHepevt().get_ID() ||
     p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(2).genHepevt().get_ID() ||
     p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(3).genHepevt().get_ID() ||
     p.relation().child(1).genHepevt().get_ID() == 
     p.relation().child(2).genHepevt().get_ID() ||
     p.relation().child(1).genHepevt().get_ID() == 
     p.relation().child(3).genHepevt().get_ID() ||
     p.relation().child(2).genHepevt().get_ID() == 
     p.relation().child(3).genHepevt().get_ID()){
    return 0;
  }
  if(p.relation().child(0).genHepevt().mother().get_ID() ==
     p.relation().child(1).genHepevt().mother().get_ID() &&
     p.relation().child(0).genHepevt().mother().get_ID() ==
     p.relation().child(2).genHepevt().mother().get_ID() &&
     p.relation().child(2).genHepevt().mother().get_ID() ==
     p.relation().child(3).genHepevt().mother().get_ID()){
    if(p.relation().child(0).genHepevt().mother().daLast() -
       p.relation().child(0).genHepevt().mother().daFirst() + 1 == 4){
      if(p.pType().lund() == p.relation().child(0).genHepevt().mother().idhep()){
	p.relation().genHepevt(p.relation().child(0).genHepevt().mother());
	return 1;
      }
    }
  }
 
  return 0;
}

unsigned
setGenHepInfoR_sub5(Particle &p)
{
  if(!p.relation().child(0).genHepevt() ||
     !p.relation().child(1).genHepevt() ||
     !p.relation().child(2).genHepevt() ||
     !p.relation().child(3).genHepevt() ||
     !p.relation().child(4).genHepevt()){
    return 0;
  }

  if(p.relation().child(0).genHepevt().idhep() == 0 ||
     p.relation().child(1).genHepevt().idhep() == 0 ||
     p.relation().child(2).genHepevt().idhep() == 0 ||
     p.relation().child(3).genHepevt().idhep() == 0 ||
     p.relation().child(4).genHepevt().idhep() == 0){
    return 0;
  }

  if(p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(1).genHepevt().get_ID() ||
     p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(2).genHepevt().get_ID() ||
     p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(3).genHepevt().get_ID() ||
     p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(4).genHepevt().get_ID() ||
     p.relation().child(1).genHepevt().get_ID() == 
     p.relation().child(2).genHepevt().get_ID() ||
     p.relation().child(1).genHepevt().get_ID() == 
     p.relation().child(3).genHepevt().get_ID() ||
     p.relation().child(1).genHepevt().get_ID() == 
     p.relation().child(4).genHepevt().get_ID() ||
     p.relation().child(2).genHepevt().get_ID() == 
     p.relation().child(3).genHepevt().get_ID() ||
     p.relation().child(2).genHepevt().get_ID() == 
     p.relation().child(4).genHepevt().get_ID() ||
     p.relation().child(3).genHepevt().get_ID() == 
     p.relation().child(4).genHepevt().get_ID()){
    return 0;
  }

  if(p.relation().child(0).genHepevt().mother().get_ID() ==
     p.relation().child(1).genHepevt().mother().get_ID() &&
     p.relation().child(2).genHepevt().mother().get_ID() ==
     p.relation().child(3).genHepevt().mother().get_ID() &&
     p.relation().child(0).genHepevt().mother().get_ID() ==
     p.relation().child(2).genHepevt().mother().get_ID() &&
     p.relation().child(0).genHepevt().mother().get_ID() ==
     p.relation().child(4).genHepevt().mother().get_ID()){
    if(p.relation().child(0).genHepevt().mother().daLast() -
       p.relation().child(0).genHepevt().mother().daFirst() + 1 == 5){
      if(p.pType().lund() == p.relation().child(0).genHepevt().mother().idhep()){
	p.relation().genHepevt(p.relation().child(0).genHepevt().mother());
	return 1;
      }
    }
  }
 
  return 0;
}

double
beamEnergyConstraint(const Particle &b, double e, double p)
{
  //b -- generally B0 or B0B
  //e -- electron beam energy
  //p -- positron beam enegry
  HepLorentzVector boost_vector(0., 0., e-p, e+p);
  HepLorentzVector b0(b.p());
  b0.boost(-(boost_vector.boostVector()));
  double mass2 = e*p - b0.vect().mag2();
  double mass  = (mass2 > 0.) ? sqrt(mass2) :  -sqrt(-mass2);
  return mass;
}

Hep3Vector
calcuThrust(vector<Particle> &list, double e, double p)
{
  vector<Hep3Vector> vec;
  HepLorentzVector boost_vector(0., 0., e-p, e+p);
  for(unsigned i=0;i<list.size();++i){
    HepLorentzVector b0(list[i].p());
    b0.boost(-(boost_vector.boostVector()));
    Hep3Vector tmp(b0.vect());
    vec.push_back(tmp);
  }
  
  return thrust(vec.begin(), vec.end(), SelfFunc(Hep3Vector()));
}

void 
makeGamma(vector<Particle> &gamma)
{
  Mdst_gamma_Manager &gamma_mag = Mdst_gamma_Manager::get_manager();
  for(vector<Mdst_gamma>::iterator i = gamma_mag.begin();
      i != gamma_mag.end(); ++i){
    gamma.push_back(Particle(*i));
  }
}

unsigned
makeMother(kvertexfitter &kv,
	   Particle &mother)
{
  unsigned n = kv.tracks();
  kmakemother kmm;
  for(unsigned i=0;i<n;++i){
    kmm.addTrack(kv.momentum(i),
		 kv.position(i),
		 kv.error(i),
		 mother.relation().child(i).pType().charge());
    kmm.errVertexTrack(kv.errVertexTrack(i));
    for(unsigned j=i+1;j<n;++j){
      kmm.correlation(kv.correlation(i,j));
      HepMatrix tmp1(kv.correlation(i,j));
      mother.relation().child(i).fittedMomentum().coMatrix(&mother.relation().child(j),tmp1);
    }
    mother.relation().child(i).fittedMomentum().momentumPosition(kv.momentum(i),
								 kv.position(i),
								 kv.error(i));
    mother.relation().child(i).fittedMomentum().vertex(kv.vertex(), kv.errVertex());
    HepMatrix tmp2(kv.errVertexTrack(i));
    mother.relation().child(i).fittedMomentum().coVertex(tmp2);
    mother.relation().child(i).fittedMomentum().cl(kv.cl());
    mother.relation().child(i).fittedMomentum().chisq(kv.chisq());
    mother.relation().child(i).fittedMomentum().dof(kv.dgf());
  }
  for(unsigned i=0;i<n;++i){
    for(unsigned j=i+1;j<n;++j){
      HepMatrix tmp1(kv.correlation(j,i));
      mother.relation().child(j).fittedMomentum().coMatrix(&mother.relation().child(i),tmp1);
    }
  }
  kmm.vertex(kv.vertex());
  kmm.errVertex(kv.errVertex());
  unsigned err = kmm.make();
  if(err != 0)return 0;
  mother.momentum().momentumPosition(kmm.momentum(),
				     kmm.position(),
				     kmm.error());
  mother.momentum().decayVertex(kv.vertex(), kv.errVertex());
  return 1;
}

unsigned
makeMother(kmassvertexfitter &kmv,
	   Particle &mother)
{
  unsigned n = kmv.tracks();
  kmakemother kmm;
  for(unsigned i=0;i<n;++i){
    kmm.addTrack(kmv.momentum(i),
		 kmv.position(i),
		 kmv.error(i),
		 mother.relation().child(i).pType().charge());
    kmm.errVertexTrack(kmv.errVertexTrack(i));
    for(unsigned j=i+1;j<n;++j){
      kmm.correlation(kmv.correlation(i,j));
      HepMatrix tmp1(kmv.correlation(i,j));
      mother.relation().child(i).fittedMomentum().coMatrix(&mother.relation().child(j),tmp1);
    }
    mother.relation().child(i).fittedMomentum().momentumPosition(kmv.momentum(i),
								 kmv.position(i),
								 kmv.error(i));
    mother.relation().child(i).fittedMomentum().vertex(kmv.vertex(), kmv.errVertex());
    HepMatrix tmp2(kmv.errVertexTrack(i));
    mother.relation().child(i).fittedMomentum().coVertex(tmp2);
    mother.relation().child(i).fittedMomentum().cl(kmv.cl());
    mother.relation().child(i).fittedMomentum().chisq(kmv.chisq());
    mother.relation().child(i).fittedMomentum().dof(kmv.dgf());
  }
  for(unsigned i=0;i<n;++i){
    for(unsigned j=i+1;j<n;++j){
      HepMatrix tmp1(kmv.correlation(j,i));
      mother.relation().child(j).fittedMomentum().coMatrix(&mother.relation().child(i),tmp1);
    }
  }
  kmm.vertex(kmv.vertex());
  kmm.errVertex(kmv.errVertex());
  unsigned err = kmm.make();
  if(err != 0)return 0;
  mother.momentum().momentumPosition(kmm.momentum(),
				     kmm.position(),
				     kmm.error());
  mother.momentum().decayVertex(kmv.vertex(), kmv.errVertex());
  return 1;
}

}
