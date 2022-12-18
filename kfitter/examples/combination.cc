#include "combination.h"

#include "panther/panther.h"
#include MDST_H

#include "belleutil/debugout.h"

namespace exkfit {

  using std::Particle;

  bool 
    checkSame(const Particle &p1,
	      const Particle &p2){
    //if(p1.mdstCharged() == p2.mdstCharged())return true;
    if(p1.mdstCharged() && p2.mdstCharged() &&
       p1.mdstCharged().get_ID() == p2.mdstCharged().get_ID())return true;
    return false;
  }
  
  /************************/
  /* 2 Bodies Combination */
  /************************/
  //...p1 != p2
  void combination_type01(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2);
  void combination_type01(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  const double &massL,
			  const double &massR);
  //...p1 == p2
  void combination_type02(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2);
  void combination_type02(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  const double &massL,
			  const double &massR);
  /************************/
  /* 3 Bodies Combination */
  /************************/
  //...p1 != p2, p1 != p3, p2 != p3
  void combination_type01(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3);
  void combination_type01(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  const double &massL,
			  const double &massR);
  //...p1 == p2, p1 != p3
  //...p1 == p3, p1 != p2
  //...p2 == p3, p2 != p1
  void combination_type02(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3);
  void combination_type02(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  const double &massL,
			  const double &massR);
  //...p1 == p2 == p3
  void combination_type03(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3);
  void combination_type03(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  const double &massL,
			  const double &massR);
  /************************/
  /* 4 Bodies Combination */
  /************************/
  //...p1 != p2, p1 != p3, p1 != p4, p2 != p3, p2 != p4, p3 != p4
  void combination_type01(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4);
  void combination_type01(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  const double &massL,
			  const double &massR);
  //...p1 == p2, p1 != p3, p1 != p4, p3 != p4
  //...p1 == p3, p1 != p2, p1 != p4, p2 != p4
  //...p1 == p4, p1 != p2, p1 != p3, p2 != p3
  void combination_type02(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4);
  void combination_type02(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  const double &massL,
			  const double &massR);
  //...p1 == p2, p3 == p4, p1 != p3
  //...p1 == p3, p2 == p4, p1 != p2
  //...p1 == p4, p2 == p3, p1 != p2
  void combination_type03(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4);
  void combination_type03(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  const double &massL,
			  const double &massR);
  //...p1 == p2 == p3, p1 != p4
  //...p1 == p2 == p4, p1 != p3
  //...p1 == p3 == p4, p1 != p2
  //...p2 == p3 == p4, p2 != p1
  void combination_type04(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4);
  void combination_type04(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  const double &massL,
			  const double &massR);
  //...p1 == p2 == p3 == p4
  void combination_type05(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4);
  void combination_type05(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  const double &massL,
			  const double &massR);
  /************************/
  /* 5 Bodies Combination */
  /************************/
  //...p1 != p2, p1 != p3, p1 != p4, p1 != p5, p2 != p3, p2 != p4, p2 != p5, p3 != p4, p3 != p5, p4 != p5
  void combination_type01(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  vector<Particle> &p5);
  void combination_type01(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  vector<Particle> &p5,
			  const double &massL,
			  const double &massR);
  //...p1 == p2, p1 != p3, p1 != p4, p1 != p5, p3 != p4, p3 != p5, p4 != p5
  //...p1 == p3, p1 != p3, p1 != p4, p1 != p5, p3 != p4, p3 != p5, p4 != p5
  //...p1 == p4, p1 != p3, p1 != p4, p1 != p5, p3 != p4, p3 != p5, p4 != p5
  //...p1 == p5, p1 != p3, p1 != p4, p1 != p5, p3 != p4, p3 != p5, p4 != p5
  //...
  void combination_type02(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  vector<Particle> &p5);
  void combination_type02(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  vector<Particle> &p5,
			  const double &massL,
			  const double &massR);
  //...p1 == p2, p3 == p4, p1 != p3, p1 != p5, p3 != p5
  //...
  void combination_type03(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  vector<Particle> &p5);
  void combination_type03(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  vector<Particle> &p5,
			  const double &massL,
			  const double &massR);
  //...p1 == p2 == p3, p4 != p5, p1 != p4
  //...
  void combination_type04(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  vector<Particle> &p5);
  void combination_type04(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  vector<Particle> &p5,
			  const double &massL,
			  const double &massR);
  //...p1 == p2 == p3, p4 == p5, p1 != p4
  //...
  void combination_type05(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  vector<Particle> &p5);
  void combination_type05(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  vector<Particle> &p5,
			  const double &massL,
			  const double &massR);
  //...p1 == p2 == p3 == p4, p1 != p5
  //...
  void combination_type06(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  vector<Particle> &p5);
  void combination_type06(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  vector<Particle> &p5,
			  const double &massL,
			  const double &massR);
  //...p1 == p2 == p3 == p4 == p5
  void combination_type07(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  vector<Particle> &p5);
  void combination_type07(vector<Particle> &new_p,
			  const Ptype &ptype,
			  vector<Particle> &p1,
			  vector<Particle> &p2,
			  vector<Particle> &p3,
			  vector<Particle> &p4,
			  vector<Particle> &p5,
			  const double &massL,
			  const double &massR);
  /*************/
  /* functions */
  /*************/
  
  
  /************/
  /* 2 bodies */
  /************/
  void 
    combination(vector<Particle> &new_p,
		const Ptype &ptype,
		vector<Particle> &p1,
		vector<Particle> &p2){
    if(&p1 != &p2){
      combination_type01(new_p, ptype, p1, p2);
    }else{
      combination_type02(new_p, ptype, p1, p2);
    }
  }

  void 
    combination(vector<Particle> &new_p,
		const Ptype &ptype,
		vector<Particle> &p1,
		vector<Particle> &p2,
		const double &mass_width){
    if(&p1 != &p2){
      combination_type01(new_p, ptype, p1, p2, mass_width, mass_width);
    }else{
      combination_type02(new_p, ptype, p1, p2, mass_width, mass_width);
    }
  }

  void 
    combination(vector<Particle> &new_p,
		const Ptype &ptype,
		vector<Particle> &p1,
		vector<Particle> &p2,
		const double &massL,
		const double &massR){
    if(&p1 != &p2){
      combination_type01(new_p, ptype, p1, p2, massL, massR);
    }else{
      combination_type02(new_p, ptype, p1, p2, massL, massR);
    }
  }
  
  void 
    combination_nocut(vector<Particle> &new_p,
		      const Ptype &ptype,
		      Particle &p1,
		      Particle &p2){
    if(checkSame(p1,p2))return;
    Particle cand(p1.momentum().p() + p2.momentum().p(), ptype);
    cand.relation().append(p1);
    cand.relation().append(p2);
    new_p.push_back(cand);
  }

void 
combination_cut(vector<Particle> &new_p,
		     const Ptype &ptype,
		     Particle &p1,
		     Particle &p2,
		     const double &massL,
		     const double &massR)
{
  if(checkSame(p1,p2))return;
  //...checks mass
  double mass = (p1.momentum().p() + p2.momentum().p()).mag();
  double new_mass = ptype.mass();
  if(new_mass - massL <= mass && mass <= new_mass + massR){
    Particle cand(p1.momentum().p() + p2.momentum().p(), ptype);
    cand.relation().append(p1);
    cand.relation().append(p2);
    new_p.push_back(cand);
  }
}

void 
combination_type01(vector<Particle> &new_p,
			const Ptype &ptype,
			vector<Particle> &p1,
			vector<Particle> &p2)
{
  for(vector<Particle>::iterator i = p1.begin();
      i != p1.end(); ++i){
    for(vector<Particle>::iterator j = p2.begin();
        j != p2.end(); ++j){
      combination_nocut(new_p,ptype,*i,*j);
    }
  }
}

void 
combination_type01(vector<Particle> &new_p,
			const Ptype &ptype,
			vector<Particle> &p1,
			vector<Particle> &p2,
			const double &massL,
			const double &massR)
{
  for(vector<Particle>::iterator i = p1.begin();
      i != p1.end(); ++i){
    for(vector<Particle>::iterator j = p2.begin();
        j != p2.end(); ++j){
      combination_cut(new_p,ptype,*i,*j,massL,massR);
    }
  }
}

void 
combination_type02(vector<Particle> &new_p,
			const Ptype &ptype,
			vector<Particle> &p1,
			vector<Particle> &p2)
{
  if(p1.size() < 2)return;
  for(vector<Particle>::iterator i = p1.begin();
      i != p1.end()-1; ++i){
    for(vector<Particle>::iterator j = i+1;
        j != p1.end(); ++j){
      combination_nocut(new_p,ptype,*i,*j);
    }
  }
}

void 
combination_type02(vector<Particle> &new_p,
			const Ptype &ptype,
			vector<Particle> &p1,
			vector<Particle> &p2,
			const double &massL,
		   const double &massR)
{
  if(p1.size() < 2)return;
  for(vector<Particle>::iterator i = p1.begin();
      i != p1.end()-1; ++i){
    for(vector<Particle>::iterator j = i+1;
        j != p1.end(); ++j){
      combination_cut(new_p,ptype,*i,*j,massL,massR);
    }
  }
}

/************/
/* 3 bodies */
/************/
void 
combination(vector<Particle> &new_p,
		 const Ptype &ptype,
		 vector<Particle> &p1,
		 vector<Particle> &p2,
		 vector<Particle> &p3)
{
  if(&p1 != &p2 &&
     &p1 != &p3 &&
     &p2 != &p3){
    combination_type01(new_p, ptype, p1, p2, p3);
  }else if((&p1 == &p2) &&
	   (&p1 != &p3)){
    combination_type02(new_p, ptype, p1, p2, p3);
  }else if((&p1 == &p3) &&
	   (&p1 != &p2)){
    combination_type02(new_p, ptype, p1, p3, p2);
  }else if((&p2 == &p3) &&
	   (&p2 != &p1)){
    combination_type02(new_p, ptype, p2, p3, p1);
  }else{
    dout(Debugout::ERR,"combination") << "NOT Support in 3 bodies" << std::endl;
  }
}

void 
combination(vector<Particle> &new_p,
		 const Ptype &ptype,
		 vector<Particle> &p1,
		 vector<Particle> &p2,
		 vector<Particle> &p3,
		 const double &mass_width)
{
  if(&p1 != &p2 &&
     &p1 != &p3 &&
     &p2 != &p3){
    combination_type01(new_p, ptype, p1, p2, p3, mass_width, mass_width);
  }else if((&p1 == &p2) &&
	   (&p1 != &p3)){
    combination_type02(new_p, ptype, p1, p2, p3, mass_width, mass_width);
  }else if((&p1 == &p3) &&
	   (&p1 != &p2)){
    combination_type02(new_p, ptype, p1, p3, p2, mass_width, mass_width);
  }else if((&p2 == &p3) &&
	   (&p2 != &p1)){
    combination_type02(new_p, ptype, p2, p3, p1, mass_width, mass_width);
  }else{
    dout(Debugout::ERR,"combination") << "NOT Support in 3 bodies" << std::endl;
  }
}

void 
combination(vector<Particle> &new_p,
		 const Ptype &ptype,
		 vector<Particle> &p1,
		 vector<Particle> &p2,
		 vector<Particle> &p3,
		 const double &massL,
		 const double &massR)
{
  if(&p1 != &p2 &&
     &p1 != &p3 &&
     &p2 != &p3){
    combination_type01(new_p, ptype, p1, p2, p3, massL, massR);
  }else if((&p1 == &p2) &&
	   (&p1 != &p3)){
    combination_type02(new_p, ptype, p1, p2, p3, massL, massR);
  }else if((&p1 == &p3) &&
	   (&p1 != &p2)){
    combination_type02(new_p, ptype, p1, p3, p2, massL, massR);
  }else if((&p2 == &p3) &&
	   (&p2 != &p1)){
    combination_type02(new_p, ptype, p2, p3, p1, massL, massR);
  }else{
    dout(Debugout::ERR,"combination") << "NOT Support in 3 bodies" << std::endl;
  }
}

void 
combination_nocut(vector<Particle> &new_p,
		       const Ptype &ptype,
		       Particle &p1,
		       Particle &p2,
		       Particle &p3)
{
  if(checkSame(p1,p2))return;
  if(checkSame(p1,p3))return;
  if(checkSame(p2,p3))return;
  Particle cand(p1.momentum().p() + p2.momentum().p() + p3.momentum().p(), ptype);
  cand.relation().append(p1);
  cand.relation().append(p2);
  cand.relation().append(p3);
  new_p.push_back(cand); 
}

void 
combination_cut(vector<Particle> &new_p,
		const Ptype &ptype,
		Particle &p1,
		Particle &p2,
		Particle &p3,
		const double &massL,
		const double &massR){
  if(checkSame(p1,p2))return;
  if(checkSame(p1,p3))return;
  if(checkSame(p2,p3))return;
  //...checks mass
  double mass = (p1.momentum().p() + p2.momentum().p() + p3.momentum().p()).mag();
  double new_mass = ptype.mass();
  if(new_mass - massL <= mass && mass <= new_mass + massR){
    Particle cand(p1.momentum().p() + p2.momentum().p() + p3.momentum().p(), ptype);
    cand.relation().append(p1);
    cand.relation().append(p2);
    cand.relation().append(p3);
    new_p.push_back(cand);
  }
}

void 
combination_type01(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   vector<Particle> &p3)
{
  for(vector<Particle>::iterator i = p1.begin();
      i != p1.end(); ++i){
    for(vector<Particle>::iterator j = p2.begin();
        j != p2.end(); ++j){
      for(vector<Particle>::iterator k = p3.begin();
	  k != p3.end(); ++k){
	combination_nocut(new_p,ptype,*i,*j,*k);
      }
    }
  }
}

void 
combination_type01(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   vector<Particle> &p3,
		   const double &massL,
		   const double &massR)
{
  for(vector<Particle>::iterator i = p1.begin();
      i != p1.end(); ++i){
    for(vector<Particle>::iterator j = p2.begin();
        j != p2.end(); ++j){
      for(vector<Particle>::iterator k = p3.begin();
	  k != p3.end(); ++k){
	combination_cut(new_p,ptype,*i,*j,*k,massL,massR);
      }
    }
  }
}

void 
combination_type02(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   vector<Particle> &p3)
{
  if(p1.size() < 2)return;
  for(vector<Particle>::iterator i = p1.begin();
      i != p1.end()-1; ++i){
    for(vector<Particle>::iterator j = i+1;
        j != p1.end(); ++j){
      for(vector<Particle>::iterator k = p3.begin();
	  k != p3.end(); ++k){
	combination_nocut(new_p,ptype,*i,*j,*k);
      }
    }
  }
}

void 
combination_type02(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   vector<Particle> &p3,
		   const double &massL,
		   const double &massR)
{
  if(p1.size() < 2)return;
  for(vector<Particle>::iterator i = p1.begin();
      i != p1.end()-1; ++i){
    for(vector<Particle>::iterator j = i+1;
        j != p1.end(); ++j){
      for(vector<Particle>::iterator k = p3.begin();
	  k != p3.end(); ++k){
	combination_cut(new_p,ptype,*i,*j,*k,massL,massR);
      }
    }
  }
}

/************/
/* 4 bodies */
/************/
void 
combination(vector<Particle> &new_p,
	    const Ptype &ptype,
	    vector<Particle> &p1,
	    vector<Particle> &p2,
	    vector<Particle> &p3,
	    vector<Particle> &p4)
{
  if(&p1 != &p2 &&
     &p1 != &p3 &&
     &p1 != &p4 &&
     &p2 != &p3 &&
     &p2 != &p4 &&
     &p3 != &p4){
    combination_type01(new_p, ptype, p1, p2, p3, p4);
  }else if((&p1 == &p2) &&
	   (&p1 != &p3) &&
	   (&p1 != &p4) &&
	   (&p3 != &p4)){
    combination_type02(new_p, ptype, p1, p2, p3, p4);
  }else if((&p1 == &p3) &&
	   (&p1 != &p2) &&
	   (&p1 != &p4) &&
	   (&p2 != &p4)){
    combination_type02(new_p, ptype, p1, p3, p2, p4);
  }else if((&p1 == &p4) &&
	   (&p1 != &p2) &&
	   (&p1 != &p3) &&
	   (&p2 != &p3)){
    combination_type02(new_p, ptype, p1, p4, p2, p3);
  }else if((&p2 == &p3) &&
	   (&p2 != &p1) &&
	   (&p2 != &p4) &&
	   (&p1 != &p4)){
    combination_type02(new_p, ptype, p2, p3, p1, p4);
  }else if((&p2 == &p4) &&
	   (&p2 != &p1) &&
	   (&p2 != &p3) &&
	   (&p1 != &p3)){
    combination_type02(new_p, ptype, p2, p4, p1, p3);
  }else if((&p3 == &p4) &&
	   (&p3 != &p1) &&
	   (&p3 != &p2) &&
	   (&p1 != &p2)){
    combination_type02(new_p, ptype, p3, p4, p1, p2);
  }else{
    dout(Debugout::INFO,"combination") << "NOT Support in 4 bodies" << std::endl;
  }
}

void 
combination(vector<Particle> &new_p,
	    const Ptype &ptype,
	    vector<Particle> &p1,
	    vector<Particle> &p2,
	    vector<Particle> &p3,
	    vector<Particle> &p4,
	    const double &mass_width)
{
  if(&p1 != &p2 &&
     &p1 != &p3 &&
     &p1 != &p4 &&
     &p2 != &p3 &&
     &p2 != &p4 &&
     &p3 != &p4){
    combination_type01(new_p, ptype, p1, p2, p3, p4, mass_width, mass_width);
  }else if((&p1 == &p2) &&
	   (&p1 != &p3) &&
	   (&p1 != &p4) &&
	   (&p3 != &p4)){
    combination_type02(new_p, ptype, p1, p2, p3, p4, mass_width, mass_width);
  }else if((&p1 == &p3) &&
	   (&p1 != &p2) &&
	   (&p1 != &p4) &&
	   (&p2 != &p4)){
    combination_type02(new_p, ptype, p1, p3, p2, p4, mass_width, mass_width);
  }else if((&p1 == &p4) &&
	   (&p1 != &p2) &&
	   (&p1 != &p3) &&
	   (&p2 != &p3)){
    combination_type02(new_p, ptype, p1, p4, p2, p3, mass_width, mass_width);
  }else if((&p2 == &p3) &&
	   (&p2 != &p1) &&
	   (&p2 != &p4) &&
	   (&p1 != &p4)){
    combination_type02(new_p, ptype, p2, p3, p1, p4, mass_width, mass_width);
  }else if((&p2 == &p4) &&
	   (&p2 != &p1) &&
	   (&p2 != &p3) &&
	   (&p1 != &p3)){
    combination_type02(new_p, ptype, p2, p4, p1, p3, mass_width, mass_width);
  }else if((&p3 == &p4) &&
	   (&p3 != &p1) &&
	   (&p3 != &p2) &&
	   (&p1 != &p2)){
    combination_type02(new_p, ptype, p3, p4, p1, p2, mass_width, mass_width);
  }else{
    dout(Debugout::INFO,"combination") << "NOT Support in 4 bodies" << std::endl;
  }
}

void 
combination(vector<Particle> &new_p,
	    const Ptype &ptype,
	    vector<Particle> &p1,
	    vector<Particle> &p2,
	    vector<Particle> &p3,
	    vector<Particle> &p4,
	    const double &massL,
	    const double &massR)
{
  if(&p1 != &p2 &&
     &p1 != &p3 &&
     &p1 != &p4 &&
     &p2 != &p3 &&
     &p2 != &p4 &&
     &p3 != &p4){
    combination_type01(new_p, ptype, p1, p2, p3, p4, massL, massR);
  }else if((&p1 == &p2) &&
	   (&p1 != &p3) &&
	   (&p1 != &p4) &&
	   (&p3 != &p4)){
    combination_type02(new_p, ptype, p1, p2, p3, p4, massL, massR);
  }else if((&p1 == &p3) &&
	   (&p1 != &p2) &&
	   (&p1 != &p4) &&
	   (&p2 != &p4)){
    combination_type02(new_p, ptype, p1, p3, p2, p4, massL, massR);
  }else if((&p1 == &p4) &&
	   (&p1 != &p2) &&
	   (&p1 != &p3) &&
	   (&p2 != &p3)){
    combination_type02(new_p, ptype, p1, p4, p2, p3, massL, massR);
  }else if((&p2 == &p3) &&
	   (&p2 != &p1) &&
	   (&p2 != &p4) &&
	   (&p1 != &p4)){
    combination_type02(new_p, ptype, p2, p3, p1, p4, massL, massR);
  }else if((&p2 == &p4) &&
	   (&p2 != &p1) &&
	   (&p2 != &p3) &&
	   (&p1 != &p3)){
    combination_type02(new_p, ptype, p2, p4, p1, p3, massL, massR);
  }else if((&p3 == &p4) &&
	   (&p3 != &p1) &&
	   (&p3 != &p2) &&
	   (&p1 != &p2)){
    combination_type02(new_p, ptype, p3, p4, p1, p2, massL, massR);
  }else{
    dout(Debugout::INFO,"combination") << "NOT Support in 4 bodies" << std::endl;
  }
}

void 
combination_nocut(vector<Particle> &new_p,
		  const Ptype &ptype,
		  Particle &p1,
		  Particle &p2,
		  Particle &p3,
		  Particle &p4)
{
  if(checkSame(p1,p2))return;
  if(checkSame(p1,p3))return;
  if(checkSame(p1,p4))return;
  if(checkSame(p2,p3))return;
  if(checkSame(p2,p4))return;
  if(checkSame(p3,p4))return;
  Particle cand(p1.momentum().p() + p2.momentum().p() + p3.momentum().p() + p4.momentum().p(), ptype);
  cand.relation().append(p1);
  cand.relation().append(p2);
  cand.relation().append(p3);
  cand.relation().append(p4);
  new_p.push_back(cand);
}

void 
combination_cut(vector<Particle> &new_p,
		const Ptype &ptype,
		Particle &p1,
		Particle &p2,
		Particle &p3,
		Particle &p4,
		const double &massL,
		const double &massR)
{
  if(checkSame(p1,p2))return;
  if(checkSame(p1,p3))return;
  if(checkSame(p1,p4))return;
  if(checkSame(p2,p3))return;
  if(checkSame(p2,p4))return;
  if(checkSame(p3,p4))return;
  //...checks mass
  double mass = (p1.momentum().p() + p2.momentum().p() + p3.momentum().p() + p4.momentum().p()).mag();
  double new_mass = ptype.mass();
  if(new_mass - massL <= mass && mass <= new_mass + massR){
    Particle cand(p1.momentum().p() + p2.momentum().p() + p3.momentum().p() + p4.momentum().p(), ptype);
    cand.relation().append(p1);
    cand.relation().append(p2);
    cand.relation().append(p3);
    cand.relation().append(p4);
    new_p.push_back(cand);
    return;
  }else{
    return;
  }
}

void 
combination_type01(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   vector<Particle> &p3,
		   vector<Particle> &p4)
{
  for(vector<Particle>::iterator i = p1.begin();
      i != p1.end(); ++i){
    for(vector<Particle>::iterator j = p2.begin();
        j != p2.end(); ++j){
      for(vector<Particle>::iterator k = p3.begin();
	  k != p3.end(); ++k){
	for(vector<Particle>::iterator l = p4.begin();
	    l != p4.end(); ++l){
	  combination_nocut(new_p, ptype, *i, *j, *k, *l);
	}
      }
    }
  }
}

void 
combination_type01(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   vector<Particle> &p3,
		   vector<Particle> &p4,
		   const double &massL,
		   const double &massR)
{
  for(vector<Particle>::iterator i = p1.begin();
      i != p1.end(); ++i){
    for(vector<Particle>::iterator j = p2.begin();
        j != p2.end(); ++j){
      for(vector<Particle>::iterator k = p3.begin();
	  k != p3.end(); ++k){
	for(vector<Particle>::iterator l = p4.begin();
	    l != p4.end(); ++l){
	  combination_cut(new_p, ptype, *i, *j, *k, *l, massL, massR);
	}
      }
    }
  }
}

void 
combination_type02(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   vector<Particle> &p3,
		   vector<Particle> &p4)
{
  if(p1.size() < 2)return;
  for(vector<Particle>::iterator i = p1.begin();
      i != p1.end()-1; ++i){
    for(vector<Particle>::iterator j = i+1;
        j != p1.end(); ++j){
      for(vector<Particle>::iterator k = p3.begin();
	  k != p3.end(); ++k){
	for(vector<Particle>::iterator l = p4.begin();
	    l != p4.end(); ++l){
	  combination_nocut(new_p, ptype, *i, *j, *k, *l);
	}
      }
    }
  }
}

void 
combination_type02(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   vector<Particle> &p3,
		   vector<Particle> &p4,
		   const double &massL,
		   const double &massR)
{
  if(p1.size() < 2)return;
  for(vector<Particle>::iterator i = p1.begin();
      i != p1.end()-1; ++i){
    for(vector<Particle>::iterator j = i+1;
        j != p1.end(); ++j){
      for(vector<Particle>::iterator k = p3.begin();
	  k != p3.end(); ++k){
	for(vector<Particle>::iterator l = p4.begin();
	    l != p4.end(); ++l){
	  combination_cut(new_p, ptype, *i, *j, *k, *l, massL, massR);
	}
      }
    }
  }
}

/************/
/* 5 Bodies */
/************/
void 
combination(vector<Particle> &new_p,
	    const Ptype &ptype,
	    vector<Particle> &p1,
	    vector<Particle> &p2,
	    vector<Particle> &p3,
	    vector<Particle> &p4,
	    vector<Particle> &p5)
{
  dout(Debugout::INFO,"combination") << "NOT Support in 5 bodies" << std::endl;
}

void 
combination(vector<Particle> &new_p,
	    const Ptype &ptype,
	    vector<Particle> &p1,
	    vector<Particle> &p2,
	    vector<Particle> &p3,
	    vector<Particle> &p4,
	    vector<Particle> &p5,
	    const double &mass_width)
{
  dout(Debugout::INFO,"combination") << "NOT Support in 5 bodies" << std::endl;
}

void 
combination(vector<Particle> &new_p,
	    const Ptype &ptype,
	    vector<Particle> &p1,
	    vector<Particle> &p2,
	    vector<Particle> &p3,
	    vector<Particle> &p4,
	    vector<Particle> &p5,
	    const double &massL,
	    const double &massR)
{
  dout(Debugout::INFO,"combination") << "NOT Support in 5 bodies" << std::endl;
}
}
