#ifndef EXKFIT_COMBINATION_FUNCTIONS_PC_H
#define EXKFIT_COMBINATION_FUNCTIONS_PC_H

#include "particle/Particle.h"

/// Utilies for Making Combination of Particle. 
///
/// Note 1: Only for particle constructed by MdstCharged. 
///
/// mass_width -->
/// nominal mass - mass_width < MASS < nominal mass + mass_width 
///
/// massL, massR -->
/// nominal mass - massL < MASS < nominal mass + massR 
///
/// nominal mass is defined by ptype.mass().
///

namespace exkfit {

  using std::Particle;

  bool checkSame(const Particle &p1,
		 const Particle &p2);
  
  /**********************/
  /* 2 Boby Combination */
  /**********************/
  void combination(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2);
  void combination(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   const double &mass_width);
  void combination(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   const double &massL,
		   const double &massR);
  /**********************/
  /* 3 Boby Combination */
  /**********************/
  void combination(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   vector<Particle> &p3);
  void combination(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   vector<Particle> &p3,
		   double mass_width);
  void combination(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   vector<Particle> &p3,
		   const double &massL,
		   const double &massR);
  /**********************/
  /* 4 Boby Combination */
  /**********************/
  void combination(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   vector<Particle> &p3,
		   vector<Particle> &p4);
  void combination(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   vector<Particle> &p3,
		   vector<Particle> &p4,
		   double mass_width);
  void combination(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   vector<Particle> &p3,
		   vector<Particle> &p4,
		   const double &massL,
		   const double &massR);
  /**********************/
  /* 5 Boby Combination */
  /**********************/
  void combination(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   vector<Particle> &p3,
		   vector<Particle> &p4,
		   vector<Particle> &p5);
  void combination(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   vector<Particle> &p3,
		   vector<Particle> &p4,
		   vector<Particle> &p5,
		   double mass_width);
  void combination(vector<Particle> &new_p,
		   const Ptype &ptype,
		   vector<Particle> &p1,
		   vector<Particle> &p2,
		   vector<Particle> &p3,
		   vector<Particle> &p4,
		   vector<Particle> &p5,
		   const double &massL,
		   const double &massR);
}

#endif /* EXKFIT_COMBINATION_FUNCTIONS_PC_H */
