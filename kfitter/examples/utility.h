#ifndef EXKFIT_FUNDAMENTAL_FUNCTIONS_PC_H
#define EXKFIT_FUNDAMENTAL_FUNCTIONS_PC_H

#include "particle/Particle.h"
#include "kfitter/kvertexfitter.h"
#include "kfitter/kmassvertexfitter.h"

namespace exkfit {
  
  const double HER = 7.996;
  const double LER = 3.5;
  
  /// makes Kaon and Pion from MdstCharged without cut.
  void makeKPi(vector<Particle> &k_p, 
	       vector<Particle> &k_m, 
	       vector<Particle> &pi_p, 
	       vector<Particle> &pi_m);
  
  /// makes Leptons from MdstCharged without cut.
  void makeLepton(vector<Particle> &e_p,
		  vector<Particle> &e_m,
		  vector<Particle> &mu_p,
		  vector<Particle> &mu_m);
  
  /// makes Gammas from MdstGamma without cut.
  void makeGamma(vector<Particle> &gamma);
  
  /// makes Kaon and Pion from GenHepEvt.
  void makeGenHepKPi(vector<Particle> &k_p, 
		     vector<Particle> &k_m, 
		     vector<Particle> &pi_p, 
		     vector<Particle> &pi_m);
  
  /// makes Leptons from GenHepEvt.
  void makeGenHepLepton(vector<Particle> &e_p,
			vector<Particle> &e_m,
			vector<Particle> &mu_p,
			vector<Particle> &mu_m);
  
  /// If each plist element is not within l <= momentum(lab) <= r,
  /// its element is removed from plist.
  void withPCut(vector<Particle> &plist, 
		const double &l,
		const double &r);
  
  /// If each plist element is not within min <= momentum(lab),
  /// its element is removed from plist.
  void withPCut(vector<Particle> &plist, 
		const double &min);
  
  /// If each plist element is not within l <= momentum(Upsilon) <= r,
  /// its element is removed from plist.
  void withPSCut(vector<Particle> &plist, 
		 const double &l,
		 const double &r);
  
  /// If each plist element is not within min <= momentum(Upsilon),
  /// its element is removed from plist.
  void withPSCut(vector<Particle> &plist, 
		 const double &min);
  
  /// If each plist element is not within l <= mass <= r,
  /// its element is removed from plist.
  void withMassCut(vector<Particle> &plist, 
		   const double &l,
		   const double &r);
  
  /// If each plist element is not within l <= mass <= r,
  /// its element is removed from plist and it is filled in nlist.
  void withMassCut(vector<Particle> &plist, 
		   vector<Particle> &nlist,
		   const double &l,
		   const double &r);
  
  /// If each plist element is within nominal-w <= mass <= nominal+w,
  /// its element is filled in new_plist.
  void massCut(const vector<Particle> &plist,
	       const vector<Particle> &new_plist,
	       const double &w);
  
  /// If each plist element is not within l <= massdif <= r,
  /// its element is removed from plist.
  /// massdif = mass - mass(child).
  void withMassDifCut(vector<Particle> &plist,
		      const double &l,
		      const double &r,
		      const unsigned &child);
  
  /// If p is within l <= massdif <= r,
  /// returns 1. If not, returns 0.
  /// massdif = mass - mass(child).
  unsigned withMassDifCut(Particle &p,
			  const double &l,
			  const double &r,
			  const unsigned &child);
  
  /// If each plist element is within nominal-w <= massdif <= nominal+w,
  /// its element is filled in new_plist.
  /// massdif = mass - mass(child).
  void massDifCut(vector<Particle> &plist,
		  vector<Particle> &new_plist, 
		  const double &w, 
		  const unsigned &child);
  
  /// If each plist element is not within muon ID >= th,
  /// its element is removed from plist.
  void withMuonId(vector<Particle> &plist, 
		  const unsigned &th = 1);
  
  /// If each plist element is not within electron ID >= th,
  /// its element is removed from plist.
  /// Recommened Value : 199901xx ... 0.998
  ///                  : 199904xx ... 0.9
  void withEId(vector<Particle> &plist,
	       const double& = 0.9);
  
  /// If each plist element is not associated with svd hits
  /// whose number is equal to or larger than nSvdHit,
  /// its element is removed from plist.
  void withSVD(vector<Particle> &list, 
	       const unsigned &nSvdHit);
  
  /// If child of each plist element is not associated with svd hits
  /// whose number is equal to or larger than nSvdHit,
  /// its element is removed from plist.
  void withSVD(vector<Particle> &list, 
	       const unsigned &nSvdHit,
	       const unsigned &child);
  
  /// sets GenHepEvt Info. to final state charged particles(K,PI,E,MU,Proton).
  /// Note: Relation is 1(Hep) to Many(Rec).
  void setGenHepInfoF(vector<Particle> &plist);
  
  /// sets GenHepEvt Info. to final state charged particles(K,PI,E,MU,Proton).
  /// Note 1: Creteria = First, Max. # of SVD Hits --> Second, Min delta P.
  ///         delta P = |rec Momentum - hep Momentum|
  /// Note 2: Called after setGenHepEvtInfoF(plist).
  /// Note 3: Relation is 1(Hep) to 1(Rec).
  void setUniqueGenHepInfoFBySvdAndDeltaP(vector<Particle> &plist);
  
  /// sets GenHepEvt Info. to reconstructed particles.
  /// Note 1: To sets this info, uses hep info of child particles. 
  /// Note 2: Not support resonauce modes.
  void setGenHepInfoR(vector<Particle> &list);
  
  /// adds tracks of plist to kvertexfitter.
  void addTrack2fit(kvertexfitter &kv, vector<Particle> &plist);
  
  /// adds a track to kvertexfitter.
  /// Note 1: Only done "add" command. 
  ///         Other infomation should be done by ourselves.
  /// Note 2: Only non-correlation tracks.
  void addTrack2fit(kvertexfitter &kv, Particle &p);
  
  /// adds tracks of plist to kmassvertexfitter.
  void addTrack2fit(kmassvertexfitter &kmv, vector<Particle> &plist);
  
  /// adds a track to kmassvertexfitter.
  /// Note 1: Only done "add" command. 
  ///         Other infomation should be done by ourselves.
  /// Note 2: Only non-correlation tracks.
  void addTrack2fit(kmassvertexfitter &kmv, Particle &p);
  
  /// removes p from plist.
  /// returns # of removed elements.
  unsigned removeParticle(vector<Particle> &plist, const Particle &p);
  
  /// returns momentum(Upsilon).
  HepLorentzVector pStar(Particle &p,
			 double e = HER, 
			 double p = LER);
  
  /// returns momentum(Upsilon).
  HepLorentzVector pStar(HepLorentzVector p,
			 double e = HER, 
			 double p = LER);
  
  /// returns Thrust.
  Hep3Vector calcuThrust(vector<Particle> &plist, 
			 double e = HER, 
			 double p = LER);
  
  /// foxWolfram[3] is FoxWolfram parameters.
  void calcuFoxWolfram(vector<Particle> &plist, 
		       double *foxWolfram,
		       double e = HER, 
		       double p = LER);

  /// new_plist is copied from plist deeply.
  void deepCopy(vector<Particle> &plist, 
		vector<Particle> &new_plist);
  
  /// deletes plist which is made using deepCopy().
  void deleteDeepCopiedObjects(vector<Particle> &plist);
  
  /// returns mass with beam energy constraint.
  double beamEnergyConstraint(const Particle &b, 
			      double e = HER,
			      double p = LER);
  
  unsigned makeMother(kvertexfitter &kv,
		      Particle &mother);
  
  unsigned makeMother(kmassvertexfitter &kmv,
		      Particle &mother);

  /// erases Particle list.
  void eraseVector(vector<Particle> &plist);
}

#endif /* EXKFIT_FUNDAMENTAL_FUNCTIONS_PC_H */
