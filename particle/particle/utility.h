//
// $Id: utility.h 10627 2008-09-08 07:42:48Z katayama $
//
#ifndef BPCU_FUNDAMENTAL_FUNCTIONS_PC_H
#define BPCU_FUNDAMENTAL_FUNCTIONS_PC_H

#include "belle.h"
#include "particle/Particle.h"
#include "kfitter/kvertexfitter.h"
#include "kfitter/kmassfitter.h"
#include "kfitter/kmassvertexfitter.h"
#include "kfitter/kmakemother.h"

#include <set>
#include <functional>
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class Helix;

/// Utilies for BELLE Particle Class. 
///

const double HER = 7.996;
const double LER = 3.5;
const double CROSS = 22.0; // mrad

const unsigned TYPE_E    = 0;
const unsigned TYPE_MU   = 1;
const unsigned TYPE_PI   = 2;
const unsigned TYPE_K    = 3;
const unsigned TYPE_P    = 4;
const unsigned TYPE_AUTO = 9;

/// makes Kaon and Pion from MdstCharged without cut.
/// 1 : with good_charged, 0 : w/o good_charged
void makeKPi(std::vector<Particle> &k_p, 
	     std::vector<Particle> &k_m, 
	     std::vector<Particle> &pi_p, 
	     std::vector<Particle> &pi_m,
	     const int = 1);

/// makes Leptons from MdstCharged without cut.
/// 1 : with good_charged, 0 : w/o good_charged
void makeLepton(std::vector<Particle> &e_p,
		std::vector<Particle> &e_m,
		std::vector<Particle> &mu_p,
		std::vector<Particle> &mu_m,
		const int = 1);

/// makes Gammas from MdstGamma without cut.
void makeGamma(std::vector<Particle> &gamma);

/// makes Gammas from MdstEcl without cut.
void makeEcl(std::vector<Particle> &ecl);

/// makes PI0s from MdstPi0 without cut.
void makePi0(std::vector<Particle> &pi0);

/// makes Kaon and Pion from GenHepEvt.
void makeGenHepKPi(std::vector<Particle> &k_p, 
		   std::vector<Particle> &k_m, 
		   std::vector<Particle> &pi_p, 
		   std::vector<Particle> &pi_m);

/// makes Leptons from GenHepEvt.
void makeGenHepLepton(std::vector<Particle> &e_p,
		      std::vector<Particle> &e_m,
		      std::vector<Particle> &mu_p,
		      std::vector<Particle> &mu_m);

/// If each plist element is not within l <= momentum(lab) <= r,
/// its element is removed from plist.
void withPCut(std::vector<Particle> &plist, 
	      const double l,
	      const double r);

/// If each plist element is not within min <= momentum(lab),
/// its element is removed from plist.
void withPCut(std::vector<Particle> &plist, 
	      const double min);

/// If each plist element is not within l <= momentum(Upsilon) <= r,
/// its element is removed from plist.
void withPSCut(std::vector<Particle> &plist, 
	       const double l,
	       const double r);
  
/// If each plist element is not within min <= momentum(Upsilon),
/// its element is removed from plist.
void withPSCut(std::vector<Particle> &plist, 
	       const double min);

/// If each plist element is not within l <= mass <= r,
/// its element is removed from plist.
void withMassCut(std::vector<Particle> &plist, 
		 const double l,
		 const double r);

/// If each plist element is not within l <= mass <= r,
/// its element is removed from plist and it is filled in nlist.
void withMassCut(std::vector<Particle> &plist, 
		 std::vector<Particle> &nlist,
		 const double l,
		 const double r);

/// If each plist element is within nominal-w <= mass <= nominal+w,
/// its element is filled in new_plist.
void massCut(const std::vector<Particle> &plist,
	     std::vector<Particle> &new_plist,
	     const double w);

/// If each plist element is not within l <= massdif <= r,
/// its element is removed from plist.
/// massdif = mass - mass(child).
void withMassDifCut(std::vector<Particle> &plist,
		    const double l,
		    const double r,
		    const unsigned child);

/// If p is within l <= massdif <= r,
/// returns 1. If not, returns 0.
/// massdif = mass - mass(child).
unsigned withMassDifCut(const Particle &p,
			const double l,
			const double r,
			const unsigned child);

/// If each plist element is within nominal-w <= massdif <= nominal+w,
/// its element is filled in new_plist.
/// massdif = mass - mass(child).
void massDifCut(const std::vector<Particle> &plist,
		std::vector<Particle> &new_plist, 
		const double w, 
		const unsigned child);

/// Kaon ID Probability
/// "atc" is used in this function.
double kaonId(const Particle &p,
	      int accq0 = 0, int tofq0 = 1, int cdcq0 = 0,
	      int ids0 = 3, int idb0 = 2);


/// If each plist element is not within atc_pID.prob >= prob,
/// its element is removed from plist.
/// "atc" is used in this function.

void withKaonId(std::vector<Particle> &plist, 
		const double prob = 0.6,
		int accq0 = 0, int tofq0 = 1, int cdcq0 = 0,
		int ids0 = 3, int idb0 = 2);

/// If each plist element is not within atc_pID.prob < prob,
/// its element is removed from plist.
/// "atc" is used in this function.

void withPionId(std::vector<Particle> &plist,
		const double prob = 0.6,
		int accq0 = 0, int tofq0 = 1, int cdcq0 = 0,
		int ids0 = 3, int idb0 = 2);

/// If each plist element is not within muon ID >= th or Likelihood >= prob,
/// its element is removed from plist.
/// withMuonId  ... Muid_mdst Class(recommend)
///             ... whichCriteria != 2 --> th( = 2) is used
///             ... whichCriteria = 2 --> prob( = 0.95) is used
/// withMuonId  ... MDST_MuId
/// withMuodId2 ... MDST_KLM
void withMuId(std::vector<Particle> &plist, 
	      const unsigned th = 2,
	      const double prob = 0.95,
	      const unsigned whichCriteria = 2);
void withMuonId(std::vector<Particle> &plist, 
		const unsigned th = 1);
void withMuonId2(std::vector<Particle> &plist, 
		 const unsigned th = 1);
  
/// If each plist element is not within electron ID >= th,
/// its element is removed from plist.
/// withEId  ... eid(recommend?)
/// withEId2 ... MDST_ElId
void withEId(std::vector<Particle> &plist,
	     const double = 0.9);
void withEId2(std::vector<Particle> &plist,
	      const double = 0.9);

/// If each plist element is not associated with z-svd hits
/// whose number is equal to or larger than nSvdHit,
/// its element is removed from plist.
void withSVD(std::vector<Particle> &list, 
	     const unsigned nZSvdHit = 1);

/// If each plist element is not associated with rphi & z-svd hits
/// whose number is equal to or larger than nRSvdHit and nZSvdHit,
/// its element is removed from plist.
void withSVD2(std::vector<Particle> &list, 
	      const unsigned nRSvdHit = 1,
	      const unsigned nZSvdHit = 1);

/// If child of each plist element is not associated with z-svd hits
/// whose number is equal to or larger than nSvdHit,
/// its element is removed from plist.
void withSVD(std::vector<Particle> &list, 
	     const unsigned nZSvdHit,
	     const unsigned child);

/// If child of each plist element is not associated with rphi & z-svd hits
/// whose number is equal to or larger than nRSvdHit and nZSvdHit,
/// its element is removed from plist.
void withSVD2(std::vector<Particle> &list, 
	      const unsigned nRSvdHit,
	      const unsigned nZSvdHit,
	      const unsigned child);

/// sets GenHepEvt Info. to final state charged particles(K,PI,E,MU,Proton).
/// Note: Relation is 1(Hep) to Many(Rec).
void setGenHepInfoF(std::vector<Particle> &plist);

/// sets GenHepEvt Info. to Gamma.
/// Note: Relation is 1(Hep) to Many(Rec).
void setGenHepInfoG(std::vector<Particle> &list);

/// sets GenHepEvt Info. to final state charged particles(K,PI,E,MU,Proton).
/// Note 1: Creteria = First, Max. # of SVD Hits --> Second, Min delta P.
///         delta P = |rec Momentum - hep Momentum|
/// Note 2: Called after setGenHepEvtInfoF(plist).
/// Note 3: Relation is 1(Hep) to 1(Rec).

void setUniqueGenHepInfoFBySvdAndDeltaP(std::vector<Particle> &plist);

/// sets GenHepEvt Info. to final state charged particles(K,PI,E,MU,Proton) and Gamma(at least?)
/// Note 1: Creteria = Min delta P.
///         delta P = |rec Momentum - hep Momentum|
/// Note 2: Called after setGenHepEvtInfoF(plist).
/// Note 3: Relation is 1(Hep) to 1(Rec).
void setUniqueGenHepInfoByDeltaP(std::vector<Particle> &list);

/// sets GenHepEvt Info. to reconstructed particles.
/// Note 1: To sets this info, uses hep info of child particles. 
/// Note 2: Not support resonauce modes.
unsigned setGenHepInfoR(Particle &);
void setGenHepInfoR(std::vector<Particle> &list);

/// sets GenHepEvt Info. to reconstructed particles.
/// Note 1: To sets this info, uses hep info of child particles. 
/// Note 2: Support resonauce modes. ... maybe ... it is test version
void setGenHepInfoR2(std::vector<Particle> &list);

/// adds tracks of plist to kvertexfitter.
void addTrack2fit(kvertexfitter &kv, const std::vector<Particle> &plist);

/// adds a track to kvertexfitter.
/// Note 1: Only done "add" command. 
///         Other infomation should be done by ourselves.
/// Note 2: Only non-correlation tracks.
void addTrack2fit(kvertexfitter &kv, const Particle &p);

/// adds a beam profile to kvertexfitter.
/// Note: Only non-correlation tracks.
void addBeam2fit(kvertexfitter &kv, 
		 const HepPoint3D &beam, const HepSymMatrix &errBeam);

/// added by T.H 2006/05/09
/// adds a default beam profile to kvertexfitter.
/// Note: Only non-correlation tracks.
void addBeam2fit(kvertexfitter &kv); 

/// added by T.H 2006/05/09
/// adds a ip tube to kvertexfitter.
void addTube2fit(kvertexfitter &kv);

/// adds tracks of plist to kmassvertexfitter.
void addTrack2fit(kmassvertexfitter &kmv, const std::vector<Particle> &plist);

/// adds a track to kmassvertexfitter.
/// Note 1: Only done "add" command. 
///         Other infomation should be done by ourselves.
/// Note 2: Only non-correlation tracks.
void addTrack2fit(kmassvertexfitter &kmv, const Particle &p);

/// adds tracks of plist to kmassfitter.
void addTrack2fit(kmassfitter &km, const std::vector<Particle> &plist);

/// adds a track to kmassfitter.
/// Note 1: Only done "add" command. 
///         Other infomation should be done by ourselves.
/// Note 2: Only non-correlation tracks.
void addTrack2fit(kmassfitter &km, const Particle &p);

/// removes p from plist.
/// returns # of removed elements.
unsigned removeParticle(std::vector<Particle> &plist, const Particle &p);

/// returns momentum(Upsilon).
HepLorentzVector pStar(const Particle &p,
		       const double e = HER, 
		       const double pos = LER,
		       const double a = CROSS);
HepLorentzVector pStar(const Particle &p,
		       const HepLorentzVector &el, 
		       const HepLorentzVector &po);

/// returns momentum(Upsilon).
HepLorentzVector pStar(HepLorentzVector p,
		       const double e = HER, 
		       const double pos = LER,
		       const double a = CROSS);
HepLorentzVector pStar(HepLorentzVector p,
		       const HepLorentzVector &el, 
		       const HepLorentzVector &po);

/// returns Thrust.
Hep3Vector calcuThrust(const std::vector<Particle> &plist, 
		       const double e = HER, 
		       const double p = LER,
		       const double a = CROSS);
Hep3Vector calcuThrust(const std::vector<Particle> &plist, 
		       const HepLorentzVector &el, 
		       const HepLorentzVector &po);

/// foxWolfram[3] is FoxWolfram parameters.
void calcuFoxWolfram(const std::vector<Particle> &plist, 
		     double *foxWolfram,
		     const double e = HER, 
		     const double p = LER,
		     const double a = CROSS);
void calcuFoxWolfram(const std::vector<Particle> &plist, 
		     double *foxWolfram,
		     const HepLorentzVector &el,
		     const HepLorentzVector &po);

/// new_plist is copied from plist deeply.
void deepCopy(std::vector<Particle> &plist, 
	      std::vector<Particle> &new_plist);

/// deletes plist which is made using deepCopy().
void deleteDeepCopiedObjects(std::vector<Particle> &plist);

/// returns mass with beam energy constraint.
double beamEnergyConstraint(const Particle &b, 
			    const double e = HER,
			    const double p = LER,
			    const double a = CROSS);
double beamEnergyConstraint(const Particle &b,
			    const HepLorentzVector &el, 
			    const HepLorentzVector &po);

/// makes mother particle from kfitter information.
unsigned makeMother(kvertexfitter &kv,
		    Particle &mother);

unsigned makeMother(kmassvertexfitter &kmv,
		    Particle &mother);

unsigned makeMother(kmassfitter &km,
		    Particle &mother);

unsigned makeMother(kmakemother &kmm,
		    kvertexfitter &kv,
		    Particle &mother,
		    const unsigned notMake = 0);

unsigned makeMother(kmakemother &kmm,
		    kmassvertexfitter &kmv,
		    Particle &mother,
		    const unsigned notMake = 0);

unsigned makeMother(kmakemother &kmm,
		    kmassfitter &km,
		    Particle &mother,
		    const unsigned notMake = 0);

/// erases Particle list.
void eraseVector(std::vector<Particle> &plist);

/// momentum and position --> helix at pivot(0,0,0)
/// Note: Only charged tracks
/// Note: If you use helix with pivot != (0,0,0),
///       please use Helix Class after the function.
bool px2helix(const Hep3Vector &p,
	      const HepPoint3D &x, 	 
	      const double charge,
	      HepVector &a,
	      const double alpha = 222.376063);

bool px2helix(const Hep3Vector &p,
	      const HepPoint3D &x,          
	      const HepSymMatrix &dpx,
	      const double charge,
	      HepVector &a,
	      HepSymMatrix &da,
	      const double alpha = 222.376063);

/// returns Ecm.
double calEcm(const double elec = HER, 
	      const double posi = LER,
	      const double angle = CROSS);

/// returns "HadronA" flag.
bool isHadronA(void);

/// returns "HadronC" flag.
bool isHadronC(void);

/// sets GenHepEvt Info. to reconstructed particles.
/// For resonance decay particles.
/// This is the same with "setGenHepInfoR2".....
/// Maybe this is faster than "setGenHepInfoR2".
void setGenHepInfoT(Particle &p);
void setGenHepInfoT(std::vector<Particle> &plist);

/// sets GenHepEvt Info. to pi0 and its daughter two-gammas.
void setGenHepInfoP(Particle &p);
void setGenHepInfoP(std::vector<Particle> &plist);

/// sets error matrix of gamma.
/// p          = gamma
/// gVertex    = decay point of gamma.
/// errGVertex = its error matrix.
void setGammaError(Particle &p,
		   const HepPoint3D &gVertex,
		   const HepSymMatrix &errGVertex);

/// sets error matrix of two-gammas of pi0.
/// p          = pi0
/// gVertex    = decay point of gamma.
/// errGVertex = its error matrix.
void setGammasError(Particle &p,
		    const HepPoint3D &gVertex,
		    const HepSymMatrix &errGVertex);

/// fits particle with simple mass-constraint.
/// Particularly pi0 after "setGammasError".
void doMassFit(Particle &p);
void doMassFit(std::vector<Particle> &plist);

/// gets expNo runNo evtNo(16BitMask) McFlag
void getEventInfo(int &expNo,
		  int &runNo,
		  int &evtNo,
		  bool &McFlag);

/// gets SVD hit number.
unsigned getNRSvd(const Particle &p, const unsigned id=TYPE_AUTO);
unsigned getNZSvd(const Particle &p, const unsigned id=TYPE_AUTO);
// for V particles (id = 0:plus, 1:minus)
unsigned getNRSvdVee(const Particle &p, const unsigned id);
unsigned getNZSvdVee(const Particle &p, const unsigned id);
unsigned getNRSvdVeeP(const Particle &p);
unsigned getNRSvdVeeM(const Particle &p);
unsigned getNZSvdVeeP(const Particle &p);
unsigned getNZSvdVeeM(const Particle &p);

//gets SVD hit pattern
unsigned getHitSvd(const Particle &p, const unsigned id=TYPE_AUTO);
// for V particles (id = 0:plus, 1:minus)
unsigned getHitSvdVee(const Particle &p, const unsigned id);
unsigned getHitSvdVeeP(const Particle &p);
unsigned getHitSvdVeeM(const Particle &p);

/// gets Helix Class from Mdst Charged Bank.
Helix calMdstChargedHelix(const Particle &p, const unsigned id=TYPE_AUTO);

/// sets GenHepEvt Info. to Ks and its daughter two-gammas.
void setGenHepInfoKs(Particle &p);
void setGenHepInfoKs(std::vector<Particle> &plist);

/// make Ks from MdstVee2(type=2) or MdstVee(type!=2)
void makeKs(std::vector<Particle> &ks, unsigned type = 2);

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* BPCU_FUNDAMENTAL_FUNCTIONS_PC_H */
