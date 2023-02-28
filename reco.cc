#include "particle/Particle.h"
#include "particle/PID.h"
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"

#include "belle.h"
#include MDST_H
#include HEPEVT_H

#include "reco.h"
#include "helix/Helix.h"
#include "ddk_skim/findKs.h"

#include "particle/combination.h"
#include "particle/utility.h"
#include "tables/belletdf.h"
#include "tables/brecon.h"
#include "tables/fullrecon.h"
#include "tables/evtcls.h"
#include "mdst/mdst.h"
#include "tables/evtvtx.h"
#include "fullrecon/frec_util.h"
#include "benergy/BeamEnergy.h"
#include "kfitter/kvertexfitter.h"
#include "kfitter/kmakemother.h"

#include "belleutil/debugout.h"

#include "userinfo.h"
#include "ip/IpProfile.h"

//#if defined(BELLE_NAMESPACE)
namespace Belle {
//#endif

    //***********************************************************
std::string pclType(int id) {
    std::string sid;
    if      (id ==   0) sid = "--";
    else if (id ==   1) sid = "d";
    else if (id ==  -1) sid = "anti-d";
    else if (id ==   2) sid = "u";
    else if (id ==  -2) sid = "anti-u";
    else if (id ==   3) sid = "s";
    else if (id ==  -3) sid = "anti-s";
    else if (id ==   4) sid = "c";
    else if (id ==  -4) sid = "anti-c";
    else if (id ==  11) sid = "el+";
    else if (id == -11) sid = "el-";
    else if (id ==  12) sid = "nu_e";
    else if (id == -12) sid = "anti-nu_e";
    else if (id ==  13) sid = "mu+";
    else if (id == -13) sid = "mu-";
    else if (id ==  14) sid = "nu_mu";
    else if (id == -14) sid = "anti-nu_mu";
    else if (id ==  15) sid = "tau+";
    else if (id == -15) sid = "tau-";
    else if (id ==  21) sid = "g";
    else if (id ==  22) sid = "gamma";
    else if (id ==  23) sid = "Z0";
    else if (id ==  24) sid = "W+";
    else if (id == -24) sid = "W-";
    else if (id ==  25) sid = "Higgs0";
    else if (id ==  32) sid = "Z'0";
    else if (id ==  33) sid = "Z''0";
    else if (id ==  34) sid = "W'+";
    else if (id == -34) sid = "W'-";
    else if (id ==  35) sid = "Higgs'0";
    else if (id ==  36) sid = "A0";
    else if (id ==  37) sid = "Higgs+";
    else if (id == -37) sid = "Higgs-";
    else if (id == 111) sid = "pi0";
    else if (id == 113) sid = "RHO0";
    else if (id == 130) sid = "K_L0";
    else if (id == 211) sid = "pi+";
    else if (id ==-211) sid = "pi-";
    else if (id == 213) sid = "RHO+";
    else if (id ==-213) sid = "RHO-";
    else if (id == 221) sid = "eta";
    else if (id == 223) sid = "omega";
    else if (id == 310) sid = "K_S0";
    else if (id ==-311) sid = "anti-K0";
    else if (id == 311) sid = "K0";
    else if (id == 313) sid = "K*0";
    else if (id ==-313) sid = "anti-K*0";
    else if (id == 321) sid = "K+";
    else if (id ==-321) sid = "K-";
    else if (id == 331) sid = "eta'";
    else if (id == 323) sid = "K*+";
    else if (id ==-323) sid = "K*-";
    else if (id == 333) sid = "phi";
    else if (id ==-411) sid = "D-";
    else if (id == 411) sid = "D+";
    else if (id ==-413) sid = "D*-";
    else if (id == 413) sid = "D*+";
    else if (id ==-415) sid = "D_2*-";
    else if (id == 415) sid = "D_2*+";
    else if (id ==-421) sid = "anti-D0";
    else if (id == 421) sid = "D0";
    else if (id ==-431) sid = "DS-";
    else if (id == 431) sid = "DS+";
    else if (id ==-423) sid = "anti-D*0";
    else if (id == 423) sid = "D*0";
    else if (id ==-425) sid = "anti-D_2*0";
    else if (id == 425) sid = "D_2*0";
    else if (id ==-433) sid = "D*S-";
    else if (id == 433) sid = "D*S+";
    else if (id == 441) sid = "eta_c";
    else if (id == 443) sid = "J/psi";
    else if (id ==-511) sid = "anti-B0";
    else if (id == 511) sid = "B0";
    else if (id ==-521) sid = "B-";
    else if (id == 521) sid = "B+";
    else if (id == 553) sid = "Upsilon";
    else if (id == 2112) sid = "n0";
    else if (id ==-2112) sid = "anti-n0";
    else if (id == 2212) sid = "p+";
    else if (id ==-2212) sid = "anti-p-";
    else if (id == 2224) sid = "Delta++";
    else if (id ==-2224) sid = "anti-Delta--";
    else if (id == 3122) sid = "Lam0";
    else if (id ==-3122) sid = "Lam0bar";
    else if (id ==-4122) sid = "LamC-";
    else if (id == 4122) sid = "LamC+";
    else if (id ==-4214) sid = "anti-Sigma_c*-";
    else if (id == 4214) sid = "Sigma_c*+";
    else if (id == 10022) sid = "vpho";
    else if (id == 10221) sid = "f_0";
    else if (id ==-10413) sid = "D_1-";
    else if (id == 10413) sid = "D_1+";
    else if (id ==-10423) sid = "anti-D_10";
    else if (id == 10423) sid = "D_10";
    else if (id ==-10431) sid = "D_s0*-";
    else if (id == 10431) sid = "D_s0*+";
    else if (id ==-10433) sid = "D_s1-";
    else if (id == 10433) sid = "D_s1+";
    else if (id ==-20213) sid = "a_1-";
    else if (id == 20213) sid = "a_1+";
    else if (id ==-20413) sid = "D'_1-";
    else if (id == 20413) sid = "D'_1+";
    else if (id ==-20423) sid = "anti-D'_10";
    else if (id == 20423) sid = "D'_10";
    else if (id ==-20433) sid = "D'_s1-";
    else if (id == 20433) sid = "D'_s1+";
    else if (id ==  300553) sid = "Upsilon(4S)";
    else if (id == 9000553) sid = "Upsilon(5S)";
    else sid = "pcl";
    return sid;
}

//======
// Set proper error matrix for gamma. Note errCart length is 4.
// Copied from FindGamma in icpv_skim package.
// From Tagir Aushev for pi0 make
//======
void setGammaError(Particle& gamma) {

  HepSymMatrix errEcl(3, 0); // 3x3 initialize to zero
  errEcl[0][0] = gamma.mdstGamma().ecl().error(0); // Energy
  errEcl[1][0] = gamma.mdstGamma().ecl().error(1);
  errEcl[1][1] = gamma.mdstGamma().ecl().error(2); // Phi
  errEcl[2][0] = gamma.mdstGamma().ecl().error(3);
  errEcl[2][1] = gamma.mdstGamma().ecl().error(4);
  errEcl[2][2] = gamma.mdstGamma().ecl().error(5); // Theta

  double  cp = cos(gamma.mdstGamma().ecl().phi());
  double  sp = sin(gamma.mdstGamma().ecl().phi());
  double  ct = cos(gamma.mdstGamma().ecl().theta());
  double  st = sin(gamma.mdstGamma().ecl().theta());
  double  E  =     gamma.mdstGamma().ecl().energy();

  HepMatrix jacobian(4, 3, 0);
  jacobian[0][0] =       cp * st;
  jacobian[0][1] =  -E * sp * st;
  jacobian[0][2] =   E * cp * ct;
  jacobian[1][0] =       sp * st;
  jacobian[1][1] =   E * cp * st;
  jacobian[1][2] =   E * sp * ct;
  jacobian[2][0] =            ct;
  jacobian[2][1] =           0.0;
  jacobian[2][2] =  -E      * st;
  jacobian[3][0] =           1.0;
  jacobian[3][1] =           0.0;
  jacobian[3][2] =           0.0;

  HepSymMatrix errCart = errEcl.similarity(jacobian);

  const HepPoint3D origin; // Default constructor sets (0,0,0).
  static double largeError = 1.0; // 1.0*1.0cm^2.
  HepSymMatrix dx(3, 0);
  dx[0][0] = largeError;
  dx[1][1] = largeError;
  dx[2][2] = largeError;
  // Convert Mdst_gamma into Particle object.

  gamma.momentum().momentum(gamma.p(), errCart);
  gamma.momentum().position(origin, dx);
  return;
}

void setGammaError(std::vector<Particle>& p) {
  for(std::vector<Particle>::iterator i = p.begin(); i != p.end(); ++i)
    setGammaError(*i);
}

// ************ Set error matrix (dpx) for pi0 ******************/
void setPi0Error(Particle& particle) {
    if (particle.nChildren() != 2) return;
    HepSymMatrix tmpErr(3, 0);
    HepPoint3D org(0., 0., 0.);
    if (!particle.child(0).mdstGamma() || !particle.child(1).mdstGamma()) return;
    for (unsigned i = 0; i < 2; ++i)
        setGammaError(particle.child(i));
    // setGammaError( p.child(i), org, tmpErr);

    // kmassfitter kmv;
    kmassvertexfitter kmv;
    if (useBF) kmv.magneticField(BF);
    kmv.invariantMass(particle.pType().mass());
    for (unsigned i = 0; i < 2; ++i)
        addTrack2fit(kmv, particle.child(i));
    // kmv.vertex(org);
    // kmv.atDecayPoint();
    int err = kmv.fit();
    if (!err) {
        kmakemother kmm2;
        if (useBF) kmm2.magneticField(BF);
        makeMother(kmm2, kmv, particle, 0);
        // makeMother(kmv, particle);
        particle.momentum().vertex(kmv.vertex(), kmv.errVertex());
    }
}

void setPi0Error(std::vector<Particle>& p_list) {
  for (std::vector<Particle>::iterator itr = p_list.begin(); itr != p_list.end(); ++itr)
      setPi0Error(*itr);
}
// ***********************************************************
int IDhep(Particle& particle) {
    if(!particle.genHepevt()) return 0;
    return particle.genHepevt().idhep();
}
// ***********************************************************
bool isLikeTrk(int lund, const std::string& chrgType = "Charged_pi0") {
    int ln = abs(lund);
    bool isTrk = (ln == 11) || (ln == 13) || (ln == 211) ||(ln == 321) || (ln == 11) || (ln == 2212); // el, mu, pi, K, p
    if ( chrgType != "onlyCharged" )
        isTrk = isTrk || (ln == 111) || (ln == 22) ;  // pi0 & gamma
    return isTrk;
}
// ***********************************************************
void printPclDebug(Particle& p, const std::string& comment = "", const std::string& comment2 = "") {
    int expNo, runNo, evtNo;
    bool McFlag;
    getEventInfo( expNo, runNo, evtNo, McFlag); // utility.cc

    int lund = (int)p.lund();
    std::string trkType = pclType(lund);
    int ln = abs(lund);
    bool isTrk = isLikeTrk( lund );
    Hep3Vector pcl3D( p.px(), p.py(), p.pz() );
    printf("\n----- Particle [%s]  %9s %s ---\n", trkType.c_str(), comment.c_str(), comment2.c_str() );
    printf(" massRECO:%8.5f,   ", p.p().m());
    // if (!isTrk)
    // if ( p.userInfo() ) printUserInfo( p );
    printf("px,py,pz: (%8.5f,%8.5f,%8.5f)  ", p.px(), p.py(), p.pz());
    printf("  pT,phi,theta: [%6.4f,%6.4f,%6.4f] \n", pcl3D.perp(), pcl3D.phi(), pcl3D.theta());
    printf("-- Ptype info.  lund:%i,  name:%s,  massPDG:%8.5f,  cTau:%9.6f, charge:%3.1f, spin:%3.1f,  nDecay:%i, stable:%i \n", 
           (int)Ptype(lund).lund(), Ptype(lund).name().c_str(), Ptype(lund).mass(), Ptype(lund).cTau(), 
           Ptype(lund).charge(), Ptype(lund).spin(), Ptype(lund).nDecay(), (int)Ptype(lund).stable() );
    if (!isTrk) {
        double vx = p.momentum().decayVertex().x();
        double vy = p.momentum().decayVertex().y();
        double vz = p.momentum().decayVertex().z();
        printf("              vtx:[%8.5f,%8.5f,%8.5f] ", vx, vy, vz );
        // printf( "      distanceToIP:%8.5f \n", distanceToIP(p) );
        printf("Children ID: ");
        for (int jM =0 ; jM < p.nChildren(); ++jM) {
            Particle& Child = p.child(jM);
            int lundChild = (int)Child.lund() ;
            std::string childType = pclType(lundChild);
            bool isTrkChild = isLikeTrk( lundChild, "onlyCharged" );
            if (isTrkChild) printf("( %s : %i )  ", childType.c_str(), (int)Child.mdstCharged().get_ID() );
            else printf("( %s : -- )  ", childType.c_str() );
        }
        printf("\n");
    }
    if (McFlag) {
        int ID = p.genHepevt() ? p.genHepevt().get_ID() : -999;
        printf("-- MC info. ID:%3i  IDhep:%i,  type:%s \n\n", 
               ID, IDhep(p), pclType( IDhep(p) ).c_str() );
    }
}
// ***********************************************************
void printPclDebug(vector<Particle>& p_list, const std::string& comment = "", const std::string& comment2 = "") {
    size_t np = p_list.size();
    printf("\n----- List of Particles (%i) %9s %s ---\n", np, comment.c_str(), comment2.c_str());
    for (size_t i = 0; i < np; ++i) {
        printPclDebug( p_list[i], comment, comment);
    }
}
// ***********************************************************
void printUserInfo(Particle& p) {
    UserInfo& info = dynamic_cast<UserInfo&>(p.userInfo());
    printf("----- printUserInfo -----  msComb:%8.5f, msKvf:%8.5f, chisq:%8.3f, chisqKvf:%8.3f, cl:%8.6f, clKvf:%8.6f, \n              dist2IP:%8.5f, dist2IPmvf:%8.5f, useTube:%i, useKmvf:%i,  isAdoptCut:%i, wMass:%8.5f, maxChi2:%8.3f,  helicity:%6.4f \n", 
           info.msComb(), info.msKvf(), info.chisq(), info.chisqKvf(), info.cl(), info.clKvf(), info.dist2IP(), info.dist2IPKmvf(),
           (int)info.useTube(), (int)info.useKmvf(), (int)info.isAdoptCut(), info.wMass(), info.maxChi2(), info.helicity() );
}
// ***********************************************************
void dRdZ(const Mdst_charged& charged, int massHyp, const HepPoint3D& ip_position,
                                                            double dR, double dZ) {

    if (charged.trk()) {
        Mdst_trk& trk = charged.trk();
        if (trk.mhyp(massHyp)) {
            Mdst_trk_fit& trkFit = trk.mhyp(massHyp);

            HepPoint3D pivot;
            pivot.setX((double) trkFit.pivot(0));
            pivot.setY((double) trkFit.pivot(1));
            pivot.setZ((double) trkFit.pivot(2));

            HepVector a(5);
            a[0] = (double) trkFit.helix(0);
            a[1] = (double) trkFit.helix(1);
            a[2] = (double) trkFit.helix(2);
            a[3] = (double) trkFit.helix(3);
            a[4] = (double) trkFit.helix(4);

            Helix h(pivot, a);
            h.pivot(ip_position);
            dR = fabs(h.dr());
            dZ = fabs(h.dz());

            return;
        }
    }

    dR = 100.;
    dZ = 100.;
}
// ***********************************************************
void withdRdZcut(std::vector<Particle>& p_list, const HepPoint3D& ip_position, double dRcut, double dZcut) {
    for (size_t i = 0; i < p_list.size(); ++i) {
        int id = abs(p_list[i].pType().lund());

        int mhyp;
        if (id == 11) mhyp = 0;
        else if (id == 13) mhyp = 1;
        // else if (id == 211) mhyp = 2; // repeatition
        else if (id == 321) mhyp = 3;
        else if (id == 2212) mhyp = 4;
        else mhyp = 2;

        double dr, dz;
        dRdZ(p_list[i].mdstCharged(), mhyp, ip_position, dr, dz);

        if (dr > dRcut || dz > dZcut) {
            p_list.erase(p_list.begin() + i);
            --i;
        }
    }
}
// **********************************************************
bool withPi0MassGamGamCut(Particle& pi0, double M_min, double M_max) {
    Particle& g1 = pi0.child(0);
    Particle& g2 = pi0.child(1);
    double msPi0_gg = (g1.p() + g2.p()).m();
    if (msPi0_gg < M_min) return false;
    if (msPi0_gg > M_max) return false;
    return true;
}
// **********************************************************
bool withPi0GammPCut(Particle& p, double p_min) {
    Particle& g1 = p.child(0);
    Particle& g2 = p.child(1);
    if(p.child(0).ptot() < p_min) return false;
    if(p.child(1).ptot() < p_min) return false;
    return true;
}
// **********************************************************
void withPi0MassGamGamCut(std::vector<Particle>& pi0, double dM_max) {
    double massPi0_PDG = Ptype(111).mass();
    for (std::vector<Particle>::iterator i = pi0.begin(); i != pi0.end();) {
        Particle& g1 = i->child(0);
        Particle& g2 = i->child(1);
        double msPi0_gg = (g1.p() + g2.p()).m();
            
        if (abs(msPi0_gg - massPi0_PDG) > dM_max) pi0.erase(i);
        else
            ++i;
    }
}
/* **********************************************************
 * performs a neutral pion's (pi0) gammas' energies in order
 * to select better candidates. Operates with gammas that
 * belong to a pi0
*/
// ***********************************************************
void withPi0GammPCut(std::vector<Particle>& pi0, double p_min) {
    for (std::vector<Particle>::iterator i = pi0.begin(); i != pi0.end();)
        if (i->child(0).ptot() < p_min || i->child(1).ptot() < p_min)
            pi0.erase(i);
        else
            ++i;
}
// ***********************************************************
void createUserInfo(Particle& particle) {
    particle.userInfo(*(new UserInfo(particle)));
    UserInfo& info = dynamic_cast<UserInfo&>(particle.userInfo());
    int lund = (int)particle.lund();
    double cTau = Ptype(lund).cTau();
    bool useTube = false;
    double wMass = dM_Dgr;

    // !!! used for this reconstruction
    if ((abs(lund) > 500) && (abs(lund) < 600)) {                         // B0, Bc, Bs
        useTube = true;
        wMass = wB;
    }
    // !!! used for this reconstruction (only Phi0)
    else if ((lund == 310) || (abs(lund) == 3122) || (lund == 333)) { // K0s, Lam0, Phi0
        wMass = dM_V0;
    }
    // !!! used for this reconstruction
    else if (abs(lund)  == 313) {                                          // K*0
        wMass = dM_Ksr0;
    }
    else if ((abs(lund) == 113) || (abs(lund) == 213)) {                // RHO0,+.-
        wMass = dM_Rho;
    }
    // !!! used for this reconstruction
    else if ((abs(lund) == 413) || (abs(lund) == 423)) {                // D*+, D*0feve
        useTube = true; 
        wMass = wDst;
    }
    else if (abs(lund)  == 431) {                                          // DS+
        wMass = dM_Dss;
    }
    else if (abs(lund)  == 433) {                                          // D*S+
        useTube = true;
        wMass = dM_Dsst;
    }
    else if (abs(lund)  == 10431) {                                        // D**S+(D_sJ(2317))
        useTube = true;
        wMass = dM_2317;
    }
    else if (abs(lund)  == 443) {                                          // J/psi
        useTube = true; 
    }
    
    info.msComb(particle.p().m());
    info.maxChi2(-1.);
    info.wMass(wMass);
    info.useTube(useTube);
    info.useKmvf(cTau > 1.e-5);
    info.isAdoptCut(true);
    info.chisqKvf(-1.);
    info.chisqKmvf(-1.);
    info.chisqKvfNdf(-1.);
    info.chisqKmvfNdf(-1.);
    info.helicity(-1.);
}
// ***********************************************************
void createUserInfo(std::vector<Particle>& p_list) {
    std::vector<Particle>::iterator particle;
    for (particle = p_list.begin(); particle != p_list.end(); ++particle) {
        if (!&(*particle).userInfo()) createUserInfo(*particle);
    }
}
// ***********************************************************
void setBadVtx(Particle& particle) {
    HepPoint3D vtx(999., 999., 999.);
    HepSymMatrix errVtx(3, 0);
    particle.momentum().decayVertex(vtx, errVtx);
    UserInfo& info = dynamic_cast<UserInfo&>(particle.userInfo());
    info.isAdoptCut(false);
}
// ***********************************************************
double distanceToIP(Particle& particle) {
    const HepPoint3D& ip_position = IpProfile::position();
    Hep3Vector vtxIP(IpProfile::position().x(), IpProfile::position().y(), 0.0);
    Hep3Vector vtxD(particle.momentum().decayVertex().x(), particle.momentum().decayVertex().y(), 0.0);
    Hep3Vector vtx_D_IP = vtxD - vtxIP;
    return vtx_D_IP.mag();
}
// ***********************************************************
void makeRecursiveVertexFit(Particle& Mother, bool debugDump = false, bool useKmvf = false,
                   bool addBeam = false) {
    /*
     * Fitting DecayMother -> DecayChild1 + DecayChild2 + ... + trk1 + trk2 + ...
     * makeRecursiveVertexFit(Particle&, bool, bool, bool)
     */

    UserInfo& infoMother = dynamic_cast<UserInfo&>(Mother.userInfo());


    int lundMother = (int)Mother.lund();
    std::string motherType = pclType(lundMother);

    if (debugDump) {
        printf("\n ========  makeRecursiveVertexFit ==========   %s [%i] --> ",
            motherType.c_str(), Mother.nChildren());
        for (int jChild = 0; jChild < Mother.nChildren(); ++jChild) {
            std::string dghtType = pclType((int)Mother.child(jChild).lund());
            printf("%s ", dghtType.c_str());
        }
        printf("\n");
    }


    kvertexfitter kvfMother;
    if (useBF) kvfMother.magneticField(BF);
    int err = -2;
    for (int jChild = 0; jChild < Mother.nChildren(); ++jChild) {
        Particle& Child = Mother.child(jChild);
        int lundChild   = (int)Child.lund();
        bool isTrkChild = isLikeTrk(lundChild);
        if (!isTrkChild) {
            if (!&Child.userInfo()) createUserInfo(Child);
            UserInfo& infoChild = dynamic_cast<UserInfo&>(Child.userInfo());
            double chisqChild = infoChild.chisqKvf();
            if (chisqChild < 0.) {
                bool idUseTube = false;
                bool idAddMF   = true;
                makeRecursiveVertexFit(Child, debugDump);
            }
        }
        addTrack2fit(kvfMother, Child);
        if (addBeam && IpProfile::usable()) addBeam2fit(kvfMother);
    }

    if (infoMother.useTube()) addTube2fit(kvfMother);
    err = kvfMother.fit();

    // If OK, 0 is returned
    if (err) {
        setBadVtx(Mother);
        if (debugDump) {
            printf( "-----  bad vertexing!!! ---- \n\n" );
        }
        return;
    }

    makeMother(kvfMother, Mother);
    infoMother.msKvf(Mother.p().m());
    infoMother.chisq(kvfMother.chisq());
    infoMother.chisqKvf(kvfMother.chisq());
    infoMother.chisqKvfNdf( kvfMother.chisq() / kvfMother.dgf());
    infoMother.cl(kvfMother.cl());
    infoMother.clKvf(kvfMother.cl());
    infoMother.dist2IP(distanceToIP(Mother));

    for (int jChild = 0; jChild < Mother.nChildren(); ++jChild) {
        Particle& Child = Mother.child(jChild);
        int lundChild = (int)Child.lund();
        bool isTrkChild = isLikeTrk(lundChild);
        if (!isTrkChild) {
            if (!&Child.userInfo()) createUserInfo(Child);
            UserInfo& infoChild = dynamic_cast<UserInfo&>(Child.userInfo());
            double dx = Mother.momentum().decayVertex().x() - Child.momentum().decayVertex().x();
            double dy = Mother.momentum().decayVertex().y() - Child.momentum().decayVertex().y();
            double dz = Mother.momentum().decayVertex().z() - Child.momentum().decayVertex().z();
            double dist2Mother = sqrt(dx * dx + dy * dy + dz * dz);
            infoChild.dist2Mother(dist2Mother);
        }
    }
                
    // Checking for long-living particles as a criteria of mass-constraint fit
    double cTauMother = Ptype(lundMother).cTau();
    if (cTauMother > 1.e-5) useKmvf =  true;
    // making it unenabled for Ds2317. PDG id: 10431 of abs. units
    if (abs(lundMother) == 10431) useKmvf = false;

    infoMother.useKmvf(useKmvf);
        
    if (infoMother.useKmvf()) {
        // use additional mass-vertex fitter to correct position (?) and LV
        kmassvertexfitter kmvMother;
        if (useBF) kmvMother.magneticField(BF);
        kmvMother.invariantMass(Mother.pType().mass());
        for (int jChild = 0; jChild < Mother.nChildren(); ++jChild) addTrack2fit(kmvMother, Mother.child(jChild));
        err = kmvMother.fit();
        if (err) {
            setBadVtx(Mother);
            return;
        }
        makeMother(kmvMother, Mother);
        Mother.momentum().decayVertex(kmvMother.vertex(), kmvMother.errVertex());

        infoMother.clKmvf(kmvMother.cl());
        infoMother.dist2IPKmvf(distanceToIP(Mother));
        infoMother.msKmvf(Mother.p().m());
        infoMother.chisqKmvfNdf(kmvMother.chisq() / kmvMother.dgf());
        infoMother.chisqKmvf(kmvMother.chisq());
    }

    if (debugDump) {
        printPclDebug( Mother, "(Mother)" );
        printf("Mother  --- vtx:[%8.5f,%8.5f,%8.5f], err[00,11,22]:[%8.5f,%8.5f,%8.5f]\n",
            Mother.momentum().decayVertex().x(),
            Mother.momentum().decayVertex().y(),
            Mother.momentum().decayVertex().z(),
            sqrt(kvfMother.errVertex()[0][0]),
            sqrt(kvfMother.errVertex()[1][1]),
            sqrt(kvfMother.errVertex()[2][2])
        );
        std::string sind[5] = { "(pcl_1)", "(pcl_2)", "(pcl_3)", "(pcl_4)" , "(pcl_5)" };
        for (int jChild = 0; jChild<Mother.nChildren(); ++jChild) {
            printPclDebug( Mother.child(jChild), sind[jChild] );
        }
        printf("\n");
    }

}
// ***********************************************************
void makeRecursiveVertexFit(vector<Particle>& p_list, bool debugDump = false, bool useKmvf = false) {
    // fit DecayMother -> DecayChild1 + DecayChild2 + ... + trk1 + trk2 + ...
    for (size_t i = 0; i < p_list.size(); ++i) {
        makeRecursiveVertexFit(p_list[i], debugDump, useKmvf);
    }
}
// ***********************************************************
void printPi0(vector<Particle>& pi0, const std::string& comment = "") {
    printf("------  %s Pi0 (%i) -------\n", comment.c_str(), pi0.size());
    
    for (size_t iPi0 = 0; iPi0 < pi0.size(); ++iPi0) {
        Particle& p0 = pi0[iPi0];
        Particle& g1 = pi0[iPi0].child(0);
        Particle& g2 = pi0[iPi0].child(1);
        double msPi0_gg = (g1.p() + g2.p()).m();
        double psrPi0 = pStar(p0.p(), BeamEnergy::E_HER(), BeamEnergy::E_LER(), BeamEnergy::Cross_angle()).vect().mag();
        printf(" Pi0 (%i)  mass:%7.5f, [px,py,pz]:[%7.4f, %7.4f, %7.4f], p:%6.4f,  p_cm:%6.4f,  Eg(1,2): [%6.4f, %6.4f], ms_gg:%7.5f\n",
        iPi0, p0.mass(), p0.px(), p0.py(), p0.pz(), p0.ptot(), psrPi0, g1.ptot(), g2.ptot(), msPi0_gg);
    }
}
// ***********************************************************
void printTrkPID(vector<Particle>& trkList, const std::string& pType, const std::string& comment = "") {
    int ip, ibg1, ibg2;
    if ((pType == "pi+") || (pType == "pi-")) {
        ip   = 2;
        ibg1 = 3;
        ibg2 = 4;
    }
    else if ((pType == "K+") || (pType == "K-")) {
        ip   = 3;
        ibg1 = 2;
        ibg2 = 4;
    }
    else if ((pType == "p+") || (pType == "p-")) {
        ip   = 4;
        ibg1 = 2;
        ibg2 = 3;
    }
    else return;
    printf("--- %3s (%i) ---%6s  trkID(eid,bg1,bg2): ", pType.c_str(), trkList.size(), comment.c_str());
    for (int i = 0; i < trkList.size(); i++)
        printf("(%2i | %5.3f,%5.3f,%5.3f) ", (int) trkList[i].mdstCharged().get_ID(),
                (double) eid(trkList[i].mdstCharged()).prob(3, -1, 5),
                (double) atc_pid(3, 1, 5, ip, ibg1).prob(trkList[i].mdstCharged()),
                (double) atc_pid(3, 1, 5, ip, ibg2).prob(trkList[i].mdstCharged()));
    printf("\n");
}
// ***********************************************************
void printVtxDebug( std::vector<Particle>& plist, std::string comment = "", std::string comment2 = "") {
    size_t np = plist.size();
    printf("\n----- Vtx (%i) %9s %s ---\n", np, comment.c_str(), comment2.c_str() );
    for (size_t i = 0; i < np; ++i) {
        UserInfo& info = dynamic_cast<UserInfo&>(plist[i].userInfo());
        double vx = plist[i].momentum().decayVertex().x();
        double vy = plist[i].momentum().decayVertex().y();
        double vz = plist[i].momentum().decayVertex().z();
        printf("[ %i ]  mass,chisq,chisqKvf: (%8.5f,%7.3f,%7.3f)  ", 
            i, plist[i].p().m(), info.chisq(), info.chisqKvf());
        printf("vx,vy,vz: (%8.5f,%8.5f,%8.5f)  ", vx, vy, vz);
        printf("px,py,pz: (%8.5f,%8.5f,%8.5f)\n", 
            plist[i].px(), plist[i].py(), plist[i].pz());
    }
    printf("\n");
}

// ***********************************************************
void clearVectors() {
    // clear all vectors
    for (int i = 0; i < nTrk; ++i) trkV[i].clear();
    gammaV.clear();
    pi0.clear();
    phi0.clear();
    Ksr0.clear();
    Ksr0bar.clear();
    Dss_p.clear();
    Dss_m.clear();
    Dss_m_2317.clear();
    Dss_p_2317.clear();
    Bs0.clear();
    Bs0bar.clear();
    // BsStar0.clear();
    // BsStar0bar.clear();
    // Upsilon_5S.clear();
}
// ***********************************************************
double getHelicity(VectorL& mother, VectorL& dght, VectorL& grndDght) {
    Hep3Vector V1_boost = -dght.boostVector();
    VectorL mother_boosted = mother;
    VectorL grndDght_boosted = grndDght;
    mother_boosted.boost(V1_boost);
    grndDght_boosted.boost(V1_boost);
    double helic = -mother_boosted.vect().unit().dot(grndDght_boosted.vect().unit());

    if (debugHel) {
        printf("\n======================== getHelicity ====================\n");
        cout << "  VectorL   mother:" << mother << endl;
        cout << "  VectorL     dght:" << dght << endl;
        cout << "  VectorL grndDght:" << grndDght << endl;
        printf("         dght helic: %7.5f\n\n",helic);
    }

    
    return helic;
}
// ***********************************************************

double getHelicity(Particle& p, int indDough=0) {
    if (p.nChildren() == 0) return -1;
    
    VectorL mother   = p.p(); 
    VectorL dght     = p.child(indDough).p();
    VectorL grndDght = p.child(indDough).child(0).p();
    double helic = getHelicity(mother, dght, grndDght);

    if (debugHel) {
        printf("\n======================== getHelicity ====================");
        printPclDebug( p, "(Mother)" );
        printf("         dght helic: %7.5f\n\n",helic);
    }

    
    return helic;
}
// **********************************************************
void withKaonId(std::vector<Particle>& p_list, const double prob, int accq0, int tofq0, int cdcq0, int ids0, int idb0) {
    atc_pid kid(accq0, tofq0, cdcq0, ids0, idb0);
    for (size_t i = 0; i < p_list.size(); ++i) {
        Particle& P = p_list[i];
        double probParticle = kid.prob(&(P.mdstCharged()));

        if (P.mdstCharged() && probParticle >= prob) {

            if (!&P.userInfo()) createUserInfo(P);
            UserInfo& info = dynamic_cast<UserInfo&>(P.userInfo());
            info.isAdoptCut(true);
        }
        else {
            p_list.erase(p_list.begin() + i);
            --i;
        }
  }
}

void withPionId(std::vector<Particle>& p_list,
                const double prob, int accq0, int tofq0, int cdcq0,
                int ids0, int idb0) {
    /* ************* THE FUNCTION'S ANNOTATION **************
     * Combine informations from ACC, TOF and CDC (dE/dx)
     * in order to separate e+-, mu+-, K+-, pi+-, p+-
     * atc_pid returns a likelihood ratio which is called
     * probParticle that compares two particles i and j:
     * Prob(i : j) = L(i) / (L(i) + L(j))
     * For example, Prob(K : pi) tends to be 1 if the
     * charged paritcle is K-like,
     * and tends to be 0 if it is pi-like.
     * Apparently, Prob(K : pi) = 1 - Prob(pi : K).
     * Prob(K : pi) is NOT a probability to be a kaon.
     * accq0, tofq0, cdcq0: flags for each detector
     * accq0 = 3, tofq0 = 1, cdcq0 = 5 are the
     * default and recommended values;
     * ids0, idb0 -- signal and background species.
     * 0 for e,
     * 1 for mu,
     * 2 for pi,
     * 3 for K and
     * 4 for p.
     * By default, ids0 = 3 and idb0 = 2 (K agaist pi).
     * */

    atc_pid kid(accq0, tofq0, cdcq0, ids0, idb0);

    for (size_t i = 0; i < p_list.size(); ++i) {
        Particle& P = p_list[i];
        double probParticle = kid.prob(&(P.mdstCharged()));
        if (P.mdstCharged() && probParticle < prob) { // a big difference compated to withKaonId

            if (!&P.userInfo()) createUserInfo(P);
            UserInfo& info = dynamic_cast<UserInfo&>(P.userInfo());
            info.isAdoptCut(true);
        }
        else {
            p_list.erase(p_list.begin() + i);
            --i;
    }
  }
}
// #=============== ================ ================#

// **********************************************************
// Modify this function
// Make sure that this function returns a correct value
int getEvtGenType() {
    // 0:"Data", 1:"evtgen-charged", 2:"evtgen-mixed", 3:"evtgen-charm", 4:"evtgen-uds", 5:"evtgen-bsbs", 6:"evtgen-nonbsbs"
    Gen_hepevt_Manager& genMgr  = Gen_hepevt_Manager::get_manager();
    int status = 0;
    for (std::vector<Gen_hepevt>::iterator i = genMgr.begin(); i != genMgr.end(); ++i) {
        int idhep = i->idhep();
        switch (abs(idhep)) {
            case 533:
                status = 1; // Bs0* (vector)
            case 531:
                status = 2; // Bs0  (scalar)
            case 431:
                status = 3; // Ds+-
            case 10431:
                status = 4; // Ds1*(2317)
            case 4:
                status = 5; // ccbar
            case 3:
                status = 6; // uds
            case 2:
                status = 6; // uds
            case 1:
                status = 6; // uds

        }
    }
    return status;
}
// ***********************************************************
void dumpEventInfo(BelleTuple* tt, bool debugDump = false) {
    // Event Information
    int expNo, runNo, evtNo;
    getEventInfo(expNo, runNo, evtNo, McFlag); // utility.cc
    const HepPoint3D& ip_position = IpProfile::position(1);
    const HepSymMatrix& ip_error = IpProfile::position_err(1);
    
    int idGenType = getEvtGenType();

    // Event Shape 
    double  r2 = -1.;
    Evtcls_hadron_info_Manager& clsMgr = Evtcls_hadron_info_Manager::get_manager();
    if (clsMgr.count()) r2 = clsMgr[0].R2();

    if (debugDump)
        printf("\n        ---- EvtInfo_DUMP --- EvtGenType[%i],  exp:%2i,  run:%2i, evt:%i, ip_position:[%f, %f, %f] ----\n",
            idGenType, expNo, runNo, evtNo,
            ip_position.x(), ip_position.y(), ip_position.z());

    tt->column("expn",  expNo); // Exp #
    tt->column("runn",  runNo); // Run #
    tt->column("evtn",  evtNo);
    tt->column("ipx",   ip_position.x());
    tt->column("ipy",   ip_position.y());
    tt->column("ipz",   ip_position.z());
    //************* Signal shape ********
    tt->column("r2", r2);
    //*************  MC  ****************
    tt->column("evtgen", idGenType);
}
// ***********************************************************
VectorL getGenVectorL(int idhPcl) {
    VectorL pclL;
    int ID_Pcl = -1; 
    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    for (std::vector<Gen_hepevt>::iterator itr = genMgr.begin(); itr != genMgr.end(); ++itr) {
        int idh = (*itr).idhep();
        if (idh == idhPcl) { // Ds+/-
            ID_Pcl = (*itr).get_ID();
            pclL.setPx(genMgr[ID_Pcl - 1].PX());
            pclL.setPy(genMgr[ID_Pcl - 1].PY());
            pclL.setPz(genMgr[ID_Pcl - 1].PZ());
            pclL.setE (genMgr[ID_Pcl - 1].E());
        }
    }
    return pclL;
}
// ***********************************************************
void dumpPi0(BelleTuple* tt, Particle& p0, const std::string& sfx, bool debugDump) {

    Particle& g1 = p0.child(0);
    Particle& g2 = p0.child(1);

    double msPi0_gg = (g1.p() + g2.p()).m();
    double psrPi0   = pStar(p0.p(), BeamEnergy::E_HER(), BeamEnergy::E_LER(), BeamEnergy::Cross_angle()).vect().mag();
    double pPi0     = p0.momentum().p().mag();

    bool gen_pi0 = IDhep(p0) == 0 ? false : true;
    bool gen_g1  = IDhep(g1) == 0 ? false : true;
    bool gen_g2  = IDhep(g2) == 0 ? false : true;

    const int nValI = 3;
    const int nValD = 4;

    std::string pclTitI[nValI] = {"gen", "gg1", "gg2"};
    int valPclI[nValI]    = {gen_pi0, gen_g1, gen_g2};
    std::string pclTitD[nValD] = {"eg1", "eg2", "psr", "mgg"};
    double valPclD[nValD] = {g1.ptot(), g2.ptot(), psrPi0, msPi0_gg};

    if (debugDump) {
        printf("  ==== dumpValues ==== pi0 (%s): ", sfx.c_str());
        for (int iVal = 0; iVal < nValI; iVal++) 
            printf("(%s:%6i) ", (pclTitI[iVal] + sfx ).c_str(), valPclI[iVal]);
        for (int iVal = 0; iVal < nValD; iVal++) 
            printf("(%s:%6.3f) ", (pclTitD[iVal] + sfx ).c_str(), valPclD[iVal]);
        printf("\n");
    }

    for (int iVal = 0; iVal < nValI; ++iVal) tt->column(pclTitI[iVal] + sfx, valPclI[iVal]);
    for (int iVal = 0; iVal < nValD; ++iVal) tt->column(pclTitD[iVal] + sfx, valPclD[iVal]);
}

// ***********************************************************
void dumpValues(BelleTuple* tt, int nValI, int nValD, int* valPclI, double* valPclD, std::string* pclTitI, std::string* pclTitD,
                                                                                            const std::string& sfx, bool debugDump) {
    if (debugDump) {
        printf("  ==== dumpValues ==== dec:(%s): ", sfx.c_str());
        for (int iVal = 0; iVal < nValI; ++iVal)
            printf("(%s:%6i) ", (pclTitI[iVal] + sfx).c_str(), valPclI[iVal]);
        for (int iVal = 0; iVal < nValD; ++iVal)
            printf("(%s:%6.3f) ", (pclTitD[iVal] + sfx).c_str(), valPclD[iVal]);
        printf("\n");
    }

    for (int iVal = 0; iVal < nValI; ++iVal)
        tt->column(pclTitI[iVal] + sfx, valPclI[iVal]);
    for (int iVal = 0; iVal < nValD; ++iVal)
        tt->column(pclTitD[iVal] + sfx, valPclD[iVal]);
}
// ***********************************************************
void dumpGenValues(BelleTuple* tt, bool gen_pcl, VectorL pclL, const std::string& sfx, bool debugDump) {
    const int nValD = 9;
    std::string pclTitD[nValD] = {"ms", "px", "py", "pz", "e", "th", "ct", "ph", "rh"};
    double valPclD[nValD] = {pclL.m(), pclL.px(), pclL.py(), pclL.pz(), pclL.e(), pclL.theta(), pclL.cosTheta(), pclL.phi(), pclL.rho()};
    if (!gen_pcl) {
        for (int iVal = 0; iVal < nValD; ++iVal)
            valPclD[iVal] = -.99;
    }
    if (debugDump) {
        printf("  ==== dumpGenValues ==== dec:(%s): ", sfx.c_str());
        for (int iVal = 0; iVal < nValD; ++iVal)
            printf("(%s:%6.3f) ", (pclTitD[iVal] + sfx ).c_str(), valPclD[iVal]);
        printf("\n");
    }

    for (int iVal = 0; iVal < nValD; ++iVal) {
        tt->column(pclTitD[iVal] + sfx, valPclD[iVal]);
    }
}
// ***********************************************************
void dumpDsChild(BelleTuple* tt, Particle& P, const std::string& sfxDs = "",
                 bool evtInfoDump = false,
                 bool stDump = true, bool debugDump = false) {


    if (evtInfoDump) dumpEventInfo(tt, debugDump);

    int lund = (int)P.lund();
    
    if (!&P.userInfo()) createUserInfo(P);
    UserInfo& info = dynamic_cast<UserInfo&>(P.userInfo());
    if (info.chisqKvf() < 0.) {
        makeRecursiveVertexFit(P, debugDump);
    }



    double chisqKvf      = info.chisqKvf();
    double chisqKmvf     = info.chisqKmvf();
    double chisqKvfNdf  = info.chisqKvfNdf();
    double chisqKmvfNdf = info.chisqKmvfNdf();
    double msComb        = info.msComb();
    double msKvf         = info.msKvf();
    double msKmvf        = info.msKmvf();
    double cl            = info.cl();
    double clKvf         = info.clKvf();
    double clKmvf        = info.clKmvf();
    double helic         = info.helicity();
    
    int ind_ch1 = 0;
    int lundChild = (int)P.lund();
    std::string ChildPcl;
    if (abs(lundChild) == 333) {
        ChildPcl ="phipi";
        ind_ch1  = 1;
    }
    if (abs(lundChild) == 313) {
        ChildPcl = "KsrK";
        ind_ch1  = 2;
    }

    double psr_ds_child           = pStar(P.p(), BeamEnergy::E_HER(), BeamEnergy::E_LER(), BeamEnergy::Cross_angle()).vect().mag();
    double px_ds_child            = P.px();
    double py_ds_child            = P.py();
    double pz_ds_child            = P.pz();
    double p_ds_child             = P.momentum().p().mag();
    double decay_vx_ds_child      = P.momentum().decayVertex().x();
    double decay_vy_ds_child      = P.momentum().decayVertex().y();
    double decay_vz_ds_child      = P.momentum().decayVertex().z();

    Hep3Vector P3D                (px_ds_child, py_ds_child, pz_ds_child);
    double perp_ds_child          = P3D.perp();
    double phi_ds_child           = P3D.phi();
    double theta_ds_child         = P3D.theta();


    std::string dgrSuff = "_ch" + sfxDs;
    int signPcl = (int)P.lund() > 0 ? 1 : -1;
    bool gen_pcl = IDhep(P) == 0 ? false : true;

    const int nValI = 3;
    const int nValD = 18;
    int valPclI[nValI] = {signPcl, gen_pcl, ind_ch1};
    double valPclD[nValD] = {msKvf, msKmvf, chisqKvf, chisqKmvf, chisqKvfNdf, chisqKmvfNdf,
                             clKvf, clKmvf, decay_vx_ds_child, decay_vy_ds_child, decay_vz_ds_child,
                             perp_ds_child, phi_ds_child, theta_ds_child, psr_ds_child, px_ds_child, py_ds_child, pz_ds_child};

    std::string pclTitI[nValI] = {"chg", "gen", "ind"};
    std::string pclTitD[nValD] = {"ms", "msm", "chi", "chim", "chn", "chmn",
                                  "cl", "clm", "dvx", "dvy", "dvz",
                                  "pt", "ph", "th", "psr", "px", "py", "pz"};

    if (debugDump) {
        printf("\n======== dumpDsChild  ========= (%s) chg_dss_child:%i, ms_dss_child:%7.3f ) \n",
               ChildPcl.c_str(), (int)P.lund(), msKvf);
        printUserInfo(P);
    }

    dumpValues(tt, nValI, nValD, valPclI, valPclD, pclTitI, pclTitD, dgrSuff, debugDump);
    
    if (stDump) tt->dumpData();
}
// ***********************************************************
std::pair<int, double> getDsChildID(Particle& DsChild,
                                    const int accq0 = 3,
                                    const int tofq0 = 1,
                                    const int cdcq0 = 5) {

    double ids0 = 3,
            idb0 = 2;
    int lund_dschild = DsChild.lund();
    switch (std::abs(lund_dschild)) {
        case 321:         // K
            break;
        case 211:         // pi
            ids0 = 2,
                    idb0 = 3;
            break;
    }
    atc_pid kid(accq0, tofq0, cdcq0, ids0, idb0);
    return std::make_pair(DsChild.lund(), kid.prob(&(DsChild.mdstCharged())));
}
// ***********************************************************
std::vector<std::pair<int, double> > getDsChildrenIdentification(Particle& DsMeson) {
    std::vector<std::pair<int, double> > values;
    std::vector<Particle> ds_children = {
            DsMeson.child(0).child(0),
            DsMeson.child(0).child(1),
            DsMeson.child(1)
    };
    for (std::size_t i = 0; i < ds_children.size(); ++i)
        values.push_back(getDsChildID(ds_children[i]));

    return values;
}
// ***********************************************************
void dumpDs(BelleTuple* tt, Particle& P, std::string sfxDs = "", bool evtInfoDump = false, bool stDump = true, bool debugDump = false) {

    if (evtInfoDump) dumpEventInfo(tt, debugDump);
    int lundDs = (int)P.lund();

    if (!&P.userInfo()) createUserInfo(P);
    UserInfo& info = dynamic_cast<UserInfo&>(P.userInfo());
    if (info.chisqKvf() < 0.) {
        makeRecursiveVertexFit(P, debugDump, true);
    }

    double chisqKvf      = info.chisqKvf();
    double chisqKmvf     = info.chisqKmvf();
    double chisqKvfNdf  = info.chisqKvfNdf();
    double chisqKmvfNdf = info.chisqKmvfNdf();
    double msKvf         = info.msKvf();
    double msKmvf        = info.msKmvf();
    double clKvf         = info.clKvf();
    double clKmvf        = info.clKmvf();


    double psr_ds           = pStar(P.p(), BeamEnergy::E_HER(), BeamEnergy::E_LER(), BeamEnergy::Cross_angle()).vect().mag();
    double px_ds            = P.px();
    double py_ds            = P.py();
    double pz_ds            = P.pz();
    Hep3Vector P3D          (px_ds, py_ds, pz_ds);
    double perp_ds          = P3D.perp();
    double phi_ds           = P3D.phi();
    double theta_ds         = P3D.theta();
    int    signPcl          = lundDs > 0 ? 1 : -1;
    bool   gen_pcl          = IDhep(P) == 0 ? false : true;
    int    c00              = 0,
           c01              = 0,
           c10              = 0;
    double pid00            = -1,
           pid01            = -1,
           pid10            = -1;
    double decay_vx_ds      = P.momentum().decayVertex().x();
    double decay_vy_ds      = P.momentum().decayVertex().y();
    double decay_vz_ds      = P.momentum().decayVertex().z();
    double helicChild1      = getHelicity(P);

    std::vector<std::pair<int, double> > ds_child_ids = getDsChildrenIdentification(P);
    c00              = ds_child_ids[0].first;
    c01              = ds_child_ids[1].first;
    c10              = ds_child_ids[2].first;
    pid00            = ds_child_ids[0].second;
    pid01            = ds_child_ids[1].second;
    pid10            = ds_child_ids[2].second;

    Particle& Child0        = P.child(0); // phi0 (--> K+ K-) or K* (--> K+ pi-), or K*bar (--> K- pi+)

    std::string dgrSuff = "_ds" + sfxDs, genDgrSuff = "_ds_t" + sfxDs;

    const int nValI = 5;
    const int nValD = 24;
    int valPclI[nValI]    = {signPcl, gen_pcl, c00, c01, c10};
    double valPclD[nValD] = {msKvf, msKmvf, chisqKvf, chisqKmvf, chisqKvfNdf, chisqKmvfNdf,
                             clKvf, clKmvf, decay_vx_ds, decay_vy_ds, decay_vz_ds,
                             psr_ds, px_ds, py_ds, pz_ds, perp_ds,
                             phi_ds, theta_ds, pid00, pid01, pid10};
    std::string pclTitI[nValI] = {"chg", "gen", "c00", "c01", "c10"};
    std::string pclTitD[nValD] = {
            "ms", "msm", "chi", "chim", "chn", "chmn",
            "cl", "clm", "dvx", "dvy", "dvz",
            "psr", "px", "py", "pz", "pt",
            "ph", "th", "p00", "p01", "p10"};


    VectorL dssL = getGenVectorL(IDhep(P));

    if (debugDump) {
        printf("\n======== Dss  ========= sfx:%s, chg_ds:%i, gen_ds:%i, ms_ds:%7.3f , child ( ms:%7.3f ) \n", 
               genDgrSuff.c_str(), signPcl, gen_pcl, msKvf, Child0.p().m());
        printUserInfo(P);
    }

    tt->column("hel_ch" + sfxDs, helicChild1);
    dumpValues(tt, nValI, nValD, valPclI, valPclD, pclTitI, pclTitD, dgrSuff, debugDump);
    dumpGenValues(tt, gen_pcl, dssL, genDgrSuff, debugDump);
    dumpDsChild(tt, Child0, sfxDs, false, false, debugDump);
    if (stDump) tt->dumpData();
}
// ***********************************************************
void dumpDs2317(BelleTuple* tt, Particle& P, std::string sfxDs = "", bool evtInfoDump = false, bool stDump = true, bool debugDump = false) {

    if (evtInfoDump) dumpEventInfo(tt, debugDump);

    if (!&P.userInfo()) createUserInfo(P);
    UserInfo& info = dynamic_cast<UserInfo&>(P.userInfo());
    if (info.chisqKvf() < 0.) {
        makeRecursiveVertexFit(P, debugDump);
    }
    
    double chisqKvf      = info.chisqKvf();
    double chisqKmvf     = info.chisqKmvf();
    double chisqKvfNdf  = info.chisqKvfNdf();
    double chisqKmvfNdf = info.chisqKmvfNdf();
    double msKvf         = info.msKvf();
    double msKmvf        = info.msKmvf();
    double cl            = info.cl();
    double clKvf         = info.clKvf();
    double clKmvf        = info.clKmvf();
    double helic_2317    = -1.;
    

    Particle& Child           = P.child(0);
    Particle& pi0_2317        = P.child(1);
    UserInfo& infoChild       = dynamic_cast<UserInfo&>(Child.userInfo());
    UserInfo& infoPi0_2317    = dynamic_cast<UserInfo&>(pi0_2317.userInfo());

    double msKvfChild         = infoChild.msKvf();
    double psr_d17            = pStar(P.p(), BeamEnergy::E_HER(), BeamEnergy::E_LER(), BeamEnergy::Cross_angle()).vect().mag();
    double px_d17             = P.px();
    double py_d17             = P.py();
    double pz_d17             = P.pz();
    double p_d17              = P.momentum().p().mag();
    double decay_vx_d17       = P.momentum().decayVertex().x();
    double decay_vy_d17       = P.momentum().decayVertex().y();
    double decay_vz_d17       = P.momentum().decayVertex().z();


    Hep3Vector P3D            (px_d17, py_d17, pz_d17);
    double perp_d17           = P3D.perp();
    double phi_d17            = P3D.phi();
    double theta_d17          = P3D.theta();
    int chg_d17               = (int)P.lund() >  0 ? 1 : -1;
    int gen_d17               = (int)IDhep(P) == 0 ? 0 :  1;

    const int nValI = 2;
    const int nValD = 20;
    int valPclI[nValI] = {chg_d17, gen_d17};
    double valPclD[nValD] = {msKvf, msKmvf, chisqKvf, chisqKmvf, chisqKvfNdf, chisqKmvfNdf,
                             clKvf, clKmvf, decay_vx_d17, decay_vy_d17, decay_vz_d17, psr_d17,
                             px_d17, py_d17, pz_d17, p_d17, perp_d17, phi_d17, theta_d17, helic_2317};
    std::string pclTitI[nValI] = {"chg", "gen"};
    std::string pclTitD[nValD] = {"ms","msm","chi", "chim", "chn", "chmn",
                                  "cl", "clm", "dvx", "dvy", "dvz", "psr",
                                  "px", "py", "pz", "p", "pt", "ph", "th", "hel"};

    std::string dgrSuff = "_d17", genDgrSuff = "_d17_t";

    if (debugDump) {
        printf("\n======== Ds(2317)  ========= chg_2317:%i, gen_2317:%i, ms_2317:%7.3f , child ( ms:%7.3f ) \n", 
               chg_d17, gen_d17, msKvf, msKvfChild);
        printUserInfo(P);
    }

    VectorL ds17L = getGenVectorL(IDhep(P));

    dumpValues(tt, nValI, nValD, valPclI, valPclD, pclTitI, pclTitD, dgrSuff, debugDump);
    dumpGenValues(tt, gen_d17, ds17L, genDgrSuff, debugDump);
    dumpPi0(tt, pi0_2317, "_p0_d", debugDump);
    dumpDs(tt, Child, sfxDs, false, false, debugDump);

    if (stDump) tt->dumpData();
}   
// ***********************************************************
void dumpBs0(BelleTuple* tt, Particle& P, bool evtInfoDump = false,
                 bool stDump = true, bool debugDump = true) {

    if (evtInfoDump) dumpEventInfo(tt, debugDump);

    if (!&P.userInfo()) createUserInfo(P);
    UserInfo& info = dynamic_cast<UserInfo&>(P.userInfo());

    if (info.chisqKvf() < 0.) {
        makeRecursiveVertexFit(P, debugDump);
    }
    double chisqKvf      = info.chisqKvf();
    double chisqKmvf     = info.chisqKmvf();
    double chisqKvfNdf   = info.chisqKvfNdf();
    double chisqKmvfNdf  = info.chisqKmvfNdf();
    double msKvf         = info.msKvf();
    double msKmvf        = info.msKmvf();
    double clKvf         = info.clKvf();
    double clKmvf        = info.clKmvf();

    Particle& Dss_Bs0 = P.child(0);
    Particle& Dss2317_Bs0 = P.child(1);
    Particle& pi0_Bs0 = P.child(2);

    double px_bs = P.px();
    double py_bs = P.py();
    double pz_bs = P.pz();
    double p_bs  = P.momentum().p().mag();

    double xx_bs = P.momentum().x().x();
    double xy_bs = P.momentum().x().y();
    double xz_bs = P.momentum().x().z();

    Hep3Vector P3D(px_bs, py_bs, pz_bs);

    double perp_bs = P3D.perp();
    double phi_bs  = P3D.phi();
    double theta_bs = P3D.theta();

    double decay_vx = P.momentum().decayVertex().x();
    double decay_vy = P.momentum().decayVertex().y();
    double decay_vz = P.momentum().decayVertex().z();
    int gen_bs = (int)IDhep(P) == 0 ? 0 : 1;
    int chg_bs = (int)P.lund() > 0 ? 1 : -1;
    

    VectorL pB = pStar(P.p());
    VectorL pB2 = pStar(P.p(), BeamEnergy::E_HER(), BeamEnergy::E_LER(), BeamEnergy::Cross_angle());
    double de_bs_old = pB.e() - Benergy();
    double de_bs = pB.e() - BeamEnergy::E_beam_corr();
    double de2_bs = pB2.e() - BeamEnergy::E_beam_corr();
    double mbc_bs_old = beamEnergyConstraint(P);
    double mbc_bs    = beamEnergyConstraint(P, BeamEnergy::E_HER(), BeamEnergy::E_LER(), BeamEnergy::Cross_angle() * 1.e3);

    const int nValI = 2;
    const int nValD = 21;
    int valPclI[nValI] = {chg_bs, gen_bs};
    double valPclD[nValD] = {msKvf, msKmvf, chisqKvf, chisqKmvf, chisqKvfNdf, chisqKmvfNdf,
                             clKvf, clKmvf, decay_vx, decay_vy, decay_vz,
                             px_bs, py_bs, pz_bs, p_bs, perp_bs,
                             phi_bs, theta_bs, mbc_bs, de_bs, de2_bs};
    std::string pclTitI[nValI] = {"chg", "gen"};
    std::string pclTitD[nValD] = {"ms", "msm", "chi", "chim", "chn", "chmn",
                                  "cl", "clm", "dvx", "dvy", "dvz",
                                  "px", "py", "pz", "p", "pt",
                                  "ph", "th", "mbc", "de", "de2"};

    if (debugDump) {
        printf("\n\n======== Bs0  ========= chg_bs:%i,  gen_bs:%i, ms_bs:%7.3f , de_bs: %7.3f , de_bs_old: %7.3f, mbc_bs: %7.3f , mbc_bs_old: %7.3f \n", 
               chg_bs, gen_bs, msKvf, de_bs, de_bs_old, mbc_bs, mbc_bs_old);
        printUserInfo(P);
    }

    VectorL bs0L = getGenVectorL(IDhep(P));
    std::string genSfx = "_bs_t";

    dumpValues   (tt,  nValI, nValD, valPclI, valPclD, pclTitI, pclTitD, "_bs", debugDump);
    dumpGenValues(tt, gen_bs, bs0L, genSfx, debugDump);
    dumpPi0      (tt, pi0_Bs0,     "_p0_b", debugDump);
    dumpDs       (tt, Dss_Bs0,     "1", false, false, debugDump);
    dumpDs2317   (tt, Dss2317_Bs0, "2", false, false, debugDump);

    if (stDump) tt->dumpData();
}
// ***********************************************************
void printVectPclWithChildren(std::vector<Particle>& pcl, const std::string& tit = "") {
    if (pcl.size() > 0) {
        printf("    ---- %s [%i] ----- \n", tit.c_str(), pcl.size());
        for (int i = 0; i < pcl.size(); ++i) {
            printPclDebug(pcl[i]);
            for (int ich = 0; ich < pcl[i].nChildren(); ++ich)
                printPclDebug(pcl[i].child(ich));
        }
    }
    printf("\n");
}
// ***********************************************************

/*
 * The main function. Events processing takes place here
*/
void Reco::event(BelleEvent* evptr, int* status) {
    *status = 0;
    bool debugTrkPID   = false,
         debugVtx      = false,
         debugPi0      = false,
         debugPhiKsr   = false,
         debugDss      = false,
         debugDss_2317 = false,
         debugBs0      = false,
         debugDumpDss  = false,
         debugDump2317 = false,
         debugDumpBs0  = false;

    // IP information
    const HepPoint3D& ip_position = IpProfile::position();
    const HepSymMatrix& ip_error  = IpProfile::position_err();
    // Gen_hepevt_Manager& hepevt = Gen_hepevt_Manager::get_manager();

    // Event Information
    int expNo = 0, runNo = 0, evtNo = 0;
    getEventInfo(expNo, runNo, evtNo, McFlag); // utility.cc
    std::cout << "\tExpNo: " << expNo << "\tRunNo: " << runNo << "\tEventNo: " << evtNo << "\tMcFlag: " << McFlag << std::endl;
  
    // Event Shape
    double r2 = -1.;

    Evtcls_hadron_info_Manager& clsMgr = Evtcls_hadron_info_Manager::get_manager();
    if (clsMgr.count()) r2 = clsMgr[0].R2();

    // //////////////////  make charged particles - tracks //////////////////
    // makes Kaon and Pion from MdstCharged w/o cut. 1 : w/ good_charged, 0 : w/ol
    // ////////////////// ORDER: k_p, k_m, pi_p, pi_m ///////////////////////
    makeKPi(trkV[2], // k_p
            trkV[3], // k_m
            trkV[0], // pi_p
            trkV[1], // pi_m
            1);

    if (debugTrkPID) {
        printf("\n");
        printf("---- exp:%2i,  run:%2i, evt:%i ----\n", expNo, runNo, evtNo);
        for (int itr = 0; itr < nTrk; ++itr)
            printTrkPID(trkV[itr], trkTit[itr], "before");
    }

    // If each plist element is not within atc_pID.prob >= prob,
    // its element is removed from plist.
    // ----------------- pos kaons ---- defalut values - kaons - pions
    withKaonId(trkV[2], minProbPID_Kp,  3,    1,     5,  3,      2);      // K+ vs bg pi
    // ----------------- neg kaons ---- defalut values - kaons - pions
    withKaonId(trkV[3], minProbPID_Kn,  3,    1,     5,  3,      2);      // K- vs bg pi


    // If each plist element is not within atc_pID.prob < prob,
    // its element is removed from plist.
    // ----------------- pos pions ---- defalut values - kaons - pions
    withPionId(trkV[0], maxProbPion,    3,    1,     5,  3,      2);        // pi+ vs bg K
    // ----------------- neg pions ---- defalut values - kaons - pions
    withPionId(trkV[1], maxProbPion,    3,    1,     5,  3,      2);        // pi- vs bg K

    // If each plist element is not associated with rphi & z-svd hits
    // whose number is equal to or larger than nRSvdHit and nZSvdHit,
    // its element is removed from plist.
    for (int itr = 0; itr < nTrk; ++itr) {
        withSVD2(trkV[itr], 1, 1); // nRSvdHit, nZSvdHit
        withdRdZcut(trkV[itr], ip_position, dRcut, dZcut);
    }

    if (debugTrkPID) {
        for (int itr = 0; itr < nTrk; ++itr)
            printTrkPID(trkV[itr], trkTit[itr], "after");
    }

    // ================= WORKING WITH Gamma CANDIDATES ================= //
    makeGamma(gammaV);
    withPCut(gammaV, eGammaMin);

    if(useVTX) {
        for (std::vector<Particle>::iterator itr = gammaV.begin(); itr != gammaV.end(); ++itr) {
            setGammaError(*itr, ip_position, ip_error); // changed from setGammasError
        }
    }

    // Match gamma candidates with their genhep info
    if (McFlag) {
        setGenHepInfoG(gammaV);
    }

    createUserInfo(gammaV);


    // =================   WORKING WITH Pi0 CANDIDATES ================= //
    makePi0(pi0);
    // Checking gammas' energies for all the pi0 daughters
    withPi0GammPCut(pi0, minPi0GammP);
    // !!!! PAY ATTENTION. Start point

    // Creating UserInfo objects in memory storage for the selected pi0's
    createUserInfo(pi0);

    // Setting error matrices for both pi0 daughters gamma
    // it's necessary for vertex fitting
    for (std::vector<Particle>::iterator itr = pi0.begin(); itr != pi0.end(); ++itr) {
        setGammasError(*itr, ip_position, ip_error);
    }

    // Making a simple Kvf fit
    makeRecursiveVertexFit(pi0, false, false);

    // Selecting pi0 candidates considering their reconstructed masses (mass of 2 gammas)
    withPi0MassGamGamCut(pi0, wMassPi0GG);

    // Making a mass-constraint fit for a pi0 candidate, which passed through all cuts (including gammas)
    makeRecursiveVertexFit(pi0, false, true);


    if (debugPi0) {
        printf(" pi0[%i]  \n", pi0.size());
        printPi0(pi0);
    }

    // Match candidates with genhep info 
    if (McFlag) {
        for (int itr = 0; itr < nTrk; ++itr) setGenHepInfoF(trkV[itr]);
        setGenHepInfoP(pi0);
    }    
    
    combination(phi0,      Ptype("PHI"),  trkV[2],   trkV[3], dM_V0);   // k_p, k_m
    combination(Ksr0,      Ptype("K*0"),  trkV[2],   trkV[1], dM_Ksr0); // k_p, pi_m
    combination(Ksr0bar,   Ptype("K*B"),  trkV[3],   trkV[0], dM_Ksr0); // k_m, pi_p
    setGenHepInfoT(phi0);
    setGenHepInfoT(Ksr0);
    setGenHepInfoT(Ksr0bar);

    
    if (debugPhiKsr) {
        printf("         phi0[%i]  \n", phi0.size());
        printf("         Ksr0[%i]  \n", Ksr0.size());
        printf("   Ksr0bar[%i]  \n", Ksr0bar.size());
    }

    combination(Dss_p, Ptype("DS+"),    phi0,     trkV[0], dM_Dss); // phi0, pi_p
    combination(Dss_m, Ptype("DS-"),    phi0,     trkV[1], dM_Dss); // phi0, pi_m
    combination(Dss_p, Ptype("DS+"),    Ksr0bar,  trkV[2], dM_Dss); // K*0bar (K-pi+), k_p
    combination(Dss_m, Ptype("DS-"),    Ksr0,     trkV[3], dM_Dss); // K*0 (K+pi-), k_m
    setGenHepInfoT(Dss_p);
    setGenHepInfoT(Dss_m);

    
    if (debugDss) {
        if (Dss_p.size() + Dss_m.size() > 0) {
            printf("\n **** debugDss ****** exp:%2i,  run:%2i, evt:%i,  *******\n", expNo, runNo, evtNo);
            printVectPclWithChildren(Dss_p, "Ds+");
            printVectPclWithChildren(Dss_m, "Ds-");
        }
    }


    combination(Dss_m_2317, Ptype(-10431), Dss_m, pi0, dM_2317);
    combination(Dss_p_2317, Ptype( 10431), Dss_p, pi0, dM_2317);
    setGenHepInfoT(Dss_m_2317);
    setGenHepInfoT(Dss_p_2317);

    
    if (debugDss_2317) {
        if (Dss_p_2317.size() + Dss_m_2317.size() > 0) {
            printf("\n **** debugDss_2317 ****** exp:%2i,  run:%2i, evt:%i,  *******\n", expNo, runNo, evtNo);
            printVectPclWithChildren(Dss_p_2317, "DsJ(2317)+");
            printVectPclWithChildren(Dss_m_2317, "DsJ(2317)-");
        }
    }

    combination(Bs0, Ptype(531), Dss_m, Dss_p_2317, pi0, dM_Bs0);
    combination(Bs0, Ptype(531), Dss_p, Dss_m_2317, pi0, dM_Bs0);
     combination(Bs0bar, Ptype(-531), Dss_m, Dss_p_2317, pi0, dM_Bs0);
     combination(Bs0bar, Ptype(-531), Dss_p, Dss_m_2317, pi0, dM_Bs0);
    setGenHepInfoT(Bs0);
    setGenHepInfoT(Bs0bar);
    
    // ----------------------------  Dumping  ---------------------------
    // Ds
    /*
    std::string sfxDs = "";
    if (stDumpDss) {
        for (int iEvt = 0; iEvt < Dss_p.size(); ++iEvt)
            dumpDs(TP_Dss, Dss_p[iEvt], sfxDs, true, stDumpDss, debugDumpDss);
        for (int iEvt = 0; iEvt < Dss_m.size(); ++iEvt)
            dumpDs(TP_Dss, Dss_m[iEvt], sfxDs, true, stDumpDss, debugDumpDss);
    }    
    // Ds(2317)
    if (stDump2317) {
        for (int iEvt = 0; iEvt < Dss_p_2317.size(); ++iEvt)
            dumpDs2317(TP_Dss_2317, Dss_p_2317[iEvt], sfxDs, true, stDump2317, debugDump2317);
        for (int iEvt = 0; iEvt < Dss_m_2317.size(); ++iEvt)
            dumpDs2317(TP_Dss_2317, Dss_m_2317[iEvt], sfxDs, true, stDump2317, debugDump2317);
    }
    */

    // Bs0
    if (stDumpBs0) {
        for (int iEvt = 0; iEvt < Bs0.size(); ++iEvt)
            dumpBs0(TP_Bs0, Bs0[iEvt], true, stDumpBs0, debugDumpBs0);
        for (int iEvt = 0; iEvt < Bs0bar.size(); ++iEvt)
            dumpBs0(TP_Bs0, Bs0bar[iEvt], true, stDumpBs0, debugDumpBs0);
    }
    printf("---------  clearVectors (final) --------- \n");
    clearVectors();
    printf("------------------  Reco end -------------------- \n");
    *status = 1;

}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

