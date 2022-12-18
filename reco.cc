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
string pclType(int id) {
    string sid;
    if      (id ==   0) sid = "--";
    else if (id ==   1) sid = "d";
    else if (id==  -1) sid = "anti-d";
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
  double   E =     gamma.mdstGamma().ecl().energy();

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

/************ Set error matrix (dpx) for pi0 ******************/
void setPi0Error(Particle& p) {
    if (p.nChildren() != 2) return;
    HepSymMatrix tmpErr(3, 0);
    HepPoint3D org(0., 0., 0.);
    if (!p.child(0).mdstGamma() || !p.child(1).mdstGamma()) return;
    for (unsigned i = 0; i < 2; ++i) 
        setGammaError(p.child(i));
//         setGammaError( p.child(i), org, tmpErr );

//     kmassfitter kmv;
    kmassvertexfitter kmv;
    if (useBF) kmv.magneticField(BF);
    kmv.invariantMass(p.pType().mass());
    for (unsigned i = 0; i < 2; ++i)
        addTrack2fit(kmv, p.child(i));
//     kmv.vertex(org);
//     kmv.atDecayPoint();
    int err=kmv.fit();
    if (!err) {
        kmakemother kmm2;
        if (useBF) kmm2.magneticField(BF);
        makeMother(kmm2, kmv, p, 0);
//         makeMother(kmv,p);
        p.momentum().vertex(kmv.vertex(), kmv.errVertex());
    }
}

void setPi0Error(std::vector<Particle>& p) {
  for (std::vector<Particle>::iterator i = p.begin(); i!=p.end(); ++i)
      setPi0Error(*i);
}
//***********************************************************
double getMsGammaInPi0(Particle& gam) {
//     std::vector<Particle> pi0;
//     makePi0(pi0);
//     withPi0GammPCut( pi0,  0.030 );
//     withPi0pStarCut( pi0,  0.050 );

    double msGG = 0.25;
    for (int iPi0 = 0; iPi0 < pi0.size(); ++iPi0) {
        Particle& g1 = pi0[iPi0].child(0);
        Particle& g2 = pi0[iPi0].child(1);
        bool idGamPi0 = false;
        if (abs(g1.ptot() - gam.ptot()) < 1.e-3) idGamPi0 = true;
        if (abs(g2.ptot() - gam.ptot()) < 1.e-3) idGamPi0 = true;
        if (idGamPi0) {
            double ms_gg = (g1.p() + g2.p()).m();
            if (abs(ms_gg - Ptype("PI0").mass()) < msGG) 
                msGG = ms_gg;
        }
    }
//     pi0.clear();
    return msGG;
}
//***********************************************************
void withGammaInPi0(std::vector<Particle>& Gamma, std::vector<Particle>& pi0, 
            double minMsPi0=0.118, double maxMsPi0=0.150) {
    for(int i = 0; i < (int)Gamma.size(); ++i) {
        bool idErase = false;
        for (int iPi0 = 0; iPi0 < pi0.size(); ++iPi0) {
            Particle& g1 = pi0[iPi0].child(0);
            Particle& g2 = pi0[iPi0].child(1);
            double ms_gg = (g1.p() + g2.p()).m();
            if (ms_gg > minMsPi0 && ms_gg < maxMsPi0) {
                if (abs(g1.ptot() - Gamma[i].ptot()) < 1.e-3) idErase = true;
                if (abs(g2.ptot() - Gamma[i].ptot()) < 1.e-3) idErase = true;
            }
        }

        if(idErase) {
            Gamma.erase(Gamma.begin() + i);
            --i;
        }
    }
}
//***********************************************************
int IDhep(Particle& part) {
    if(!part.genHepevt()) return 0;
    return part.genHepevt().idhep();
}

//***********************************************************
bool isLikeTrk(int lund, string chrgType="Charged_pi0") {
    int ln = abs(lund);
    bool isTrk = (ln==11) || (ln==13) || (ln==211) || 
            (ln==321) || (ln==11) || (ln==2212); // el, mu, pi, K, p
    if ( chrgType!="onlyCharged" )
        isTrk = isTrk || (ln==111) || (ln==22) ;  // pi0 & gamma
    return isTrk;
}

//***********************************************************
void printPclDebug(Particle& p, string comment="", string comment2=""){
    int expNo, runNo, evtNo;
    bool McFlag;
    getEventInfo( expNo, runNo, evtNo, McFlag); // utility.cc

    int lund = (int)p.lund();
    string trkType = pclType(lund);
    int ln = abs(lund);
    bool isTrk = isLikeTrk( lund );
    Hep3Vector pcl3D( p.px(), p.py(), p.pz() );
    printf("\n----- Particle [%s]  %9s %s ---\n", trkType.c_str(), comment.c_str(), comment2.c_str() );
    printf(" massRECO:%8.5f,   ", p.p().m());
//     if (!isTrk)     
//         if ( p.userInfo() ) printUserInfo( p );
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
//         printf( "      distanceToIP:%8.5f \n", distanceToIP(p) );
        printf("Children  ID: ");
        for (int jM =0 ; jM < p.nChildren(); ++jM) {
            Particle& Child = p.child(jM);
            int lundChild = (int)Child.lund() ;
            string childType = pclType(lundChild);
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
//***********************************************************
void printPclDebug(vector<Particle>& plist, string comment="", string comment2=""){
    int np = plist.size();
    printf("\n----- List of Particles (%i) %9s %s ---\n", np, comment.c_str(), comment2.c_str() );
    for (int i=0; i<np; i++) {
        printPclDebug( plist[i], comment, comment);
    }
}
//***********************************************************
void printUserInfo(Particle& p) {
    UserInfo &info = dynamic_cast<UserInfo&>(p.userInfo());
    printf("----- printUserInfo -----  msComb:%8.5f, msKvf:%8.5f, chisq:%8.3f, chisqKvf:%8.3f, cl:%8.6f, clKvf:%8.6f, \n              dist2IP:%8.5f, dist2IPmvf:%8.5f, useTube:%i, useKmvf:%i,  isAdoptCut:%i, wMass:%8.5f, maxChi2:%8.3f,  helicity:%6.4f \n", 
           info.msComb(), info.msKvf(), info.chisq(), info.chisqKvf(), info.cl(), info.clKvf(), info.dist2IP(), info.dist2IPmvf(), 
           (int)info.useTube(), (int)info.useKmvf(), (int)info.isAdoptCut(), info.wMass(), info.maxChi2(), info.helicity() );
}
//***********************************************************

void dRdZ(const Mdst_charged& charged, int massHyp, const HepPoint3D& ip_position,
                                                            double& dR, double& dZ) {

    if (charged.trk()) {
        Mdst_trk& trk = charged.trk();
        if (trk.mhyp(massHyp)) {
            Mdst_trk_fit &trkFit = trk.mhyp(massHyp);

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


void withdRdZcut(std::vector<Particle>& list, const HepPoint3D& ip_position, double dRcut, double dZcut) {
    for (int i = 0; i < (int) list.size(); ++i) {
        int id = abs(list[i].pType().lund());

        int mhyp;
        if (id == 11) mhyp = 0;
        else if (id == 13) mhyp = 1;
        else if (id == 211) mhyp = 2;
        else if (id == 321) mhyp = 3;
        else if (id == 2212) mhyp = 4;
        else mhyp = 2;

        double dr, dz;
        dRdZ(list[i].mdstCharged(), mhyp, ip_position, dr, dz);

        if (dr > dRcut || dz > dZcut) {
            list.erase(list.begin() + i);
            --i;
        }
    }
}
//**********************************************************
bool withPi0MassGamGamCut(Particle& pi0, double M_min, double M_max) {
    Particle& g1 = pi0.child(0);
    Particle& g2 = pi0.child(1);
    double msPi0_gg = (g1.p() + g2.p()).m();
    if (msPi0_gg < M_min) return false;
    if (msPi0_gg > M_max) return false;
    return true;
}
//**********************************************************
bool withPi0GammPCut(Particle& p, double p_min) {
    Particle& g1 = p.child(0);
    Particle& g2 = p.child(1);
    if(p.child(0).ptot() < p_min) return false;
    if(p.child(1).ptot() < p_min) return false;
    return true;
}
//**********************************************************
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

//**********************************************************
void withPi0GammPCut(std::vector<Particle>& pi0, double p_min) {
    for (std::vector<Particle>::iterator i = pi0.begin(); i != pi0.end();)
        if (i->child(0).ptot() < p_min || i->child(1).ptot() < p_min)
            pi0.erase(i);
        else
            ++i;
}
//***********************************************************
void createUserInfo(Particle& p) {
    p.userInfo(*(new UserInfo(p)));
    UserInfo& info = dynamic_cast<UserInfo&>(p.userInfo());
    int lund = (int)p.lund();
    double cTau = Ptype(lund).cTau();
    bool useTube = false;
    double wMass = dM_Dgr;

    if ((abs(lund) > 500) && (abs(lund) < 600)) { // B0, Bc
        useTube = true;
        wMass = wB;
    } else if ((lund      == 310) || (abs(lund) == 3122) || (lund == 333)) { // K0s, Lam0, Phi0
        wMass = dM_V0;
    } else if (abs(lund)  == 313) { // K*0
        wMass = dM_Ksr0;
    } else if ((abs(lund) == 113) || (abs(lund) == 213)) { // RHO0,+.-
        wMass = dM_Rho;
    } else if ((abs(lund) == 413) || (abs(lund) == 423)) { // D*+, D*0feve
        useTube = true; 
        wMass = wDst;
    } else if (abs(lund)  == 431) { // DS+
        wMass = dM_Dss;
    } else if (abs(lund)  == 433) { // D*S+
        useTube = true; 
        wMass = dM_Dsst;
    } else if (abs(lund)  == 10431) { // D**S+(D_sJ(2317))
        useTube = true;
        wMass = dM_2317;
    } else if (abs(lund)  == 443) { // J/psi
        useTube = true; 
    }
    
    info.msComb(p.p().m());
    info.maxChi2(maxChisq);
    info.wMass(wMass);
    info.useTube(useTube);
    info.useKmvf(cTau > 1.e-5 ? true : false);
    info.isAdoptCut(true);
    info.chisqKvf(-1.);
    info.helicity(-1.);
}

void createUserInfo(std::vector<Particle>& p_list) {
    std::vector<Particle>::iterator particle;
    for (particle = p_list.begin(); particle != p_list.end(); ++particle) {
        if (!&(*particle).userInfo()) createUserInfo(*particle);
    }
}
//***********************************************************
void setBadVtx(Particle& p) {
    HepPoint3D vtx(999., 999., 999.);
    HepSymMatrix errVtx(3, 0);
    p.momentum().decayVertex(vtx, errVtx);
    UserInfo &info = dynamic_cast<UserInfo&>(p.userInfo());
    info.isAdoptCut(false);
}
//***********************************************************

double distanceToIP(Particle& p) {
    const HepPoint3D& ip_position = IpProfile::position();
    Hep3Vector vtxIP(IpProfile::position().x(), IpProfile::position().y(), 0.0);
    Hep3Vector vtxD(p.momentum().decayVertex().x(), p.momentum().decayVertex().y(), 0.0);
    Hep3Vector vtx_D_IP = vtxD - vtxIP;
    return vtx_D_IP.mag();
}
//***********************************************************
void vertex_staged(Particle& Mother, bool debugDump=false, bool useKmvf=false) {
    // fit DecayMother -> DecayChild1 + DecayChild2 + ... + trk1 + trk2 + ...
//         Particle &Child ( Mother.child(0) );  // Daughter Decaying particle
    
        if (!&Mother.userInfo()) createUserInfo(Mother);
        UserInfo& infoMother = dynamic_cast<UserInfo&>(Mother.userInfo());

        int lundMother = (int)Mother.lund();
        string motherType = pclType(lundMother);
        if (debugDump) {
            printf("\n ========  vertex_staged ==========   %s [%i] --> ", 
                   motherType.c_str(), Mother.nChildren() );
            for (int j = 0; j < Mother.nChildren(); ++j) {
                string dghtType = pclType((int)Mother.child(j).lund());
                printf("%s ", dghtType.c_str());
            }
            printf("\n");
        }

        kvertexfitter kvfMother;
        if (useBF) kvfMother.magneticField(BF);
        int err;
        for (int jM = 0; jM < Mother.nChildren(); ++jM) {
            Particle& Child = Mother.child(jM);
            int lundChild = (int)Child.lund();
            bool isTrkChild = isLikeTrk(lundChild);
            if (!isTrkChild) {
                if (!&Child.userInfo()) createUserInfo(Child);
                UserInfo& infoChild = dynamic_cast<UserInfo&>(Child.userInfo());
                double chisqChild = infoChild.chisqKvf(); // -1.; //
                if (chisqChild < 0.) {
                    // particle not vertexed yet
                    bool idUseTube = false;
                    bool idAddMF = true; /*false*/  // magneticField (needed?)
                    vertex_staged(Child, debugDump);
                }
            }
            addTrack2fit(kvfMother, Child);

        }
        if (infoMother.useTube()) addTube2fit(kvfMother);
        err = kvfMother.fit();
//         double chisq = 1.e3;
        if(err) { 
            setBadVtx(Mother); 
            if (debugDump) {
                printf( "-----  bad vertexing!!! ---- \n\n" );
            }
            return;
        }
//         chisq = kvfMother.chisq();
        makeMother(kvfMother, Mother);
        infoMother.msKvf(Mother.p().m());
        infoMother.chisq(kvfMother.chisq());
        infoMother.chisqKvf(kvfMother.chisq());
        infoMother.cl(kvfMother.cl());
        infoMother.clKvf(kvfMother.cl());
        infoMother.dist2IP(distanceToIP(Mother));
        for (int jM = 0; jM < Mother.nChildren(); ++jM) {
            Particle& Child = Mother.child(jM);
            int lundChild = (int)Child.lund() ;
            bool isTrkChild = isLikeTrk(lundChild);
            if (!isTrkChild) {
                if (!&Child.userInfo()) createUserInfo(Child);
                UserInfo& infoChild = dynamic_cast<UserInfo&>(Child.userInfo());
                double dx = Mother.momentum().decayVertex().x() - Child.momentum().decayVertex().x();
                double dy = Mother.momentum().decayVertex().y() - Child.momentum().decayVertex().y();
                double dz = Mother.momentum().decayVertex().z() - Child.momentum().decayVertex().z();
                double dist2Mother = sqrt(dx*dx + dy*dy + dz*dz);
                infoChild.dist2Mother(dist2Mother);
            }
        }
                
        // check on long-lived particles as criteria to use mass-vertex
        double cTauMother = Ptype(lundMother).cTau();
        // bool useKmvf = false;
//         if (cTauMother>1.e-5) useKmvf =  true;
        if (abs(lundMother)==10431) useKmvf = false;

        infoMother.useKmvf(useKmvf);
        
        if (infoMother.useKmvf()) { 
            // use additional mass-vertex fitter to correct position (?) and LV
            kmassvertexfitter kmvMother;
            if (useBF) kmvMother.magneticField(BF);
            kmvMother.invariantMass(Mother.pType().mass());
            for (int j = 0; j < Mother.nChildren(); ++j)
                addTrack2fit(kmvMother, Mother.child(j));
            err = kmvMother.fit();
            if (err) {
                setBadVtx(Mother);
                return;
            }
            makeMother(kmvMother, Mother);
//             chisq = kmvMother.chisq();
            Mother.momentum().decayVertex(kmvMother.vertex(), kmvMother.errVertex());
            infoMother.chisq(kmvMother.chisq());
            infoMother.cl(kmvMother.cl());
            infoMother.dist2IPmvf(distanceToIP(Mother));
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
            string sind[5] = { "(pcl_1)", "(pcl_2)", "(pcl_3)", "(pcl_4)" , "(pcl_5)" };
            for(int jM = 0; jM<Mother.nChildren(); ++jM) {
                printPclDebug( Mother.child(jM), sind[jM] );
            }
            printf("\n");
        }
}
//***********************************************************
void vertex_staged(vector<Particle>& plist, bool debugDump=false) {
    // fit DecayMother -> DecayChild1 + DecayChild2 + ... + trk1 + trk2 + ...
    for (int i = 0;i < plist.size(); ++i) {
        vertex_staged(plist[i], debugDump);
    }
}
//***********************************************************
void printPi0(vector<Particle> &pi0, string comment="") {
    printf("------  %s Pi0 (%i) -------\n", comment.c_str(), pi0.size());
    double E_HER = BeamEnergy::E_HER();
    double E_LER = BeamEnergy::E_LER();
    
    for (int iPi0=0; iPi0<pi0.size(); iPi0++) {
        Particle &p0 = pi0[iPi0];
        Particle &g1 = pi0[iPi0].child(0);
        Particle &g2 = pi0[iPi0].child(1);
        double msPi0_gg = ( g1.p()+g2.p() ).m();
        double ptot_cm = pStar(p0.p(), E_HER, E_LER).vect().mag();
        printf(" Pi0 (%i)  mass:%7.5f, [px,py,pz]:[%7.4f, %7.4f, %7.4f], p:%6.4f,  p_cm:%6.4f,  Eg(1,2): [%6.4f, %6.4f], ms_gg:%7.5f\n",
        iPi0, p0.mass(), p0.px(), p0.py(), p0.pz(), p0.ptot(), ptot_cm, g1.ptot(), g2.ptot(), msPi0_gg);
    }
}
//***********************************************************
void printTrkPID(vector<Particle> &trkList, string pType, string comment = "") {
    int ip, ibg1, ibg2;
    if ((pType == "pi+") || (pType == "pi-")) {
        ip = 2;
        ibg1 = 3;
        ibg2 = 4;
    }
    else if ((pType == "K+") || (pType == "K-")) {
        ip = 3;
        ibg1 = 2;
        ibg2 = 4;
    }
    else if ((pType == "p+") || (pType == "p-")) {
        ip = 4;
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
//***********************************************************
void printVtxDebug( std::vector<Particle> &plist, string comment="", string comment2="" ){
    int np = plist.size();
    printf("\n----- Vtx (%i) %9s %s ---\n", np, comment.c_str(), comment2.c_str() );
    for (int i = 0; i < np; ++i) {
        UserInfo &info = dynamic_cast<UserInfo&>(plist[i].userInfo());
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

double get_pcm(Particle &particle) {
    //return particle.p().vect().mag();
    //printf("pcm: %2f\n", BeamEnergy::p_cm(particle.p()).vect().mag());
    return BeamEnergy::p_cm(particle.p()).vect().mag();
}

double get_delta_e(Particle &particle) {
    //printf("Beam energy2: %2f\n", BeamEnergy::E_beam2());
    int n_childs = particle.nChildren();
    double summ_e_i = 0;
    for (int i = 0; i < n_childs; ++i) {
        //printf("particle.child(n_childs).p().e(): %2f\n", particle.child(i).p().e());
        summ_e_i += particle.child(i).p().e();
    }
    //printf("summ_e_i: %2f\n", summ_e_i);
    double mbc = summ_e_i - BeamEnergy::E_beam2();//(const int version, const int expno, const int runno);
//     double mbc = summ_e_i - BeamEnergy::E_beam_corr();
    return mbc;
}

double get_m_bc(Particle &particle) {
    int n_childs = particle.nChildren();
    double summ_pi = 0;
    for (int i = 0; i < n_childs; ++i) {
        summ_pi += get_pcm(particle.child(i));
    }
    //printf("Beam energy2: %2f\n", BeamEnergy::E_beam2());
    //printf("Beam energy: %2f\n", BeamEnergy::E_beam_corr());
//     double e_beam = BeamEnergy::E_beam2();
    double e_beam = BeamEnergy::E_beam_corr();
    return sqrt(e_beam * e_beam - summ_pi * summ_pi);
}
//***********************************************************
void clearVectors() {
    // clear all vectors
    for (int i=0; i<nTrk; i++) trkV[i].clear();
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
//     BsStar0.clear();
//     BsStar0bar.clear();
//     Upsilon_5S.clear();
}
//***********************************************************
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
//***********************************************************
double getHelicity(Particle& p, int indDough=0) {
    if (p.nChildren() == 0) return -1;
    
    VectorL mother   = p.p(); 
    VectorL dght     = p.child(indDough).p();
    VectorL grndDght = p.child(indDough).child(0).p();
    double helic = getHelicity(mother, dght, grndDght);
    
// 
    if (debugHel) {
        printf("\n======================== getHelicity ====================");
        printPclDebug( p, "(Mother)" );
        printf("         dght helic: %7.5f\n\n",helic);
    }
    
    return helic;
}

//**********************************************************
void checkAdoptCutMassChisqVtx(Particle& p, double pL, double pR, double maxChisq=1.e3, string status="", int iChild=0) {
    // if maxChisq<0 vtx cuts not used
    // status: "", "massdif"
    UserInfo& info = dynamic_cast<UserInfo&>(p.userInfo());
    double ms = info.msKvf();
//     double ms = p.p().m();
    if (status == "massdif") {
        Particle& pChild = p.child(iChild);
        UserInfo& infoChild = dynamic_cast<UserInfo&>(pChild.userInfo());
//         ms -= p.relation().child(iChild).p().m();
        ms = info.msComb() - infoChild.msComb();
    }
    bool idErase = false;
    double chisqCand = info.chisqKvf();
    double vx = p.momentum().decayVertex().x();
    double vy = p.momentum().decayVertex().y();
    double vz = p.momentum().decayVertex().z();
    if (maxChisq > 0.) {
        idErase = (chisqCand > maxChisq || chisqCand < 0. || abs(vz) > 50. || sqrt(vx*vx + vy*vy) > 30. || ms < pL || ms > pR);
    } else {
        idErase = ms < pL || ms > pR;
    }
    if (idErase) {
        info.isAdoptCut(false);
    }
// printf("----checkAdoptCutMassChisqVtx [%s] --- ms:%f, pL:%f, pR:%f, maxChi2:%f, chisqCand:%f, vx:%f, vy:%f, vz:%f, --- idErase:%i \n",
// pclType((int)p.lund()).c_str(), ms, pL, pR, maxChisq, chisqCand, vx,vy,vz, (int)idErase);
}

void withKaonId(std::vector<Particle>& list,
               const double prob,
               int accq0,
               int tofq0,
               int cdcq0,
	           int ids0,
               int idb0) {

  atc_pid kid(accq0, tofq0, cdcq0, ids0, idb0);
  for (size_t i = 0; i < list.size(); ++i) {
      Particle& P = list[i];
      double probParticle = kid.prob(&(P.mdstCharged()));

	  if (P.mdstCharged() && probParticle >= prob) {

         if (!&P.userInfo()) createUserInfo(P);
         UserInfo& info = dynamic_cast<UserInfo&>(P.userInfo());
         info.probPID(probParticle);
        }
    else {
      list.erase(list.begin() + i);
      --i;
    }
  }
}

void withPionId(std::vector<Particle>& list,
	           const double prob,
	           int accq0,
               int tofq0,
               int cdcq0,
	           int ids0,
               int idb0) {

  atc_pid kid(accq0, tofq0, cdcq0, ids0, idb0);
  for (size_t i = 0; i < list.size(); ++i) {
      Particle& P = list[i];
      double probParticle = kid.prob(&(P.mdstCharged()));
      if (P.mdstCharged() && probParticle < prob) {

         if (!&P.userInfo()) createUserInfo(P);
         UserInfo& info = dynamic_cast<UserInfo&>(P.userInfo());
         info.probPID(probParticle);
        }
    else {
      list.erase(list.begin() + i);
      --i;
    }
  }
}

void checkKaonPionPID(std::vector<Particle>& DssList, double k1MinProb=0.6, double k2MinProb=0.2, double piMaxProb = 0.9) {
    for (size_t iDss = 0; iDss < DssList.size(); iDss++) {
        // First Ds child
        Particle& DssChild1 = DssList[iDss].child(0);
        // Its children
        Particle& DssChild11 = DssChild1.child(0);
        Particle& DssChild12 = DssChild1.child(1);
        // Second Ds child
        Particle& DssChild2 = DssList[iDss].child(1);
        int child1Lund = DssChild1.lund();
        int child2Lund = DssChild2.lund();
        std::string child1Name = pclType(child1Lund);
        std::string child2Name = pclType(child2Lund);
        std::cout << "First Ds child: " << child1Name << " , second Ds child: " << child2Name << std::endl;
        // First Dss child children
        int child11Lund = DssChild11.lund();
        int child12Lund = DssChild12.lund();
        std::string child11Name = pclType(child11Lund);
        std::string child12Name = pclType(child12Lund);
        std::cout << "--- First child of the 1st Ds child: " << child11Name << " , second child of the 1st Ds child: " << child12Name << std::endl;
        
        if (child1Lund == 333) {          // phi -> K+ K-, and the second Dss child is pi
            Particle& kaon1 = DssChild11; // positive -- >= 0.6
            Particle& kaon2 = DssChild12; // negative -- >= 0.2
            Particle& pion  = DssChild2;  // pion     -- <  0.9
            // PID kaon1PID_KPI(kaon1, 3, 1, 5, 3, 2);
            // PID kaon2PID_KPI(kaon2, 3, 1, 5, 3, 2);
            // double kaon1Prob_KPI = kaon1PID_KPI.kIdProb();
            // double kaon2Prob_KPI = kaon2PID_KPI.kIdProb();

            if (!&kaon1.userInfo()) createUserInfo(kaon1);
            UserInfo& info_kaon1 = dynamic_cast<UserInfo&>(kaon1.userInfo());

            if (!&kaon2.userInfo()) createUserInfo(kaon2);
            UserInfo& info_kaon2 = dynamic_cast<UserInfo&>(kaon2.userInfo());

            if (!&pion.userInfo()) createUserInfo(pion);
            UserInfo& info_pion = dynamic_cast<UserInfo&>(pion.userInfo());

            double kaon1Prob_KPI = info_kaon1.probPID();
            double kaon2Prob_KPI = info_kaon2.probPID();
            double pionProb_KPI = info_pion.probPID();
            std::cout << "--- --- Compared to pion, the first kaon prob: " << kaon1Prob_KPI << " ,the second kaon prob: " << kaon2Prob_KPI << std::endl;
            if (kaon1Prob_KPI >= k1MinProb && pionProb_KPI < piMaxProb) {
                if (kaon2Prob_KPI < k2MinProb) {
                    DssList.erase(DssList.begin() + iDss);
                    std::cout << "--- --- --- Delete" << std::endl;
                }
            }
            else if (kaon1Prob_KPI < k1MinProb || pionProb_KPI >= piMaxProb) {
                    DssList.erase(DssList.begin() + iDss);
                    std::cout << "--- --- --- Delete" << std::endl;
            }   
        }
        else if (child1Lund == 313) {     // K*0 -> K+ pi-, and another Dss child is K
            Particle& kaon1 = DssChild11; // positive kaon -- >= 0.6
            Particle& kaon2 = DssChild2;  // negative kaon -- >= 0.2
            Particle& pion = DssChild12;  // nagative pion  -- <  0.9

            // PID kaon1PID_KPI(kaon1, 3, 1, 5, 3, 2);
            // PID kaon2PID_KPI(kaon2, 3, 1, 5, 3, 2);
            // double kaon1Prob_KPI = kaon1PID_KPI.kIdProb();
            // double kaon2Prob_KPI = kaon2PID_KPI.kIdProb();

            if (!&kaon1.userInfo()) createUserInfo(kaon1);
            UserInfo& info_kaon1 = dynamic_cast<UserInfo&>(kaon1.userInfo());

            if (!&kaon2.userInfo()) createUserInfo(kaon2);
            UserInfo& info_kaon2 = dynamic_cast<UserInfo&>(kaon2.userInfo());

            if (!&pion.userInfo()) createUserInfo(pion);
            UserInfo& info_pion = dynamic_cast<UserInfo&>(pion.userInfo());

            double kaon1Prob_KPI = info_kaon1.probPID();
            double kaon2Prob_KPI = info_kaon2.probPID();
            double pionProb_KPI = info_pion.probPID();

            std::cout << "--- --- Compared to pion, first kaon prob: " << kaon1Prob_KPI << ", second kaon prob: " << kaon2Prob_KPI << std::endl;
            if (kaon1Prob_KPI >= k1MinProb && pionProb_KPI < piMaxProb) {
                if (kaon2Prob_KPI < k2MinProb) {
                    DssList.erase(DssList.begin() + iDss);
                    std::cout << "--- --- --- Delete" << std::endl;
                }
            }
            else if (kaon1Prob_KPI < k1MinProb || pionProb_KPI >= piMaxProb) {
                    DssList.erase(DssList.begin() + iDss);
                    std::cout << "--- --- --- Delete" << std::endl;
            }   
        }
        else if (child1Lund == -313) {    // K*0bar -> K- pi+, and another Dss child is K
            Particle& kaon2 = DssChild11; // negative kaon -- >= 0.2
            Particle& kaon1 = DssChild2;  // positive kaon -- >= 0.6
            Particle& pion = DssChild12;  // positive pion  -- <  0.9
            // PID kaon1PID_KPI(kaon1, 3, 1, 5, 3, 2);
            // PID kaon2PID_KPI(kaon2, 3, 1, 5, 3, 2);
            // double kaon1Prob_KPI = kaon1PID_KPI.kIdProb();
            // double kaon2Prob_KPI = kaon2PID_KPI.kIdProb();

            if (!&kaon1.userInfo()) createUserInfo(kaon1);
            UserInfo& info_kaon1 = dynamic_cast<UserInfo&>(kaon1.userInfo());

            if (!&kaon2.userInfo()) createUserInfo(kaon2);
            UserInfo& info_kaon2 = dynamic_cast<UserInfo&>(kaon2.userInfo());

            if (!&pion.userInfo()) createUserInfo(pion);
            UserInfo& info_pion = dynamic_cast<UserInfo&>(pion.userInfo());

            double kaon1Prob_KPI = info_kaon1.probPID();
            double kaon2Prob_KPI = info_kaon2.probPID();
            double pionProb_KPI = info_pion.probPID();

            std::cout << "--- --- Compared to pion, first kaon prob: " << kaon1Prob_KPI << ", second kaon prob: " << kaon2Prob_KPI << std::endl;
            if (kaon1Prob_KPI >= 0.6 && pionProb_KPI < 0.9) {
                if (kaon2Prob_KPI < 0.2) {
                    DssList.erase(DssList.begin() + iDss);
                    std::cout << "--- --- --- Delete" << std::endl;
                }
            }
            else if (kaon1Prob_KPI < 0.6 || pionProb_KPI >= 0.9) {
                    DssList.erase(DssList.begin() + iDss);
                    std::cout << "--- --- --- Delete" << std::endl;
            }   
        }
    }
}

//**********************************************************
void checkAdoptCutMassChisqVtx( std::vector<Particle>& plist, double pL, double pR, 
            double maxChisq=1.e3, string status="", int iChild=0){
    for (int i = 0; i < (int)plist.size(); ++i) {
        checkAdoptCutMassChisqVtx(plist[i], pL, pR, maxChisq, status, iChild);
    }
}
//**********************************************************
int getEvtGenType() {
//    0:"Data", 1:"evtgen-charged", 2:"evtgen-mixed", 3:"evtgen-charm", 4:"evtgen-uds", 5:"evtgen-bsbs", 6:"evtgen-nonbsbs"
    Gen_hepevt_Manager &genMgr  = Gen_hepevt_Manager::get_manager();
    for(std::vector<Gen_hepevt>::iterator i = genMgr.begin(); 
                i != genMgr.end(); ++i) {
        int idh = (*i).idhep();
        if      (abs(idh)==521) {return 1;}  // B+B-
        else if (abs(idh)==511) {return 2;}  // B0B0bar
        else if (abs(idh)==  4) {return 3;}  // ccbar
        else if (abs(idh)<=  3) {return 4;}  // uds
    }
    return 0;
}
//***********************************************************
void evtInfo_dump(BelleTuple* tt, bool debugDump=false) {
    // Event Information
    int expNo, runNo, evtNo;
    getEventInfo(expNo, runNo, evtNo, McFlag); // utility.cc
    const HepPoint3D& ip_position = IpProfile::position(1);
    const HepSymMatrix& ip_error = IpProfile::position_err(1);
    
    int idGenType = getEvtGenType();

    // Event Shape 
    double  r2 = -1.;
    Evtcls_hadron_info_Manager &clsMgr = Evtcls_hadron_info_Manager::get_manager();
    if(clsMgr.count()) r2 = clsMgr[0].R2();
    
    if (debugDump)
        printf("\n        ---- EvtInfo_DUMP --- EvtGenType[%i],  exp:%2i,  run:%2i, evt:%i, ip_position:[%f, %f, %f] ----\n", 
            idGenType, expNo, runNo, evtNo, 
            ip_position.x(), ip_position.y(), ip_position.z());

    tt->column("expn",  expNo ); // Exp #
    tt->column("runn",  runNo ); // Run #
    tt->column("evtn",  evtNo );
    tt->column("ipx",  ip_position.x() );
    tt->column("ipy",  ip_position.y() );
    tt->column("ipz",  ip_position.z() );
    //************* Signal shape ********
    tt->column("r2", r2 );
    //*************  MC  ****************
    tt->column("evtgen", idGenType );
}

VectorL getGenVectorL(int idhPcl) {
    VectorL pclL;
    int ID_Pcl = -1; 
    Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
    for (std::vector<Gen_hepevt>::iterator i = genMgr.begin(); i != genMgr.end(); ++i) {
        int idh = (*i).idhep();
        if (idh == idhPcl) { // Ds+/-
            ID_Pcl = (*i).get_ID();
            pclL.setPx(genMgr[ID_Pcl - 1].PX());
            pclL.setPy(genMgr[ID_Pcl - 1].PY());
            pclL.setPz(genMgr[ID_Pcl - 1].PZ());
            pclL.setE (genMgr[ID_Pcl - 1].E());
        }
    }
    return pclL;
}
//***********************************************************
void pi0_dump(BelleTuple* tt, Particle& p0, string sfx, bool debugDump) {
    double E_HER = BeamEnergy::E_HER();
    double E_LER = BeamEnergy::E_LER();
    // string s_pi_d17 =  " gen_p0_d eg1_p0_d  eg2_p0_d  psr_p0_d  mgg_p0_d ";
    
    Particle& g1 = p0.child(0);
    Particle& g2 = p0.child(1);
    double msPi0_gg = (g1.p() + g2.p()).m();
    double ptot_cm = pStar(p0.p(), E_HER, E_LER).vect().mag();

    bool gen_pi0 = IDhep(p0) == 0 ? false : true;
    bool gen_g1  = IDhep(g1) == 0 ? false : true;
    bool gen_g2  = IDhep(g2) == 0 ? false : true;

    const int nValI = 1; 
    const int nValD = 6; 

    string pclTitI[nValI] = {"gen"};
    int valPclI[nValI]    = {gen_pi0};
    string pclTitD[nValD] = {"eg1",     "eg2",     "psr",   "mgg",    "gg1",  "gg2"};
    double valPclD[nValD] = {g1.ptot(), g2.ptot(), ptot_cm, msPi0_gg, gen_g1, gen_g2};
    
    if (debugDump) {
        printf("  ==== val_dump ==== pi0 (%s): ", sfx.c_str());
        for (int iVal = 0; iVal < nValI; iVal++) 
            printf("(%s:%6i) ", (pclTitI[iVal] + sfx ).c_str(), valPclI[iVal]);
        for (int iVal = 0; iVal < nValD; iVal++) 
            printf("(%s:%6.3f) ", (pclTitD[iVal] + sfx ).c_str(), valPclD[iVal]);
        printf("\n");
    }
    for (int iVal = 0; iVal < nValI; iVal++) 
        tt->column(pclTitI[iVal] + sfx, valPclI[iVal]);
    for (int iVal = 0; iVal < nValD; iVal++) 
        tt->column(pclTitD[iVal] + sfx, valPclD[iVal]);
}

//***********************************************************
void val_dump(BelleTuple* tt, int nValI, int nValD, int* valPclI, double* valPclD, string* pclTitI, string* pclTitD, 
                                                                                            string sfx, bool debugDump) {
    if (debugDump) {
        printf("  ==== val_dump ==== dec:(%s): ", sfx.c_str() );
        for (int iVal=0; iVal<nValI; iVal++) 
            printf("(%s:%6i) ", ( pclTitI[iVal]+sfx ).c_str(), valPclI[iVal]);
        for (int iVal=0; iVal<nValD; iVal++) 
            printf("(%s:%6.3f) ", ( pclTitD[iVal]+sfx ).c_str(), valPclD[iVal]);
        printf("\n");
    }
    for (int iVal = 0; iVal < nValI; iVal++) 
        tt->column(pclTitI[iVal] + sfx, valPclI[iVal]);
    for (int iVal = 0; iVal < nValD; iVal++) 
        tt->column(pclTitD[iVal] + sfx, valPclD[iVal]);
}
// ----------------------------------------------------------
void gen_val_dump(BelleTuple* tt, bool gen_pcl, VectorL pclL, string sfx, bool debugDump) {
    const int nValD = 5; 
    string pclTitD[nValD] = {"ms", "px", "py", "pz", "e"};
    double valPclD[nValD] = {pclL.m(), pclL.px(), pclL.py(), pclL.pz(), pclL.e()};
    if (!gen_pcl) 
        for (int iVal = 0; iVal < nValD; iVal++)
            valPclD[iVal] = -.99;

    if (debugDump) {
        printf("  ==== gen_val_dump ==== dec:(%s): ", sfx.c_str() );
        for (int iVal = 0; iVal < nValD; iVal++) 
            printf("(%s:%6.3f) ", ( pclTitD[iVal]+sfx ).c_str(), valPclD[iVal]);
        printf("\n");
    }
    for (int iVal = 0; iVal < nValD; iVal++) {
//             printf("(%s:%6.3f) ", ( pclTitD[iVal]+sfx ).c_str(), valPclD[iVal]);
        tt->column(pclTitD[iVal]+sfx, valPclD[iVal]);
    }
}
//***********************************************************
void dumpDssChild(BelleTuple* tt, Particle& P, string sfxDs="", bool evtInfoDump=false, 
             bool stDump=true, bool debugDump=false) {
// printf("\n======== dumpDssChild  ========= chg_dss_child:%i, ms_dss_child:%7.3f ) \n",
//                (int)P.lund(), P.p().m() );
    if (evtInfoDump) evtInfo_dump(tt, debugDump);
// printf("---- 1 ----- \n");
    
    int lund = (int)P.lund();
    
    if (!&P.userInfo()) createUserInfo(P);
    UserInfo& info = dynamic_cast<UserInfo&>(P.userInfo());
    if (info.chisqKvf() < 0.) {
        // particle not vertexed yet
        vertex_staged(P, debugDump);
    }

    double massPDG = Ptype(lund).mass();
    double msLimLeft  = massPDG - info.wMass();
    double msLimRight = massPDG + info.wMass();

    checkAdoptCutMassChisqVtx(P, msLimLeft, msLimRight, info.maxChi2());
    if (!info.isAdoptCut()) {
        if (debugDump) {
            printf("\n       *****  DssChild is not adopted ******\n" );
            printPclDebug( P );
            double vx = P.momentum().decayVertex().x();
            double vy = P.momentum().decayVertex().y();
            double vz = P.momentum().decayVertex().z();
            printf("msKvf:%f, msLimLeft:%f, msLimRight:%f, maxChi2:%f, vx:%f, vy:%f, vz:%f \n",
                   info.msKvf(), msLimLeft, msLimRight, info.maxChi2(), vx, vy, vz);
        }
        return;
    }
    
    double chisq  = info.chisqKvf(); // -1.; //
    double cl     = info.clKvf(); // -1.; // 
    double d2m    = info.dist2Mother();
    double msKvf  = info.msKvf(); // dgr.p().m() ; //
    double msComb = info.msComb();
//     double helic = getHelicity(dgr);
    double helic = info.helicity();
    double pidprob = info.probPID();
    
    int ind_ch1;
    int lundChild = P.lund();
    string ChildPcl;
    if (abs(lundChild) == 333) {
        ChildPcl ="phipi";
        ind_ch1  = 1;
    }
    if (abs(lundChild) == 313) {
        ChildPcl = "KsrK";
        ind_ch1  = 2;
    }

    Hep3Vector P3D(P.px(), P.py(), P.pz());
// printf("---- 2 ----- \n");
    
    string dgrSuff = "_ch" + sfxDs;
    int signPcl = (int)P.lund() > 0 ? 1 : -1;
    bool gen_pcl = IDhep(P) == 0 ? false : true;
    double vx = P.momentum().decayVertex().x();
    double vy = P.momentum().decayVertex().y();
    double vz = P.momentum().decayVertex().z();

    
//     double helicChild1 = getHelicity( P );
// printf("---- 3 ----- \n");
    
    const int nValI = 3; 
    const int nValD = 5; 
    int valPclI[nValI] = {signPcl, gen_pcl, ind_ch1};
    double valPclD[nValD] = {msKvf, chisq, P3D.perp(), P3D.phi(), P3D.theta()};
    string pclTitI[nValI] = {"chg", "gen", "ind"};
    string pclTitD[nValD] = {"ms", "chi", "pt", "ph", "th"};

    if (debugDump) {
        printf("\n======== dumpDssChild  ========= (%s) chg_dss_child:%i, ms_dss_child:%7.3f ) \n", 
               ChildPcl.c_str(), (int)P.lund(), msKvf);
        printUserInfo(P);
    }

    val_dump(tt, nValI, nValD, valPclI, valPclD, pclTitI, pclTitD, dgrSuff, debugDump);
    
    if (stDump) tt->dumpData();
}
//-----------------------------------------------------------
void dumpDss(BelleTuple* tt, Particle& P, string sfxDs="", bool evtInfoDump=false, bool stDump=true, bool debugDump=false) {

    if (evtInfoDump) evtInfo_dump(tt, debugDump);
    int lund = (int)P.lund();

    if (!&P.userInfo()) createUserInfo(P);
    UserInfo& info = dynamic_cast<UserInfo&>(P.userInfo());

    if (info.chisqKvf() < 0.) {
        // particle not vertexed yet
        vertex_staged(P, debugDump, true);
    }
    
    double massPDG = Ptype(lund).mass();
    double msLimLeft  = massPDG - info.wMass();
    double msLimRight = massPDG + info.wMass();

    // validation of chi2 value
    checkAdoptCutMassChisqVtx(P, msLimLeft, msLimRight, info.maxChi2());
    if (!info.isAdoptCut()) {
        if (debugDump) {
            printf("\n       *****  Dss is not adopted ******\n" );
            printPclDebug(P);
            double vx = P.momentum().decayVertex().x();
            double vy = P.momentum().decayVertex().y();
            double vz = P.momentum().decayVertex().z();
            printf("msKvf:%f, msLimLeft:%f, msLimRight:%f, maxChi2:%f, vx:%f, vy:%f, vz:%f \n",
                   info.msKvf(), msLimLeft, msLimRight, info.maxChi2(), vx, vy, vz);
        }
        return;
    }
    double chisq  = info.chisqKvf(); // -1.; // ???
    double msKvf  = info.msKvf();    // dgr.p().m() ; // 
    double helic = info.helicity();


    
    Particle& Child = P.child(0);    // phi0 (--> K+ K-) or K* (--> K+ pi-), or K*bar (--> K- pi+)
    // Particle& Child = P.child(1);    // pi or K
    Hep3Vector P3D(P.px(), P.py(), P.pz());
    Hep3Vector Child3D(Child.px(), Child.py(), Child.pz());
    
    string dgrSuff = "_ds" + sfxDs, genDgrSuff = "_ds_t" + sfxDs;
    int signPcl = (int)P.lund() > 0 ? 1 : -1;
    bool gen_pcl = IDhep(P) == 0 ? false : true;
    double vx = P.momentum().decayVertex().x();
    double vy = P.momentum().decayVertex().y();
    double vz = P.momentum().decayVertex().z();
//     double chisq = -1.;
    
    double helicChild1 = getHelicity(P);
    
    const int nValI = 2; 
    const int nValD = 6; 
    int valPclI[nValI]    = {signPcl, gen_pcl};
    double valPclD[nValD] = {msKvf, chisq, P3D.perp(), P3D.phi(), P3D.theta(), helicChild1};
    string pclTitI[nValI] = {"chg", "gen"};
    string pclTitD[nValD] = {"ms", "chi", "pt", "ph", "th", "hel"};

    VectorL dssL = getGenVectorL(IDhep(P));

    if (debugDump) {
        printf("\n======== Dss  ========= sfx:%s, chg_ds:%i, gen_ds:%i, ms_ds:%7.3f , child ( ms:%7.3f ) \n", 
               genDgrSuff.c_str(), signPcl, gen_pcl, msKvf, Child.p().m() );
        printUserInfo(P);
    }
    tt->column("hel_ch" + sfxDs, helicChild1);
    val_dump(tt, nValI, nValD, valPclI, valPclD, pclTitI, pclTitD, dgrSuff, debugDump);
    gen_val_dump(tt, gen_pcl, dssL, genDgrSuff, debugDump);
    dumpDssChild(tt, Child, sfxDs, false, false, debugDump);
    if (stDump) tt->dumpData();
}
//***********************************************************
void dumpDs2317(BelleTuple* tt, Particle& P, string sfxDs="", bool evtInfoDump=false, bool stDump=true, bool debugDump=false) {

    if (evtInfoDump) evtInfo_dump(tt,debugDump);

    if (!&P.userInfo()) createUserInfo(P);
    UserInfo& info = dynamic_cast<UserInfo&>(P.userInfo());
    if (info.chisqKvf() < 0.) {
        // particle not vertexed yet
        vertex_staged(P, debugDump);
    }
    
    double chisq        = info.chisqKvf(); // -1.; //
    double msKvf        = info.msKvf(); // dgr.p().m() ; // 
    double helic_2317   = -1;
    

    Particle& Child     = P.child(0);
    Particle& pi0_2317  = P.child(1);
    UserInfo& infoChild = dynamic_cast<UserInfo&>(Child.userInfo());
    double msKvfChild   = infoChild.msKvf();
    double E_HER        = BeamEnergy::E_HER();
    double E_LER        = BeamEnergy::E_LER();
    double psr_d17      = pStar(P.p(), E_HER, E_LER).vect().mag();
    double px_d17       = P.px();
    double py_d17       = P.py();
    double pz_d17       = P.pz();
    double vx_d17       = P.momentum().decayVertex().x();
    double vy_d17       = P.momentum().decayVertex().y();
    double vz_d17       = P.momentum().decayVertex().z();


    Hep3Vector P3D(px_d17, py_d17, pz_d17);
    int chg_d17 = (int)P.lund() >  0 ? 1 : -1;
    int gen_d17 = (int)IDhep(P) == 0 ? 0 :  1;
    
    const int nValI = 2;
    const int nValD = 13; 
    int valPclI[nValI] = {chg_d17, gen_d17};
    double valPclD[nValD] = {msKvf, chisq, P3D.perp(), psr_d17, P3D.phi(), P3D.theta(), helic_2317, px_d17, py_d17, pz_d17, vx_d17, vy_d17, vz_d17};
    string pclTitI[nValI] = {"chg", "gen"};
    string pclTitD[nValD] = {"ms", "chi", "pt", "psr", "ph", "th", "hel", "px", "py", "pz", "vx", "vy", "vz"};

    string dgrSuff="_d17", genDgrSuff="_d17_t";
    
    if (debugDump) {
        printf("\n======== Ds(2317)  ========= chg_2317:%i, gen_2317:%i, ms_2317:%7.3f , child ( ms:%7.3f ) \n", 
               chg_d17, gen_d17, msKvf, msKvfChild);
        printUserInfo(P);
    }

// string s_2317_gen = " ms_d17_t px_d17_t py_d17_t pz_d17_t e_d17_t ";
    VectorL ds17L = getGenVectorL(IDhep(P));

    val_dump( tt, nValI, nValD, valPclI, valPclD, pclTitI, pclTitD, dgrSuff, debugDump);
    gen_val_dump(tt, gen_d17, ds17L, genDgrSuff, debugDump);
    pi0_dump(tt, pi0_2317, "_p0_d", debugDump);
    
    dumpDss(tt, Child, sfxDs, false, false, debugDump);

    if (stDump) tt->dumpData();
}   
//***********************************************************
void dumpBs0(BelleTuple* tt, Particle& P, bool evtInfoDump=false, 
                 bool stDump=true, bool debugDump=true) {

    if (evtInfoDump) evtInfo_dump(tt, debugDump);

    if (!&P.userInfo()) createUserInfo(P);
    UserInfo& info = dynamic_cast<UserInfo&>(P.userInfo());
    if (info.chisqKvf() < 0.) {
        // particle not vertexed yet
        vertex_staged(P, debugDump);
    }
    double chisq  = info.chisqKvf(); // -1.; // 
    double msKvf  = info.msKvf(); // dgr.p().m() ; // 

    Particle& Dss_Bs0 = P.child(0);
    Particle& Dss2317_Bs0 = P.child(1);
    Particle& pi0_Bs0 = P.child(2);
    Particle& pi0_Ds2317 = Dss2317_Bs0.child(1);
    
    Hep3Vector P3D(P.px(), P.py(), P.pz());
    
    double vx = P.momentum().decayVertex().x();
    double vy = P.momentum().decayVertex().y();
    double vz = P.momentum().decayVertex().z();
    int gen_bs = (int)IDhep(P) == 0 ? 0 : 1;
    int chg_bs = (int)P.lund() > 0 ? 1 : -1;
    

    VectorL pB  = pStar(P.p());
    double de_bs_old = pB.e() - Benergy();
    double de_bs = pB.e() - BeamEnergy::E_beam_corr();
    double mbc_bs_old = beamEnergyConstraint(P);
    
    double energyEl = BeamEnergy::E_HER();
    double energyPos = BeamEnergy::E_LER();
    double angle = BeamEnergy::Cross_angle();
    double mbc_bs = beamEnergyConstraint(P, energyEl, energyPos, angle);
//     double mbc_bs = get_m_bc( P );
    
    const int nValI = 2; 
    const int nValD = 10; 
    int valPclI[nValI] = {chg_bs, gen_bs};
    double valPclD[nValD] = {msKvf, chisq, vx, vy, vz, P3D.perp(), P3D.phi(), P3D.theta(), mbc_bs, de_bs};
    string pclTitI[nValI] = {"chg", "gen"};
    string pclTitD[nValD] = {"ms", "chi", "vx", "vy", "vz", "pt", "ph", "th", "mbc", "de"};
    
    
    if (debugDump) {
        printf("\n\n======== Bs0  ========= chg_bs:%i,  gen_bs:%i, ms_bs:%7.3f , de_bs: %7.3f , de_bs_old: %7.3f, mbc_bs: %7.3f , mbc_bs_old: %7.3f \n", 
               chg_bs, gen_bs, msKvf, de_bs, de_bs_old, mbc_bs, mbc_bs_old );
        printUserInfo( P );
    }
    
    val_dump(  tt, nValI, nValD, valPclI, valPclD, pclTitI, pclTitD, "_bs", debugDump );
    pi0_dump(  tt, pi0_Bs0,    "_p0_b",            debugDump);
    pi0_dump(  tt, pi0_Ds2317, "_p0_d",            debugDump);
    dumpDss(   tt, Dss_Bs0,     "1", false, false, debugDump);
    dumpDs2317(tt, Dss2317_Bs0, "2", false, false, debugDump);

    if (stDump) tt->dumpData();
}
//***********************************************************
void printVectPclWithChildren(std::vector<Particle>& pcl, string tit="") {
    if (pcl.size() > 0) {
        printf("    ---- %s [%i] ----- \n", tit.c_str(), pcl.size() );
        for (int i = 0; i < pcl.size(); i++) {
            printPclDebug(pcl[i]);
            for (int ich = 0; ich < pcl[i].nChildren(); ich++) 
                printPclDebug(pcl[i].child(ich));
        }
    }
    printf("\n");
}
//***********************************************************
void Reco::event(BelleEvent *evptr, int *status) {
    *status = 0;
    bool debugTrkPID = false;
    bool debugVtx = false;
    bool debugPi0 = false;
    bool debugPhiKsr = false;
    bool debugDss = false;
    bool debugDss_2317 = false;
    bool debugBs0 = false;
    
    bool debugDumpDss=false;
    bool debugDump2317=false;
    bool debugDumpBs0=false;

    const HepPoint3D &ip_position = IpProfile::position(1);
    const HepSymMatrix &ip_error = IpProfile::position_err(1);
//  Gen_hepevt_Manager& hepevt = Gen_hepevt_Manager::get_manager();

    // Event Information
    int expNo = 0, runNo = 0, evtNo = 0;
    getEventInfo(expNo, runNo, evtNo, McFlag); // utility.cc
//     printf("\n---- exp:%2i,  run:%2i, evt:%i, ip_position:[%f, %f, %f] ----\n", 
//     expNo,runNo,evtNo, ip_position.x(), ip_position.y(), ip_position.z() );
    printf("\n\n***************** exp:%2i,  run:%2i, evt:%i *********************\n", expNo, runNo, evtNo);
  
    // Event Shape
    double r2 = -1.;

    Evtcls_hadron_info_Manager &clsMgr = Evtcls_hadron_info_Manager::get_manager();
    if (clsMgr.count()) r2 = clsMgr[0].R2();

    ////////////////////  make charged particles - tracks //////////////////
    // makes Kaon and Pion from MdstCharged w/o cut. 1 : w/ good_charged, 0 : w/o
    makeKPi(trkV[2], trkV[3], trkV[0], trkV[1], 0); //k_p, k_m, pi_p, pi_m

    if (debugTrkPID) {
        printf("\n");
        printf("---- exp:%2i,  run:%2i, evt:%i ----\n", expNo, runNo, evtNo);
        for (int itr = 0; itr < nTrk; itr++)
            printTrkPID(trkV[itr], trkTit[itr], "before");
    }

    double minProbPID_Kn = 0.2, minProbPID_Kp = 0.6,
    minProbProtPID = 0.0, maxProbEl = 1.0, maxProbPion = 0.9;  // preselected
// If each plist element is not within atc_pID.prob >= prob,
// its element is removed from plist.
    withKaonId(trkV[2], minProbPID_Kp, 3, 1, 5, 3, 2);      // K+ vs bg pi
    withKaonId(trkV[3], minProbPID_Kn, 3, 1, 5, 3, 2);      // K- vs bg pi
    withKaonId(trkV[2], minProbPID_Kp, 3, 1, 5, 3, 4);       // K+ vs bg p
    withKaonId(trkV[3], minProbPID_Kn, 3, 1, 5, 3, 4);       // K- vs bg p

// If each plist element is not within atc_pID.prob < prob,
// its element is removed from plist.
    withPionId(trkV[0], maxProbPion, 3, 1, 5, 2, 3); // pi+ vs bg K
    withPionId(trkV[1], maxProbPion, 3, 1, 5, 2, 3); // pi- vs bg K
    withPionId(trkV[0], maxProbPion, 3, 1, 5, 2, 4);  // pi+ vs bg p
    withPionId(trkV[1], maxProbPion, 3, 1, 5, 2, 4);  // pi- vs bg p


// If each plist element is not associated with rphi & z-svd hits
// whose number is equal to or larger than nRSvdHit and nZSvdHit,
// its element is removed from plist.
    for (int itr = 0; itr < nTrk; itr++) {
        withSVD2(trkV[itr], 1, 1); // nRSvdHit, nZSvdHit
        withdRdZcut(trkV[itr], ip_position, 0.5, 3.0);
    }

    if (debugTrkPID) {
        for (int itr = 0; itr < nTrk; itr++)
            printTrkPID(trkV[itr], trkTit[itr], "after");
    }

    // =================   WORKING WITH Gamma CANDIDATES ================= //
    makeGamma(gammaV);
    withPCut(gammaV, eGammaMin);

    if(useVTX) {
        for (std::vector<Particle>::iterator it = gammaV.begin(); it != gammaV.end(); ++it) {
            setGammaError(*it, ip_position, ip_error);
        }
    }

    // Match gamma candidates with their genhep info
    if (McFlag) {
        setGenHepInfoG(gammaV);
    }

    createUserInfo(gammaV);


    // =================   WORKING WITH PI0 CANDIDATES ================= // 
    makePi0(pi0);
    withPi0GammPCut(pi0, minPi0GammP);
    withPi0MassGamGamCut(pi0, wMassPi0GG);
    setPi0Error(pi0);

    if (debugPi0) {
        printf(" pi0[%i]  \n", pi0.size() );
//         printPi0( pi0 );
    }

    // Match candidates with genhep info 
    if (McFlag) {
        for (int itr = 0; itr < nTrk; ++itr) setGenHepInfoF(trkV[itr]);
        setGenHepInfoP(pi0);
    }    
    
    combination(phi0,    Ptype("PHI"),  trkV[2],   trkV[3], dM_V0); // k_p, k_m
    combination(Ksr0,    Ptype("K*0"),  trkV[2],   trkV[1], dM_Ksr0); // k_p, pi_m
    combination(Ksr0bar, Ptype("K*B"),  trkV[3],   trkV[0], dM_Ksr0); // k_m, pi_p
    setGenHepInfoT(phi0);
    setGenHepInfoT(Ksr0);
    setGenHepInfoT(Ksr0bar);
    
//     if(useVTX) {
//         vertex_staged( phi0, debugVtx );
//         vertex_staged( Ksr0, debugVtx );
//         vertex_staged( Ksr0bar, debugVtx );
//     }
    
    if (debugPhiKsr) {
        printf("      phi0[%i]  \n", phi0.size());
        printf("      Ksr0[%i]  \n", Ksr0.size());
        printf("   Ksr0bar[%i]  \n", Ksr0bar.size());
    }

    combination(Dss_p, Ptype("DS+"),    phi0,     trkV[0], dM_Dss); // phi0, pi_p
    combination(Dss_m, Ptype("DS-"),    phi0,     trkV[1], dM_Dss); // phi0, pi_m
    combination(Dss_p, Ptype("DS+"),    Ksr0bar,  trkV[2], dM_Dss); // K*0bar (K-pi+), k_p
    combination(Dss_m, Ptype("DS-"),    Ksr0,     trkV[3], dM_Dss);  // K*0 (K+pi-), k_m
    setGenHepInfoT(Dss_p);
    setGenHepInfoT(Dss_m);

    // checkKaonPionPID(Dss_m);
    // checkKaonPionPID(Dss_p);
    
//     if(useVTX) {
//         vertex_staged( Dss_p, debugVtx );
//         vertex_staged( Dss_m, debugVtx );
//     }
    
    if (debugDss) {
        if (Dss_p.size() + Dss_m.size() > 0) {
            printf("\n **** debugDss ****** exp:%2i,  run:%2i, evt:%i,  *******\n", expNo, runNo, evtNo);
            printVectPclWithChildren(Dss_p, "Ds+");
            printVectPclWithChildren(Dss_m, "Ds-");
        }
    }


    combination(Dss_m_2317, Ptype(-10431), Dss_m, pi0, dM_2317);
    combination(Dss_p_2317, Ptype( 10431), Dss_p, pi0, dM_2317);
//     combination( Dss_m_2317, Ptype(-10431), Dss_m, pi0, M_2317_min, M_2317_max );
//     combination( Dss_p_2317, Ptype( 10431), Dss_p, pi0, M_2317_min, M_2317_max );
    setGenHepInfoT(Dss_m_2317);
    setGenHepInfoT(Dss_p_2317);
    
//     if(useVTX) {
//         vertex_staged( Dss_p_2317, debugVtx );
//         vertex_staged( Dss_m_2317, debugVtx );
//     }
    
    if (debugDss_2317) {
        if (Dss_p_2317.size() + Dss_m_2317.size() > 0) {
            printf("\n **** debugDss_2317 ****** exp:%2i,  run:%2i, evt:%i,  *******\n", expNo, runNo, evtNo);
            printVectPclWithChildren(Dss_p_2317, "DsJ(2317)+");
            printVectPclWithChildren(Dss_m_2317, "DsJ(2317)-");
        }
    }

    combination(Bs0, Ptype(531), Dss_m, Dss_p_2317, pi0, dM_Bs0);
    combination(Bs0, Ptype(531), Dss_p, Dss_m_2317, pi0, dM_Bs0);
    // combination(Bs0bar, Ptype(-531), Dss_m, Dss_p_2317, pi0, dM_Bs0);
    // combination(Bs0bar, Ptype(-531), Dss_p, Dss_m_2317, pi0, dM_Bs0);
    setGenHepInfoT(Bs0);
    setGenHepInfoT(Bs0bar);
    
//     if(useVTX) {
//         vertex_staged( Bs0, debugVtx );
//         vertex_staged( Bs0bar, debugVtx );
//     }
    
//     combination( BsStar0,    Ptype( 533), Bs0,    gammaV, dM_Bs0 );
//     combination( BsStar0bar, Ptype(-533), Bs0bar, gammaV, dM_Bs0 );
//     setGenHepInfoT(BsStar0);
//     setGenHepInfoT(BsStar0bar);

//     combination( Upsilon_5S, Ptype(9000553),  BsStar0, BsStar0bar, dM_Bs0 );
    
    // ----------------------------  Dump  ---------------------------
    //   Dss
    string sfxDs = "";
    if (stDumpDss) {
        for (int iEvt=0; iEvt < Dss_p.size(); iEvt++) 
                dumpDss(TP_Dss, Dss_p[iEvt], sfxDs, true, stDumpDss, debugDumpDss);
        for (int iEvt=0; iEvt < Dss_m.size(); iEvt++) 
                dumpDss(TP_Dss, Dss_m[iEvt], sfxDs, true, stDumpDss, debugDumpDss);
    }    
    //  Dss(2317)
    if (stDump2317) {
        for (int iEvt=0; iEvt < Dss_p_2317.size(); iEvt++) 
                dumpDs2317(TP_Dss_2317, Dss_p_2317[iEvt], sfxDs, true, stDump2317, debugDump2317);
        for (int iEvt=0; iEvt < Dss_m_2317.size(); iEvt++) 
                dumpDs2317(TP_Dss_2317, Dss_m_2317[iEvt], sfxDs, true, stDump2317, debugDump2317);
    }
    //  Bs0
    if (stDumpBs0) {
        for (int iEvt=0; iEvt < Bs0.size(); iEvt++) 
                dumpBs0(TP_Bs0, Bs0[iEvt], true, stDumpBs0, debugDumpBs0);
        for (int iEvt=0; iEvt < Bs0bar.size(); iEvt++)
                dumpBs0(TP_Bs0, Bs0bar[iEvt], true, stDumpBs0, debugDumpBs0);
    }    
//     printf("---------  clearVectors (final)   ---------------\n");
    clearVectors();
//     printf("------------------  Reco end --------------------  \n");

}

#if defined(BELLE_NAMESPACE)
}//namespace Belle
#endif

