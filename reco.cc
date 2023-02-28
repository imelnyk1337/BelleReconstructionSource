#include <unordered_map>
#include <utility>
#include <map>
#include <cmath>
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
#include "tables/fullrecon.h"re
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

namespace Belle {

    //***********************************************************
    std::string getParticleType(int id) {
        const std::unordered_map<int, std::string> particle_ids {
            {   0, "--"           }, {   1, "d"             }, {  -1, "anti-d"         },
            {   2, "u"            }, {  -2, "anti-u"        }, {   3, "s"              },
            {  -3, "anti-s"       }, {   4, "c"             }, {  -4, "anti-c"         },
            {  11, "el+"          }, { -11, "el-"           }, {  12, "nu_e"           },
            { -12, "anti-nu_e"    }, {  13, "mu+"           }, { -13, "mu-"            },
            {  14, "nu_mu"        }, { -14, "anti-nu_mu"    }, { 15, "tau+"            },
            { -15, "tau-"         }, {  21, "g"             }, {  22, "gamma"          },
            {  23, "Z0"           }, {  24, "W+"            }, { -24, "W-"             },
            {  25, "Higgs0"       }, {  32, "Z'0"           }, {  33, "Z''0"           },
            {  34, "W'+",         }, { -34, "W'-"           }, {  35, "Higgs'0"        },
            {  36, "A0"           }, {  37, "Higgs+"        }, { -37, "Higgs-"         },
            { 111, "pi0"          }, { 113, "RHO0"          }, { 130, "K_L0"           },
            { 211, "pi+"          }, { -211, "pi-"          }, { 213, "RHO+"           },
            { -213, "RHO-"        }, { 221, "eta"           }, { 223, "omega"          },
            { 310, "K_S0"         }, { -311, "anti-K0"      }, { 311, "K0"             },
            { 313, "K*0"          }, { -313, "anti-K*0"     }, { 321, "K+"             },
            { -321, "K-"          }, { 331, "eta'"          }, { 323, "K*+"            },
            { -323, "K*-"         }, { 333, "phi"           }, { -411, "D-"            },
            { 411, "D+"           }, { -413, "D*-"          }, { 413, "D*+"            },
            { -415, "D_2*-"       }, { 415, "D_2*+"         }, { -421, "anti-D0"       },
            { 421, "D0"           }, { -431, "DS-"          }, { 431, "DS+"            },
            { -423, "anti-D*0"    }, { 423, "D*0"           }, { -425, "anti-D_2*0"    },
            { 425, "D_2*0"        }, { -433, "D*S-"         }, { 433, "D*S+"           },
            { 441, "eta_c"        }, { 443, "J/psi"         }, {-511, "anti-B0"        },
            { 511, "B0"           }, {-521, "B-"            }, { 521, "B+"             },
            { 553, "Upsilon"      }, {2112, "n0"            }, {-2112, "anti-n0"       },
            {2212, "p+"           }, {-2212, "anti-p-"      }, {2224, "Delta++"        },
            {-2224, "anti-Delta--"}, {3122, "Lam0"          }, {-3122, "Lam0bar"       },
            {-4122, "LamC-"       }, { 4122, "LamC+"        }, {-4214, "anti-Sigma_c*-"},
            { 4214, "Sigma_c*+"   }, {10022, "vpho"         }, {10221, "f_0"           },
            {-10413, "D_1-"       }, { 10413, "D_1+"        }, {-10423, "anti-D_10"    },
            { 10423, "D_10"       }, {-10431, "D_s0*-"      }, { 10431, "D_s0*+"       },
            {-10433, "D_s1-"      }, { 10433, "D_s1+"       }, {-20213, "a_1-"         },
            { 20213, "a_1+"       }, {-20413, "D'_1-"       }, { 20413, "D'_1+"        },
            {-20423, "anti-D'_10" }, { 20423, "D'_10"       }, {-20433, "D'_s1-"       },
            { 20433, "D'_s1+"     }, { 300553, "Upsilon(4S)"}, { 9000553, "Upsilon(5S)"}};

        auto it = particle_ids.find(id);
        return (it != particle_ids.end() ? it->second : "pcl");
    }
    // ***********************************************************
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
    }

    void setGammaError(std::vector<Particle>& p) {
      for(std::vector<Particle>::iterator i = p.begin(); i != p.end(); ++i)
        setGammaError(*i);
    }
    // ***********************************************************
    int getIdHep(Particle& particle) {
        if (!particle.genHepevt()) return 0;
        return particle.genHepevt().idhep();
    }
    // ***********************************************************
    bool isLikeTrk(int lund, const std::string& chrgType = "Charged_pi0") {
        const int absLund = std::abs(lund);
        const bool isCharged = (absLund == 11) || (absLund == 13) || (absLund == 211) || (absLund == 321) || (absLund == 2212);
        const bool isNeutral = !isCharged && (absLund == 111 || absLund == 22);  // pi0 & gamma
        const bool isOnlyCharged = chrgType == "onlyCharged";

        return (isOnlyCharged && isCharged) || (!isOnlyCharged && (isCharged || isNeutral));
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
                dR = std::fabs(h.dr());
                dZ = std::fabs(h.dz());

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

        int id_lund = std::abs(particle.lund());
        bool use_tube = false;
        double w_mass = dM_Bs0;
        bool use_kmvf = false;
        std::string vertex_mode = "vertex-fit";

        switch (id_lund) {
            case 531:  // Bs
                use_tube = true;
                w_mass = dM_Bs0;
                break;

            case 333:  // Phi0
                w_mass = dM_V0;
                break;

            case 313:  // K*0
                w_mass = dM_Ksr0;
                break;

            case 431:  // Ds
                use_kmvf = true;
                w_mass = dM_Dss;
                vertex_mode = "mass-constraint-fit";
                break;

            case 10431:  // D_sJ(2317)
                w_mass = dM_2317;
                break;

            case 111:  // pi0
                use_kmvf = true;
                w_mass = wMassPi0GG;
                vertex_mode = "mass-constraint-fit";
                break;

            default:
                break;

        }

        info.msComb(particle.p().m());
        info.wMass(w_mass);
        info.useTube(use_tube);
        info.useKmvf(use_kmvf);
        info.isAdoptCut(true);
        info.chisqKvf(-1.);
        info.chisqKmvf(-1.);
        info.chisqKvfNdf(-1.);
        info.chisqKmvfNdf(-1.);
        info.helicity(-1.);
        info.vertexMode(vertex_mode);

        particle.userInfo(info);
    }




    // ***********************************************************
    void createUserInfo(std::vector<Particle>& p_list) {
        std::vector<Particle>::iterator particle;
        for (particle = p_list.begin(); particle != p_list.end(); ++particle) {
            if (!&(*particle).userInfo()) createUserInfo(*particle);
        }
    }
    // ***********************************************************
    void setIncorrectVertexFitResult(Particle& Particle) {
        HepPoint3D vtx(999., 999., 999.);
        HepSymMatrix errVtx(3, 0);
        Particle.momentum().decayVertex(vtx, errVtx);
        UserInfo& info = dynamic_cast<UserInfo&>(Particle.userInfo());
        info.isAdoptCut(false);
        Particle.userInfo(info);
    }
    // ***********************************************************
    double distanceToIP(Particle& Particle) {
        const HepPoint3D& ip_position = IpProfile::position(1);
        Hep3Vector vtxIP(IpProfile::position(1).x(), IpProfile::position(1).y(), 0.0);
        Hep3Vector vtxD(Particle.momentum().decayVertex().x(), Particle.momentum().decayVertex().y(), 0.0);
        Hep3Vector vtx_D_IP = vtxD - vtxIP;
        return vtx_D_IP.mag();
    }
    // ***********************************************************
    void makeRecursiveVertexFit(Particle& Mother, bool debugDump = false, bool useKmvf = false,
                       bool addBeam = true) {
        /*
         * Fitting DecayMother -> DecayChild1 + DecayChild2 + ... + trk1 + trk2 + ...
         */
        int expNo, runNo, evtNo;
        getEventInfo(expNo, runNo, evtNo, McFlag); // utility.cc
        int lund_mother = int(Mother.lund());
        if (evtNo == 470) {
            std::cout << "**1, " << lund_mother  << std::endl;
        }
        if (!&Mother.userInfo()) createUserInfo(Mother);
        UserInfo& info_mother = dynamic_cast<UserInfo&>(Mother.userInfo());
        if (evtNo == 470) {
            std::cout << "**2, " << lund_mother << std::endl;
        }
        kvertexfitter kvf_mother;
        if (useBF) kvf_mother.magneticField(BF);
        std::size_t kvf_error = 0;
        std::size_t kmm_error = 1;

        for (int j_child = 0; j_child < Mother.nChildren(); ++j_child) {
            Particle& child = Mother.child(j_child);
            int lund_child   = int(child.lund());
            bool is_child_trk = isLikeTrk(lund_child);
            if (!is_child_trk) {
                if (!&child.userInfo()) createUserInfo(child);
                UserInfo& info_child = dynamic_cast<UserInfo&>(child.userInfo());
                double chi2_child = info_child.chisqKvf();
                if (chi2_child < 0.) {
                    makeRecursiveVertexFit(child, debugDump);
                    }
                }
            addTrack2fit(kvf_mother, child);
            }
        if (addBeam && Belle::IpProfile::usable()) addBeam2fit(kvf_mother);
    //    if (info_mother.useTube()) addTube2fit(kvf_mother);
        kvf_error = kvf_mother.fit();

        // If OK, 0 is returned
        if (kvf_error) {
            setIncorrectVertexFitResult(Mother);
//            return;
        }

        kmm_error = makeMother(kvf_mother, Mother);
        // if OK, 1 is returned
        if (!kmm_error) {
            setIncorrectVertexFitResult(Mother);
//            return;
        }
        Mother.momentum().decayVertex(kvf_mother.vertex(), kvf_mother.errVertex());
        info_mother.msKvf(Mother.p().m());
        info_mother.chisqKvf(kvf_mother.chisq());
        info_mother.chisqKvfNdf(kvf_mother.chisq() / kvf_mother.dgf());
        info_mother.clKvf(kvf_mother.cl());
        info_mother.dist2IP(distanceToIP(Mother));

        for (int j_child = 0; j_child < Mother.nChildren(); ++j_child) {
            Particle& child = Mother.child(j_child);
            int lund_child = int(child.lund());
            bool is_child_trk = isLikeTrk(lund_child);
            if (!is_child_trk) {
                UserInfo& infoChild = dynamic_cast<UserInfo&>(child.userInfo());
                double dx = Mother.momentum().decayVertex().x() - child.momentum().decayVertex().x();
                double dy = Mother.momentum().decayVertex().y() - child.momentum().decayVertex().y();
                double dz = Mother.momentum().decayVertex().z() - child.momentum().decayVertex().z();
                double dist2Mother = std::sqrt(dx * dx + dy * dy + dz * dz);
                infoChild.dist2Mother(dist2Mother);
                }
            }

        if (evtNo == 470) {
            std::cout << "**3, "  << lund_mother << std::endl;
            }
        // Checking for long-living particles as a criteria of mass-constraint fit
        if ((info_mother.useKmvf() || useKmvf) && info_mother.chisqKmvf() < 0.) {
            // use additional mass-vertex fitter to correct the position (?) and LV
            kmassvertexfitter kmvf_mother;
            if (useBF) kmvf_mother.magneticField(BF);
            int kmvf_error = 0;
            if (useBF) kmvf_mother.magneticField(BF);
            kmvf_mother.invariantMass(Mother.pType().mass());
            for (int j_child = 0; j_child < Mother.nChildren(); ++j_child)
                addTrack2fit(kmvf_mother, Mother.child(j_child));

            kmvf_error = kmvf_mother.fit();
            if (kmvf_error) {
                setIncorrectVertexFitResult(Mother);
            }
            kmm_error = makeMother(kmvf_mother, Mother);
            if (!kmm_error) {
                setIncorrectVertexFitResult(Mother);
            }
            if (evtNo == 470) {
                std::cout << "**4, "  << lund_mother << std::endl;
            }
            Mother.momentum().decayVertex(kmvf_mother.vertex(), kmvf_mother.errVertex());
            info_mother.clKmvf(kmvf_mother.cl());
            info_mother.dist2IPKmvf(distanceToIP(Mother));
            info_mother.msKmvf(Mother.p().m());
            info_mother.chisqKmvf(kmvf_mother.chisq());
            info_mother.chisqKmvfNdf(kmvf_mother.chisq() / kmvf_mother.dgf());

        }
        if (evtNo == 470) {
            std::cout << "**5, "  << lund_mother << std::endl;
        }
        Mother.userInfo(info_mother);
    }
    // ***********************************************************
    void makeRecursiveVertexFit(vector<Particle>& p_list, bool debugDump = false, bool useKmvf = false) {
        // fit DecayMother -> DecayChild1 + DecayChild2 + ... + trk1 + trk2 + ...
        for (size_t i = 0; i < p_list.size(); ++i) {
            makeRecursiveVertexFit(p_list[i], debugDump, useKmvf);
        }
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
        return helic;
    }
    // ***********************************************************

    double getHelicity(Particle& p, int indDough = 0) {
        if (p.nChildren() == 0) return -1;

        VectorL mother   = p.p();
        VectorL dght     = p.child(indDough).p();
        VectorL grndDght = p.child(indDough).child(0).p();
        double helic = getHelicity(mother, dght, grndDght);
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
                info.probPid(probParticle);
                P.userInfo(info);
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
         * By default, ids0 = 3 and idb0 = 2 (K against pi).
         * */

        atc_pid kid(accq0, tofq0, cdcq0, ids0, idb0);

        for (size_t i = 0; i < p_list.size(); ++i) {
            Particle& P = p_list[i];
            double probParticle = kid.prob(&(P.mdstCharged()));
            if (P.mdstCharged() && probParticle < prob) { // a big difference compated to withKaonId

                if (!&P.userInfo()) createUserInfo(P);
                UserInfo& info = dynamic_cast<UserInfo&>(P.userInfo());
                info.probPid(probParticle);
                P.userInfo(info);
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
        Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
        int status = 0;

        for (std::vector<Gen_hepevt>::iterator i = genMgr.begin(); i != genMgr.end(); ++i) {
            int idhep = i->idhep();
            switch (std::abs(idhep)) {
                case 533:
                    status = 1; // Bs0* (vector)
                    break;
                case 531:
                    status = 2; // Bs0  (scalar)
                    break;
                case 431:
                    status = 3; // Ds+-
                    break;
                case 10431:
                    status = 4; // Ds1*(2317)
                    break;
                case 4:
                    status = 5; // ccbar
                    break;
                case 3:
                case 2:
                case 1:
                    status = 6; // uds
                    break;
            }
        }
        return status;
    }
    // ***********************************************************
    void addEventValues2Columns(BelleTuple* tt, bool debugDump = false) {
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

        tt->column("expn",  expNo); // Exp #
        tt->column("runn",  runNo); // Run #
        tt->column("evtn",  evtNo);
        tt->column("ipx",   ip_position.x());
        tt->column("ipy",   ip_position.y());
        tt->column("ipz",   ip_position.z());
        // ************* Signal shape ********
        tt->column("r2", r2);
        // *************  MC  ****************
        tt->column("evtgen", idGenType);
    }
    // ***********************************************************
    VectorL getGenVectorL(int idhPcl) {
        VectorL pclL;
        int ID_Pcl = -1;
        Gen_hepevt_Manager& genMgr = Gen_hepevt_Manager::get_manager();
        for (std::vector<Gen_hepevt>::iterator itr = genMgr.begin(); itr != genMgr.end(); ++itr) {
            int idh = itr->idhep();
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
    void addPi0Values2Columns(BelleTuple* tt, Particle& p0, const std::string& sfx, bool debugDump) {

        if (!&p0.userInfo()) createUserInfo(p0);
        UserInfo& p0_info = dynamic_cast<UserInfo&>(p0.userInfo());

        Particle& g1 = p0.child(0);
        Particle& g2 = p0.child(1);

        double msPi0_gg = (g1.p() + g2.p()).m();
        double psrPi0   = pStar(p0.p(), BeamEnergy::E_HER(), BeamEnergy::E_LER(), BeamEnergy::Cross_angle()).vect().mag();
        double pPi0     = p0.momentum().p().mag();

        bool gen_pi0 = getIdHep(p0) == 0 ? false : true;
        bool gen_g1  = getIdHep(g1) == 0 ? false : true;
        bool gen_g2  = getIdHep(g2) == 0 ? false : true;

        const int nValI = 3;
        const int nValD = 13;

        std::string pclTitI[nValI] = {"gen",   "gg1",  "gg2"};
        int valPclI[nValI]         = {gen_pi0, gen_g1, gen_g2};
        std::string pclTitD[nValD] = {"eg1", "eg2", "psr", "p", "mgg", "ms", "msm", "chi", "chm", "chn", "cmn", "cl", "clm"};
        double valPclD[nValD]      = {g1.ptot(), g2.ptot(), psrPi0, pPi0, msPi0_gg, p0_info.msKvf(), p0_info.msKmvf(),
                                      p0_info.chisqKvf(), p0_info.chisqKmvf(), p0_info.chisqKvfNdf(),
                                      p0_info.chisqKmvfNdf(), p0_info.clKvf(), p0_info.clKmvf()};

        for (int iVal = 0; iVal < nValI; ++iVal) tt->column(pclTitI[iVal] + sfx, valPclI[iVal]);
        for (int iVal = 0; iVal < nValD; ++iVal) tt->column(pclTitD[iVal] + sfx, valPclD[iVal]);
    }

    // ***********************************************************
    void addValues2Columns(BelleTuple* tt, int nValI, int nValD, int* valPclI, double* valPclD, std::string* pclTitI, std::string* pclTitD,
                           const std::string& sfx, bool debugDump) {
        for (int iVal = 0; iVal < nValI; ++iVal)
            tt->column(pclTitI[iVal] + sfx, valPclI[iVal]);
        for (int iVal = 0; iVal < nValD; ++iVal)
            tt->column(pclTitD[iVal] + sfx, valPclD[iVal]);
    }
    // ***********************************************************
    void addGenValues2Columns(BelleTuple* tt, bool gen_pcl, VectorL pclL, const std::string& sfx, bool debugDump) {
        const std::size_t nValD = 9;
        const std::string pclTitD[nValD] = {"ms", "px", "py", "pz", "e", "th", "ct", "ph", "rh"};
        double valPclD[nValD];

        if (gen_pcl) {
            valPclD[0] = pclL.m();
            valPclD[1] = pclL.px();
            valPclD[2] = pclL.py();
            valPclD[3] = pclL.pz();
            valPclD[4] = pclL.e();
            valPclD[5] = pclL.theta();
            valPclD[6] = pclL.cosTheta();
            valPclD[7] = pclL.phi();
            valPclD[8] = pclL.rho();
        } else {
            std::fill_n(valPclD, nValD, -0.99);
        }

        for (std::size_t iVal = 0; iVal < nValD; ++iVal) {
            tt->column(pclTitD[iVal] + sfx, valPclD[iVal]);
        }
    }
    // ***********************************************************
    void dumpDsChild(BelleTuple* tt, Particle& P, const std::string& sfxDs = "",
                     bool evtInfoDump = false,
                     bool stDump = true, bool debugDump = false) {

        if (evtInfoDump) addEventValues2Columns(tt, debugDump);

        int lund = int(P.lund());

        if (!&P.userInfo()) createUserInfo(P);
        UserInfo& info = dynamic_cast<UserInfo&>(P.userInfo());

        double chisqKvf      = info.chisqKvf();
        double chisqKmvf     = info.chisqKmvf();
        double chisqKvfNdf   = info.chisqKvfNdf();
        double chisqKmvfNdf  = info.chisqKmvfNdf();
        double msKvf         = info.msKvf();
        double msKmvf        = info.msKmvf();
        double clKvf         = info.clKvf();
        double clKmvf        = info.clKmvf();
        double helic         = info.helicity();

        int ind_ch1 = 0;
        int lundChild = int(P.lund());
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
        double decay_vx_ds_child      = P.momentum().decayVertex().x();
        double decay_vy_ds_child      = P.momentum().decayVertex().y();
        double decay_vz_ds_child      = P.momentum().decayVertex().z();

        Hep3Vector P3D                (px_ds_child, py_ds_child, pz_ds_child);
        double perp_ds_child          = P3D.perp();
        double phi_ds_child           = P3D.phi();
        double theta_ds_child         = P3D.theta();


        std::string dgrSuff = "_ch" + sfxDs;
        int signPcl = int(P.lund() > 0 ? 1 : -1);
        bool gen_pcl = getIdHep(P) == 0 ? false : true;

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

        addValues2Columns(tt, nValI, nValD, valPclI, valPclD, pclTitI, pclTitD, dgrSuff, debugDump);

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

        if (evtInfoDump) addEventValues2Columns(tt, debugDump);
        int lundDs = int(P.lund());

        // Creating a userInfo object for the Ds meson candidate
        if (!&P.userInfo()) createUserInfo(P);
        UserInfo& info = dynamic_cast<UserInfo&>(P.userInfo());

        std::vector<std::pair<int, double> > ds_child_ids = getDsChildrenIdentification(P);



        double chisqKvf      = info.chisqKvf();
        double chisqKmvf     = info.chisqKmvf();
        double chisqKvfNdf   = info.chisqKvfNdf();
        double chisqKmvfNdf  = info.chisqKmvfNdf();
        double msComb        = info.msComb();
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
        bool   gen_pcl          = getIdHep(P) == 0 ? false : true;
        double decay_vx_ds      = P.momentum().decayVertex().x();
        double decay_vy_ds      = P.momentum().decayVertex().y();
        double decay_vz_ds      = P.momentum().decayVertex().z();
        double helicChild1      = getHelicity(P);
        double pid00            = ds_child_ids[0].second,
               pid01            = ds_child_ids[1].second,
               pid10            = ds_child_ids[2].second;
        int    c00              = ds_child_ids[0].first,
               c01              = ds_child_ids[1].first,
               c10              = ds_child_ids[2].first;


        // Working with Ds+(-) children (K+(-), K-(+) and pi+(-))
        Particle& Child0        = P.child(0);             // phi0 (--> K+ K-) or K* (--> K+ pi-), or K*bar (--> K- pi+)

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

        VectorL dssL = getGenVectorL(getIdHep(P));


        tt->column("hel_ch" + sfxDs, helicChild1);
        addValues2Columns(tt, nValI, nValD, valPclI, valPclD, pclTitI, pclTitD, dgrSuff, debugDump);
        addGenValues2Columns(tt, gen_pcl, dssL, genDgrSuff, debugDump);
        dumpDsChild(tt, Child0, sfxDs, false, false, debugDump);
        if (stDump) tt->dumpData();
    }
    // ***********************************************************
    void dumpDs2317(BelleTuple* tt, Particle& P, std::string sfxDs = "", bool evtInfoDump = false, bool stDump = true, bool debugDump = false) {

        if (evtInfoDump) addEventValues2Columns(tt, debugDump);

        if (!&P.userInfo()) createUserInfo(P);
        UserInfo& info = dynamic_cast<UserInfo&>(P.userInfo());


        double chisqKvf      = info.chisqKvf();
        double chisqKmvf     = info.chisqKmvf();
        double chisqKvfNdf   = info.chisqKvfNdf();
        double chisqKmvfNdf  = info.chisqKmvfNdf();
        double msKvf         = info.msKvf();
        double msKmvf        = info.msKmvf();
        double clKvf         = info.clKvf();
        double clKmvf        = info.clKmvf();
        double helic_2317    = -1.;


        Particle& Child           = P.child(0);
        Particle& pi0_2317        = P.child(1);
        UserInfo& infoChild           = dynamic_cast<UserInfo&>(Child.userInfo());
        UserInfo& infoPi0_2317        = dynamic_cast<UserInfo&>(pi0_2317.userInfo());

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
        int chg_d17               = int(P.lund() >  0 ? 1 : -1);
        int gen_d17               = int(getIdHep(P) == 0 ? 0 : 1);

        const int nValI = 2;
        const int nValD = 20;
        int valPclI[nValI] = {chg_d17, gen_d17};
        double valPclD[nValD] = {msKvf, msKmvf, chisqKvf, chisqKmvf, chisqKvfNdf, chisqKmvfNdf,
                                 clKvf, clKmvf, decay_vx_d17, decay_vy_d17, decay_vz_d17, psr_d17,
                                 px_d17, py_d17, pz_d17, p_d17, perp_d17, phi_d17, theta_d17, helic_2317};
        std::string pclTitI[nValI] = {"chg", "gen"};
        std::string pclTitD[nValD] = {"ms","msm","chi", "chim", "chn", "chmn",
                                 "cl", "clm", "dvx", "dvy", "dvz", "psr",
                                 "px", "py", "pz", "p", "pt", "ph", "th","hel"};

        std::string dgrSuff = "_d17", genDgrSuff = "_d17_t";

        VectorL ds17L = getGenVectorL(getIdHep(P));

        addValues2Columns(tt, nValI, nValD, valPclI, valPclD, pclTitI, pclTitD, dgrSuff, debugDump);
        addGenValues2Columns(tt, gen_d17, ds17L, genDgrSuff, debugDump);
        addPi0Values2Columns(tt, pi0_2317, "_p0_d", debugDump);
        dumpDs(tt, Child, sfxDs, false, false, debugDump);

        if (stDump) tt->dumpData();
    }
    // ***********************************************************
    void dumpBs0(BelleTuple* tt, Particle& P, bool evtInfoDump = false,
                     bool stDump = true, bool debugDump = true) {

        if (evtInfoDump) addEventValues2Columns(tt, debugDump);

        if (!&P.userInfo()) createUserInfo(P);
        UserInfo &info = dynamic_cast<UserInfo &>(P.userInfo());

        double chisqKvf = info.chisqKvf();
        double chisqKmvf = info.chisqKmvf();
        double chisqKvfNdf = info.chisqKvfNdf();
        double chisqKmvfNdf = info.chisqKmvfNdf();
        double msKvf = info.msKvf();
        double msKmvf = info.msKmvf();
        double clKvf = info.clKvf();
        double clKmvf = info.clKmvf();

        Particle &Dss_Bs0 = P.child(0);
        Particle &Dss2317_Bs0 = P.child(1);
        Particle &pi0_Bs0 = P.child(2);

        double px_bs = P.px();
        double py_bs = P.py();
        double pz_bs = P.pz();
        double p_bs = P.momentum().p().mag();

        Hep3Vector P3D(px_bs, py_bs, pz_bs);

        double perp_bs = P3D.perp();
        double phi_bs = P3D.phi();
        double theta_bs = P3D.theta();

        double decay_vx = P.momentum().decayVertex().x();
        double decay_vy = P.momentum().decayVertex().y();
        double decay_vz = P.momentum().decayVertex().z();
        int gen_bs = int(getIdHep(P) == 0 ? 0 : 1);
        int chg_bs = int(P.lund() > 0 ? 1 : -1);


        VectorL pB = pStar(P.p());
        VectorL pB2 = pStar(P.p(), BeamEnergy::E_HER(), BeamEnergy::E_LER(), BeamEnergy::Cross_angle());
        double de_bs_old = pB.e() - Benergy();
        double de_bs = pB.e() - BeamEnergy::E_beam_corr();
        double de2_bs = pB2.e() - BeamEnergy::E_beam_corr();
        double mbc_bs_old = beamEnergyConstraint(P);
        double mbc_bs = beamEnergyConstraint(P, BeamEnergy::E_HER(), BeamEnergy::E_LER(),
                                             BeamEnergy::Cross_angle() * 1.e3);

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

        VectorL bs0L = getGenVectorL(getIdHep(P));
        std::string genSfx = "_bs_t";

        addValues2Columns(tt, nValI, nValD, valPclI, valPclD, pclTitI, pclTitD, "_bs", debugDump);
        addGenValues2Columns(tt, gen_bs, bs0L, genSfx, debugDump);
        addPi0Values2Columns(tt, pi0_Bs0, "_p0_b", debugDump);
        dumpDs(tt, Dss_Bs0, "1", false, false, debugDump);
        dumpDs2317(tt, Dss2317_Bs0, "2", false, false, debugDump);

        if (stDump) tt->dumpData();
    }

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
        const HepPoint3D& ip_position = IpProfile::position(1);
        const HepSymMatrix& ip_error  = IpProfile::position_err(1);

        // Event Information
        int expNo = 0, runNo = 0, evtNo = 0;
        getEventInfo(expNo, runNo, evtNo, McFlag); // utility.cc
        std::cout << "\tExpNo: " << expNo << "\tRunNo: " << runNo << "\tEventNo: " << evtNo << "\tMcFlag: " << McFlag << std::endl;

        if (evtNo == 470) {
            std::cout << "--1" << std::endl;
        }

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


        // If each plist element is not within atc_pID.prob >= prob,
        // its element is removed from plist.
        // ----------------- pos kaons ---- default values - kaons - pions
        withKaonId(trkV[2], minProbPID_Kp,  3,    1,     5,  3,      2);      // K+ vs bg pi

        // ----------------- neg kaons ---- default values - kaons - pions
        withKaonId(trkV[3], minProbPID_Kn,  3,    1,     5,  3,      2);      // K- vs bg pi


        // If each plist element is not within atc_pID.prob < prob,
        // its element is removed from plist.
        // ----------------- pos pions ---- default values - kaons - pions
        withPionId(trkV[0], maxProbPion,    3,    1,     5,  2,      3);        // pi+ vs bg K
        // ----------------- neg pions ---- default values - kaons - pions
        withPionId(trkV[1], maxProbPion,    3,    1,     5,  2,      3);        // pi- vs bg K


        // If each plist element is not associated with rphi & z-svd hits
        // whose number is equal to or larger than nRSvdHit and nZSvdHit,
        // its element is removed from plist.
        for (int itr = 0; itr < nTrk; ++itr) {
            withSVD2(trkV[itr], 1, 1); // nRSvdHit, nZSvdHit
            withdRdZcut(trkV[itr], ip_position, dRcut, dZcut);
        }


        // ================= WORKING WITH Gamma CANDIDATES ================= //
        makeGamma(gammaV);
        withPCut(gammaV, eGammaMin);

        if(useVTX) {
            for (std::vector<Particle>::iterator itr = gammaV.begin(); itr != gammaV.end(); ++itr) {
                setGammaError(*itr, ip_position, ip_error); // changed from setGammasError
            }
        }

        // Match gamma candidates with their gen-hep info
        if (McFlag) {
            setGenHepInfoG(gammaV);
        }

        if (evtNo == 470) {
            std::cout << "--2" << std::endl;
        }
        createUserInfo(gammaV);


        // =================   WORKING WITH Pi0 CANDIDATES ================= //
        makePi0(pi0);
        // Checking gammas' energies for all the pi0 daughters
        withPi0GammPCut(pi0, minPi0GammP);
        // !!!! PAY ATTENTION. Start point

        // Creating UserInfo objects in memory storage for the selected pi0's
        createUserInfo(pi0);
        if (evtNo == 470) {
            std::cout << "--3" << std::endl;
        }

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


        // Match candidates with gen-hep info
        if (McFlag) {
            for (int itr = 0; itr < nTrk; ++itr) setGenHepInfoF(trkV[itr]);
            setGenHepInfoP(pi0);
        }
        if (evtNo == 470) {
            std::cout << "--4" << std::endl;
        }
        combination(phi0,      Ptype("PHI"),  trkV[2],   trkV[3], dM_V0);   // k_p, k_m
        combination(Ksr0,      Ptype("K*0"),  trkV[2],   trkV[1], dM_Ksr0); // k_p, pi_m
        combination(Ksr0bar,   Ptype("K*B"),  trkV[3],   trkV[0], dM_Ksr0); // k_m, pi_p
        setGenHepInfoT(phi0);
        setGenHepInfoT(Ksr0);
        setGenHepInfoT(Ksr0bar);

        combination(Dss_p, Ptype("DS+"),    phi0,     trkV[0], dM_Dss); // phi0, pi_p
        combination(Dss_m, Ptype("DS-"),    phi0,     trkV[1], dM_Dss); // phi0, pi_m
        combination(Dss_p, Ptype("DS+"),    Ksr0bar,  trkV[2], dM_Dss); // K*0bar (K-pi+), k_p
        combination(Dss_m, Ptype("DS-"),    Ksr0,     trkV[3], dM_Dss); // K*0 (K+pi-), k_m
        setGenHepInfoT(Dss_p);
        setGenHepInfoT(Dss_m);

        combination(Dss_m_2317, Ptype(-10431), Dss_m, pi0, dM_2317);
        combination(Dss_p_2317, Ptype( 10431), Dss_p, pi0, dM_2317);
        setGenHepInfoT(Dss_m_2317);
        setGenHepInfoT(Dss_p_2317);

        combination(Bs0, Ptype(531), Dss_m, Dss_p_2317, pi0, dM_Bs0);
        combination(Bs0, Ptype(531), Dss_p, Dss_m_2317, pi0, dM_Bs0);
        setGenHepInfoT(Bs0);
        combination(Bs0bar, Ptype(-531), Dss_m, Dss_p_2317, pi0, dM_Bs0);
        combination(Bs0bar, Ptype(-531), Dss_p, Dss_m_2317, pi0, dM_Bs0);
        setGenHepInfoT(Bs0bar);

        /* combination(BsStar0, Ptype(533), Bs0, gammaV, dM_Bs0);
         * combination(BsStar0bar, Ptype(-533), Bs0bar, gammaV, dM_Bs0);
         * setGenHepInfoT(BsStar0);
         * setGenHepInfoT(BsStar0bar);
         * combination(Upsilon_5S, Ptype(9000553),  BsStar0, BsStar0bar, dM_Bs0);
         */

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
        if (evtNo == 470) {
            std::cout << "--5" << std::endl;
        }
        // Bs0
        if (stDumpBs0) {
            for (int iEvt = 0; iEvt < Bs0.size(); ++iEvt) {
                if (evtNo == 470) {
                    std::cout << "--6" << std::endl;
                }
                makeRecursiveVertexFit(Bs0[iEvt]);
                if (evtNo == 470) {
                    std::cout << "--7" << std::endl;
                }
                dumpBs0(TP_Bs0, Bs0[iEvt], true, stDumpBs0, debugDumpBs0);
                }

            for (int iEvt = 0; iEvt < Bs0bar.size(); ++iEvt) {
                makeRecursiveVertexFit(Bs0bar[iEvt]);
                dumpBs0(TP_Bs0, Bs0bar[iEvt], true, stDumpBs0, debugDumpBs0);
                }

    //        for (int iEvt = 0; iEvt < Bs0bar.size(); ++iEvt)
    //            dumpBs0(TP_Bs0, Bs0bar[iEvt], true, stDumpBs0, debugDumpBs0);
                }
        std::cout << "---------  clearVectors (final) ----------" << std::endl;
        clearVectors();
        if (evtNo == 470) {
            std::cout << "--8" << std::endl;
        }
        std::cout << "------------------  Reco end ---------------------" << std::endl;
        *status = 1;

    }

} // namespace Belle

