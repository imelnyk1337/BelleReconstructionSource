### Selection Criterion
**This is a short overview of the cuts, which were applied on different stages of the reconstruction/analysis**

#### Pre-selection
The main cuts are considered to be applied either on the reconstruction stage of the analysis or after all candidates are recorded into a ```ROOT```-based (or converted from ```hbook``` files).

##### BASF reconstruction

Charged particles' identification;
since we have two $D_{s}^{\pm}$ meson in the decay mode, we expect two deal with 6 charged tracks (3 from the 1st (ground) $D_{s}$ and other 3 from the secondary one). These 3 tracks are $K^{\pm}$, $K^{\mp}$, $\pi^{\pm}$.

- Likelihood ratio: $P(K/\pi) = L(K)/(L(K) + L(\pi))$:
  - greater than $0.6$ for the 1st $K$;
  - greater than $0.2$ for the 2nd $K$;
  - less than $0.9$ for $\pi^{0}$ candidates.

- Distance to the POCA (the point of the closest approach):
  - $d_{0} < 0.5$ cm;
  - $d_{z} < 3$ cm.
- SVD (Silicon Vertex Detector) hits requirements:
  - ```withSVD2(trkV[itr], 1, 1); // nRSvdHit, nZSvdHit```
- $D_{s}^{\pm}$ pre-selection criteria:
  - $\phi^{0}(1020)$ mass window: $[1009.461 ,1029.461]$ MeV (the nominal mass is 1019.461)



