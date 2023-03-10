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
  - $\phi^{0}(1020) \rightarrow K^{+}K^{-}$ mass window: $[1009.461 ,1029.461]$ MeV (the nominal mass is $1019.461$ MeV)
  - $K^{*}(892) \rightarrow K^{+}\pi^{-}$ mass window: $[845.81, 945.81]$ MeV (the nominal mass is $895.81$ MeV)
  - $D_{s}^{\pm}$ mass window: $[1.955, 1.979]$ MeV (the nominal mass is $1968.30$ MeV)
- $D_{s_{0}}^{*}(2317)^{\pm} \rightarrow D_{s}^{\pm} \pi^{0}$ pre-selection criteria:
  - $E(\gamma) > 100$ MeV;
  - Reconstructed $\pi^{0}$ invariant mass: $122 < m_{reco}(\pi^{0}) = m(\gamma\gamma) < 148$ MeV;
  - $\pi^{0}$ momentum: $p(\pi^{0}) > 150$ MeV;
  - Vertex fit (mass-constraint fit to $\pi^{0}$ nominal mass): $\chi^{2}_{vertexing} < 200$

#### Post-selection criteria

- A second Fox-Wolfram momentum $R_{2}$:
  - should be greater than 0.8 to get rid of $B\bar{B}$ contribution ?;
- $\phi^{0}$ helicity angle:
  - $|\cos{\theta_{H}}| > 0.42$ 
