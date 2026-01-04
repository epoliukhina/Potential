## Extraction of KBI from cryo-ET tomograms of proteins

The quantitative description of protein–protein interactions (PPIs) can be given in the form of an interaction potential, defined as the work required to bring two protein macromolecules from an infinite distance to a finite distance. There exist two forms of interaction potential: the **potential of mean force** (PMF, *W(r)*) and the **pair interaction potential** (PIP, *U(r)*), which are linked through:

\[
W(r)= U(r)+\Delta w(r),
\tag{1}
\]

where \(\Delta w(r)\) is a correction accounting for the effect of surrounding particles, and \(r\) is the interparticle centroid-to-centroid distance. In other words, the concentration-dependent \(W(r)\) converges to the concentration-independent \(U(r)\) in the infinite dilution limit.

In statistical-mechanical theory, the PMF can be obtained via the reversible work theorem:

\[
W(r)=-k_B T \ln[g(r)],
\tag{2}
\]

where \(g(r)\) is the radial distribution function (RDF), which can be directly calculated from the spatial distribution of particles.

Another thermodynamic quantity providing insight into PPIs is the second virial coefficient \(B_{22}\) and the Kirkwood–Buff integral (KBI, \(G_{22}^{\infty}\)):

\[
B_{22}=2\pi\int_{0}^{\infty}\left(e^{-U(r)/(k_B T)}-1\right) r^2\,dr,
\tag{3}
\]

\[
G_{22}^{\infty}=4\pi\int_{0}^{\infty}\left(e^{-W(r)/(k_B T)}-1\right) r^2\,dr,
\tag{4}
\]

where \(k_B\) is the Boltzmann constant and \(T\) is the absolute temperature. From Eqs. 1 and 3–4, in the dilute regime typically reached for proteins, \(B_{22}\) and KBI are proportional by a factor of \(-2\):

\[
G_{22,0}^{\infty} = -2B_{22}.
\tag{5}
\]

The particle spatial distribution readily available from cryo-ET tomograms suggests that statistical-mechanical treatments can be applied to extract useful thermodynamic parameters. Here, we use cryo-ET to determine the KBI of protein systems in several ways.

### Direct integration of the RDF

By definition, KBI can be obtained by direct integration of the RDF (Fig. SI5a):

\[
G_{22}^{\infty}=4\pi\int_{0}^{\infty}\left(g(r)-1\right) r^2\,dr.
\tag{6}
\]

The cumulative \(G_{22}(r)\) for BSA, computed by replacing the upper limit with \(r\) in Eq. 6, is shown in Fig. SI5b. Because \(g(r)\) starts fluctuating around 1 at a distance of approximately 18 nm, \(G_{22}(r)\) is expected to reach its asymptote \(G_{22}^{\infty}\) near that distance. This, however, is only partially observed, revealing two main disadvantages of this approach:

1. The experimentally determined \(g(r)\) is noisy; since the noise is multiplied by \(r^2\), it strongly amplifies the integration noise.
2. In a finite-size box system, the RDF may not converge to 1.

The second issue is largely addressed using advanced RDF calculation methods. However, the noise problem is clearly present, as evidenced by the growing deviation between replicates at distances of 15 nm and above.

### Sub-box (Schnell’s) method

To circumvent these drawbacks, we apply the **sub-box (Schnell’s) method**, a statistical-mechanics approach developed for finite-volume samples. This method relies on density fluctuation theory, which states that in the grand-canonical ensemble (\(\mu,V,T = \mathrm{const}\)), the KBI can be rewritten as:

\[
G_{22}= \frac{V}{\langle N_2\rangle^2}\left(\langle N_2^2\rangle-\langle N_2\rangle^2-\langle N_2\rangle\right),
\tag{7}
\]

where \(V\) is the sub-box volume, \(\langle N_2\rangle\) is the ensemble-averaged number of particles in the sub-box, and \(\langle N_2^2\rangle\) is the ensemble average of its square.

Practically, this is performed by dividing the cryo-ET tomogram (“bath” in statistical-mechanical terms) into sub-boxes—which can exchange particles and energy with the bath—of a chosen geometric shape (Fig. SI5c). In our case, the sub-boxes are parallelepipeds with a square base of side \(L\) and fixed height equal to the tomogram thickness \(H\). Sub-boxes overlap with a user-chosen step. Decreasing the step smooths the resulting curve without changing its trend (Fig. SI5d); accordingly, we fix the step at 3 nm for all analyses.

By performing this division for sub-boxes of various sizes, one can determine the appropriate value of \(G_{22}^{\infty}\). This is achieved when the volume of the sub-box is sufficiently large such that \(G_{22}\) becomes independent of box size. In contrast to gold nanoparticles, which exhibit a clear plateau with increasing sub-box size, proteins do not: the effective KBI monotonically decreases with \(L\) (Fig. SI5e).

We found that a robust approach in this case is to determine KBI as the intercept of the dependence of effective \(G_{22}\) on the surface area-to-volume ratio \(S/V\) of a sub-box (Fig. SI5f):

\[
G_{22}(S/V)=G_{22}^{\infty}+\frac{S}{V}G_{22}^{\prime},
\tag{8}
\]

where \(S\) is the surface area of the sub-box and \(G_{22}^{\prime}\) is the surface contribution to \(G_{22}\). The linear regime of \(G_{22}(S/V)\) corresponds to \(L\) values below 50 nm.
