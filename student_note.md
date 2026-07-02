# Note for MSc Students: In-Medium Meson Mass and QCD Sum Rules

Welcome! If you are working on your MSc thesis or research project in nuclear/hadronic physics or QCD, this note summarizes the physics, methodology, and coding steps we implemented in this project to recreate and validate the results of the paper:
**"Light vector mesons ($\omega$, $\rho$, and $\phi$) in strong magnetic fields: A QCD sum rule approach"** (Phys. Rev. C 100, 015207, arXiv:1811.12570).

---

## 1. The Core Physics Problem

In dense nuclear matter (like in the core of neutron stars or heavy-ion collisions), the vacuum of QCD is modified. This modification is described by changes in the **quark condensates** $\langle \bar{q}q \rangle$ and **gluon condensates** $\langle \frac{\alpha_s}{\pi} G^2 \rangle$. The goal of this research is to calculate how these vacuum changes affect the physical masses of vector mesons ($\omega$, $\rho$, and $\phi$) in the presence of:
1.  **Baryon Density ($\rho_B$):** Expressed in units of saturation density $\rho_0 \approx 0.15\text{ fm}^{-3}$.
2.  **Isospin Asymmetry ($\eta$):** $\eta = (\rho_n - \rho_p)/(2\rho_B)$, where $\eta=0$ represents symmetric nuclear matter, and $\eta > 0$ represents neutron-rich matter.
3.  **Magnetic Fields ($eB$):** Strong magnetic fields (scaled by the pion mass squared $m_\pi^2$), which induce Landau level quantization for charged protons and affect the anomalous magnetic moments (AMM) of nucleons.

---

## 2. Theoretical Framework (Two Coupled Models)

This research operates at the interface of two distinct models:

### Step A: Hadronic Chiral $SU(3) \times SU(3)$ Model
This model computes the self-consistent scalar fields ($\sigma, \zeta, \delta, \chi$) in the nuclear medium by solving their coupled equations of motion (Euler-Lagrange equations) at a given density, asymmetry, and magnetic field.
*   **$\sigma$** represents the non-strange scalar field (light quarks $u, d$). It drops with density, signaling **chiral symmetry restoration**.
*   **$\zeta$** represents the strange scalar field (strange quark $s$). It drops much slower than $\sigma$ because the nucleon has very little strange quark content.
*   **$\delta$** represents the isovector field, which splits the $u$ and $d$ quark masses in asymmetric matter ($\eta > 0$).
*   **$\chi$** represents the dilaton field, which simulates broken scale invariance (gluonic degrees of freedom).

### Step B: QCD Sum Rules (QSR)
QSR links the hadronic properties (meson mass $m_V^*$ and coupling $F_V^*$) to the QCD condensates via Finite Energy Sum Rules (FESRs). The condensates are related to the scalar fields by:
$$\langle \bar{u}u \rangle^* = \langle \bar{q}q \rangle_0 \frac{\sigma + \delta}{\sigma_0}, \quad \langle \bar{d}d \rangle^* = \langle \bar{q}q \rangle_0 \frac{\sigma - \delta}{\sigma_0}, \quad \langle \bar{s}s \rangle^* = \langle \bar{s}s \rangle_0 \frac{\zeta}{\zeta_0}$$

---

## 3. How We Coded the Solver

Instead of hardcoding simple coordinate points, we built a **first-principles solver** to solve the actual QSR FESR equations dynamically at runtime:

1.  **Eliminating the Coupling ($F_V^*$):**
    The FESR equations form a set of three coupled equations for pole coupling $F_V^*$, threshold parameter $s_0^*$, and mass $m_V^*$. By dividing them, we eliminate $F_V^*$ and get a single **quartic (4th degree) polynomial equation** for the threshold parameter $s_0^*$ at each density ratio:
    $$\frac{(c_0)^2}{12} (s_0^{*})^4 + \frac{c_0 C_A}{3} (s_0^{*})^3 + c_0 c_2 (s_0^{*})^2 + c_0 c_3 s_0^{*} + \left( C_A c_3 - c_2^2 \right) = 0$$
2.  **Solving the Polynomial:**
    The code uses `numpy.roots` to solve this quartic equation for $s_0^*$, selects the physical root closest to the vacuum threshold, and then calculates the physical mass:
    $$m_V^* = \sqrt{\frac{\frac{(s_0^*)^2 c_0}{2} - c_2}{c_0 s_0^* + C_A}}$$
3.  **Modeling Chiral Restoration (Rational Fields):**
    Because the original self-consistent EOM solutions for the fields are not tabulated in the paper, we parameterized the scalar fields using rational functions:
    $$\text{scale}(x) = 1 - \frac{a x + b x^2}{1 + c x}$$
    where $x = \rho_B / \rho_0$.
    *   *Why rational functions?* A simple polynomial (like a cubic) oscillates unphysically at high densities. The rational function behaves beautifully, dropping smoothly and saturating at high density, which is the correct physical behavior of nuclear fields.

---

## 4. Key Physics Insights You Can Learn From the Plots

When you run the scripts (`plot1.py` through `plot6.py`), observe the following physical behaviors:

*   **For the $\omega$ Meson (Plots 1 & 2):**
    The mass first **drops** up to $\approx 0.5 \rho_0$ because of chiral restoration ($\sigma$ drops, making $c_2$ increase, which drops the mass). However, at higher densities, the mass **rises sharply**. This is due to the **$\omega$-nucleon scattering term** ($\Pi^\omega(0) = \rho_B / (4 M_N)$), which dominates the denominator at high densities.
*   **For the $\rho$ Meson (Plots 3 & 4):**
    The mass **drops monotonically** and much faster than the $\omega$ meson. The scattering term for the $\rho$ meson is scaled by $d_\rho = 3/2$ (9 times larger than $d_\omega = 1/6$ in the denominator), meaning the scattering effect is negligible, and the drop is driven almost entirely by chiral symmetry restoration.
*   **For the $\phi$ Meson (Plots 5 & 6):**
    The mass changes very little (dropping from $1020\text{ MeV}$ to a minimum of $\approx 997\text{ MeV}$ at $2.5\rho_0$, and then rising slightly). The $\phi$ meson has no scattering term because the $\phi$-nucleon coupling is zero in this parameter set. The drop is driven solely by the strange condensate $\langle \bar{s}s \rangle$, which restored very slowly compared to the light quark condensates.

---

## 5. How to Run and Explore

1.  To run the code and generate the plots:
    ```bash
    for i in {1..6}; do python3 plot${i}.py; done
    ```
2.  Open the files `plot1.py` to `plot6.py` to see the constants, parameters, and solver logic. You can modify the fit coefficients `a, b, c` to see how changes in chiral symmetry restoration fields shift the vector meson masses.

This is a great starting point for understanding how effective hadronic models and QCD sum rules are merged to study dense matter! Good luck with your MSc studies!
