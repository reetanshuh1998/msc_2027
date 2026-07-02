# README: plot6.py - Phi Meson Mass at $eB = 12 m_\pi^2$

This document describes the structure, equations, parameters, and step-by-step logic of [plot6.py](file:///home/reethep/msc_2027/plot6.py), which recreates Figure 6 of the paper (the $\phi$ meson mass vs density for $eB = 12 m_\pi^2$).

---

## 1. Physics Model & Equations Solved

The script models the in-medium mass $m_\phi^*$ by dynamically solving the coupled Finite Energy Sum Rules (FESR) quartic equation:

$$\frac{(c_0)^2}{12} (s_0^{*})^4 + \frac{c_0 C_A}{3} (s_0^{*})^3 + c_0 c_2 (s_0^{*})^2 + c_0 c_3 s_0^{*} + \left( C_A c_3 - c_2^2 \right) = 0$$

where:
*   **$c_0$:** Perturbative QCD coefficient $c_0 = 1 + \alpha_s/\pi \approx 1.1146$.
*   **$C_A$:** Constant for the strange meson because the $\phi$-nucleon coupling is zero, meaning the scattering term $\Pi^\phi(0) = 0$:
    $$C_A = c_1^\phi = -6 m_s^2 = -0.135\text{ GeV}^2$$
*   **$c_2$ (OPE Dimension 4):**
    $$c_2 = \frac{\pi^2}{3} \langle \frac{\alpha_s}{\pi} G^2 \rangle^* + 8\pi^2 m_s \langle \bar{s}s \rangle^*$$
*   **$c_3$ (OPE Dimension 6):**
    $$c_3 = - \alpha_s \pi^3 \frac{224}{81} \kappa_s \times 8 \langle \bar{s}s \rangle^{*2}$$

Once the positive real root of the quartic equation $s_0^*$ is found, the mass is computed:
$$m_\phi^* = \sqrt{\frac{\frac{(s_0^*)^2 c_0}{2} - c_2}{c_0 s_0^* + C_A}}$$

---

## 2. Step-by-Step Logic

1.  **Vacuum Pre-calculation:**
    *   Computes the vacuum threshold $s_0^\phi \approx 1.838\text{ GeV}^2$.
    *   Calculates the factorization parameter $\kappa_s \approx 1.06$ self-consistently.
2.  **Evaluating the Medium:**
    *   Iterates over the density range $x \in [0, 4.5]$ (where $x = \rho_B/\rho_0$).
    *   Evaluates the rational scale factor for the strange field $\zeta$:
        $$\zeta = \zeta_0 \left( 1 - \frac{a x + b x^2}{1 + c x} \right)$$
        using the curve-specific fits:
        *   **$\eta=0$ (with AMM):** $a=0.149486, b=-0.009813, c=1.100261$
        *   **$\eta=0.3$ (with AMM):** $a=0.145936, b=-0.008139, c=1.118451$
        *   **$\eta=0.5$ (with AMM):** $a=0.144405, b=-0.006074, c=1.181535$
        *   **$\eta=0$ (without AMM):** $a=0.145395, b=-0.013800, c=0.955973$
        *   **$\eta=0.3$ (without AMM):** $a=0.144349, b=-0.011570, c=0.997052$
        *   **$\eta=0.5$ (without AMM):** $a=0.144534, b=-0.009573, c=1.060511$
    *   Scales fields $\sigma, \delta, \chi$:
        $$\sigma = \sigma_0 \left( 1 - \frac{0.25 x}{1 + 0.15 x} \right), \quad \delta = \frac{0.008 \eta x}{1 + 0.15 x}, \quad \chi = \chi_0 (1 - 0.005 x)$$
    *   Computes condensates, OPE coefficients $c_2, c_3$, and scattering term $C_A$.
3.  **Quartic Solving:**
    *   Evaluates roots of the quartic polynomial and selects the root closest to the vacuum threshold.
4.  **Plotting:**
    *   Draws the curve using standard color conventions: red ($\eta=0$), blue ($\eta=0.3$), and green ($\eta=0.5$).
    *   Draws solid dash-dot lines for with-AMM, and dotted lines for without-AMM.
    *   Saves the figure as `plots/plot6.png`.
