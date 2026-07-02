# README: plot1.py - Omega Meson Mass at $eB = 4 m_\pi^2$

This document describes the structure, equations, parameters, and step-by-step logic of [plot1.py](file:///home/reethep/msc_2027/plot1.py), which recreates Figure 1 of the paper (the $\omega$ meson mass vs density for $eB = 4 m_\pi^2$).

---

## 1. Physics Model & Equations Solved

The script models the in-medium mass $m_\omega^*$ by dynamically solving the coupled Finite Energy Sum Rules (FESR) quartic equation:

$$\frac{(c_0)^2}{12} (s_0^{*})^4 + \frac{c_0 C_A}{3} (s_0^{*})^3 + c_0 c_2 (s_0^{*})^2 + c_0 c_3 s_0^{*} + \left( C_A c_3 - c_2^2 \right) = 0$$

where:
*   **$c_0$:** Perturbative QCD coefficient $c_0 = 1 + \alpha_s/\pi \approx 1.1146$.
*   **$C_A$:** Incorporates the mass correction $c_1$ and the $\omega$-nucleon scattering term $\Pi^\omega(0) = \rho_B / (4 M_N)$ (scaled by current coupling $d_\omega = 1/6$):
    $$C_A = c_1 - \frac{12\pi^2}{d_\omega} \Pi^\omega(0) \approx -0.000195 - 710.616 \times \left( \frac{\rho_B}{3.752} \right)$$
*   **$c_2$ (OPE Dimension 4):**
    $$c_2 = \frac{\pi^2}{3} \langle \frac{\alpha_s}{\pi} G^2 \rangle^* + 4\pi^2 (m_u \langle \bar{u}u \rangle^* + m_d \langle \bar{d}d \rangle^*)$$
*   **$c_3$ (OPE Dimension 6):**
    $$c_3 = - \alpha_s \pi^3 \frac{448}{81} \kappa_q \left( \langle \bar{u}u \rangle^{*2} + \langle \bar{d}d \rangle^{*2} \right)$$

Once the positive real root of the quartic equation $s_0^*$ is found, the mass is computed:
$$m_\omega^* = \sqrt{\frac{\frac{(s_0^*)^2 c_0}{2} - c_2}{c_0 s_0^* + C_A}}$$

---

## 2. Step-by-Step Logic

1.  **Vacuum Pre-calculation:**
    *   Computes the vacuum threshold $s_0^\omega \approx 1.272\text{ GeV}^2$.
    *   Calculates the factorization parameter $\kappa_q \approx 7.695$ self-consistently.
2.  **Evaluating the Medium:**
    *   Iterates over the density range $x \in [0, 2.0]$ (where $x = \rho_B/\rho_0$).
    *   Evaluates the rational scale factor for the nonstrange field $\sigma$:
        $$\sigma = \sigma_0 \left( 1 - \frac{a x + b x^2}{1 + c x} \right)$$
        using the curve-specific fits:
        *   **$\eta=0$ (with AMM):** $a=0.486572, b=-0.077940, c=-0.057508$
        *   **$\eta=0.3$ (with AMM):** $a=0.470063, b=-0.064288, c=-0.017268$
        *   **$\eta=0.5$ (with AMM):** $a=0.671448, b=0.544399, c=2.097651$
        *   **$\eta=0$ (without AMM):** $a=0.487945, b=-0.074682, c=-0.058666$
        *   **$\eta=0.3$ (without AMM):** $a=0.470189, b=-0.071593, c=-0.051749$
        *   **$\eta=0.5$ (without AMM):** $a=0.468673, b=-0.076350, c=-0.040097$
    *   Scales fields $\zeta, \delta, \chi$:
        $$\zeta = \zeta_0 \left( 1 - \frac{0.02 x}{1 + 0.1 x} \right), \quad \delta = \frac{0.008 \eta x}{1 + 0.1 x}, \quad \chi = \chi_0 (1 - 0.005 x)$$
    *   Computes condensates, OPE coefficients $c_2, c_3$, and scattering term $C_A$.
3.  **Quartic Solving:**
    *   Evaluates roots of the quartic polynomial and selects the root closest to the vacuum threshold.
4.  **Plotting:**
    *   Draws the curve using standard color conventions: red ($\eta=0$), blue ($\eta=0.3$), and green ($\eta=0.5$).
    *   Draws solid dash-dot lines for with-AMM, and dotted lines for without-AMM.
    *   Saves the figure as `plots/plot1.png`.
