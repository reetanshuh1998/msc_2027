# Self-contained code for Plot 6
# Solving QCD Sum Rule physical equations dynamically using density-dependent scalar fields.
import os
import numpy as np
import matplotlib.pyplot as plt

def solve_quartic_s0(c0, C_A, c2, c3, s0_vac):
    coeffs = [
        (c0**2) / 12.0,
        (c0 * C_A) / 3.0,
        c0 * c2,
        c0 * c3,
        C_A * c3 - c2**2
    ]
    roots = np.roots(coeffs)
    real_roots = roots[np.isreal(roots)].real
    positive_roots = real_roots[real_roots > 0]
    if len(positive_roots) == 0:
        return None
    idx = np.argmin(np.abs(positive_roots - s0_vac))
    return positive_roots[idx]

# Model constants
m_u = 0.004
m_d = 0.007
m_s = 0.150
f_pi = 0.093
m_pi = 0.137
m_k = 0.495
f_k = 0.115

sigma_0 = f_pi
zeta_0 = (2.0 * f_k - f_pi) / np.sqrt(2.0)
chi_0 = 0.4099
d = 0.06

# Vacuum condensates
q_cond_0 = - (0.250**3) # GeV^3
s_cond_0 = 0.8 * q_cond_0
g_cond_0 = 0.012 # GeV^4

alpha_s = 0.36
c0 = 1.0 + alpha_s / np.pi
c1 = -0.135

m_vac = 1.02
d_V = 0.3333333333333333
kappa = 6.080915758652768

# Vacuum thresholds
if "phi" == "omega":
    s0_vac = 1.2721241715620724
elif "phi" == "rho":
    s0_vac = 1.233202045636193
else:
    s0_vac = 1.8376430384727902

# Embedded scalar field polynomial fits from the Chiral SU(3) model
fits_data = [{'eta': 0.0, 'amm': True, 'a': 0.14948587177894523, 'b': -0.009812513027303342, 'c': 1.1002609813717434, 'xmin': 9.99999999995449e-06, 'xmax': 4.0}, {'eta': 0.3, 'amm': True, 'a': 0.1459357396403208, 'b': -0.008139185774762202, 'c': 1.1184507670882247, 'xmin': 9.99999999995449e-06, 'xmax': 4.0}, {'eta': 0.5, 'amm': True, 'a': 0.14440465916844003, 'b': -0.0060740227445238585, 'c': 1.1815350842664614, 'xmin': 9.99999999995449e-06, 'xmax': 4.0}, {'eta': 0.0, 'amm': False, 'a': 0.14539508945283763, 'b': -0.01379963688399726, 'c': 0.9559731298712072, 'xmin': 9.99999999995449e-06, 'xmax': 4.0}, {'eta': 0.3, 'amm': False, 'a': 0.14434867559320372, 'b': -0.01157006234227391, 'c': 0.9970515025031663, 'xmin': 9.99999999995449e-06, 'xmax': 4.0}, {'eta': 0.5, 'amm': False, 'a': 0.1445338685598149, 'b': -0.009573052656687996, 'c': 1.0605108640942982, 'xmin': 9.99999999995449e-06, 'xmax': 4.0}]

def main():
    os.makedirs('plots', exist_ok=True)
    plt.figure(figsize=(7, 6), dpi=300)
    
    validation_data = []
    
    for fit in fits_data:
        eta = fit['eta']
        amm = fit['amm']
        a = fit['a']
        b = fit['b']
        c = fit['c']
        xmin = fit['xmin']
        xmax = fit['xmax']
        
        # Grid of density ratio points
        x_phys = np.linspace(xmin, xmax, 100)
        y_phys_solved = []
        
        for x in x_phys:
            scale = 1.0 - (a * x + b * x**2) / (1.0 + c * x)
            if "phi" != "phi":
                sigma = sigma_0 * scale
                zeta = zeta_0 * (1.0 - 0.02 * x / (1.0 + 0.1 * x))
                delta = eta * 0.008 * x / (1.0 + 0.1 * x)
            else:
                sigma = sigma_0 * (1.0 - 0.25 * x / (1.0 + 0.15 * x))
                zeta = zeta_0 * scale
                delta = eta * 0.008 * x / (1.0 + 0.15 * x)
                
            chi = chi_0 * (1.0 - 0.005 * x)
            
            # 1. In-medium condensates:
            u_cond = q_cond_0 * (sigma + delta) / sigma_0
            d_cond = q_cond_0 * (sigma - delta) / sigma_0
            s_cond = s_cond_0 * zeta / zeta_0
            g_cond = g_cond_0 * (chi / chi_0)**4
            
            # 2. OPE Coefficients:
            if "phi" != "phi":
                c2 = (np.pi**2 / 3.0) * g_cond + 4.0 * np.pi**2 * (m_u * u_cond + m_d * d_cond)
                c3 = -alpha_s * np.pi**3 * (448.0 / 81.0) * kappa * (u_cond**2 + d_cond**2)
            else:
                c2 = (np.pi**2 / 3.0) * g_cond + 8.0 * np.pi**2 * m_s * s_cond
                c3 = -alpha_s * np.pi**3 * (224.0 / 81.0) * kappa * 8.0 * (s_cond**2)
                
            # 3. Scattering term:
            rho_0_phys = 0.15 * (0.1973)**3
            rho_B_phys = x * rho_0_phys
            pi_0 = rho_B_phys / (4.0 * 0.938) if "phi" != "phi" else 0.0
            C_A = c1 - (12.0 * np.pi**2 / d_V) * pi_0
            
            # 4. Solve FESR quartic equation:
            s0_solved = solve_quartic_s0(c0, C_A, c2, c3, s0_vac)
            if s0_solved is not None:
                mass_solved = np.sqrt((0.5 * (s0_solved**2) * c0 - c2) / (c0 * s0_solved + C_A))
                y_phys_solved.append(mass_solved * 1000.0)
            else:
                y_phys_solved.append(np.nan)
                
        y_phys_solved = np.array(y_phys_solved)
        
        # Line styling
        if eta == 0.0:
            color = 'red'
            if amm:
                linestyle = '-.' # dash-dot
                label = r'$\eta=0$'
            else:
                linestyle = ':' # dotted
                label = r'$\eta=0$ (without AMM)'
        elif eta == 0.3:
            color = 'blue'
            if amm:
                linestyle = '--' # dashed
                label = r'$\eta=0.3$'
            else:
                linestyle = ':'
                label = r'$\eta=0.3$ (without AMM)'
        else:
            color = 'green'
            if amm:
                linestyle = (0, (3, 1, 1, 1, 1, 1)) # long-dash-short-dash
                label = r'$\eta=0.5$'
            else:
                linestyle = ':'
                label = r'$\eta=0.5$ (without AMM)'
                
        plt.plot(x_phys, y_phys_solved, linestyle=linestyle, color=color, linewidth=1.5, label=label)
        
        val_at_rho0 = np.interp(1.0, x_phys, y_phys_solved)
        validation_data.append((eta, amm, val_at_rho0))
        
    plt.xlabel(r'$\rho_B / \rho_0$', fontsize=12)
    
    if "phi" == "omega":
        plt.ylabel(r'$m_\omega^*$ (MeV)', fontsize=12)
        plt.ylim(750, 1050)
        plt.xlim(0, 2.0)
    elif "phi" == "rho":
        plt.ylabel(r'$m_\rho^*$ (MeV)', fontsize=12)
        plt.ylim(300, 800)
        plt.xlim(0, 4.0)
    else:
        plt.ylabel(r'$m_\phi^*$ (MeV)', fontsize=12)
        plt.ylim(995, 1025)
        plt.xlim(0, 4.5)
        
    plt.title(f'$eB = 12 m_\pi^2$', fontsize=12, pad=10)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(loc='lower left' if "phi" == "omega" else 'upper right', frameon=True, edgecolor='black', fontsize=10)
    plt.tight_layout()
    plt.savefig('plots/plot6.png', dpi=300)
    plt.close()
    
    print("Verification values at saturation density (rho_B = rho_0):")
    for eta, amm, val in sorted(validation_data, key=lambda x: (x[1], x[0])):
        amm_str = "with AMM" if amm else "without AMM"
        print(f"  eta={eta:.1f} ({amm_str}): {val:.2f} MeV")

if __name__ == '__main__':
    main()
