import numpy as np
from scipy.integrate import quad
from scipy.special import gamma
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Define the s function
def s(N, t, y, h, delta, nu):
    # t is in the units of 1/(2 alpha)
    lmbda = np.exp(h*t)/(4 * N * delta**2 * np.cos(np.pi * nu / 2))
    def integrand(yl):
        return (yl**(2 * delta - 1)) * np.exp(-lmbda * (y * yl)**h - yl)
    
    integral, _ = quad(integrand, 0, np.inf)
    return 0.5 * (1 - ((np.cos(np.pi * nu / 2)**(2 * delta)) / gamma(2 * delta)) * integral)

# Define the ArgQ function
def ArgQ(y, delta, nu):
    return np.sin(np.pi * nu / 2) * y - np.pi * nu * delta

def ArgQ_analytic(s, N, t, h, delta, nu):
    lmbda = np.exp(h*t)/(4 * N * delta**2 * np.cos(np.pi * nu / 2))
    s0 = (1 - np.cos(np.pi * nu / 2)**(2*delta))/2
    y = (gamma(2*delta+1)/(np.cos(np.pi * nu / 2)**(2*delta) * gamma(2*delta + h) * delta))**(1/h)*(s-s0)**(1/h)/lmbda**(1/h)
    return np.sin(np.pi * nu / 2) * y - np.pi * nu * delta

# Define the P function
def P(N, t, y, h, delta, nu):
    # t is in the units of 1/(2 alpha)
    lmbda = np.exp(h*t)/(4 * N * delta**2 * np.cos(np.pi * nu / 2))
    def integrand(yl):
        return ((yl**(2 * delta + h - 1)) * np.exp(-lmbda * (y * yl)**h - yl))
    
    integral, _ = quad(integrand, 0, np.inf)
    numerator = (2 * y**(2 * delta - h)) * np.exp(-y * np.cos(np.pi * nu / 2))
    denominator = lmbda * h * integral
    
    return numerator / denominator

def P_analytic(s, N, t, h, delta, nu):
    lmbda = np.exp(h*t)/(4 * N * delta**2 * np.cos(np.pi * nu / 2))
    s0 = (1 - np.cos(np.pi * nu / 2)**(2*delta))/2
    y = (gamma(2*delta+1)/(np.cos(np.pi * nu / 2)**(2*delta) * gamma(2*delta + h) * delta))**(1/h)*(s-s0)**(1/h)/lmbda**(1/h)
    return np.sin(np.pi * nu / 2) * y - np.pi * nu * delta

def plot_P(lmbda, delta, nu,logY):
    plt.rcParams.update({'font.size': 12, 'figure.figsize': (4.5,4.5)})
    # Generate values for y
    y_values = np.linspace(0, 15, 10000)
    h_arr = [0.4, 0.5, 0.8, 1]
    num_hs = len(h_arr)
    norm = Normalize(vmin=0, vmax=num_hs - 1) 


    for (i,h) in enumerate(h_arr):
        P_values = [P(lmbda, y, h, delta, nu) for y in y_values]
        s_values = [s(lmbda, y, h, delta, nu) - s(lmbda, 0, h, delta, nu) for y in y_values]
        plt.plot(s_values, P_values, label=f"h = {h}", color = plt.get_cmap("cool")(norm(i)))

    plt.xlabel(r"$s-s_0$", fontsize=14)
    plt.ylabel(r"$P(s)$", fontsize=14)
    plt.legend(fontsize=12)
    plt.xlim(0, 0.015)
    if logY:
        plt.yscale('log')
        plt.ylim(10**(-2),10**(4))
        plt.savefig(f'P_logY_lmbda_{lmbda}_delta_{delta}_nu_{nu}.pdf', bbox_inches='tight')
    else: 
        plt.ylim(0, 25)
        plt.savefig(f'P_lmbda_{lmbda}_delta_{delta}_nu_{nu}.pdf', bbox_inches='tight')
    plt.close()


def plot_PQ(N, t, delta, nu):
    gr = (5**.5 - 1)/2
    plt.rcParams.update({'font.size': 8, 'figure.figsize': (3+3/8, (3+3/8))})
    fig, (ax2,ax1) = plt.subplots(2,1, sharex= True)

    # Generate values for y
    y_values1 = np.linspace(0, 15, 1000)
    h_arr = [0.4, 0.5, 0.7, 1]
    num_hs = len(h_arr)
    # norm = Normalize(vmin=0, vmax=num_hs - 1) 
    norm = [0.3 + i*0.6/(num_hs-1) for i in range(num_hs)]

    # plot_cmap = "cool"
    plot_cmap = "BuGn"

    for (i,h) in enumerate(h_arr):
        ArgQ_values = [ArgQ(y, delta, nu) - ArgQ(0, delta, nu) for y in y_values1]
        s_values = [s(N, t, y, h, delta, nu) - s(N, t, 0, h, delta, nu) for y in y_values1]
        s0 = (1 - np.cos(np.pi * nu / 2)**(2*delta))/2
        ArgQ_analyticValues = [ArgQ_analytic(s+s0, N, t, h, delta, nu) - ArgQ_analytic(s0, N, t, h, delta, nu) for s in s_values]
        plt.plot(s_values, ArgQ_values, label=f"h = {h}", color = plt.get_cmap("BuGn")(norm[i]))
        # plt.plot(s_values, ArgQ_analyticValues, 'x', color = plt.get_cmap("BuGn")(norm[i]), markersize = 0.1)

    ax1.set_xlabel(r"$s-s_0$")
    ax1.set_ylabel(r"$\mathrm{Arg}~q(s, t)-\mathrm{Arg}~q(s_0, t)$")
    ax1.legend(fontsize=7)
    ax1.set_xlim(0, 0.015)
    ax1.set_ylim(0,10)
    ax1.text(-0.004, 10, '(b)')
    # ax1.grid()
    
    # Generate values for y
    y_values2 = np.linspace(0, 15, 10000)

    for (i,h) in enumerate(h_arr):
        P_values = [P(N, t, y, h, delta, nu) for y in y_values2]
        s_values = [s(N, t, y, h, delta, nu) - s(N, t, 0, h, delta, nu) for y in y_values2]
        ax2.plot(s_values, P_values, label=f"h = {h}", color = plt.get_cmap(plot_cmap)(norm[i]))

    # ax2.set_xlabel(r"$s-s_0$")
    # ax2.grid()
    ax2.set_ylabel(r"$p(s,t)$")
    ax2.legend(fontsize=7)
    # ax2.set_xlim(0, 0.015)
    ax2.set_yscale('log')
    ax2.set_ylim(10**(-2),10**(4))
    ax2.text(-0.004, 10**4, '(a)')
    fig.tight_layout()
    plt.savefig(f'PQ_N_{N}_t_{t}_delta_{delta}_nu_{nu}.pdf', bbox_inches='tight')
    plt.close()
    
# Parameters
N = 3000 
t = 0.9
delta = 1/6
nu = 0.5
plot_PQ(N, t, delta, nu)