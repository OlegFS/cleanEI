# MF-equtions. Following Brunel et al. 2000, Tartaglia Brunel
from scipy.integrate import odeint,solve_ivp,ode
from numpy import arange
from scipy.optimize import fsolve as fsolve
import numpy as np
from math import erf as erfp
from math import exp as exp
from scipy.optimize import root as root
import scipy
from scipy.signal import find_peaks
na = np.array
from math import erf
from scipy.integrate import quad
from scipy.special import erfcx,ndtr,erf
pi_sqrt = np.sqrt(np.pi);
def F(mu,sigma,params):
    """ 
    FI for LIF neuron under white noise
    Args:
        mu (float): mean input
        std (float): variance of the input
        params (dict): paramters of the network

    Returns:
        rate(float): output rate
    """
    theta = params['theta']
    Vr = params['V_res']#/theta
    tau_m = params['tauMem']#/1000
#     mu/=theta
#     sigma/=theta
#     r= tau_m*pi_sqrt
    def integrand(x):
        return erfcx(-x)#exp((x**2))*(1+erfp(x))#
    upper=(theta-mu)/sigma
    lower=(Vr-mu)/sigma
#     print(upper,lower)
    I= quad(integrand, lower, upper)[0]
    # add tau_rp??
    tau_rp = params['t_ref']#/1000
    nu0 = tau_rp +tau_m*pi_sqrt*I#
    return 1/nu0

def mu_(nu0, params):
#     print(nu0)
    nu_ext=params['p_rate']
    Ke = params['Ke']
    tau_m = params['tauMem']#/1000
    tau_w = params['tau_w']
    g = params['g']
    gamma = params['gamma']
    J = params['J']
    J_ext = params['J_ext']
    b = params['b']
    mu_ext = nu_ext*J_ext*tau_m
    mu_rec =tau_m*Ke*J*nu0*(1-(g*gamma))#tau_m*b
#     print(mu_ext,mu_rec)
    return  mu_ext+mu_rec-(b*nu0)#*100

def mu_w(nu0,w, params):
#     print(nu0)
    nu_ext=params['p_rate']
    Ke = params['Ke']
    tau_m = params['tauMem']#/1000
    tau_w = params['tau_w']
    g = params['g']
    gamma = params['gamma']
    J = params['J']
    J_ext = params['J_ext']
    mu_ext = nu_ext*J_ext*tau_m
    mu_rec =tau_m*Ke*J*nu0*(1-(g*gamma))
#     print(mu_ext+mu_rec)
    return  mu_ext+mu_rec-w #*100    

def sigma_sq(nu0,params):
    nu_ext=params['p_rate']
    Ke = params['Ke']
    tau_m = params['tauMem']#/1000
    g = params['g']
    gamma = params['epsilon']
    J = params['J']
    J_ext = params['J_ext']
    sigma_ex =tau_m*nu_ext*(J_ext**2)
    sigma_rec = (tau_m*Ke*(J**2)*nu0)*(1+((g**2)*gamma))
    return np.sqrt((sigma_ex+sigma_rec))
    
# def adap(a,nu,params):
#     b= params['b']
#     return -(a/0.5)+(b*nu)
    
def ps_dyn(nu,t,params):
    """ 
    Pseudo dynamics of the FI 
    nu' = -nu + F(mu,sigma)
    Args:
            path (str): Directory of the simulation with "/".

    Returns:
            list: Burst times 
    """
    mu = mu_(nu,params)
    sigma = sigma_sq(nu,params)
    return -nu+F(mu,sigma,params)

