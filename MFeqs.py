import matplotlib.pyplot as plt
from scipy.integrate import odeint,solve_ivp,ode
from numpy import arange
from scipy.optimize import fsolve as fsolve
import numpy as np
from math import erf as erfp
from math import exp as exp
import seaborn as sns
from scipy.optimize import root as root
import scipy
# from func.intersect  import *
from scipy.signal import find_peaks
sns.set_style('whitegrid')
na = np.array
from math import erf
from scipy.integrate import quad
from scipy.special import erfcx,ndtr,erf
sns.set_style('ticks')
sns.set_context('poster')
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

import matplotlib
# matplotlib.use('PDF')
import matplotlib.pylab as plt
from matplotlib import rc
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
pi_sqrt = np.sqrt(np.pi);
from func.eq_helpers.FI import *


import numpy as np
from scipy import stats
import sys
sys.path.insert(0,'/home/ovinogradov/Projects/EI_project/adLIF_NEST')
from func.helpers import *
from func.network import *
na = np.array
# %

j= 2.;
N = 1000#40000
epsilon = 0.2
g_mod =4.
NI = N*epsilon
NE = N-NI
g = NE/NI
p = 0.1
nu_thr = 20./(NE*p*0.2*0.01)

aparams = {'J':j,#(j*np.sqrt(5))/np.sqrt((N-(epsilon*N))*0.1),
        'g':4.,#g*g_mod,
        'N':N,
        'epsilon':epsilon,
        'gamma':NI/NE,
        'Ke':N*p,
        'p_rate':700.,#2*nu_thr,#((ex_rate/80)*((N-(epsilon*N))*0.1)),
        'J_ext':1.0,
        'tauMem' :0.02,
        'CMem':1.,
        'theta':20.,
        'V_res':10.,
        'V_m':0.0,
        'b' :0.3,#2./40,
        'a':0.,
        'tau_w':4.,#0.500,
        'p':p,
        't_ref':0.002,
        'd':0.0035,
        }


t =np.arange(0,1000,0.01)    
conds = np.arange(0,100.,1)#[0.,1.,5.,10.,100.,]
sol = [odeint(ps_dyn, cond, t, args =(aparams,)) for cond in conds]
[plt.plot(t[:1000],s[:1000]) for s in sol]
#print([s[-1] for s in sol])
plt.show()


sns.set_context('talk')
sns.set_style('ticks')
mus  = np.arange(0,60)
plt.plot([F(i,1,params) for i in mus ],label ='$\\sigma=1 $')
plt.plot([F(i,3,params) for i in mus],label ='$\\sigma=3 $')
plt.plot([F(i,10,params) for i in mus],label ='$\\sigma=10 $')
plt.plot([0,50],[0,50],'-',color = 'gray')
plt.xlim([-5,55])
plt.ylim([-5,65])
plt.xlabel('input $\\mu$')
plt.ylabel('output rate $\\nu$')
plt.legend()
sns.despine(trim=1)
plt.tight_layout()
plt.show(block=False)

t =np.arange(0,1000,0.01)    
sol = odeint(ps_dyn, (0.2), t, args =(params,))
# scipy.integrate.RK45(ps_dyn, 0., [10.], t_bound= 100.)
# nu,a = sol.T
nu_overG = np.zeros((3,len(np.arange(0.,6.,0.1))))
nus = [450.,500.,600.]
for eta_i, eta in enumerate(nus):
    params['p_rate']= eta#eta*nu_thr
    for g_i,g in enumerate(np.arange(0.,6.,0.1)):
        params['g']= g
        sol = odeint(ps_dyn, (10.2), t, args =(params,))[-1]
        nu_overG[eta_i, g_i ]= sol
        
for i in range(len(nus)):    
    plt.plot(np.arange(0.,6.,0.1), nu_overG[i],label="Input per Ke= %s Hz"%(nus[i]) )
plt.legend()
plt.xlabel('g')
plt.ylabel('nu_0 (Hz)')
plt.tight_layout()
plt.show(block = False)


#%%
# %config InlineBackend.figure_format = 'png'

params= {'J': 2.0,#10.0,
 'g': 4.0,
 'N': 1000,
 'epsilon': 0.2,
 'eta': 0.0,
 'p_rate': 700.,#np.random.uniform(2070,2080,1)[0],#2077.4792278063533,
 'J_ext': 1.,#3.820498723458609,
 'tauMem': 20.0,
 'CMem': 1.0,
 'theta': 20.0,
 'V_res': 10.0,
 'Ks': [80,20],
 'V_m': 0.0,
 'b': 0.3,#1.0,
 'a': 0.0,
 'tau_w': 4000.,#10000.0,#17000.0,
 'p': 0.1,
 't_ref': 2.0,
 'd': 3.5,
 'N_rec': 1000,
 'voltage': False,
 'chunk': False,
 'chunk_size': 1000.0,
 'directory': 'sim/',
 'simulation': 'hash',
 'simtime':200003.0,
 'master_seed': 1000,
 'dt': 0.5,
 'threads': 20}


def prepare_data(params, bin_size =20):
    """ 
    minifunction to read the data, detect bursts and calculate the spike
    counts
    Args:
            params(dict): see adLIF network

    Returns:
            tuple: bursts,spike counts, spike count times
    """
    name = get_hash(params)
    print('file:', name)
    st,gid = read_gdf(params['directory'],name,(5000,params['simtime']),threads = params['threads'])
    bursts = MI_bursts(st)
    sc,times =np.histogram(st,np.arange(0,params['simtime'],bin_size))
    return bursts, sc, times[1:]

force = False
try: 
    if force:
        raise Exception('force simulate')
    bursts1,sc1,times1 = prepare_data(params)
except:
    A = adLIFNet(params)
    A.build()
    A.connect()
    A.run()
    bursts1,sc1,times1 = prepare_data(params)

# %%
bursts1,sc1,times1 = prepare_data(params)
# %%
plt.figure()
plt.plot(st,gid,'.', markersize = 0.8)
plt.show(block = False)
print(st[1:20])
print('done')
# %%
