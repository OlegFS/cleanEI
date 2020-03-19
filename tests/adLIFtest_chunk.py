#Test simulation with different chunk sizes 
# Using Prepare, Run, Cleanu speeds up simulations to the nest.Simulate speed
# 03.03.20
# %config InlineBackend.figure_format = 'png'
import os

#os.chdir('/home/ovinogradov/Documents/tests')
import numpy as np
from scipy import stats
import sys
sys.path.insert(0,'../../cleanEI')

from func.helpers import get_hash
from func.helpers import read_gdf
from func.network import *
# from func_.helpers import *
#import seaborn as sns|
#sns.set_context('talk')
na = np.array
# plt.style.use('dark_background')

# %
epsilon = 0.2
N= 1000
NI = np.int(N*epsilon)
NE = N-NI
g = 4.#NE/NI
KE = 0.5490304776651601*((NE*0.14159415833146138)+32.43607416673515)
nest.set_verbosity('M_FATAL')
params= {'J': 2.0,#10.0,
 'g': 4.0,
 'N': 1000,
 'epsilon': epsilon,
 'eta': 0.0,
 'p_rate':800.,#np.random.uniform(2070,2080,1)[0],#2077.4792278063533,
 'J_ext': 1.,#3.820498723458609,
 'tauMem': 20.0,
 'CMem': 1.0,
 'theta': 20.0,
 'V_res': 10.0,
 'Ks': [int(KE),int(KE/4)],
 'V_m': 0.0,
 'b': 0.4,#1.0,
 'a': 0.0,
 'tau_w': 5000.,#10000.0,#17000.0,
 'p': 0.1,
 't_ref': 2.0,
 'd': 3.5,
 'N_rec': 1000,
 'voltage': False,
 'chunk': True,
 'chunk_size':105000.0,
 'directory': 'sim/ABC_small/',
 'simulation': 'hash',
 'simtime': 500500.0,
 'master_seed': 1000,
 'dt': 0.5,
 'threads': 20}

 # %%
A = adLIFNet(params)
A.build()
A.connect()
A.run()
# %%
name = get_hash(params)
st,gid = read_gdf('sim/ABC_small/',name,(5000,params['simtime']),threads = params['threads'])
print('Firing Rate is: ',len(st)/500000/200)
# %%
#plt.figure()
#plt.plot(st,gid,'.', markersize = 0.8)
#plt.show()
#plt.savefig('foo.png', bbox_inches='tight')
#print(st[1:20])
#print('done')
# %%
