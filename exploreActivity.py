# %% 

%config InlineBackend.print_figure_kwargs

from simple_abc_truncatedGaus import Model, basic_abc, pmc_abc
import pylab as plt
import numpy as np
from scipy import stats
import sys
sys.path.insert(0,'../adLIF_NEST/')

from func.helpers import *
from func.network import *
import matplotlib.pyplot as plt
from func_.helpers import *
#import seaborn as sns|
#sns.set_context('talk')
na = np.array
# plt.style.use('dark_background')
import seaborn as sns
import pandas as pd

# %%
from matplotlib import rc
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
# plt.style.use('ggplot')
# from func_.helpers import *

# %%
# Start with 20% of inhibition 
params= {'J': 9.0,
 'g': 4.0,
 'N': 1000,
 'epsilon': 0.2,
 'eta': 0.0,
 'p_rate': 2077.4792278063533,
 'J_ext': 3.820498723458609,
 'tauMem': 20.0,
 'CMem': 1.0,
 'theta': 20.0,
 'V_res': 10.0,
 'Ks': [33,8 ],
 'V_m': 0.0,
 'b': 1.0,
 'a': 0.0,
 'tau_w': 17000.0,
 'p': 0.1,
 't_ref': 2.0,
 'd': 3.5,
 'N_rec': 1000,
 'voltage': False,
 'chunk': False,
 'chunk_size': 1000.0,
 'directory': 'sim/ABC_small/',
 'simulation': 'hash',
 'simtime': 500000.0,
 'master_seed': 1000,
 'dt': 1.0,
 'threads': 20}
KE = 0.5490304776651601*((params['N']*0.14159415833146138)+32.43607416673515)
params['Ks'] = [int(KE), int(KE/params['g'])]
def get_model(params,conn= False):
    """ 
    Load/Simulate adLIF model
    Args:
            [1] (dict) params
    Returns:
            st,gid (tuple) spike times, unit ids 
    """
    visual = False
    try: 
        name = get_hash(params)
        st,gid = read_gdf('sim/ABC_small/',name,(5000,params['simtime']),threads = params['threads'])
        if conn:
            with h5py.File('sim/ABC_small/'+name+'-conn','r') as f:
                conn = na(f['conn'])
                return st,gid,conn
    except:
        if visual:
            raise LoacError('wont load')
        else:
            print('Simulating...')
            A = adLIFNet(params); 
            A.build()
            A.connect()
            A.run()
        try:
            name = get_hash(params)
            st,gid = read_gdf('sim/ABC_small/',name,(5000,params['simtime']),threads =  params['threads'])
        except:
            st,gid = [0],[0]
        if conn:
            conn =A.get_connectivity()
            return st,gid,conn
    return st,gid,[]


#%%
params['g'] = 0.
params['p_rate'] = 135.#480.
params['J'] =1.5# #1.21
params['J_ext'] =1. #params['J']#3.5
params['tau_w'] =4000.
params['b'] = .4#.4
params['simtime']= 500000
params['i_noise_off'] = 1
# %%

# %time 
st,gid,conn = get_model(params,conn=False)
bin_size = 20
sc,_ = np.histogram(st,np.arange(0,params['simtime'],bin_size))
# params['g'] = 4.
# st1,gid1,conn = get_model(params)
# sc1,_ = np.histogram(st1,np.arange(0,params['simtime'],bin_size))

# plt.plot(sc)
# plt.plot(sc1)

#%%

A = np.load('/home/ovinogradov/Documents/ABC_bursting/ABC_results/N=1000_balance_0.2fit_15_multi_fit_KEmodel_Bw0.25_bFit_25_240220meta.npy')

#%%
data_path = 'sim/ABC_small'
simulation = get_hash(params)

filenames = [i for i in listdir(data_path) if simulation in i and 'gdf' in i]

print(filenames)

with open(data_path+'/'+simulation+'-all.gdf', 'w') as outfile:
     for fname in filenames:
        with open(data_path+'/'+fname) as infile:
            for line in infile:
                outfile.write(line)
#%% check if can read from one file

st,gid = read_gdf('sim/ABC_small/',simulation,(5000,params['simtime']),threads = params['threads'])
from_file_pandas(data_path+simulation,params['simtime'])

def read_gdf(mypath,simulation,t,threads=4):
    """Read spikes from native NEST recordings
    Args:
            mypath (str): Directory of the simulation with "/".
            simulation (str): Simulation name
            t (tuple: int,int): time of the simulation to load (ms)
            
    Returns:
            arr: spikes timestamp 
            arr: corresponding unit ids
    
    """
    t1,t2 = t
    # list all files (for multithreading)
    files = [mypath +i for i in listdir(mypath) if simulation in i and 'gdf' in i]#
    concFile = [i for i,f in enumerate(files) if '-all' in f]
    if len(concFile)==1:
        #reading from one file
        print('read from all')
        data = from_file_pandas(files[concFile[0]],t2)
    else:

        #print('Simulations: %s'% simulation)
        #print('Reading from %s files'%len(files))
        if len(files)>threads:
            print(len(files))
            raise SizeError()
        data = from_file_pandas(files,t2)

    ts,gids = data[:,1],data[:,0]
    
    ts1 = ts[(ts>t1)*(ts<t2)]
    gids1 = gids[(ts>t1)*(ts<t2)]
    
    return ts1,gids1

# %%
def plot_st(st,gid):
    plt.figure()
    # plt.subplot(211)
    st_plot = st[gid<50]
    gid_plot = gid[gid<50]
    plt.plot(st_plot,gid_plot,'|',markersize = 0.5)
    # plt.xlim([0,250000])
    # plt.subplot(212)
    # plt.plot(st1,gid1,'|',markersize = 0.5)

# %%
def print_features(st,gid,st1,gid1):
    bursts = MI_bursts(st, 
        maxISIstart=4.5,
        maxISIb=4.5,
        minBdur=10,
        minIBI=20,
        minSburst=100)
    ibis =na(bursts)[1:,0]-na(bursts)[:-1,1]
    dur = np.mean(na(bursts)[:,1]- na(bursts)[:,0])
    mibis = np.mean(ibis)
    cvs = np.std(ibis)/mibis
    print('E only')
    print('Duration',dur)
    print('CV',cvs)
    print('mIBI',mibis)

    bursts1 = MI_bursts(st1, 
        maxISIstart=4.5,
        maxISIb=4.5,
        minBdur=10,
        minIBI=20,
        minSburst=100)
    ibis =na(bursts1)[1:,0]-na(bursts1)[:-1,1]
    dur = np.mean(na(bursts1)[:,1]- na(bursts1)[:,0])
    mibis = np.mean(ibis)
    cvs = np.std(ibis)/mibis
    print('E+I')
    print('Duration',dur)
    print('CV',cvs)
    print('mIBI',mibis)


# %%
def filter_bursts(st,gid):
    bursts1 = MI_bursts(st, 
        maxISIstart=4.5,
        maxISIb=4.5,
        minBdur=10,
        minIBI=20,
        minSburst=100)

    for burst in bursts1:
        mask = (st<(burst[1]+100)*(st>(burst[0])-100))==False
        st = st[mask]
        gid = gid[mask]

    return st,gid
def get_input(st,gid,conn,params):
    bin_size = 100
    input_spikes = np.zeros([len(np.arange(0,params['simtime'],bin_size))-1,2])#[]
    NE = params['N']-(params['N']*params['epsilon'])
    nid = 4
    for nid in range(100):
        input_neurons = conn[conn[:,1]==nid,0]

        for inn in input_neurons:
            sc_,_ = np.histogram(st[gid==inn],np.arange(0,params['simtime'],bin_size))
            if inn<NE:
                input_spikes[:,0] +=sc_*params['J']
    #         elif inn==params['N']+1:
    #             input_spikes[:,2] +=sc_*(params['J_ext'])
            else:
                input_spikes[:,1] +=sc_*(-params['J']*params['g'])
    input_spikes/=100
    return input_spikes


# %%
st_,gid_ = filter_bursts(st,gid)
st1_,gid1_ = filter_bursts(st1,gid1)
params['g'] = 0.
input1 = get_input(st_,gid_,conn,params)
params['g'] = 4.
input2 = get_input(st1_,gid1_,conn,params)
# %%

plt.plot(input1)

# %%
input_spikes_= np.hstack(input1)
input_spikes_ =input_spikes_[input_spikes_!=0]
# input_spikes_ = input_spikes_[np.abs(input_spikes_)<5]

plt.hist(input_spikes_,bins = 20,density = True)
# input_spikes_ = input_spikes[:,1]
# %%
with h5py.File('sim/ABC_small/'+name+'-conn','r') as f:
    conn = na(f['conn'])
# plt.hist(-input_spikes_[(input_spikes_<15)*(input_spikes_>0)])
input_spikes_= np.hstack(input2)
input_spikes_ =input_spikes_[input_spikes_!=0]
# input_spikes_ = input_spikes_[np.abs(input_spikes_)<5]
plt.hist(input_spikes_,bins = 20, alpha =0.4,density = True)

# %%

