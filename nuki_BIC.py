# Check posteriors of ABC fitting
#
import numpy as np
from scipy import stats
import sys
sys.path.insert(0,'/home/ovinogradov/Documents/adLIF_NEST/')
sys.path.insert(0,'/home/ovinogradov/Documents/ABC_bursting/')
from func.helpers import *
from func.network import *
from func_.model_setup_experimental import *
na = np.array
import pandas as pd
na = np.array
import seaborn as sns
import pandas as pd


## Load Experimental Data
origEps2 = [0.05, 0.1, 0.2, 0.25, 0.3, 0.5, 0.7, 0.8, 0.95]
IBIs = pd.read_excel('../ABC_bursting/data/meanFeatures.xlsx','meanIBI')
sdIBI = pd.read_excel('../ABC_bursting/data/meanFeatures.xlsx','SD')
CVs = pd.read_excel('../ABC_bursting/data/meanFeatures.xlsx','meanCV')
DatamIBI =np.nanmean(IBIs,0)*1000
DatamCV =np.nanmean(CVs,0)
#Load bic data
meanBic = []

Eps = [0.2,0.8]#[0.2,0.5,0.8]#[0.1,0.2,0.3,0.5,0.7,0.8,0.95]
mean_bic=[ ]
std_bic = []
for i,epsilon in enumerate(Eps):
    bic_ = na(pd.read_excel('data/bic_IBI_raw_corr.xlsx',str(epsilon)))[:7,1:]
    data_len = bic_.shape[1]
    bic_dat = bic_/bic_[0,:]
    mean_bic.append(np.nanmean(bic_dat,1))
    std_bic.append(np.nanstd(na(bic_dat,dtype = 'float'),1))
popMean_bic = np.nanmean(mean_bic,0)


# set BIC concentrations
c = [0,0.5,1,1.5,3.,10.,40.]
Kd =3
bic_conc =(1/(1+(na(c)/Kd)))
G_mod = np.round(bic_conc,2) # decrease in synaptic strength

# Get the model
# Load Posterioes
Post = []
Eps = [0.2,0.8]#,0.5,0.8]#[0.05,0.1,0.2,0.3,0.5,0.65,0.8,0.95]
for epsilon in Eps:
    name =  'ABC_results/nuki_fit_%s.npy'%epsilon
    posterior_= np.load(name, allow_pickle = True)
    posterior_ = [posterior_[na([p[4] for p in posterior_])>0]]
    Post.append(posterior_[0])

def get_model(eps,sim_type = 'ML', plot=True,
                index=False):
    """ 
    return the model with fiitted parameters
    Args:
            eps (float): Percentage of inhibitory neurons
            type (string: 'ML' or 'MAP'): if ML try to retrieve the best model from the accepted paramters, MAP - simulate a model with MAP parameters
            plot(bool): plot raster
            index(int): which simulations out of accepted to plot, if 
    Returns:
        ST (tuple of 2 lists): spike times,gids
        sc (array): spike counts
        signal (array): Ca-activity
        bursts (list): burst times (see MI_burst)
    """
    # TODO: read paramters from a a simulation file
    #PAramters of the simualtions
    params = load_obj('sim/nuki_fit/%s/params'%eps)
    # get the data
    sim_index= np.where(na(Eps)==eps)[0][0]
    print(sim_index)
    data_index= np.where(na(origEps2)==eps)[0][0]
    vs= True
    if sim_type=='MAP':
        print('estimating MAP...')
        best_params =find_MAP_2D(X[sim_index])#Best_params[i]#
        best_params = np.round(best_params,2)
        vs = False# need to run simualtions
    elif index:
        print('accepted simulation', index)
        best_params =Post[sim_index][-1][0][:,index]
    else:
        min_dist = np.argmin(Post[sim_index][-1][1])
        best_params =Post[sim_index][-1][0][:,min_dist]

    model =nuki_fit(params,orig_data=[DatamIBI[data_index],DatamCV[data_index]],
            alpha=1,beta =1,visual =False)
    #model.path = 'sim/nuki_fit_old/%s/'%eps
#    model.params['directory'] = 'sim/nuki_fit_old/%s/'%eps
    print('parameters:',best_params)
    params['epsilon']=eps
    # Generate data
    dat = model.generate_data(na([best_params[0],best_params[1]]))
    print('IBI,CV:',dat[0],dat[1])
    return model,best_params

eps = 0.2
model,best_theta = get_model(eps,sim_type='MAP')#index=3
dat = model.generate_data(na([best_theta[0],best_theta[1]]))

name = get_hash(model.params)
st,gid = read_gdf('sim/nuki_fit/%s/'%eps,name,(0,model.sim_time),threads = model.threads)

# simulate BIC for posterior samples
nest.set_verbosity('M_FATAL')
t0 = time.time()
params_list = []
all_data = np.zeros([2,len(G_mod)])
for g_i,g_mod in enumerate(G_mod):
    model.params['g'] = 4*g_mod
    model.params['directory'] = 'sim/nuki_fit/0.2/BIC/'#model.params['directory']+'BIC/'
    model.path ='sim/nuki_fit/%s/BIC/'%eps
    params_list.append(model.params.copy())
    dat = model.generate_data(na([best_theta[0],best_theta[1]]))
    all_data[:,g_i] = dat
t1 = time.time()
print('time elapsed=%s'%(t1-t0))
name = get_hash(model.params)
st,gid = read_gdf(model.params['directory'],name,(0,model.sim_time),threads = model.threads)

plt.plot(st,gid,'.',markersize = 0.5)
plt.show()
#np.save('ABC_results/BIC/5%_5x10samples_deg_set1_beta0.2',[all_data,params_list])


c[0]=0.0
bic_data = all_data[0,:]/all_data[0,0]
plt.plot(c,bic_data)
plt.errorbar(c,mean_bic[0],std_bic[0])
plt.xscale('symlog')
plt.show()

