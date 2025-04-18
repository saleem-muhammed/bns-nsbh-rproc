import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def wtquantile(xs,qs,wts=[]):
    
    nan_idxs = np.isnan(xs)
    xs = np.array(xs[~nan_idxs])
    
    num_xs = len(xs)
    qs = np.array(qs, ndmin=1)
    if len(xs) < 1: return np.full(2*len(qs),np.nan)
    elif len(xs) == 1: return np.full(2*len(qs),xs[0])
    if len(wts) < 1: wts = np.full(num_xs, 1.)
    else: wts = np.array(wts[~nan_idxs])
    
    ps = wts/np.sum(wts)
    xs_sorted,ps_sorted = zip(*sorted(list(zip(xs,ps)),reverse=False))

    Ps = np.cumsum(ps_sorted)

    idxs_lb = np.array([np.where(Ps >= (1.-q)/2.)[0][0] for q in qs])
    idxs_ub = np.array([np.where(Ps >= 1.-(1.-q)/2.)[0][0] for q in qs])
    xs_sorted = np.array(xs_sorted)

    return list(xs_sorted[idxs_lb])+list(xs_sorted[idxs_ub])


def plot_abundance_history_single (data, show_error_bar):
    
    if show_error_bar:
        plt.fill_between(data['FeH_grid'],data['qs'][:,1],data['qs'][:,3],
                     facecolor=data['color'],edgecolor=None,alpha=0.25, zorder=10) # 90% CI
    
        plt.fill_between(data['FeH_grid'],data['qs'][:,0],data['qs'][:,2],
                     facecolor=data['color'],edgecolor=None,alpha=0.5,zorder=10) # 68% CI
    
    plt.plot(data['FeH_grid'],np.array(data['md'])[:,0],c=data['color'], label=data['label'], zorder=10) # median

def plot_abundance_history (list_of_data_dicts, show_error_bar = True, show_observed = True, save=False, save_label=None):
    plt.figure(figsize=(6.4,4.8))

    if show_observed:

        # load disk star observations and make gaussian likelihood model for each datapoint
        OBSPATH = '../etc/Battistini16_disk.csv' # path to disk star observations
        FeHs, EuFes, FeH_errs, EuFe_errs = np.loadtxt(OBSPATH, unpack=True, delimiter=',', skiprows=1)

        # load disk+halo star observations for plotting
        OBSPATH2 = '../etc/SAGA_MP.csv' # path to disk+halo star observations
        FeHs2, EuFes2, FeH_errs2, EuFe_errs2 = np.loadtxt(OBSPATH2, unpack=True, delimiter=',', skiprows=1)
        
        #plot them on the background
        plt.errorbar(FeHs2, EuFes2, xerr=[FeH_errs2,FeH_errs2], yerr=[EuFe_errs2,EuFe_errs2], c='g', fmt=',', lw=1, label='SAGA')
        plt.scatter(FeHs, EuFes,marker='D',facecolor='dodgerblue',edgecolor='navy', s=16, lw=0.5, label='Battistini & Bensby (2016)')
    
    # plot the abundance history
    for data_dict in list_of_data_dicts:
        plot_abundance_history_single (data_dict, show_error_bar)
     
    plt.xlim(-3.,0.5)
    plt.ylim(-1.,1.5)
    plt.xlabel('[Fe/H]')
    plt.ylabel('[Eu/Fe]')
    plt.legend(frameon=True,loc='upper right')
    if save:
        plt.savefig(save_label)
    plt.show()

    


