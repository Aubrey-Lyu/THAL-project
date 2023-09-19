import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import os
import mne
from mne.time_frequency import tfr_morlet
import scipy.io as io

def tscore_1sample(X):
    '''
    assume X is a 2D matrix (sample size * n-datapoints of 1 sample)
    '''
    n = X.shape[0]
    t = np.mean(X,axis=0)/(np.std(X,axis=0)/np.sqrt(n))
    return t

def tscore_2sample(x1,x2):
    n1 = x1.shape[0]
    n2 = x2.shape[0]
    s1 = np.std(x1,axis=0)
    s2 = np.std(x2,axis=0)
    sp = np.sqrt(((n1-1)*np.square(s1)+(n2-1)*np.square(s2))/(n1+n2-2)) # pooled standard deviation
    t = (np.mean(x1,axis=0)-np.mean(x2,axis=0))/(sp*np.sqrt(1/n1 + 1/n2))
    return t

def contrast_2sample(x1,x2):
    return np.mean(x1,axis=0)-np.mean(x2,axis=0)
    

def zscore_1sample(X):
    '''
    assume X is a 2D matrix (sample size * n-datapoints of 1 sample)
    '''
    n = X.shape[0]
    z = np.mean(X,axis=0)/(np.std(X,axis=0))
    return z

def MNE1SampleSpectraSigTest(d1, threshold='auto', outname= 'tmp.mat', ifplot=1,  n_permutations = 5000, ifsave=0):
    # convert vector to time-freq matrix
    n1=1380#mat['Vrpw'].shape[1] # 1380
    n2=1200#mat['Vrpc'].shape[1] # 1200
    ntime = 60
    #
    pw1 = d1[:, 0:n1].reshape(d1.shape[0],ntime, -1) 
    pc1 = d1[:, n1:(n1+n2)].reshape(d1.shape[0],ntime, -1)
  
   # COMPUTE ADJACENCY MATRIX
    adjacency = mne.stats.combine_adjacency(pw1.shape[1], pw1.shape[2]) #  times, freqs

    # We want a one-tailed test
    tail = 1
    degrees_of_freedom = d1.shape[0] - 1
    if threshold=='auto':
        t_thresh = scipy.stats.t.ppf(1 - 0.001/2, df=degrees_of_freedom)
    else:
        t_thresh = threshold
    print('Initial cluster threshold is %s'%(t_thresh))
    print('Degrees of freedom is %d' % (degrees_of_freedom))


    # Run the analysis
    T_obs, clusters, cluster_p_values, H0 = mne.stats.permutation_cluster_1samp_test(
        pw1, n_jobs=-1,
        n_permutations=n_permutations,
        threshold=t_thresh,
        tail=tail,
        adjacency=adjacency,
        out_type="mask",
        verbose=True
    )
    T_obs_plot = np.nan * np.ones_like(T_obs)
    
    
    for c, p_val in zip(clusters, cluster_p_values):
        if p_val <= 0.01:
            T_obs_plot[c] = T_obs[c]
            
    t_idx, f_idx = np.unravel_index(
            np.nanargmax(np.abs(T_obs_plot)), pw1.shape[1:]
        )

    if ifplot==1:
        plt.close()
        plt.figure()
        plt.subplots_adjust(0.12, 0.08, 0.96, 0.94, 0.2, 0.43)

        vmax = np.max(np.abs(T_obs))
        vmin = -vmax
        plt.close()
        plt.figure()
        plt.imshow(
            np.transpose(T_obs),
            cmap=plt.cm.gray,
           # extent=[times[0], times[-1], freqs[0], freqs[-1]],
            aspect="auto",
            origin="lower",
            vmin=vmin,
            vmax=vmax,
        )

        plt.imshow(
            np.transpose(T_obs_plot),
            cmap=plt.cm.RdBu_r,
         #   extent=[times[0], times[-1], freqs[0], freqs[-1]],
            aspect="auto",
            origin="lower",
            vmin=vmin,
            vmax=vmax,
        )
        plt.colorbar()
        plt.xlabel("Time (custome unit)")
        plt.ylabel("Frequency (custome unit)")
        #plt.title(f"Induced power ({tfr_epochs.ch_names[ch_idx]})")
        plt.show()
    
    sigClsts = {"T_obs":T_obs, "clusters":np.array(clusters), "cluster_p_values":cluster_p_values, "H0":H0,
                   "T_obs_thr":T_obs_plot, 't_idx':t_idx, 'f_idx':f_idx}
        
    if ifsave==1:
        # save significant clusters
          io.savemat("GroupLevel_SigTest/"+outname, sigClsts)
        
    return sigClsts

def MNE2SampleSpectraSigTest(d1, d2, threshold, outname= 'tmp.mat', ifplot=1,  n_permutations = 5000, ifsave=0):
    # convert vector to time-freq matrix
    n1=mat['Vrpw'].shape[1] # 1380
    n2=mat['Vrpc'].shape[1] # 1200
    ntime = 60
    #
    pw1 = d1[:, 0:n1].reshape(d1.shape[0],ntime, -1) 
    pc1 = d1[:, n1:(n1+n2)].reshape(d1.shape[0],ntime, -1)
    #
    pw2 = d2[:, 0:n1].reshape(d2.shape[0],ntime, -1) 
    pc2 = d2[:, n1:(n1+n2)].reshape(d2.shape[0],ntime, -1)
    
    # plot
    # plot to check data
    if ifplot==1:
        plotSpectra2Groups(np.rot90(np.nanmean(pw1, axis=0)),
                           np.rot90(np.nanmean(pw2, axis=0)))
        
   # COMPUTE ADJACENCY MATRIX
    adjacency = mne.stats.combine_adjacency(pw1.shape[1], pw1.shape[2]) #  times, freqs

    # We want a two-tailed test
    tail = 0

    # Set the number of permutations to run.
    # Warning: 50 is way too small for a real-world analysis (where values of 5000
    # or higher are used), but here we use it to increase the computation speed.
   

    # Run the analysis
    F_obs, clusters, cluster_p_values, H0 = mne.stats.permutation_cluster_test(
        [pw1,pw2], n_jobs=-1,
        n_permutations=n_permutations,
        threshold=threshold,
        tail=tail,
        adjacency=adjacency,
        out_type="mask",
        verbose=True
    )
    F_obs_plot = np.nan * np.ones_like(F_obs)
    
    
    for c, p_val in zip(clusters, cluster_p_values):
        if p_val <= 0.01:
            F_obs_plot[c] = F_obs[c]
            
    t_idx, f_idx = np.unravel_index(
            np.nanargmax(np.abs(F_obs_plot)), pw1.shape[1:]
        )

    if ifplot==1:
        plt.figure()
        plt.subplots_adjust(0.12, 0.08, 0.96, 0.94, 0.2, 0.43)

        vmax = np.max(np.abs(F_obs))
        vmin = -vmax
        plt.close()
        plt.figure()
        plt.imshow(
            np.transpose(F_obs),
            cmap=plt.cm.gray,
           # extent=[times[0], times[-1], freqs[0], freqs[-1]],
            aspect="auto",
            origin="lower",
            vmin=vmin,
            vmax=vmax,
        )

        plt.imshow(
            np.transpose(F_obs_plot),
            cmap=plt.cm.RdBu_r,
         #   extent=[times[0], times[-1], freqs[0], freqs[-1]],
            aspect="auto",
            origin="lower",
            vmin=vmin,
            vmax=vmax,
        )
        plt.colorbar()
        plt.xlabel("Time (custome unit)")
        plt.ylabel("Frequency (custome unit)")
        #plt.title(f"Induced power ({tfr_epochs.ch_names[ch_idx]})")
        plt.show()
    
    sigClsts = {"F_obs":F_obs, "clusters":np.array(clusters), "cluster_p_values":cluster_p_values, "H0":H0,
                   "F_obs_thr":F_obs_plot, 't_idx':t_idx, 'f_idx':f_idx}
        
    if ifsave==1:
        # save significant clusters
        io.savemat("GroupLevel_SigTest/"+outname, sigClsts)
        
    return sigClsts
    
def plotSpectra2Groups(dshow1,dshow2):
    plt.close()
    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(15,5))

    sp1 = ax1.imshow(dshow1, aspect=2.1)
    ax1.set_title('mean power (group1)', fontsize=10)
    #ax1.set(xticks=np.arange(0, 65)[::20], xticklabels=np.arange(-25, 650)[::200]);
    ax1.tick_params( labelsize=10)
    plt.colorbar(sp1,shrink=0.6)

    sp2 = ax2.imshow(dshow2, aspect=2.1)
    ax2.set_title('mean power (group2)', fontsize=10)
    #ax2.set(xticks=np.arange(0, 65)[::20], xticklabels=np.arange(-25, 650)[::200]);
    ax2.tick_params( labelsize=10)
    plt.colorbar(sp2,shrink=0.6)

    sp3 = ax3.imshow(dshow1-dshow2,
                     aspect=2.1, cmap='RdBu_r')
    ax3.set_title('mean power differences', fontsize=10)
    #ax2.set(xticks=np.arange(0, 65)[::20], xticklabels=np.arange(-25, 650)[::200]);
    ax3.tick_params(labelsize=10)
    plt.colorbar(sp3,shrink=0.6)
    plt.show()
    