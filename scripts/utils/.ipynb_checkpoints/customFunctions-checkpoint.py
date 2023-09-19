import math
import numpy as np
import scipy.io as io
import matplotlib.pyplot as plt
from sklearn import metrics
import seaborn as sns
import mat73
import sklearn.cluster as cluster
import time
sns.set_context('poster')
sns.set_color_codes()
plot_kwds = {'alpha' : 0.35, 's' : 1.5, 'linewidths':0}

def rotate(matrix):
    new_matrix = []
    for i in range(len(matrix[0]), 0, -1):
        new_matrix.append(list(map(lambda x: x[i-1], matrix)))

    return new_matrix
    
def cortexLab(anatList):
    catList = []
    nonCort = ['antTH', 'midTH', 'pstTH', 'BG', 'AMY', 'CLT']
    for i in anatList:
        if i not in nonCort:
            catList.append('COR')
        elif i in ['antTH', 'midTH', 'pstTH']:
            catList.append('THAL')
        else:
            catList.append('others')
    return catList

def thalLab(anatList):
    catList = []
    subCor = ['BG', 'AMY', 'CLT']
    thal = ['antTH', 'midTH', 'pstTH']
    for i in anatList:
        if i in thal:
            catList.append(i)
        elif i in subCor:
            catList.append('others')
        else:
            catList.append('Cor')
    return catList
 
def detHemi(x1, y1):
    '''
    x1, y1 are the first coordinate of the MNI coordinates of the out and in electrodes.
    '''
    if x1*y1<0:
        decision = 'contr' # contralateral hemisphere
    else:
        decision = 'ipsi'
    return decision

def mergeMat(sblist, keys=['filteridx_metaT','Vrpw','Vrpc'],
             inputDir = '.'):
    mat = {}
    for i,sb in enumerate(sblist):
        
        fn = "%s/SpecReduceCollapse_%s.mat" % (inputDir,sb)
        current = io.loadmat(fn)
        # account for the -1 indexing in python
        current['filteridx_metaT'] = current['filteridx_metaT']-1
        
        for key in keys:
            if i==0:
                mat[key] = current[key]
            else:
                mat[key]=np.concatenate((mat[key], current[key]),0)
    return mat
                
def makeLabel(variable,filterIdx, df):
   
    originalLabels = df[variable][filterIdx]
    labels = originalLabels.unique().tolist()
    target = []
    for i in originalLabels.tolist():
        target.append(labels.index(i))
    
    return target, labels, originalLabels
    
def anotEmbedding(embedding, center, threshold, shape='circle',angle=0):
    '''
    set boundary of each embedded space in the 2 dimensions
    '''
    diff = embedding - np.array(center)
    cl_anot = np.zeros((embedding.shape[0],) )
    
    if shape=='circle':
        dist = np.sqrt(np.sum(diff**2, axis=1))
        cl_anot[dist<=threshold]=1
    elif shape == 'rectangle':
        # threshold should have 2 dimensions: x, y
        x_thr = threshold[0]
        y_thr = threshold[1]
        cl_anot[(np.abs(diff[:,0])<=x_thr) & (np.abs(diff[:,1])<=y_thr)]=1
    elif shape == 'ellipse':
        angle = np.radians(angle)
        dist1 = ((np.cos(angle)*diff[:,0] + np.sin(angle)*diff[:,1])**2)/(threshold[0]**2)
        dist2 = ((np.sin(angle)*diff[:,0] - np.cos(angle)*diff[:,1])**2)/(threshold[1]**2)
        dist = dist1+dist2
        cl_anot[dist<=1]=1
        
    return cl_anot

def rgb_to_hex(rgb):
    r,g,b = np.round(np.dot(255,rgb))
    return '#{:02x}{:02x}{:02x}'.format(int(r), int(g),int(b))

def cm_mpl(cmap_name, n, start=0.1, end = 1, scale='equal'):
    cmap = plt.get_cmap(cmap_name)
    if scale=='equal':
        scale = np.arange(start, end+(end-start)/(n-1), (end-start)/(n-1))
    else:
        scale = (scale - np.min(scale)) / (np.max(scale) - np.min(scale))

    custom_palette = [cmap(x)[0:3] for x in scale]
    custom_palette = [rgb_to_hex(color) for color in custom_palette]
    return custom_palette

def plot_vecSpec(featureVec, n1,n2, ntime=60, 
                 title1 = 'power', title2='phase coherence',
                 vmin=-1, vmax=1):
    # (60,1)-reducer3,(65,-1)-reducer1, (45,-1)-reducer2
    pw = featureVec[0:n1].reshape(ntime, -1) 
    pc = featureVec[n1:(n1+n2)].reshape(ntime, -1)
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10,5))

    sp1 = ax1.imshow(rotate(pw), aspect=2)
    ax1.set_title(title1, fontsize=10)
    #ax1.set(xticks=np.arange(0, 65)[::20], xticklabels=np.arange(-25, 650)[::200]);
    ax1.tick_params( labelsize=10)
    
    sp2 = ax2.imshow(rotate(pc), aspect=2.1)
    ax2.set_title(title2, fontsize=10)
    #ax2.set(xticks=np.arange(0, 65)[::20], xticklabels=np.arange(-25, 650)[::200]);
    ax2.tick_params( labelsize=10)

def evalClst(clusterDB, X, labels_true):
    # Number of clusters in labels, ignoring noise if present.
    labels = clusterDB.labels_
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    print(' ')
    print('###################### metric report ##########################')
    print('Estimated number of clusters: %d' % n_clusters)
    print('Homogeneity: %0.3f' % metrics.homogeneity_score(labels_true, labels))
    print('Completeness: %0.3f' % metrics.completeness_score(labels_true, labels))
    print('V-measure: %0.3f' % metrics.v_measure_score(labels_true, labels))
    print('Adjusted Rand Index: %0.3f'
        % metrics.adjusted_rand_score(labels_true, labels))
    print('Adjusted Mutual Information: %0.3f'
        % metrics.adjusted_mutual_info_score(labels_true, labels))
    print('Silhouette Coefficient: %0.3f'
        % metrics.silhouette_score(X, labels))
    
def plot_clusters(data, algorithm, args, kwds):
   # start_time = time.time()
    labels = algorithm(*args, **kwds).fit_predict(data)
    #end_time = time.time()
    palette = sns.color_palette('deep', np.unique(labels).max() + 1)
    colors = [palette[x] if x >= 0 else (0.0, 0.0, 0.0) for x in labels]
    plt.scatter(data.T[0], data.T[1], c=colors, **plot_kwds)
    frame = plt.gca()
    frame.axes.get_xaxis().set_visible(False)
    frame.axes.get_yaxis().set_visible(False)
    plt.title('Clusters found by {}'.format(str(algorithm.__name__)), fontsize=12)
    #plt.text(-0.5, 0.7, 'Clustering took {:.2f} s'.format(end_time - start_time), fontsize=10)
    
def mergeSpectCCEP_test(sblist, 
                   cropdim=[[],[],[]], # filteridx_metaT is in the matlab counting scheme
                   keys=['CCEP_power', 'freqs', 'time', 'idx_in_metaT'], 
                   inputDir = '/data/dian/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked/spectCCEP',
                   data_type = "power" # or itpc
                  ):
    '''
    # merge subject level spectCCEP
    # input mat file's key ['idx_in_metaT'] is in the matlab counting scheme
    # crop dimension: [filteridx_py,[t1, t2), [f1, f2)]
    # input cropping parameter metaT index: filteridx must be in the python counting scheme, as a numpy array
    # input cropping parameter time t1, t2 are in the unit of seconds, and f1, f2 in Hz
    '''
    mat = {}
    for i,sb in enumerate(sblist):
        fn = "%s/spectCCEP_%s_%s.mat" % (inputDir, data_type ,sb)
        current = mat73.loadmat(fn)
        
        if (len(keys) == 0) & (i==0):
            keys = list(current.keys())
        
        # cropping parameters
        given_idx_metaT = cropdim[0] # in python counting scheme >=0
        t1 = cropdim[1][0]
        t2 = cropdim[1][1]
        f1 = cropdim[2][0]
        f2 = cropdim[2][1]
        # account for the -1 indexing in python
        current['idx_in_metaT'] = current['idx_in_metaT']-1
        boolean2select = np.isin(current['idx_in_metaT'], given_idx_metaT)
        if len(np.where(boolean2select)[0]) == 0:
            print('No data is found for %s' %(sb))
            continue
        else:
            # cropping area
            tme = current['time']
            frq = current['freqs']
            nt_2select = np.where((tme>= t1) & (tme<t2))[0]
            nf_2select = np.where((frq>= f1) & (frq<f2))[0]
            current['time']  = current['time'][nt_2select]
            current['freqs'] = current['freqs'][nf_2select]

            for key in keys:
                if current[key].shape[0] == boolean2select.shape[0]:
                    if current[key].ndim == 3: # dim = stimpair*time*freq
                        tmp = np.reshape(current[key][:,nt_2select][:,:,nf_2select],
                                                  (current[key].shape[0], len(nt_2select)*len(nf_2select))
                                                 )# crop + reshape
                        current[key] = tmp
                    if i == 0:
                        mat[key] = current[key][boolean2select]
                    else:
                        mat[key]=np.vstack((mat[key], current[key][boolean2select]))
                else:
                    if i == 0:
                        mat[key]=current[key]
                    else:
                        mat[key]=np.vstack((mat[key], current[key]))
    return mat

def stackSpectCCEP(filteridx,# numpy array
                   cropdim=[[],[]], # filteridx_metaT is in the matlab counting scheme
                   allfiles = [], # may need to be fixed when used in a parallel job
                   inputDir = '/data/dian/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked/spectCCEP/SingleTXT',
                   data_type = "power", #(pw) or itpc(pc) or phase(ph)
                   time_sample = np.arange(-0.5, 2.001, 0.005), # given sampled time points
                   freq_sample = np.array([  1.        ,   1.18920712,   1.41421356,   1.68179283,
         2.        ,   2.37841423,   2.82842712,   3.36358566,
         4.        ,   4.43827789,   4.92457765,   5.46416103,
         6.06286627,   6.72717132,   7.46426393,   8.28211939,
         9.18958684,  10.19648502,  11.3137085 ,  12.55334557,
        13.92880901,  15.45498126,  17.1483754 ,  19.02731384,
        21.11212657,  23.42537114,  25.99207668,  28.8400148 ,
        32.        ,  34.2967508 ,  36.75834736,  39.39662123,
        42.22425314,  45.254834  ,  48.50293013,  51.98415337,
        55.71523605,  59.71411146,  64.        ,  68.5935016 ,
        73.51669472,  78.79324245,  84.44850629,  90.50966799,
        97.00586026, 103.96830673, 111.4304721 , 119.42822292,
       128.        , 137.1870032 , 147.03338944, 157.58648491,
       168.89701258, 181.01933598, 194.01172051, 207.93661347,
       222.8609442 , 238.85644583, 256.        ])
                  ):
    '''
    # ----- stack trial level/stim-pair level spectCCEP -----
    # filteridx must be in the python counting scheme, as a numpy array - refer back to the metaT
    # crop dimension: [[t1, t2), [f1, f2)]
    # input cropping parameter time t1, t2 are in the unit of seconds, and f1, f2 in Hz
    '''
    mat = {} # output file
    
    if isinstance(filteridx, list):
        filteridx = np.array(filteridx)
        
    if data_type =="power":
        data_type = 'pw'
    elif data_type == "phase":
        data_type = 'pw'
    elif data_type == 'itpc':
        data_type = 'pc'
    
    if len(allfiles) == 0:
        allfiles = fnmatch.filter(os.listdir(inputDir), (data_type+'*.txt'))
        IDX = [int(i.split('IDX')[1].split('.txt')[0]) for i in allfiles]
     
    inputfiles = allfiles[np.isin(filteridx, IDX)]
    filteridx_out = filteridx[np.isin(filteridx, IDX)]
    mat['filteridx'] = filteridx_out
    
    if filteridx_out.shape[0] < filteridx.shape[0]:
        print('Given files do not have all the requested stim-pair instances.')
    
     # cropping area from the dimension of 59 (freq)*501 (time)
    t1 = cropdim[0][0]
    t2 = cropdim[0][1]
    f1 = cropdim[1][0]
    f2 = cropdim[1][1]
    nt_2select = ((time_sample>= t1) & (time_sample<t2))
    nf_2select = ((freq_sample>= f1) & (freq_sample<f2))
    # put in the output dict
    mat['time']  = time_sample[nt_2select]
    mat['freqs'] = freq_sample[nf_2select]

    dat = np.zeros((filteridx_out.shape[0],
                    mat['freqs'].shape[0], 
                    mat['time'].shape[0]))
    
    for i,fn in enumerate(inputfiles):
        
        d = np.reshape(np.loadtxt(inputDir+'/'+fn), 
                       (mat['freqs'].shape[0],  mat['time'].shape[0]))
        # crop d
        dat[i] = d[nf_2select,:][:,nt_2select]
        
    mat[data_type] = dat
    
    return mat
