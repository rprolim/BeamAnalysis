def HIDEmaps2multichannel(dict_hidemap=None):
    import os
    import numpy as np
    import astropy.io.fits as fits
    from copy import deepcopy as dcopy
    dirpath  = dict_hidemap['dirpath']
    files    = np.sort(os.listdir(dirpath))
    filenames_ = np.array([])
    indch      = np.array([])
    for i,ifile in enumerate(files): #limited to the mapname files
        if len(ifile.split(dict_hidemap['mapname']))==2:
            ich = int(ifile.split('_ch')[1].split('_')[0]) #takes the channel number
            indch      = np.hstack((indch,int(ich)))
            filenames_ = np.hstack((filenames_,ifile))
    indch       = np.int8(indch)   
    if dict_hidemap['stack_channel_sorting']:
        indch=np.sort(indch)
    else:
        indch=np.flip(np.sort(indch))
    nfilenames_ = np.array([])
    for i, ich in enumerate(np.arange(dict_hidemap['first_channel'], dict_hidemap['first_channel'] + dict_hidemap['nch'])):
        if ich in indch:
            ifilename   = filenames_[np.where(indch==ich)][0]
            nfilenames_ = np.hstack((nfilenames_,ifilename))
            ipath       = os.path.join(dirpath, ifilename)
            with fits.open(ipath) as h:
                    imap = dcopy(h[int(dict_hidemap['hdu_map'])].data)
            try:
                jmap = np.vstack((jmap,imap)) 
            except:
                jmap = imap
        else:
            print("There is no file corresponding in {} channel".format(ich))
    return jmap
        
        
def load_multFRemoval(paramsCS=None, paramsWT=None, paramsMaps=None, paramsPath=None,
                      map_input=None, ns_vec=[1,2,3,4,5],
                      mask_dir='/media/BINGODATA1/ComponentSeparation/building_dataset/dataset/M256',
                      mask_name=None, mask_apply=True):
    import tnumpy as np
    ns_vec = np.atleast_1d(ns_vec)
    M = 1 
    if mask_apply: 
        M = hdata.getmap(dirpath_ = mask_dir,filename_= mask_name, healpix_readingformat=False, hdu=1)      
    for i,ni in enumerate(ns_vec):
        params_CS['ns']=int(ni)
        mmap = dcopy(map_input*M)
        params_cs, params_wt = cs.load(paramsCS,paramsWT)
        X = cs.adaptation_maps(mmap, paramsMaps, paramsPath)
        X = cs.maps2CSmaps(    X, params_wt  , params_cs)
        if not i:
            dict = {'ns1':X}
        else:
            dict['ns{:d}'.format(i+1)]=X
        del X
    return dict  
