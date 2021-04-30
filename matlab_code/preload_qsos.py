#using python to do preprocessing
#The order of this output is not the same as in catalog.mat and needs to be matched using preload_qsos.mat
import numpy as np
from dla_cnn.desi.DesiMock import DesiMock
import scipy.io as scio
import os
all_wavelengths    =  []
all_flux           =  []
all_noise_variance =  []
all_pixel_mask     =  []
all_normalizers    = []
sightline_ids=[]
normalization_min_lambda = 1310
normalization_max_lambda = 1325
min_lambda         =  911.75
max_lambda         = 1215.75
min_num_pixels = 200
loading_min_lambda = 910
loading_max_lambda = 1217
path='desiY1-0.2-DLA/spectra-16'
item = os.listdir(path)
     
    
def make_mat(k):
    item1 = os.listdir(path+'/'+str(k))
    for j in item1:
        spectra=path+'/'+str(k)+'/'+str(j)+'/spectra-16-%s.fits'%(j)
        truth=path+'/'+str(k)+'/'+str(j)+'/truth-16-%s.fits'%(j)
        zbest=path+'/'+str(k)+'/'+str(j)+'/zbest-16-%s.fits'%(j)
        specs = DesiMock()
        specs.read_fits_file(spectra,truth,zbest)
        keys = list(specs.data.keys())
        for jj in keys:
            sightline = specs.get_sightline(jj,camera = 'all', rebin=False, normalize=False)
            this_wavelength=10**sightline.loglam
            rest_wavelength=this_wavelength/(1+sightline.z_qso)
            ind = (rest_wavelength >= normalization_min_lambda) & (rest_wavelength <= normalization_max_lambda)
            norm=np.nanmedian(sightline.flux[ind])
            ind = (rest_wavelength >= min_lambda) & (rest_wavelength<= max_lambda)
            if len(np.nonzero(ind)[0])>=min_num_pixels:
                sightline_ids.append(sightline.id)
                flux=sightline.flux/norm
                var=sightline.error**2/norm**2
                ind = (rest_wavelength >= loading_min_lambda) & (rest_wavelength <= loading_max_lambda)
                ind[max(0,np.nonzero(ind)[0][0]-1)]=True
                ind[min(np.nonzero(ind)[0][-1]+1,len(np.nonzero(ind)[0]-1))]=True
                all_wavelengths.append(list(this_wavelength[ind]))
                all_flux.append(list(flux[ind]))
                all_noise_variance.append(list(var[ind]))
                all_normalizers.append(norm)
        print(j)
    dataNew = 'desiY1-0.2-DLA/%s-sightlines.mat'%(k)
    scio.savemat(dataNew, {'sightline_ids':sightline_ids, 'wavelengths':all_wavelengths,'flux':all_flux,'noise_variance':all_noise_variance,'normalizers':all_normalizers})

#import multiprocessing
#for k in [14,16,4]:
    #make_datasets(k)
    #p=multiprocessing.Process(target=make_mat,args=(int(k),))
    #p.start()
    #print(k)
for k in ['20', '18', '12', '16', '5', '30', '23', '11', '2', '13', '19', '0', '14', '7', '22', '3', '10', '9', '4', '17']:
    datafile2='desiY1-0.2-DLA/%s-sightlines.mat'%(k)
    data2 = scio.loadmat(datafile2)
    sightline_ids=np.hstack((sightline_ids,data2['sightline_ids']))
    wavelengths=np.hstack((wavelengths,data2['wavelengths']))
    flux=np.hstack((flux,data2['flux']))
    noise_variance=np.hstack((noise_variance,data2['noise_variance']))
    normalizers=np.hstack((normalizers,data2['normalizers']))

dataall = 'desiY1-0.2-DLA/sightlines.mat'
scio.savemat(dataall, {'sightline_ids':list(sightline_ids), 'wavelengths':wavelengths,'flux':flux,'noise_variance':noise_variance,'normalizers':normalizers})
