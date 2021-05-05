#using python to do preprocessing
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
     

idlist=[80605, 80607, 80609, 80620, 80622, 80669, 80673, 80674, 80675, 80676, 80677,  80678, 80679, 80680, 80681, 80682, 
80683, 80684, 80685, 80686, 80688, 80690, 80692, 80693, 80694, 80699, 80700, 80707, 80711, 80712]
#id_with3=[80605, 80607, 80609, 80620, 80622,80712]
sightlines=[]
for i in idlist:
    loc=path+'/'+str(i)+'/deep'
    for j in [0,1,2,4,5,6,7,8,9]:
        spectra=join(loc,'coadd-%s-%s-deep.fits'%(j,i))
        zbest=join(loc,'zbest-%s-%s-deep.fits'%(j,i))
        specs = DesiMock()
        specs.read_fits_file(spectra,[],zbest)
        keys = list(specs.data.keys())
        for kk in keys:
            if (specs.data[kk]['spectype'] == 'QSO')&(specs.data[kk]['z_qso'] <5.81)&(specs.data[kk]['z_qso'] >2.33)&(len(np.nonzero(specs.data[kk]['FLUX'])[0])!=0):
                sightline = specs.get_sightline(kk,camera = 'all', rebin=False, normalize=False)
                sightline.s2n=estimate_s2n(sightline)
                if np.isnan(sightline.s2n):
                    print(i,j,kk)
                sightlines.append(sightline)
for sightline in sightlines:
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
    else:
        print(sightline.id)
   
   
dataNew = 'zoujq/cascades/sightlines/preload_qsos.mat'
scio.savemat(dataNew, {'loading_min_lambda': loading_min_lambda,'loading_max_lambda': loading_max_lambda,'normalization_min_lambda':normalization_min_lambda,'normalization_max_lambda':normalization_max_lambda,'min_num_pixels':min_num_pixels,'sightline_ids':sightline_ids, 'all_wavelengths':all_wavelengths,'all_flux':all_flux,'all_noise_variance':all_noise_variance,'all_normalizers':all_normalizers})

