#using python to build catalogs
import numpy as np
import scipy.io as scio
from astropy.table import Table, vstack

idlist=[80605, 80607, 80609, 80620, 80622, 80669, 80673, 80674, 80675, 80676, 80677,  80678, 80679, 80680, 80681, 80682, 
80683, 80684, 80685, 80686, 80688, 80690, 80692, 80693, 80694, 80699, 80700, 80707, 80711, 80712]
sightlines=[]
zwarning=[]
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
                zwarning.append(specs.data[kk]['ZWARN'])
                sightlines.append(sightline)
#set bal flags=0
bal_visual_flags = np.zeros(2766)
filter_flags = np.zeros(2766)
for i in np.nonzero(zwarning)[0]:
    filter_flags[i]=5
qso_tbl = Table(names=('Plate','FiberID','MJD','TARGET_RA','TARGET_DEC', 'ZQSO','TARGETID','S/N','bal_visual_flags','filter_flags','zwarning'),dtype=('int','int','int','float','float','float','int','float','int','int','int'),meta={'EXTNAME': 'QSOCAT'})
for i in range(0,len(sightlines)):
    sightline = sightlines[i]
    qso_tbl.add_row((sightline.id,sightline.id,sightline.id,sightline.ra,sightline.dec,sightline.z_qso,sightline.id,sightline.s2n,0,filter_flags[i],zwarning[i]))
data = 'zoujq/cascades/sightlines/catalog.mat'
scio.savemat(data, {'ras':list(qso_tbl['TARGET_RA']),'decs':list(qso_tbl['TARGET_DEC']),'target_ids':list(qso_tbl['TARGETID']),'z_qsos':list(qso_tbl['ZQSO']),'snrs':list(qso_tbl['S/N']),'bal_visual_flags':list(qso_tbl['bal_visual_flags']),'filter_flags':list(qso_tbl['filter_flags']),'zwarning':list(qso_tbl['zwarning'])})
