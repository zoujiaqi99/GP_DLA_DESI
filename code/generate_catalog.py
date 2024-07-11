import scipy.io as scio
import numpy as np
import h5py
import json
from astropy.table import Table,vstack
from tqdm import tqdm
import os
import scipy.io as scio
#this module can directly derive the dla catalog for gp model
#quite simple to extract dla catalog
def pred_gp_sightline(preloaded_file,catalogue_file, 
                      processed_file,  
                      dlacatalog_file,
                      sub_dla=1):
   
    #read files
    preloaded_file=scio.loadmat(preloaded_file)
    catalogue_file=scio.loadmat(catalogue_file)
    #learned_file   = h5py.File(learned_file,   'r')
    processed_file = h5py.File(processed_file, 'r')
    #test ind
    test_ind = processed_file['test_ind'][:,0].astype(np.bool) #size: (num_qsos, )sv:[:, 0]
    test_real_index = np.nonzero(test_ind )[0] 
    model_posteriors = processed_file['model_posteriors'][()].T
    p_dlas           = processed_file['p_dlas'][0, :]#at least one dla prob
    p_no_dlas        = processed_file['p_no_dlas'][0, :]#no dla, no sub-dla prob
    if sub_dla:
        p_no_dlas += model_posteriors[:, 1]#add sub-dla prob to p_no_dlas
        log_priors_dla   = processed_file['log_priors_dla'][0, :]
        min_z_dlas       = processed_file['min_z_dlas'][0, :]#search range
        max_z_dlas       = processed_file['max_z_dlas'][0, :]#search range
    if 'MAP_log_nhis' in processed_file.keys():
        map_log_nhis = processed_file['MAP_log_nhis'][()].T#(n_qso,4,4),contains the DLA(1,2,3,4) NHI value with highest prob
        map_z_dlas   = processed_file['MAP_z_dlas'][()].T#(n_qso,4,4),contains the DLA(1,2,3,4) Z value with highest prob
    target_ids = catalogue_file['target_ids'][0].astype(np.int)[test_ind]
    z_qsos     = catalogue_file['z_qsos'][0,:][test_ind]
    snrs   = catalogue_file['snrs'][0, :][test_ind]
    try:
        ras    = catalogue_file['ras'][0, :][test_ind]
        decs   = catalogue_file['decs'][0, :][test_ind] 
    except:
        ras    = np.zeros(len(z_qsos)) 
        decs    = np.zeros(len(z_qsos))    
    predictions_DLAs = []

    # query the maximum values of the model posteriors per spectrum first
    model_index = model_posteriors.argmax(axis=1)#find the best model with highest prob
    max_model_posteriors = model_posteriors.max(axis=1)

    if sub_dla:
            # now we need to combine sub-DLAs with null model posterior
            # to get the mosterior of no-DLAs
        inds = model_index < (1 + sub_dla)
        max_model_posteriors[inds] = p_no_dlas[inds]

        # prepare to use model_index as num_dlas
        num_dlas = model_index - sub_dla #get the best model's dla number
        num_dlas[num_dlas < 0] = 0 # num_dlas should be positive

    assert len(max_model_posteriors) == len(target_ids)

    for i,target_id in enumerate(target_ids):
        spec = dict()
        spec["p_dla"]               = p_dlas[i].item()
        spec["p_no_dla"]            = p_no_dlas[i].item()
        spec["max_model_posterior"] = max_model_posteriors[i].item()
        spec["num_dlas"]            = num_dlas[i].item()
        spec["min_z_dla"]           = min_z_dlas[i].item()
        spec["max_z_dla"]           = max_z_dlas[i].item()
        spec["snr"]                 = snrs[i].item()
        spec["ra"]                  = ras[i].item()
        spec["dec"]                 = decs[i].item()
        spec["target_id"]            = target_id.item()
        spec["z_qso"]               = z_qsos[i].item()

        dlas = []

        if num_dlas[i] > 0:
            this_map_z_dlas   = map_z_dlas[i, (num_dlas[i] - 1), :num_dlas[i]]#extract the best Z value
            this_map_log_nhis = map_log_nhis[i, (num_dlas[i] - 1), :num_dlas[i]]

            assert len(this_map_z_dlas) == num_dlas[i]

            for j in range(num_dlas[i]):
                dlas.append([
                    this_map_log_nhis[j],
                    this_map_z_dlas[j]])#change this to predict lyb
            #这里这么写下一步找lyb很费力（不知道怎么同时提取所有dla的信息）所以下次做的时候在这里改一下把nhi z分成两个dict
            
        spec["dlas"] = dlas

        predictions_DLAs.append(spec)
    dla_tbl = Table(names=('TARGET_RA','TARGET_DEC', 'Z_QSO','Z_DLA','TARGETID','S2N','DLAID','NHI','DLA_CONFIDENCE','NHI_STD','ABSORBER_TYPE'),dtype=('float','float','float','float','int','float','str','float','float','float','str'),meta={'EXTNAME': 'DLACAT'})
    for ii in range(0,len(predictions_DLAs)):
        for jj in range(0,len(predictions_DLAs[ii]['dlas'])):
            dla=predictions_DLAs[ii]['dlas'][jj]
            absorber_type =  "DLA" if dla[0] >= 20.3 else "SUBDLA"
            #add lyb prediction:
            lambda_higher = (dla[1]+1) *1215.67/ (1025.722/1215.67)
            peak_difference_spectrum = np.abs((np.array(predictions_DLAs[ii]['dlas'])[:,1]+1)*1215.67 - lambda_higher)
            nearest_peak_ix = np.argmin(peak_difference_spectrum)
            if (peak_difference_spectrum[nearest_peak_ix]<=15)&(dla[0]<predictions_DLAs[ii]['dlas'][nearest_peak_ix][0]-0.3):
                absorber_type='LYB'
           
            dla_tbl.add_row((predictions_DLAs[ii]['ra'],predictions_DLAs[ii]['dec'],predictions_DLAs[ii]['z_qso'],dla[1],predictions_DLAs[ii]['target_id'],predictions_DLAs[ii]['snr'],str(predictions_DLAs[ii]['target_id'])+'00'+str(jj),float(dla[0]),predictions_DLAs[ii]['p_dla'],0.0,absorber_type))
            
    dla_tbl.write(dlacatalog_file,overwrite=True)
    #dla_subtbl=dla_tbl[(dla_tbl['ABSORBER_TYPE']!='LYB')&(1216*(1+dla_tbl['Z_DLA'])>912*(1+dla_tbl['Z_QSO']))&(1216*(1+dla_tbl['Z_DLA'])<1250*(1+dla_tbl['Z_QSO']))&(dla_tbl['DLA_CONFIDENCE']>0.9)]
    #dla_subtbl.write(dlacatalog_sub_file,overwrite=True)

    
    #with open(outfile, "w") as json_file:
        #json.dump(predictions_DLAs, json_file, indent=2)

    #return predictions_DLAs

def rename_cats(dlacatalog_file):
    dlas=Table.read(dlacatalog_file)
    dlas.rename_column('Z','Z_DLA')
    dlas.rename_column('ZQSO','Z_QSO')
    dlas.rename_column('S/N','S2N')
    dlas.write(dlacatalog_file,overwrite=True)

def combine_cat(server='/home/zjqi/data/desi',output='main-dark-dlacat-gp-final.fits',sub_output='main-dark-dlacat-gp-final-sub.fits'):
    wrong_healpix=[]
    dlacat=Table()
    drs=os.listdir(server+'/%s/sightlines_%s_%s'%(release,survey,program))
    for d in tqdm(drs): 
        ds=os.listdir(server+'/%s/sightlines_%s_%s/%s'%(release,survey,program,d))
        if len(ds)>1:
            
            dlacat_loc=server+'/%s/new_dlacat_gp/%s-dlacat_gp.fits'%(release,d)
            try:#combine cats
                dlacat_pix=Table.read(dlacat_loc)
                    #if len(dlacat_pix)>0:
                        #dlacat_pix['release']=release
                        #dlacat_pix['survey']=survey
                        #dlacat_pix['program']=program
                        #dlacat_pix['healpix']=d
                        #dlacat=vstack((dlacat,dlacat_pix))
                dlacat_pix['healpix']=d
                dlacat=vstack((dlacat,dlacat_pix))
            except:
                print('Wrong running: %s %s %s %s'%(release,survey,program,d))
                wrong_healpix.append(d)
    dlacat.write('/home/zjqi/data/desi/iron/'+str(output),overwrite=True)
    #dla_subtbl=dlacat[(dlacat['ABSORBER_TYPE']!='LYB')&(1216*(1+dlacat['Z_DLA'])>912*(1+dlacat['Z_QSO']))&(1216*(1+dlacat['Z_DLA'])<1250*(1+dlacat['Z_QSO']))&(dlacat['DLA_CONFIDENCE']>0.9)&(dlacat['NHI']>=19)&(dlacat['NHI']<=23)]
    #dla_subtbl.write('/home/zjqi/data/desi/iron/'+str(sub_output),overwrite=True)
    print('%s do not have dlacat'%len(wrong_healpix))
    return dlacat

#*****For Iron release: make individual dla cat and then check error healpixs, finally combine catalogs
release='iron'
program='dark'
survey='sv3'
server='/home/zjqi/data/desi'
drs=os.listdir(server+'/%s/sightlines_%s_%s'%(release,survey,program))
lose_healpix=[]
for d in tqdm(drs):
    if len(os.listdir(server+'/%s/sightlines_%s_%s/%s'%(release,survey,program,d)))>1:
        dr=str(d)[0:-2]
        if dr =='':
            dr=0
        try:
            preloaded_file=server+'/%s/sightlines_%s_%s/%s/%s-%s-%s-%s-%s-preload_qsos.mat'%(release,survey,program,d,release,survey,program,dr,d)
            catalogue_file=server+'/%s/sightlines_%s_%s/%s/%s-%s-%s-%s-%s-catalog.mat'%(release,survey,program,d,release,survey,program,dr,d)
            processed_file=server+'/%s/%s_newmat/%s-process_qsos.mat'%(release,survey,d)
            dlacatalog_file=server+'/%s/%s_dlacat_gp/%s-dlacat_gp.fits'%(release,survey,d)
            pred_gp_sightline(preloaded_file,catalogue_file, processed_file, dlacatalog_file)
        except:
            qsos=Table.read('/home/zjqi/data/desi/%s/sightlines_%s_%s/%s/%s-%s-%s-%s-%s-catalog.fits'%(release,survey,program,d,release,survey,program,dr,d))
            if len(qsos)>1:
                print(d,len(qsos))
            lose_healpix.append(d)
            
print('%s healpixs do not have cats'%len(lose_healpix))
np.save(server+'/%s/new_lose_heapix_gpcat.npy'%release,lose_healpix)
scio.savemat(server+'/%s/new_lose_heapix_gpcat.mat'%release, {'los_pix':lose_healpix})
'''
combine_cat()

#****check except error: no qsos or only 1 qso
lose_healpix=np.load('/home/zjqi/data/desi/iron/new_lose_heapix_gpcat.npy',allow_pickle=True)
mat=os.listdir(server+'/%s/newmat'%release)
for d in lose_healpix:
    dr= d[0:-2]
    if dr=='':
        dr='0'
    preloaded_file=server+'/%s/sightlines_%s_%s/%s/%s-%s-%s-%s-%s-preload_qsos.mat'%(release,survey,program,d,release,survey,program,dr,d)
    catalogue_file=server+'/%s/sightlines_%s_%s/%s/%s-%s-%s-%s-%s-catalog.mat'%(release,survey,program,d,release,survey,program,dr,d)
    processed_file=server+'/%s/newmat/%s-process_qsos.mat'%(release,d)
    dlacatalog_file=server+'/%s/new_dlacat_gp/%s-dlacat_gp.fits'%(release,d)
    qsos=Table.read(server+'/%s/sightlines_%s_%s/%s/%s-%s-%s-%s-%s-catalog.fits'%(release,survey,program,d,release,survey,program,dr,d))
    if len(qsos)==0:
        print('%s no qsos'%d)
    else:
        if np.isin('%s-process_qsos.mat'%d, mat):
            pred_gp_sightline(preloaded_file,catalogue_file, processed_file, dlacatalog_file)
            print('redone %s'%d)
        else:
            print('not processed by gp %s'%d)

'''
'''
release='fuji'
pix={'special':['dark'],'sv1':['backup','bright','dark','other'],
'sv2':['dark','bright'],
'sv3':['backup','bright','dark'],'cmx':['other']}#,'bright'
pix={'main':['bright','dark'],'special':['dark','bright']}
release='guadalupe'
preloaded_file='/home/zjqi/data/highz/preload_qsos.mat'
catalogue_file='/home/zjqi/data/highz/catalog.mat'
processed_file='/home/zjqi/data/highz/process_qsos.mat'
outfile='/home/zjqi/data/highz/dlacat_gp.json'
dlacatalog_file='/home/zjqi/data/highz/dlacat_gp.fits'
dlacatalog_sub_file='/home/zjqi/data/highz/dlacat_gp_pre.fits'
pred_gp_sightline(preloaded_file,catalogue_file, 
                            processed_file, outfile,dlacatalog_file,dlacatalog_sub_file)
'''
'''
for survey in tqdm(list(pix.keys())):
    programs=pix[survey]
    for program in programs:
        print(release,survey,program)
        preloaded_file='/data-192T/home/zjqi/desi/afterburner_fuji_guadalupe/gp/preload/afterburner-%s-%s-%s-preload_qsos.mat'%(release,survey,program)
        catalogue_file='/data-192T/home/zjqi/desi/afterburner_fuji_guadalupe/gp/%s-%s-%s-catalog.mat'%(release,survey,program)
        processed_file='/data-192T/home/zjqi/desi/afterburner_fuji_guadalupe/gp/process/afterburner-%s-%s-%s-process_qsos.mat'%(release,survey,program)
        outfile='/data-192T/home/zjqi/desi/afterburner_fuji_guadalupe/gp/dlacatalog/%s-%s-%s-GP-dlacatalog.json'%(release,survey,program)
        dlacatalog_file='/data-192T/home/zjqi/desi/afterburner_fuji_guadalupe/gp/dlacatalog/%s-%s-%s-GP-dlacatalog.fits'%(release,survey,program)
        #pred_gp_sightline(preloaded_file=preloaded_file,catalogue_file=catalogue_file, 
                            #processed_file=processed_file, sub_dla=True, outfile=outfile,dlacatalog_file=dlacatalog_file)
        
        rename_cats(dlacatalog_file)
'''
'''
qsos_y1 = QSOLoader(
    preloaded_file="/home/zjqi/gp_dla_detection/desimock/data/Y1/processed/preloaded_qsos.mat",
    catalogue_file="/home/zjqi/gp_dla_detection/desimock/data/Y1/processed/catalog.mat",
    processed_file="/home/zjqi/gp_dla_detection/desimock/data/Y1/processed/processed_qsos_multi_meanfluxY1_v3.mat",
    learned_file="/home/zjqi/gp_dla_detection/desimock/data/Y1/processed/learned_qso_model_lyseries_variance_kim_Y1_v2.mat",
    dla_concordance="/home/zjqi/gp_dla_detection/desimock/data/dla_catalogs/real/processed/dla_catalog",
    los_concordance="/home/zjqi/gp_dla_detection/desimock/data/dla_catalogs/real/processed/los_catalog",
    snrs_file="/home/zjqi/gp_dla_detection/desimock/data/Y1/processed/snrs_qsos_multi_Y1.mat",
    sub_dla=True,
    sample_file="/home/zjqi/gp_dla_detection/desimock/data/Y1/processed/dla_samples_a03.mat",
    occams_razor=1)

dla_json = qsos_y1.generate_json_catalogue(outfile="/home/zjqi/gp_dla_detection/desimock/data/Y1/processed/predictions_DLA_Y1_248025.json")

dla_tbl = Table(names=('TARGET_RA','TARGET_DEC', 'ZQSO','Z','TARGETID','S/N','DLAID','NHI','DLA_CONFIDENCE','NHI_STD','ABSORBER_TYPE'),dtype=('float','float','float','float','int','float','str','float','float','float','str'),meta={'EXTNAME': 'DLACAT'})
for ii in range(0,len(dla_json)):
        for jj in range(0,len(dla_json[ii]['dlas'])):
            dla=dla_json[ii]['dlas'][jj]
            absorber_type =  "DLA" if dla['log_nhi'] >= 20.3 else "SUBDLA"
            dla_tbl.add_row((dla_json[ii]['ra'],dla_json[ii]['dec'],dla_json[ii]['z_qso'],dla['z_dla'],dla_json[ii]['thing_id'],dla_json[ii]['snr'],str(dla_json[ii]['thing_id'])+'00'+str(jj),float(dla['log_nhi']),dla_json[ii]['p_dla'],0.0,absorber_type))
dla_tbl.write('/home/zjqi/gp_dla_detection/desimock/data/Y1/processed/dla_catalog_248025.fits',overwrite=True)

qsos_dr12q = QSOLoader(
    preloaded_file="/home/zjqi/gp_dla_detection/data/dr12q/processed/preloaded_qsos_multi.mat",
    catalogue_file="/home/zjqi/gp_dla_detection/data/dr12q/processed/catalog.mat",
    processed_file="/home/zjqi/gp_dla_detection/data/dr12q/processed/processed_qsos_multi_meanfluxdr12q.mat",
    learned_file="/home/zjqi/gp_dla_detection/data/dr12q/processed/learned_qso_model_lyseries_variance_kim_dr9q_minus_concordance.mat",
    dla_concordance="/home/zjqi/gp_dla_detection/data/dla_catalogs/dr9q_concordance/processed/dla_catalog",
    los_concordance="/home/zjqi/gp_dla_detection/data/dla_catalogs/dr9q_concordance/processed/los_catalog",
    snrs_file="/home/zjqi/gp_dla_detection/data/dr12q/processed/snrs_qsos_multi_dr12q.mat",
    sub_dla=True,
    sample_file="/home/zjqi/gp_dla_detection/data/dr12q/processed/dla_samples_a03.mat",
    occams_razor=1)

#all_nspecs = np.where(qsos_dr12q.p_dlas > 0.90)
for nspec in [1,2,3,4,5,6,7]:

    #print(nspec)
    qsos_y1.plot_this_mu(nspec)


qsos_sv = QSOLoader(
    preloaded_file="/home/zjqi/gp_dla_detection/SV/data/cascades/processed/preloaded_qsos.mat",
    catalogue_file="/home/zjqi/gp_dla_detection/SV/data/cascades/processed/catalog.mat",
    processed_file="/home/zjqi/gp_dla_detection/SV/data/cascades/processed/processed_qsos_multi_meanfluxcascades_v2.mat",
    learned_file="/home/zjqi/gp_dla_detection/desimock/data/Y1/processed/learned_qso_model_lyseries_variance_kim_Y1_v2.mat",
    dla_concordance=[],
    los_concordance=[],
    snrs_file="/home/zjqi/gp_dla_detection/SV/data/cascades/processed/snrs_qsos_multi_cascades.mat",
    sub_dla=True,
    sample_file="/home/zjqi/gp_dla_detection/desimock/data/Y1/processed/dla_samples_a03.mat",
    occams_razor=1)

#all_specs = np.where(qsos_sv.p_dlas > 0.90)
#for nspec in all_specs[0:100]:
    #qsos_sv.plot_this_mu(nspec)

dla_json = qsos_sv.generate_json_catalogue(outfile="/home/zjqi/gp_dla_detection/SV/data/cascades/processed/predictions_DLA_Y1_gp_v2.json")

dla_tbl = Table(names=('TARGET_RA','TARGET_DEC', 'ZQSO','Z','TARGETID','S/N','DLAID','NHI','DLA_CONFIDENCE','NHI_STD','ABSORBER_TYPE'),dtype=('float','float','float','float','int','float','str','float','float','float','str'),meta={'EXTNAME': 'DLACAT'})
for ii in range(0,len(dla_json)):
        for jj in range(0,len(dla_json[ii]['dlas'])):
            dla=dla_json[ii]['dlas'][jj]
            absorber_type =  "DLA" if dla['log_nhi'] >= 20.3 else "SUBDLA"
            dla_tbl.add_row((dla_json[ii]['ra'],dla_json[ii]['dec'],dla_json[ii]['z_qso'],dla['z_dla'],dla_json[ii]['thing_id'],dla_json[ii]['snr'],str(dla_json[ii]['thing_id'])+'00'+str(jj),float(dla['log_nhi']),dla_json[ii]['p_dla'],0.0,absorber_type))
dla_tbl.write('/home/zjqi/gp_dla_detection/desimock/data/cascades/processed/dla_catalog_gp_487.fits',overwrite=True)

min_lambda = 911.75#qsos_dr16q.GP.min_lambda
max_lambda = 1215.75#qsos_dr16q.GP.max_lambda
scale = np.shape(qsos_y1.GP.C)[0] / (max_lambda - min_lambda)

# some wavelegnths for emission lines to plot
lya_wavelength=1215.67
lyb_wavelength=1025.72
lyman_limit=911.76
lyg_wavelength = 972
nv_wavelength = 1240.81
oi_wavelength = 1305.53
cii_wavelength = 1335.31
siv_wavelength = 1399.8
fig, ax = plt.subplots(figsize=(8, 8))
im = ax.imshow(qsos_y1.GP.C, origin="lower")
ax.set_xticks(
    [
        (lyman_limit - min_lambda) * scale,
        (lyg_wavelength - min_lambda) * scale,
        (lyb_wavelength - min_lambda) * scale,
        (lya_wavelength - min_lambda) * scale,
        (nv_wavelength - min_lambda) * scale,
        (oi_wavelength - min_lambda) * scale,
        (cii_wavelength - min_lambda) * scale,
        (siv_wavelength - min_lambda) * scale,
    ],
)
ax.set_xticklabels(
    [
        r"Ly $\infty$",
        r"Ly $\gamma$",
        r"Ly $\beta$",
        r"Ly $\alpha$",
        r"NV",
        r"OI",
        r"CII",
        r"SIV",
    ],
    rotation=45,
)
ax.set_yticks(
    [
        (lyman_limit - min_lambda) * scale,
        (lyg_wavelength - min_lambda) * scale,
        (lyb_wavelength - min_lambda) * scale,
        (lya_wavelength - min_lambda) * scale,
        (nv_wavelength - min_lambda) * scale,
        (oi_wavelength - min_lambda) * scale,
        (cii_wavelength - min_lambda) * scale,
        (siv_wavelength - min_lambda) * scale,
    ]
)
ax.set_yticklabels(
    [
        r"Ly $\infty$",
        r"Ly $\gamma$",
        r"Ly $\beta$",
        r"Ly $\alpha$",
        r"NV",
        r"OI",
        r"CII",
        r"SIV",
    ]
)

fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

plt.tight_layout()

plt.savefig('/home/zjqi/gp_dla_detection/desimock/data/Y1/processed/covariance_mat_487.png')

'''
