cd /home/zjqi/gp_dla_detection/Everest/code
set_parameters_multi
mex voigt.c -lcerf
% do not need to change lls parameters
set_lls_parameters
release = 'jura';
survey = 'main';
program = 'dark';
path = sprintf('/home/zjqi/data/desi/%s/mat_%s_%s',release,survey,program)
dirs=dir(path);
for j = 4:length(dirs)
  
  pix=dirs(j).name;
  d=pix(1:end-2);
  if isempty(d)
      d = '0';
  end
  catalog_path = sprintf('/home/zjqi/data/desi/%s/mat_%s_%s/%s/%s-%s-%s-%s-%s-catalog',release,survey,program,pix,release,survey,program,d,pix);
  preload_path = sprintf('/home/zjqi/data/desi/%s/mat_%s_%s/%s/%s-%s-%s-%s-%s-preload_qsos',release,survey,program,pix,release,survey,program,d,pix);
  filename = sprintf('/home/zjqi/data/desi/%s/mat_%s_%s/%s/%s-%s-%s-%s-%s-process_qsos',release,survey,program,pix,release,survey,program,d,pix);
  try
      process_qsos_multiple_dlas_meanflux
  catch 
      warning('Problem')   
  end
end
%for j = 1:length(healpixs)
%    pix = int2str(healpixs(j));
%    d = pix(1:end-2);
 %   if isempty(d)
  %    d = '0';c
   %   catalog_path = sprintf('/home/zjqi/data/desi/iron/sightlines/main/dark/%s/%s/catalog',d,pix);
    %  preload_path = sprintf('/home/zjqi/data/desi/iron/sightlines/main/dark/%s/%s/preload_qsos',d,pix);
     % filename = sprintf('/home/zjqi/data/desi/iron/mat/%s-%s-process_qsos',d,pix);
      %try 
       % process_qsos_multiple_dlas_meanflux
      %catch err
        %sprintf('%s:%s',pix,err.message)
      %end
    %end
%end
        
%catalog_path = '/home/zjqi/data/highz/catalog'
%preload_path = '//home/zjqi/data/highz/preload_qsos'
%filename = '/home/zjqi/data/highz/process_qsos'
%process_qsos_multiple_dlas_meanflux