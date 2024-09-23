% script_for_gp_minibatch.m
% This script is used to run the GP in minibatches
% 
% Parameters:
% ------------
%  num_quasars    : the size of the batch / saved array
%  qsos_num_offset: the starting quasar_ind for the processed quasars
% 
% Example:
% ------------
%  Run 1st to 1000th quasars in the catalog
%  num_quasars = 1000;
%  qsos_num_offset = 0;
%  script_for_gp_minibatch
%
%  Run 1001st to 2000th quasars in the catalog
%  num_quasars = 1000;
%  qsos_num_offset = 1000;
%  script_for_gp_minibatch


cd ..
% Want to use the files in embarrassingly_parallel folder
addpath embarrassingly_parallel


% set parameters for the GP
set_parameters_multi

% compile the voigt profile
mex voigt.c -lcerf

% do not need to change lls parameters
set_lls_parameters

release = 'jura';
survey  = 'main';
program = 'dark';

path    = sprintf('.../%s/mat_%s_%s',release,survey,program) %this is the dir where you save preprocess mat files.
dirs    = dir(path);

for j = 4:length(dirs)
  
  pix=dirs(j).name; %this is the healpix number
  d=pix(1:end-2);
  if isempty(d)
      d = '0';
  end
  catalog_path = sprintf('.../%s/mat_%s_%s/%s/%s-%s-%s-%s-%s-catalog',release,survey,program,pix,release,survey,program,d,pix);
  preload_path = sprintf('.../%s/mat_%s_%s/%s/%s-%s-%s-%s-%s-preload_qsos',release,survey,program,pix,release,survey,program,d,pix);

  % Set the filename with the quasar indices, e.g. jura-main-dark-0-1000-process_qsos_ind_1-1000
  % so output files are in minibatches, and we can run them in parallel
  % Example filename: .../jura/mat_main_dark/{pix}/jura-main-dark-{d}-{pix}-process_qsos_ind_1-1000
  filename = sprintf('.../%s/mat_%s_%s/%s/%s-%s-%s-%s-%s-process_qsos_ind_%d-%d', ...
      release, ...
      survey, program, ...
      pix, ... 
      release, survey, program, d, pix, ...
      quasar_ind_offset + 1, quasar_ind_offset + num_quasars);
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
