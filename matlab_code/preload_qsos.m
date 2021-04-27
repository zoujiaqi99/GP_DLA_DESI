% preload_qsos: loads spectra from SDSS FITS files, applies further
% filters, and applies some basic preprocessing such as normalization
% and truncation to the region of interest
% Modification: matching ids in catalog.mat and sightlines.mat
% load QSO catalog
release='Y1';
variables_to_load = {'z_qsos', 'target_ids', 'filter_flags'};
load(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_load{:});
load(sprintf('%s/sightlines', processed_directory(release)));

num_quasars = numel(z_qsos);

all_wavelengths    =  cell(num_quasars, 1);
all_flux           =  cell(num_quasars, 1);
all_noise_variance =  cell(num_quasars, 1);
all_pixel_mask     =  cell(num_quasars, 1);
all_normalizers    = zeros(num_quasars, 1);


for i = 1:num_quasars
  if (filter_flags(i) > 0)
    continue;
  end
  
  [sightline_ind,ind]=ismember(target_ids(i),sightline_ids);
  if sightline_ind == 0
    filter_flags(i) = bitset(filter_flags(i), 4, true);
    continue;
  end
  
  this_wavelengths=wavelengths{ind}(:);
  this_flux=flux{ind}(:);
  this_noise_variance=noise_variance{ind}(:);
  this_median=normalizers(ind);

  % bit 2: cannot normalize (all normalizing pixels are masked)
  if (isnan(this_median))
    filter_flags(i) = bitset(filter_flags(i), 3, true);
    continue;
  end

  all_wavelengths{i}    =    this_wavelengths;
  all_flux{i}           =           this_flux;
  all_noise_variance{i} = this_noise_variance;
  all_normalizers(i)     =     this_median;

  fprintf('loaded quasar %i of %i, ind=%i \n', ...
          i, num_quasars,ind);
end

variables_to_save = {'loading_min_lambda', 'loading_max_lambda', ...
                     'normalization_min_lambda', 'normalization_max_lambda', ...
                     'min_num_pixels', 'all_wavelengths', 'all_flux', ...
                     'all_noise_variance',  ...
                     'all_normalizers'};
save(sprintf('%s/preloaded_qsos', processed_directory(release)), ...
     variables_to_save{:}, '-v7.3');

% write new filter flags to catalog
%save(sprintf('%s/catalog', processed_directory(release)), ...
     %'filter_flags', '-append');
