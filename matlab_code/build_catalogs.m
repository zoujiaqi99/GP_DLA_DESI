% build_catalogs: loads existing QSO and DLA catalogs, applies some
% initial filters, and creates a list of spectra to download from SDSS
%
% ZWARNING: ensure we exclude those spectra with bad redshift status reported

% load QSO catalogs
release='Y1';
qso_catalog = ...
    fitsread(sprintf('%s/sum-qsocatalog.fits',distfiles_directory(release)), ...
             'binarytable');

% extract basic QSO information from DR12Q catalog
%sdss_names       =  qso_catalog{1};
plates           =  qso_catalog{1};
fiber_ids        =  qso_catalog{2};
mjds             =  qso_catalog{3};
ras              =  qso_catalog{4};
decs             =  qso_catalog{5};
z_qsos           =  qso_catalog{6};
target_ids        =  qso_catalog{7};
snrs             =  qso_catalog{8};

num_quasars      = numel(z_qsos);
bal_visual_flags = zeros(num_quasars, 1, 'uint8');
zwarning         =  zeros(num_quasars, 1, 'uint8');





% determine which objects in DR12Q are in DR10Q and DR9Q, using SDSS
% thing IDs
%in_dr9  = ismember(thing_ids,  dr9_catalog{4});
%in_dr10 = ismember(thing_ids, dr10_catalog{4});

% to track reasons for filtering out QSOs
filter_flags = zeros(num_quasars, 1, 'uint8');

% filtering bit 0: z_QSO < 2.15
ind = (z_qsos < z_qso_cut);
filter_flags(ind) = bitset(filter_flags(ind), 1, true);

% filtering bit 1: BAL
ind = (bal_visual_flags>0);
filter_flags(ind) = bitset(filter_flags(ind), 2, true);

% filtering bit 4: ZWARNING
ind = (zwarning > 0);
%% but include `MANY_OUTLIERS` in our samples (bit: 1000)
%ind_many_outliers      = (zwarning == bin2dec('10000'));
%ind(ind_many_outliers) = 0;
filter_flags(ind) = bitset(filter_flags(ind), 5, true);

%los_inds = containers.Map();
%dla_inds = containers.Map();
%z_dlas   = containers.Map();
%log_nhis = containers.Map();

% load available DLA catalogs

  % determine lines of sight searched in this catalog
name='real';
los_catalog = ...
      load(sprintf('%s/los_catalog', dla_catalog_directory(name)));
los_inds = ismember(target_ids, los_catalog);

dla_catalog = ...
      load(sprintf('%s/dla_catalog', dla_catalog_directory(name)));

  % determine DLAs flagged in this catalog
[dla_inds, ind] = ismember(target_ids, dla_catalog(:, 1));
ind = find(ind);

  % determine lists of DLA parameters for identified DLAs, when
  % available
z_dlas   = cell(num_quasars, 1);
log_nhis = cell(num_quasars, 1);
for i = 1:numel(ind)
  this_dla_ind = (dla_catalog(:, 1) == target_ids(ind(i)));
  z_dlas{ind(i)}   = dla_catalog(this_dla_ind, 2);
  log_nhis{ind(i)} = dla_catalog(this_dla_ind, 3);
end




% save catalog
variables_to_save = { 'ras', 'decs', 'target_ids', 'plates', ...
                     'mjds', 'fiber_ids', 'z_qsos', 'snrs', ...
                     'bal_visual_flags', 'filter_flags', ...
                     'los_inds', 'dla_inds', 'z_dlas', 'log_nhis', ...
                     'zwarning'};
save(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_save{:}, '-v7.3');

% these plates use the 5.7.2 processing pipeline in SDSS DR12

% build file list for SDSS DR12Q spectra to download (i.e., the ones
% that are not yet removed from the catalog according to the filtering
% flags)
%fid = fopen(sprintf('%s/file_list', spectra_directory(release)), 'w');
%for i = 1:num_quasars
  %if (filter_flags(i) > 0)
    %continue;
  %end

  % for 5.7.2 plates, simply print greedily print both 5.7.0 and 5.7.2 paths
  %if (v_5_7_2_ind(i))
    %fprintf(fid, 'v5_7_2/spectra/lite/./%i/spec-%i-%i-%04i.fits\n', ...
            %plates(i), plates(i), mjds(i), fiber_ids(i));
  %end

  %fprintf(fid, 'v5_7_0/spectra/lite/./%i/spec-%i-%i-%04i.fits\n', ...
          %plates(i), plates(i), mjds(i), fiber_ids(i));
%end
%fclose(fid);

for i in 1:numel(all_flux)
    all_flux{i}=all_flux{i}(:);
    all_wavelengths{i}=all_wavelengths{i}(:);
    all_noise_variance{i}=all_noise_variance{i}(:);
