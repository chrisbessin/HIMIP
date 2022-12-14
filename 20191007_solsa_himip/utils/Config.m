global cfg

cfg.dir_proj = 'D:\\SOLSA\\HIMIP';
cfg.dir_path_library = [cfg.dir_proj, '\\Hyper_Spectral_Libraries'];
cfg.path_library = [cfg.dir_path_library, '\\Library_fus_VNIR_SWIR_final.csv'];
cfg.dir_data_fus = 'D:\\SOLSA\\DATA\\111_fused\\Ech_ref';
cfg.dir_data = 'D:\\SOLSA\\DATA\\111_fused\\Ech_ref';
cfg.file_data = 'ER-NC00-0050_2_SWIR_04122018_16-04-18.raw';
cfg.darkref_tag = 'DARKREF_';
cfg.whiteref_tag = 'WHITEREF_';
cfg.red_default = 1400;
cfg.green_default = 1900;
cfg.blue_default = 2310;

% Bounds spectra
cfg.min_wl_vnir = 450;
cfg.max_wl_vnir = 999;
cfg.min_wl_swir = 1071;
cfg.max_wl_swir = 2500;
cfg.lim_vnir_swir = 1000;

% Parameters ROI
cfg.average_density_samples = 2;
cfg.diam_sample = 67; % mm
cfg.mass_ROI = 0.2; % kg

% Parameters classif
cfg.spectrum_types = [6,6];
cfg.corsen_x = 8;
cfg.corsen_y = 8;
cfg.att_names = [{'VNIR'}, {'SWIR'}];
cfg.nb_cluster = 1000;
cfg.process_sample = 'lib_tree';
cfg.nb_tree = 100;

% print times
cfg.print_times = 0;

% Unused
cfg.elements = {'Si','O','H','Fe','Mg','Al','Ni','Co','Cr','Mn','Na','Ca','C','S'};
cfg.elements_masse_mol = [28.1 16.0 1.0 55.8 24.3 27.0 58.7 58.9 52.0 54.9 23.0 40.1 12.0 32.1];
cfg.elts_mass_dict = containers.Map(cfg.elements,cfg.elements_masse_mol);
cfg.res_x_swir = 169.5;
cfg.res_y_swir = 120.0;