global cfg
%run('D:\SOLSA\HIMIP\20191007_solsa_himip\utils\Config.m');
run(fullfile(pwd, 'utils','Config.m'))

% Update of the dir_data
dir = 'D:\SOLSA\DATA\111_fused\Ech_ref\';
%samples = [{'ER-NC00-0056'},{'ER-NC00-0059'},{'ER-NC00-0049_1'},{'ER-NC00-0049'},{'ER-NC00-0054'},{'ER-NC00-0051'},{'ER-NC00-0061'},{'ER-NC00-0053'},{'ER-NC00-0060'},{'GR'}];
samples = [{'ER-NC00-0056'},{'ER-NC00-0059'},{'ER-NC00-0049_1'}];

% Library
%dir_lib = 'D:\SOLSA\HIMIP\Hyper_Spectral_Libraries';
dir_lib = fullfile(pwd, '.', 'Hyper_Spectral_Libraries');
library_file = 'Library_fus_VNIR_SWIR_final.csv';

spectrum_types = [1,7];
corsen_rate = 2;
att_names = [{'VNIR'}, {'SWIR'}];

min_max_wl = [450,999,1001,2500];

min_max_xy = [];

nb_repl = [1,2,4];
nb_clus = [500, 1000, 2000, 5000];

p_val = zeros(length(samples),length(repl)*length(nb_clus));

times = zeros(length(samples),(1 + length(repl)) * length(nb_clus));

for s = 1:length(samples)
    dir_data = fullfile(dir,'XYZ_DFusion',samples{s});
    % Load the data 
    img = Multi_att_Img();
    img = load_fused_data(img,dir,dir_data);
    
    % Load the library
    path_library = fullfile(dir_lib,library_file);
    hs_lib = Multi_att_Lib();
    hs_lib = load_csv(hs_lib,path_library);
    
    lib_names_sep = cellfun(@(x) split( x , "_"), hs_lib.entry_names, 'un',0);
    lib_names = cellfun(@(x) x{1}, lib_names_sep, 'un',0);
    lib_names_uni = unique(lib_names);
    
    % Trims the spectra 
    att_names_C = att_names(strcmp(img.att_types(contains(att_names,img.att_groups)),'C'));
    for a = 1:length(att_names_C)
        min_val.(att_names_C{a}) = min_max_wl(2*a - 1);
        max_val.(att_names_C{a}) = min_max_wl(2*a);
        % Remove the unnecessary wavelengths
        img = rm_min_max_wls(img,att_names_C{a},min_val.(att_names_C{a}),max_val.(att_names_C{a}));
    end
    
    % Crop the image
    img = img_modify(img,corsen_rate,min_max_xy);
    
    time_tmp = [];
    p_val_tmp = [];
    for n = 1:length(nb_clus)
        params.nb_clust = nb_clus(n);    
        for r = 1:length(nb_repl)
            params.nb_repl = nb_repl(r);
            nb_ite = 50;
            res_prop = zeros(length(lib_names_uni),nb_ite);
            for i = 1:nb_ite
                tic
                % Application of the method 
                find_mineral = Multi_Classifier();
                find_mineral = method_2(find_mineral,img,hs_lib,att_names,spectrum_types,min_max_wl,[],params);
                
                % Creation of mineral table
                res_tmp = zeros(length(lib_names_uni),1);
                
                for name = 1:length(find_mineral.table_minerals_ind.Minerals)
                    ind_tmp = find(ismember(lib_names_uni,find_mineral.table_minerals_ind.Minerals{name}));
                    res_tmp(ind_tmp) = str2double(find_mineral.table_minerals_ind.Prop{name});
                end
                res_prop(:,i) = res_tmp;
                
                time_tmp = [time_tmp, toc]
            end
            p_val_tmp = [p_val_tmp, get_pval(res_prop)];
        end
    end
    
    p_val(s,:) = p_val_tmp;
end

p_val

function p_val = get_pval(res_prop)
    mean_prop = mean(res_prop,2);
    chi2_all = (res_prop - mean_prop).^2./mean_prop;
    chi2 = nansum(chi2_all);
    
    p_value = zeros(1,length(chi2));
    for p = 1:length(chi2)
        p_value(p) = chi2cdf(chi2(p), size(res_prop,1), 'upper');
    end
    
    p_val = mean(p_value);
end


