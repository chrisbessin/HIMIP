% Config
% Initialisations
global cfg
run('D:\SOLSA\HIMIP\20191007_solsa_himip\utils\Config.m');

% Library
dir_lib = cfg.dir_path_library;
library_file = 'Library_fus_VNIR_SWIR_final.csv';

% Load the library
path_library = fullfile(dir_lib,library_file);
library = Multi_att_Lib();
library = load_csv(library,path_library);

%
spectrum_types = [6,6];
att_names = [{'VNIR'}, {'SWIR'}];
min_max_wl = [450,999,1001,2500];
min_max_xy = [];
methods = [{'unmixing'},...
           {'dist'},...
           {'clustering'},...
           {'lib_tree'},...
           {'custom_tree'}...
           ];
samples = [{'11MC_harz'},...
           {'11MC_harz'},...
           {'11MC_kaolinite'},...
           {'11MC_sap_terreuse'},...
           {'BST_1'},...
           {'ER-NC00-0046'},...
           {'ER-NC00-0046-53-54_slice_Anas'},...
           {'ER-NC00-0047'},...
           {'ER-NC00-0048'},...
           {'ER-NC00-0049'},...
           {'ER-NC00-0051'},...
           {'ER-NC00-0053'},...
           {'ER-NC00-0054'},...
           {'ER-NC00-0055'},...
           {'ER-NC00-0055-56-60_slice_Anas'},...
           {'ER-NC00-0056'},...
           {'ER-NC00-0057'},...
           {'ER-NC00-0059'},...
           {'ER-NC00-0060'},...
           {'ER-NC00-0061'},...
           {'ER-NC00-0063'},...
           {'ER-NC00-0063-73_slice_Anas'},...
           {'ER-NC00-0064'},...
           {'ER-NC00-0067'},...
           {'ER-NC00-0068'},...
           {'ER-NC00-0070'},...
           {'ER-NC00-0071'},...
           {'ER-NC00-0072'},...
           {'ER-NC00-0073'},...
           {'GR'},...
           {'HN5'},...
           {'LJ_2'},...
           {'LR'}];


dir_report = fullfile(cfg.dir_data,'_rapport_samples');
times = cell(length(samples),length(methods)+1);
times(:) = {''};

for s = 1:length(samples)
    dir = fullfile(cfg.dir_data,samples{s});
    dir_data = fullfile(dir,'XYZ_DFusion');
    % Load the data 
    img = Multi_att_Img();
    tic
    img = load_fused_data(img,dir,dir_data);
    img = create_img(img, 1400, 1900, 2310, 2, []);
    saveas(img.img_rgb,fullfile(dir_report,[samples{s},'_rgb.png']));
    
    times(s,1) = {num2str(toc)};
    for m = 1:length(methods)
        % Application of the method 
        multi_class = Multi_Classifier();
        tic
        multi_class = method(multi_class,methods{m},img,library,att_names,spectrum_types,min_max_wl,[]);
        
        times(s,m + 1) = {num2str(toc)};
        
        writetable(multi_class.(strcat('table_minerals_ind')),fullfile(dir_report,[samples{s},'_',methods{m},'_table_minerals.csv']));
        saveas(multi_class.('ROI_minerals_ind_rgb'),fullfile(dir_report,[samples{s},'_',methods{m},'_minerals_ind.png']));
        %csvwrite(fullfile(dir_report,[samples{s},'_',methods{m},'_misc.csv']),misc);
    end
end
% Export time table
table_times = cell2table(times);
table_times.Properties.VariableNames = [{'Loading'},methods];
writetable(table_times,fullfile(dir_report,'all_times.csv'));

