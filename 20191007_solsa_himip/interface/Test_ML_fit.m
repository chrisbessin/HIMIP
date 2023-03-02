% Config
global cfg
%run('D:\SOLSA\HIMIP\20191007_solsa_himip\utils\Config.m');
run(fullfile(pwd, 'utils','Config.m'))

% Library
%dir_lib = 'D:\SOLSA\HIMIP\Hyper_Spectral_Libraries';
dir_lib = fullfile(pwd, '.', 'Hyper_Spectral_Libraries');
library_file = 'Library_fus_VNIR_SWIR_final.csv';

% param sample
spectrum_types = [3,3];
corsen_rate = 2;
att_names = [{'VNIR'}, {'SWIR'}];
min_max_wl = [450,999,1001,2500];
min_max_xy = [];
dir = fullfile(cfg.dir_data,'ER-NC00-0056');
dir_data = fullfile(dir,'XYZ_DFusion');

% Load the library
path_library = fullfile(dir_lib,library_file);
lib_l = Multi_att_Lib();
lib_l = load_csv(lib_l,path_library);
lib_l = rm_min_max_wls(lib_l,att_names{1},min_max_wl(1),min_max_wl(2));
lib_l = rm_min_max_wls(lib_l,att_names{2},min_max_wl(3),min_max_wl(4));
library = apply_change_spec(lib_l,att_names,spectrum_types,0);

% image
img = Multi_att_Img();
img = load_fused_data(img,dir,dir_data);


% Crop the image
img_tmp = img_modify(img,corsen_rate,min_max_xy);
[library,img_tmp] = consis_multi_att(library,img_tmp,att_names,1);

% Change lib into a ratios library
library = to_subs_ref(library,'all',[]);

% Get the minerals (labels)
tmp_entry_names = split(library.entry_names{1},'_');
y = cellstr(tmp_entry_names{1});
for i = 2:length(library.entry_names)
    tmp_entry_names = split(library.entry_names{i},'_');
    y = [y; tmp_entry_names{1}];
end



% test number of weak learners
nb_learners = [25,50,75,100,200,300,400,500,600,750,1000,1500,2000];
nb_learners = [2000];
ooberrors = NaN(1,length(nb_learners));
times_t = NaN(1,length(nb_learners));
times_pred = NaN(1,length(nb_learners));
for i = 1:length(nb_learners)
    %{
    nb_learners(i)
    tic
    model = TreeBagger(nb_learners(i),library.att_data,y,'OOBPrediction','on','OOBPredictorImportance','on');
    file = [cfg.dir_path_library, '\\model_ML_library_RF.mat'];
    save(file,'model');
    ooberrors(i) = mean(model.oobError)
    times_t(i) = toc
    %}
    multi_class = Multi_Classifier();
    for j = 1:2
        tic
        cfg.nb_tree = nb_learners(i);
        multi_class = method(multi_class,'lib_tree',img,lib_l,att_names,spectrum_types,min_max_wl,[]);
        toc
    end
end
ooberrors
times_t
times_pred
figure; plot(nb_learners_s,ooberrors_s)

%{
% Check outliers
min_dist = 0.9;
model = fillprox(model);
[y,cellstr(num2str(model.OutlierMeasure))]

prox = model.Proximity;
for i=1:length(prox)
    prox(i,i) = 0;
end
ind_max = max(prox > min_dist);
[{''},y(ind_max)';y(ind_max),mat2cell(prox(ind_max,ind_max),ones(1,sum(ind_max)),ones(1,sum(ind_max)))]


% print duplicates
mat_ind = prox > min_dist;
dupl_1 = [{}];
dupl_2 = [{}];
for i = 1:length(prox) - 1
    for j = i + 1:length(prox)
        if mat_ind(i,j)
            dupl_1 = [dupl_1, library.entry_names(i)];
            dupl_2 = [dupl_2, library.entry_names(j)];
        end
    end
end
[dupl_1;dupl_2]'
%}