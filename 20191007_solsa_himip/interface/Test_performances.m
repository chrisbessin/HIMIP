global cfg
%run('D:\SOLSA\HIMIP\20191007_solsa_himip\utils\Config.m');
run(fullfile(pwd, 'utils','Config.m'))

% Update of the dir_data
dir = 'D:\SOLSA\DATA\111_fused\Ech_ref\';
%samples = [{'ER-NC00-0056'},{'ER-NC00-0059'},{'ER-NC00-0049_1'},{'ER-NC00-0049'},{'ER-NC00-0054'},{'ER-NC00-0051'},{'ER-NC00-0061'},{'ER-NC00-0053'},{'ER-NC00-0060'}];
samples = [{'GR'}];

spectrum_types = [1,7];
corsen_rate = 2;
att_names = [{'VNIR'}, {'SWIR'}];

min_max_wl = [450,999,1001,2500];

min_max_xy = [];

repl = [1,2,4,10];
nb_clus = [10, 100, 500, 1000];

times = zeros(length(samples),(1 + length(repl)) * length(nb_clus));

for s = 1:length(samples)
    dir_data = fullfile(dir,'XYZ_DFusion',samples{s});
    % Load the data 
    img = Multi_att_Img();
    img = load_fused_data(img,dir,dir_data);
    
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

    % Creation of the classifier
    classif = Unit_Classifier(img);

    time_tmp = [];
    
    for n = 1:length(nb_clus)
        nb_cluster = nb_clus(n);
        
        tic
        %{
        hierC_metric = 'correlation';
        cutoff = 0.01;
        classif_cah = Classifier_UnSuprv_HierC(classif,hierC_metric,nb_cluster);
        %}
        size(classif.ma_img.att_data,1)
        time_tmp = [time_tmp, toc]
        
        for r = 1:length(repl)
            tic
            clf_metric = 'sqeuclidean'; %'sqeuclidean' 'cityblock', 'cosine', 'correlation', 'hamming'
            nb_replicate = repl(r);
            max_iter = 100;
            classif_kmeans = Classifier_UnSuprv_Kmeans(classif,clf_metric,nb_cluster,nb_replicate,max_iter);
            time_tmp = [time_tmp, toc]
        end
    end
    
    times(s,:) = time_tmp;
end

times

