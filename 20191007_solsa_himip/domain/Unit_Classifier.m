classdef Unit_Classifier
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ma_img
        prop_class
        error
        entry_names
        nb_cluster
        
        % for unsupervized classifiers
        metric
        cutoff
        
        % for supervized classifiers
        library
        
        % unmixing
        method
        mineral_maps
    end
    
    methods
        function obj = Unit_Classifier(ma_img)
            % Creation of the object
            obj.ma_img = ma_img;
        end
        
        function classif = Classifier_UnSuprv_Kmeans(obj,metric,nb_cluster,nb_replicate,max_iter,n_sub_group)
            %
            warning ('off','all');
            
            classif = obj;
            input_attrib = classif.ma_img.att_data;
            
            if size(input_attrib,1) > nb_cluster
                if nargin < 6
                    n_sub_group = 1;
                end
                
                % Subdivision of the image to improve the performance
                opt = 'sep_kmeans';
                %opt = 'slice';
                if strcmp(opt, 'sep_kmeans')
                    [class_sep, ~,dist] = kmeans(input_attrib, n_sub_group,...
                                                 'Distance', metric,...
                                                 'Replicates', nb_replicate,...
                                                 'MaxIter', max_iter);
                elseif strcmp(opt, 'slice')
                    s_mini_in_attrib = ceil(size(input_attrib,1) / n_sub_group);
                    class_sep = ceil((1:size(input_attrib,1))./s_mini_in_attrib);
                    dist = ones(n_cluster);
                    if n_sub_group > nb_cluster/5
                        fprinf('Warning: less than 5 cluster per class, they are set to 1/5 of the n_cluster')
                        n_sub_group = ceil(nb_cluster/5);
                    end
                end
                
                % Number of cluster searched per sub_group
                [~, ~, ic] = unique(class_sep(~isnan(class_sep)));
                % Number of instance of each sub_group
                a_counts = accumarray(ic, 1);  
                % multiplication by the K Mean distance 
                a_counts = a_counts.*dist;  
                % normalizing the weights to 1, and get the number of
                % classes by sub_group
                nb_clusters = ceil(a_counts./sum(a_counts)*nb_cluster);
                
                % The excess of classes are removed from the largest group
                [~,max_v] = max(nb_clusters);
                nb_clusters(max_v) = nb_clusters(max_v) - (sum(nb_clusters) - nb_cluster);
                
                % Initilisations
                att_data_tmp = [];
                kmeans_class = [];
                max_tmp = 0;
                for i = 1:n_sub_group
                    % Apply K means
                    [mini_class, min_data] = ... 
                        kmeans(input_attrib(class_sep == i,:),...
                               nb_clusters(i),...
                               'Distance',metric,...
                               'Replicates',nb_replicate,...
                               'MaxIter',max_iter);
                    
                    %  Add results 
                    att_data_tmp = [att_data_tmp; min_data];
                    kmeans_class = [kmeans_class; mini_class + max_tmp];
                    
                    % Update the max_tmp for the class labels
                    max_tmp = max(kmeans_class);
                end
                
                % Add the NaN to classes
                kmeans_class = [kmeans_class; NaN(sum(isnan(class_sep)),1)];
                
                % Reordering the classes
                a = [class_sep,(1:length(class_sep))'];
                b = sortrows(a,1);
                c = sortrows([b,kmeans_class],2);
                kmeans_class = c(:,3);
                
                % Set the att_data
                classif.ma_img.att_data = att_data_tmp;
                
                % Reconstruction of the image and Proportion of classes
                classif = compute_output_img_prop_class(classif,kmeans_class);
                
                % Other properties
                classif.nb_cluster = nb_cluster;
                classif.metric = metric;
            else
                fprintf('Not enough classes, skipping the unsupervized classification.')
            end
        end
        
        function classif = Classifier_UnSuprv_ImSegKmeans(obj,nb_cluster,nb_replicate,max_iter,n_slice)
            % Application of the Image Segment KMeans on the Multi_att_Img
            %
            
            classif = obj;
            
            input_attrib = classif.ma_img.att_data;
            
            if size(input_attrib,1) > nb_cluster
                if nargin < 5
                    n_slice = 1;
                end
                
                if n_slice > nb_cluster/5
                    fprinf('Warning: less than 5 cluster per class, they are set to 1/5 of the n_cluster')
                    n_slice = ceil(nb_cluster/5);
                end
                
                % Transform back into an image
                data_img = reshape_att_data(classif.ma_img);

                % Transform data to int (imsegkmeans requires int)
                data_int_prep = uint16(reshape(data_img,[size(classif.ma_img.ind_img),size(classif.ma_img.att_data,2)])*1000);
                
                % Initilisations
                s_mini_in_attrib = ceil(size(data_int_prep,1) / n_slice);
                mini_n_clust = ceil(nb_cluster/n_slice);
                classif.ma_img.att_data = [];
                imkmeans_class = [];
                ind_end = 0;
                mini_n_clust_tmp = 0;
                max_tmp = 0;
                
                for i = 1:n_slice
                    % Apply K means
                    
                    % Starting and ending indices of the mini input attrib
                    ind_sta = ind_end + 1;
                    ind_end = min(ind_sta - 1 + s_mini_in_attrib,size(data_int_prep,1));
                    
                    % Mini clusters
                    mini_n_clust_tmp = mini_n_clust_tmp + mini_n_clust;
                    if mini_n_clust_tmp > nb_cluster
                        mini_n_clust = mini_n_clust + nb_cluster - mini_n_clust_tmp;
                    end
                    
                    % %%%%%%% Imseg algo
                    % Application of ImSegKMeans
                    [mini_class, min_data] = imsegkmeans(data_int_prep(ind_sta:ind_end,:,:),mini_n_clust + 1,'NumAttempts',nb_replicate,'MaxIterations',max_iter);
                    
                    % Removing the NaN center
                    nan_idx = find(sum(min_data,2) == 0,1);
                    min_data = min_data(sum(min_data,2) ~= 0,:);

                    % Delete nan index from the image
                    mini_class = double(mini_class);
                    if ~isempty(nan_idx)
                        mini_class(mini_class == nan_idx) = NaN;
                        mini_class(mini_class > nan_idx) = mini_class(mini_class > nan_idx) - 1;
                    end
                    
                    % Transformation of the centers back to double
                    min_data = double(min_data)/1000;
                    
                    % %%%%%%% Imseg algo
                    
                    classif.ma_img.att_data = [classif.ma_img.att_data; min_data];
                    
                    imkmeans_class = [imkmeans_class; mini_class + max_tmp];
                    
                    max_tmp = max(imkmeans_class(:));
                end
                                
                % Reshape the kmeans_class to a vector
                [N_line,N_samp] = size(imkmeans_class);
                imkmeans_class = reshape(imkmeans_class,[N_samp*N_line,1]);
                imkmeans_class = imkmeans_class(~isnan(imkmeans_class));
                
                % Reconstruction of the image and Proportion of classes
                classif = compute_output_img_prop_class(classif,imkmeans_class);

                % Other properties
                classif.nb_cluster = nb_cluster;
            else
                fprintf('Not enough classes, skipping the unsupervized classification.\n')
            end
        end
        
        function classif = Classifier_UnSuprv_HierC(obj,metric,max_cluster_cutoff)
            %
            classif = obj;
            input_attrib = classif.ma_img.att_data;
            
            if size(input_attrib,1)<100000
                link = 'average';
                savemem = 'off';
            else
                metric = 'euclidean';
                link = 'centroid';
                savemem = 'on';
                disp('The option savememory is activated, the metric is now set to "euclidean", and the linkage "centroid".')
            end
            
            % Apply the hierchichal classification
            if max_cluster_cutoff < 1
                HierC_class = clusterdata(input_attrib,'Distance', metric,...
                                                       'Linkage', link,...
                                                       'Cutoff', max_cluster_cutoff,...
                                                       'Criterion', 'distance',...
                                                       'SaveMemory',savemem);
                classif.cutoff = max_cluster_cutoff;
            else
                HierC_class = clusterdata(input_attrib,'Distance', metric,...
                                                       'Linkage', link,...
                                                       'MaxClust', max_cluster_cutoff,...
                                                       'Criterion', 'distance',...
                                                       'SaveMemory',savemem);
                classif.nb_cluster = max_cluster_cutoff;
            end
            
            %dendrogram(HierC_class)
            
            % Building the output attributes
            [~,N_attrib] = size(input_attrib);
            class_names = unique(HierC_class);
            class_names = class_names(~isnan(class_names));
            N_class = length(class_names);
            
            classif.ma_img.att_data = zeros(N_class,N_attrib);
            for i = 1:N_class
                classif.ma_img.att_data(i, :) = mean(input_attrib(HierC_class == class_names(i),:),1);
            end
            
            % Reconstruction of the image and Proportion of classes
            classif = compute_output_img_prop_class(classif, HierC_class);
            
            % Other properties
            classif.metric = metric;
        end
        
        function classif = Classifier_Suprv_metric(obj,metric,library)
            % Application of the classification based on library
            %
            % Parameters
            % ----------
            % data : tensor dim3
            %       hyperspectra of the sample
            % HS_lib : matrix
            %       hyperspectra of the library
            % spectrum_type: int
            %       type of the spectrum to use when applying the unmixing
            % metric: string
            %       type of metric to compare the spectra
            %
            classif = obj;
            input_attrib = classif.ma_img.att_data;
            
            % Initialisations
            lib_spectra = library.att_data;
            submineral_lib = library.entry_names;  
            N_class_in = size(input_attrib,1);
            N_lib_entry = size(lib_spectra,1);
            subminerals_id = zeros(1,N_class_in);
            error_mod = zeros(N_class_in,1);
            entry_names_mod = cellstr(num2str(zeros(1,N_class_in)'));
            
            % Comparison
            for i = 1:N_class_in
                % Comparison of a given spectrum to every spectrum of
                % the library
                score = zeros(1,N_lib_entry);
                if strcmp(metric,'correlation')
                    for j = 1:N_lib_entry
                        score(j) = 1 - corr(lib_spectra(j,:)',input_attrib(i,:)');
                    end
                elseif strcmp(metric,'sqeuclidean')
                    for j = 1:N_lib_entry
                        score(j) = sqrt(sum((lib_spectra(j,:) - input_attrib(i,:)).^2));
                    end
                else
                    error('Metric unknown...')
                end

                % Results
                [error_mod(i),subminerals_id(i)] = min(score);
                entry_names_mod(i) = submineral_lib(subminerals_id(i));
            end
            
            % Output parameters
            classif.entry_names = entry_names_mod;
            classif.error = error_mod;
            
            % Parameters
            classif.metric = metric;
        end
        
        function classif = Classifier_Suprv_Unmixing(obj,library,method,prop_glob,image_idx)
            % Application of the unmixing
            %
            % Parameters
            % ----------
            % library: :obj:
            %       multi_att_Lin containing the library
            % method : int
            %       method chosen to unmix
            %           1 : FCLS, 2 : SUnSAL, 3 : CLSUnSAL
            % prop_glob: boolean
            %       option to use the global proportion from the unmixing
            %       process instead of the one from the class image
            % image_idx : array of array
            %       for each cluster, it gives the indexes of the pixels
            %       inside
            %
            % Output
            % ------
            % 
            classif = obj;
            input_attrib = classif.ma_img.att_data;
                        
            % Initialisations
            lib_spectra = library.att_data;
            submineral_lib = library.entry_names;  
            N_class_in = size(input_attrib,1);
            N_lib_entry = size(lib_spectra,1);
            entry_names_mod = cellstr(num2str(zeros(1,N_class_in)'));
            X_hat = zeros(N_lib_entry, N_class_in);
            
            if ~exist('image_idx')
                for i = 1:N_class_in
                    image_idx{i} = [i];
                end
            end
            min_number = 0.01;
            nb_idx = length(image_idx);
            
            for i = 1:nb_idx
                % Construct the appropriate data source
                input_attrib_idx = input_attrib(image_idx{i},:)';
                input_attrib_idx(isnan(input_attrib_idx)) = 0;
                input_attrib_idx(input_attrib_idx<=0) = min_number;

                % Running unmixing methods
                if (method==1)           % FCLS method
                    if (i == 1)
                        fprintf('FCLS unmixing, is running ...\n please wait ...\n')
                    end
                    [temp] = sunsal(lib_spectra',input_attrib_idx,'ADDONE','yes', 'POSITIVITY', 'yes');

                elseif(method==2)      % SUnSAL method
                    if (i == 1)
                        fprintf('SUnSAL unmixing, is running ...\n please wait ...\n')
                    end
                    [temp] = sunsal(lib_spectra',input_attrib_idx,'POSITIVITY','yes','VERBOSE','no','ADDONE','no', ...
                                    'lambda', 40e-4,'AL_ITERS',2000, 'TOL', 1e-8);        

                elseif(method==3)    % CLSUnSAL method
                    if (i==1)
                        fprintf('CLSUnSAL unmixing, is running ...\n please wait ...\n')
                    end
                    [temp] = clsunsal(lib_spectra',input_attrib_idx,'POSITIVITY','yes','VERBOSE','no','ADDONE','yes', ...
                                    'lambda', 40e-4,'AL_ITERS',2000, 'TOL', 1e-8);
                else
                    error('Chosen unmixing method not valid...');
                end
                X_hat(:,image_idx{i}) = temp;  
            end
            
            % Replace the proportion exceeding 1 by or below an epsilon by 0
            X_hat(X_hat < 0.02) = 0 ;
            X_hat(X_hat > 1) = 1 ;
            
            % Delete weak proportions and renormalise
            sum_prop = squeeze(sum(X_hat,1));
            X_hat = bsxfun(@rdivide, X_hat, sum_prop);

            % Compute the rmse between the reconstruction and reflectance
            res = lib_spectra'*X_hat;
            error_mod = sqrt(mean((res' - input_attrib).^2, 2));
            
            % Compute Proportion of spectrum
            global_prop = squeeze(mean(X_hat, 2)); % Global proportions of classes
            
            % Make mineral map from proportions
            [~, subminerals_id] = max(X_hat,[],1);  % Main class for each input
            for i = 1:N_class_in
                entry_names_mod(i) = submineral_lib(subminerals_id(i));
            end
            
            % Option of the global prop
            %classif.prop_class
            %global_prop
            if false
                classif.prop_class = global_prop;
            end
            
            % Parameters
            classif.method = method;
            classif.mineral_maps = X_hat;
            
            % Output parameters
            classif.entry_names = entry_names_mod;
            classif.error = error_mod;
        end
        
        function classif = Classifier_Suprv_Tree(obj,library,method)
            % Application of the unmixing
            %
            % Parameters
            % ----------
            % library: :obj:
            %       multi_att_Lin containing the library
            % method : int, optionnal 
            %       method chosen
            %           1: Decision Tree, 2: Random Forest (default)
            %
            % Output
            % ------
            % 
            global cfg
            
            classif = obj;
            
            if nargin <=2
                method = 2;
            end
            
            % Get the minerals (labels)
            tmp_entry_names = split(library.entry_names{1},'_');
            y = cellstr(tmp_entry_names{1});
            for i = 2:length(library.entry_names)
                tmp_entry_names = split(library.entry_names{i},'_');
                y = [y; tmp_entry_names{1}];
            end
            
            if method == 1
                files = [[cfg.dir_path_library, '\\', library.file_name,'_feats_ML_library_DT.mat'],
                         [cfg.dir_path_library, '\\', library.file_name,'_model_ML_library_DT.mat']];
                try
                    % Tries to load a saved model
                    [feats,model] = search_files_ML(files,library,y);
                catch
                    % Training the model if cannot load
                    [feats,model] = train_ML(library,y,method,files);
                end
            elseif method == 2
                files = [[cfg.dir_path_library, '\\', library.file_name,'_feats_ML_library_RF.mat'],
                         [cfg.dir_path_library, '\\', library.file_name,'_model_ML_library_RF.mat']];
                try
                    % Tries to load a saved model
                    [feats,model] = search_files_ML(files,library,y);
                catch
                    % Training the model if cannot load
                    [feats,model] = train_ML(library,y,method,files);
                end
            end
            
            % Modification of the image
            classif.ma_img = to_subs_ref(classif.ma_img,'feats',feats);
            
            % Application of the model and Parameters setting
            [classif.entry_names,err_mat] = predict(model,classif.ma_img.att_data);
            
            % Takes the right error from the error matrix
            classif.error = NaN(length(classif.entry_names),1);
            for i = 1:length(classif.entry_names)
                ind = find(ismember(unique(y),classif.entry_names(i)));
                classif.error(i) = err_mat(i,ind);
            end
        end
        
        function classif = Classifier_Suprv_Tree_custom(obj)
            %
            %
            
            classif = obj;
            img = classif.ma_img;
            input_attrib = classif.ma_img.att_data;
            
            % Initialisations
            N_class_in = size(input_attrib,1);
            entry_names_mod = cellstr(num2str(zeros(1,N_class_in)'));
            
            % Global features
            
            % 1390/1415 absorption
                % weak signal for low reflectance
                abs_1390_low = tmp_img_multi(img,[1350,1390,1430],[0.03,-0.03],[],[1,3]);
                
                % Narrow absorption
                abs_1390_narr = tmp_img_multi(img,[1365,1390,1415],[0.15,-0.15],[],[1,3]);
                abs_1390_low = abs_1390_low & ~abs_1390_narr;
                
                % 1415 Asymmetric 
                abs_1415_asym = tmp_img_multi(img,[1350,1415,1520],[0.1,-0.1],[],1);
                abs_1390_low = abs_1390_low & ~abs_1415_asym;
                
                % 1415 Asymmetric deep
                abs_1415_asym_deep = tmp_img_multi(img,[1350,1415,1520],[0.25,-0.25],[],1);
                
                % 1415 abs
                abs_1415 = tmp_img_multi(img,[1390,1415,1520],[0.15,-0.08],[],1);
                
                % Double absorption
                abs_d1390 = tmp_img_multi(img,[1350,1400,1420,1460,1610],[0.08,-0.01,0.01,-0.03],[],[1,5]);
                abs_1390_low = abs_1390_low & ~abs_d1390;
                abs_1390_narr = abs_1390_narr & ~abs_d1390;               

                
            % Application of criteria
            % Pyroxene
            abs_pyr_1 = tmp_img_multi(img,[810,905,990],[0.05,-0.03],[],[1,3]);
            abs_pyr_2 = tmp_img_multi(img,[1600,1850,2180],[0.04,-0.04],[],[1,3]);
            
            abs_pyr = abs_pyr_1 & abs_pyr_2 & ~abs_d1390;
            
            % Olivine
            abs_broad_1000_1_st = tmp_img_multi(img,[470,550],-0.1,[],2);
            abs_broad_1000_2_st = tmp_img_multi(img,[650,970],0.15,[],1);
            abs_broad_1000_3_st = tmp_img_multi(img,[1100,1350,1600],[-0.05,-0.15],[],3);
            
            abs_broad_1000_st = abs_broad_1000_1_st & abs_broad_1000_2_st & abs_broad_1000_3_st;
            
            abs_oli = abs_broad_1000_st & ~abs_d1390 & ~abs_pyr & ~abs_1415_asym_deep;
            
            % Serpentine
            abs_serp_1 = abs_1390_low | abs_1390_narr;
            abs_serp_2 = tmp_img_multi(img,[2290,2330,2360],[0.05,-0.05],[],[1,3]);
            
            abs_serp = abs_serp_1 & abs_serp_2;
            
            % Kaolinite
            abs_kaol_1 = tmp_img_multi(img,[2130,2205,2270],[0.15,-0.02],[],1);
            
            abs_kaol = abs_kaol_1 & abs_1415;
            
            % Smectite
            abs_smec = (abs_1415_asym | abs_1415_asym_deep) & ~abs_pyr & ~abs_kaol;
            
            % Iron oxides
                abs_iro_ox = tmp_img_multi(img,[450,490,600],[0.07,-0.07],[],[1,3]) & tmp_img_multi(img,[450,600],-0.01); 
                sep_hem_goe = tmp_img_multi(img,[780,886],0.03,[],1) & tmp_img_multi(img,[1060,1370],-0.25,[],2);
                %abs_pyr = abs_pyr & ~sep_hem_goe;
                
                % Goethite
                abs_goe = abs_iro_ox & (tmp_img_multi(img,[595,775],-0.45,1,2) | ~sep_hem_goe) & ~abs_d1390;
                
                % Hematite
                abs_hem = abs_iro_ox & (tmp_img_multi(img,[595,775],-0.35,0,2) & sep_hem_goe) & ~abs_d1390;
            
            % Ni
            abs_ni_1 = tmp_img_multi(img,[570,660,840],[0.01,-0.02],[],3);
            
            abs_ni = abs_ni_1;
            
            % LDH
            abs_ldh_1 = tmp_img_wl(img,1750) < 0.95;
            
            abs_ldh = abs_ldh_1 & abs_d1390;
            
            % Final
            min_name = [{'Pyroxene'},{'Olivine'},{'Serpentine'},{'Smectite'},{'LDH'},{'Kaolinite'},{'Goethite'},{'Hematite'},{'Ni'}];
            ind_names = 10^0*abs_pyr +...
                        10^1*abs_oli +...
                        10^2*abs_serp +...
                        10^3*abs_smec +...
                        10^4*abs_ldh +...
                        10^5*abs_kaol +...
                        10^6*abs_goe +...
                        10^7*abs_hem +...
                        10^8*abs_ni;
            ind_names_uni = unique(ind_names);
            
            % Naming the different combinations (10101 -> Pyroxene-Serpentine-Goethite
            entry_names_uni = [];
            for i = 1:length(ind_names_uni)
                name_tmp = '';
                for j = 1:length(min_name)
                    test = ind_names_uni(i) - 10^(j-1);  % test if one mineral is present
                    if ~contains(num2str(test),'9') && test>=0
                        name_tmp = [name_tmp, '-', min_name{j}];  % Add mineral
                    end
                end
                if isempty(name_tmp)
                    entry_names_uni = [entry_names_uni, {'X'}]; % Add if empty
                else
                    entry_names_uni = [entry_names_uni, {name_tmp(2:end)}];  % Add combination of minerals
                end
            end
            
            % Creation of the entry_names based on the ind_names
            cell_ind_names = cellstr(num2str(ind_names));
            cell_ind_names_uni = unique(cell_ind_names);
            for i = 1:length(cell_ind_names_uni)
                entry_names_mod(strcmp(cell_ind_names,cell_ind_names_uni(i))) = entry_names_uni(i);
            end
            %[classif.ma_img.entry_names,entry_names_mod]
            % Parameters setting
            classif.entry_names = entry_names_mod;
            classif.error = zeros(size(input_attrib,1),1);
        end
        
        function classif = Classifier_Suprv_submin2min(obj,prop_class_old)
            % Classification based on the mineral name. It is only based on the
            % submineral names (mineral + ID).
            %
            % Parameters
            % ----------
            % classif: :obj:
            %       classifier to be transform from submineral to minerals
            % prop_class_old: 1 x nb_input_class
            %       proportions of the previous classifier (it has to be
            %       consistent with the error of 'classif')
            %
            classif = obj;
            input_attrib = classif.ma_img.att_data;
            
            if sum(contains(classif.entry_names,'_'))>0
                sepInd = strfind(classif.entry_names, '_');
                entry_names_mod = cellfun(@(x,y) x(1:y-1), classif.entry_names, sepInd, 'un', 0);
            else
                entry_names_mod = classif.entry_names;
            end
            classif.entry_names = unique(entry_names_mod);
            class_res = zeros(1,length(entry_names_mod));
            error_mod = zeros(length(classif.entry_names),1);
            N_attrib = size(classif.ma_img.att_data,2);
            
            % New mineral class
            for i = 1:length(entry_names_mod)
                class_res(i) = find(strcmp(classif.entry_names,entry_names_mod(i)));
            end
            
            % Error
            for i = 1:length(classif.entry_names)
                idx = strcmp(entry_names_mod,classif.entry_names(i));
                error_mod(i) = sum(classif.error(idx).*prop_class_old(idx))./sum(prop_class_old(idx));
            end
            
            % Building the output attributes
            class_names = unique(class_res);
            N_class = length(class_names);
            
            classif.ma_img.att_data = zeros(N_class,N_attrib);
            for i = 1:N_class
                classif.ma_img.att_data(i, :) = mean(input_attrib(class_res == class_names(i),:),1);
            end
            
            % Reconstruction of the image and Proportion of classes
            classif = compute_output_img_prop_class(classif,class_res,1,1);           
        end

        function classif = compute_output_img_prop_class(obj,classes,error_opt,glob_opt)
            % Build the image based on an input image and a vector of new classes
            %
            % Parameters
            % ----------
            % obj: :obj:
            %       classifier
            % classes: vector 1xc
            %       cn new classes, the index of this vector match with the old
            %       classes
            % error_opt: boolean
            %       option to compute the error of each new class
            % glob_opt: boolean
            %       option to compute the proportion and the error globally from
            %       the input and not reconstructed image, this option may lead to
            %       inconsistencies between the image displayed and the proportion
            %       / errors 
            %
            % Output
            % ------
            % classif: :obj:
            %       classifier with updated ind_img, proportion and error
            %
            classif = obj;

            if nargin <= 2
                error_opt = 0;
            end
            if nargin <= 3
                glob_opt = 0;
            end
            
            % Reconstruction of the image
            classif.ma_img = update_classes(classif.ma_img,classes);
            
            % Case of NaN values
            classif.ma_img.ind_img(classif.ma_img.ind_img<0) = NaN;
            
            % Get the number of unique classes
            uni_cla = unique(classes);
            uni_cla(isnan(uni_cla)) = [];
            nb_classes = length(uni_cla);
            
            % Proportion of classes
            classif.prop_class = zeros(nb_classes,1);
            for i = 1:nb_classes
                if glob_opt
                    classif.prop_class(i) = nansum(obj.prop_class(classes==i));
                else
                    classif.prop_class(i) = nansum(nansum(obj.ma_img.surf(classif.ma_img.ind_img==i)))/nansum(obj.ma_img.surf(:));
                end
            end
            
            if error_opt
                % Error of the classes
                classif.error = zeros(nb_classes,1);
                for i = 1:nb_classes
                    if glob_opt
                        classif.error(i) = nansum(obj.error(classes==i).*obj.prop_class(classes==i))/nansum(obj.prop_class(classes==i));
                    else
                        classif.error(i) = mean(classif.error_img(classif.ma_img.ind_img==i));
                    end
                end
            end
        end
        
        function [table_out,img_rgb] = build_table_img_res(obj,color_type)
            %
            %
            % obj: :obj:
            %       
            % color_type: string
            %       type of color scale
            %
            
            % Colors
            if strcmp(color_type,'cont')
                % Display clusters
                nb = length(obj.error);
                
                % Add colors to the clusters, the colors are ranked with regards to the
                % error values
                table_temp = table(obj.error,linspace(1,nb,nb)');
                table_sorted_err = sortrows(table_temp,1);    % sort by error
                table_sorted_err.colors = parula(nb);         % set colors
                table_temp = sortrows(table_sorted_err,2);    % sort back to original order
                colors_classes = table_temp{:,'colors'};
            elseif strcmp(color_type,'min')
                colors_classes = color_classes(obj.entry_names);
            end
            
            % Output table
            list = obj.entry_names;
            errors = obj.error;
            prop = obj.prop_class;
            
            % Color of the table
            colergen = @(color) ['<html><table border=0 width=400 bgcolor=',color,'><TR><TD> </TD></TR></table>'];
            for c = 1:length(list)
                color_column(c,1) = {colergen(rgb2hex(colors_classes(c,:)))};
            end
            
            % Table
            table_out = table(list,...
                              cellstr(num2str(round(errors,2))),...
                              color_column,...
                              cellstr(num2str(round(prop*100,2))),...
                              'VariableNames',[{'Minerals'},{'Error'},{'Color'},{'Prop'}]);
            
            % Removing the line 'X' 
            table_out(strcmp(table_out{:, "Minerals"}, {'X'}),:) = [];
            
            if ~strcmp(color_type,'cont')
                % Removing low proportions 
                table_out((str2num(cell2mat(table_out{:, "Prop"})) < 0.1),:) = [];
                
                % Normalize the proportions 
                table_out{:, "Prop"} = cellstr(num2str(round(str2num(cell2mat(table_out{:, "Prop"}))/sum(str2num(cell2mat(table_out{:, "Prop"})))*100,2)));
            end
            
            % Apply color to map
            img_rgb = apply_color2map(obj.ma_img.ind_img,colors_classes);
        end
    end
end

% Functions for naming
function y = soft(x,T)

y = max(abs(x) - T, 0);
y = y./(y+T) .* x;
end

function Y = vector_soft_row(X,tau)
%
%  computes the vector soft columnwise

NU = sqrt(sum(X.^2,2));
A = max(0, NU-tau);
Y = repmat((A./(A+tau)),1,size(X,2)).* X;
end

function [x,res_p,res_d] = sunsal(M,y,varargin)

% [x] = sunsal_v2(M,y,varargin)
%
%  SUNSAL -> sparse unmixing via variable splitting and augmented
%  Lagrangian methods 
%
% --------------- Description --------------------------------------------
%
%  SUNSAL solves the following l2-l1 optimization  problem 
%  [size(M) = (L,p); size(X) = (p,N)]; size(Y) = (L,N)]
%
%         min  (1/2) ||M X-y||^2_F + lambda ||X||_1
%          X              
%
%  where ||X||_1 = sum(sum(abs(X)).
% 
%    CONSTRAINTS ACCEPTED:
%
%    1) POSITIVITY:  X >= 0;
%    2) ADDONE:  sum(X) = ones(1,N);
%
%    NOTES: 
%       1) The optimization w.r.t each column of X is decoupled. Thus, 
%          SUNSAL solves N simultaneous problems.
%
%       2) SUNSAL solves the following  problems:
%  
%          a) BPDN - Basis pursuit denoising l2-l1 
%                    (lambda > 0, POSITIVITY = 'no', ADDONE, 'no')
%
%          b) CBPDN - Constrained basis pursuit denoising l2-l1 
%                    (lambda > 0, POSITIVITY = 'yes', ADDONE, 'no')
%      
%          c) CLS   - Constrained least squares
%                     (lambda = 0, POSITIVITY = 'yes', ADDONE, 'no')
%
%          c) FCLS   - Fully constrained least squares
%                     (lambda >=0 , POSITIVITY = 'yes', ADDONE, 'yes')
%                      In this case, the regularizer ||X||_1  plays no role, 
%                      as it is constant.
%          
%
% -------------------- Line of Attack  -----------------------------------
%
%  SUNSAL solves the above optimization problem by introducing a variable
%  splitting and then solving the resulting constrained optimization with
%  the augmented Lagrangian method of multipliers (ADMM). 
% 
% 
%         min  (1/2) ||M X-y||^2_F + lambda ||Z||_1
%          X,Z              
%         subject to: sum(X) = ones(1,N)); Z >= 0; X = Z
%
%  Augmented Lagrangian (scaled version):
%
%       L(X,Z,D) = (1/2) ||M X-y||^2_F + lambda ||Z||_1 + mu/2||X-Z-D||^2_F
%       
%  where D are the scale Lagrange multipliers
%
%
%  ADMM:
%
%      do 
%        X  <-- arg min L(X,Z,D)
%                    X, s.t: sum(X) = ones(1,N));
%        Z  <-- arg min L(X,Z,D)
%                    Z, s.t: Z >= 0;
%        D  <-- D - (X-Z);
%      while ~stop_rulde
%  
%For details see
%
%
% [1] J. Bioucas-Dias and M. Figueiredo, “Alternating direction algorithms
% for constrained sparse regression: Application to hyperspectral unmixing”, 
% in 2nd  IEEE GRSS Workshop on Hyperspectral Image and Signal 
% Processing-WHISPERS'2010, Raykjavik, Iceland, 2010. 
%
%
% ------------------------------------------------------------------------
%  ===== Required inputs =============
%
%  M - [L(channels) x p(endmembers)] mixing matrix
%
%  y - matrix with  L(channels) x N(pixels).
%      each pixel is a linear mixture of p endmembers
%      signatures y = M*x + noise,
%
%      
%
%
%  ====================== Optional inputs =============================
%
%  'AL_ITERS' - Minimum number of augmented Lagrangian iterations
%               Default: 100;
%               
%  lambda - regularization parameter. lambda is either a scalar
%           or a vector with N components (one per column of x)
%           Default: 0. 
%
%
%  'POSITIVITY'  = {'yes', 'no'}; Enforces the positivity constraint: 
%                   X >= 0
%                   Default 'no'
%
%  'ADDONE'  = {'yes', 'no'}; Enforces the sum to one constraint: sum(X) = ones(1,N)
%              Default 'no'
% 
%   'TOL'    - tolerance for the primal and  dual residuals 
%              Default = 1e-4; 
%
%
%  'verbose'   = {'yes', 'no'}; 
%                 'no' - work silently
%                 'yes' - display warnings
%                  Default 'no'
%        
%  =========================== Outputs ==================================
%
% X  =  [pxN] estimated mixing matrix
%
%

%
% -------------------------------------------------------------------------
%
% Copyright (July, 2009):        José Bioucas-Dias (bioucas@lx.it.pt)
%
% SUNSAL is distributed under the terms of
% the GNU General Public License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ---------------------------------------------------------------------



%
%--------------------------------------------------------------
% test for number of required parametres
%--------------------------------------------------------------
if (nargin-length(varargin)) ~= 2
    error('Wrong number of required parameters');
end
% mixing matrixsize
[LM,p] = size(M);
% data set size
[L,N] = size(y);
if (LM ~= L)
    error('mixing matrix M and data set y are inconsistent');
end
% if (L<p)
%     error('Insufficient number of columns in y');
% end


%
%--------------------------------------------------------------
% Set the defaults for the optional parameters
%--------------------------------------------------------------
% maximum number of AL iteration
AL_iters = 2000;
% regularizatio parameter
lambda = 0.0;
% display only sunsal warnings
verbose = 'off';
% Positivity constraint
positivity = 'no';
% Sum-to-one constraint
addone = 'no';
% tolerance for the primal and dual residues
tol = 1e-4;
% initialization
x0 = 0;

%
%--------------------------------------------------------------
% Local variables
%--------------------------------------------------------------

%--------------------------------------------------------------
% Read the optional parameters
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'AL_ITERS'
                AL_iters = round(varargin{i+1});
                if (AL_iters <= 0 )
                       error('AL_iters must a positive integer');
                end
            case 'LAMBDA'
                lambda = varargin{i+1};
                if (sum(sum(lambda < 0)) >  0 )
                       error('lambda must be positive');
                end
            case 'POSITIVITY'
                positivity = varargin{i+1};
            case 'ADDONE'
                addone = varargin{i+1};
            case 'TOL'
                tol = varargin{i+1};
            case 'VERBOSE'
                verbose = varargin{i+1};
            case 'X0'
                x0 = varargin{i+1};
                if (size(x0,1) ~= p) || (size(x0,1) ~= N)
                    error('initial X is  inconsistent with M or Y');
                end
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

%---------------------------------------------
%  If lambda is scalar convert it into vector
%---------------------------------------------
Nlambda = size(lambda);
if Nlambda == 1
    % same lambda for all pixels
    lambda = lambda*ones(p,N);
elseif Nlambda ~= N
        error('Lambda size is inconsistent with the size of the data set');
else
  %each pixel has its own lambda
   lambda = repmat(lambda(:)',p,1);
end

% compute mean norm
norm_y = sqrt(mean(mean(y.^2)));
% rescale M and Y and lambda
M = M/norm_y;
y = y/norm_y;
lambda = lambda/norm_y^2;

%
%---------------------------------------------
% just least squares
%---------------------------------------------
if sum(sum(lambda == 0)) &&  strcmp(positivity,'no') && strcmp(addone,'no')
    fprintf('Just least square \n')
    x = pinv(M)*y;
    % primal and dual residues
    res_p = 0;
    res_d = 0;
    return
end
%---------------------------------------------
% constrained least squares (sum(x) = 1)
%---------------------------------------------
SMALL = 1e-12;
B = ones(1,p);
a = ones(1,N);

if  strcmp(addone,'yes') && strcmp(positivity,'no') 
    
    F = M'*M;
    % test if F is invertible
    if rcond(F) > SMALL
        % compute the solution explicitly
        IF = inv(F);
        x = IF*M'*y-IF*B'*inv(B*IF*B')*(B*IF*M'*y-a);
        fprintf('Contrained least square with F invertible \n')
        % primal and dual residues
        res_p = 0;
        res_d = 0;
        return
    end
end


%
%---------------------------------------------
%  Constants and initializations
%---------------------------------------------
mu_AL = 0.01;
mu = 10*mean(lambda(:)) + mu_AL;

%F = M'*M+mu*eye(p);
[UF,SF] = svd(M'*M);
sF = diag(SF);
IF = UF*diag(1./(sF+mu))*UF';
%IF = inv(F);
Aux = IF*B'*inv(B*IF*B');
x_aux = Aux*a;
IF1 = (IF-Aux*B*IF);

yy = M'*y;

%
%---------------------------------------------
%  Initializations
%---------------------------------------------

% no intial solution supplied
if x0 == 0
    x= IF*M'*y;
end

z = x;
% scaled Lagrange Multipliers
d  = 0*z;


%
%---------------------------------------------
%  AL iterations - main body
%---------------------------------------------
tol1 = sqrt(N*p)*tol;
tol2 = sqrt(N*p)*tol;
i=1;
res_p = inf;
res_d = inf;
%maskz = ones(size(z));
mu_changed = 0;
while (i <= AL_iters) && ((abs (res_p) > tol1) || (abs (res_d) > tol2)) % Figs. 2, 3
    % save z to be used later
    if mod(i,10) == 1
        z0 = z;
    end
    % minimize with respect to z (u in the paper)
    z =  soft(x-d,lambda/mu);
    % teste for positivity
    if strcmp(positivity,'yes')
       maskz = (z >= 0);
       z = z.*maskz; 
    end
    % teste for sum-to-one 
    if strcmp(addone,'yes')
       x = IF1*(yy+mu*(z+d))+x_aux;
    else
       x = IF*(yy+mu*(z+d));
    end
    
    % Lagrange multipliers update
    d = d -(x-z);

    % update mu so to keep primal and dual residuals within a factor of 10
    if mod(i,10) == 1
        % primal residue
        res_p = norm(x-z,'fro');
        % dual residue
        res_d = mu*norm(z-z0,'fro');
        if  strcmp(verbose,'yes')
            fprintf(' i = %f, res_p = %f, res_d = %f\n',i,res_p,res_d)
        end
        % update mu
        if res_p > 10*res_d
            mu = mu*2;
            d = d/2;
            mu_changed = 1;
        elseif res_d > 10*res_p
            mu = mu/2;
            d = d*2;
            mu_changed = 1;
        end
        if  mu_changed
           % update IF and IF1
           IF = UF*diag(1./(sF+mu))*UF';
           Aux = IF*B'*inv(B*IF*B');
           x_aux = Aux*a;
           IF1 = (IF-Aux*B*IF);
           mu_changed = 0;
            %mu
        end 
    end
    i=i+1;            
end
end

function [x,res_p,res_d] = clsunsal(M,y,varargin)

% [x] = clsunsal(M,y,varargin)
%
%  CLSUNSAL -> collaborative sparse unmixing via variable splitting and augmented
%  Lagrangian  
%
% --------------- Description --------------------------------------------
%
%  CLSUNSAL solves the following l2-l1 optimization  problem 
%  [size(M) = (L,p); size(X) = (p,N)]; size(Y) = (L,N)]
%
%         min  (1/2) ||M X-Y||^2_F + lambda ||X||_{2,1}
%          X              
%
%  where ||X||_{2,1} = sum(sqrt(sum(X.^2),2))
% 
%    CONSTRAINTS ACCEPTED:
%
%    1) POSITIVITY:  X >= 0;
%    2) ADDONE:  sum(X) = ones(1,N);
%
%          
%
% -------------------- Line of Attack  -----------------------------------
%
%  CLSUNSAL solves the above optimization problem by introducing a variable
%  splitting and then solving the resulting constrained optimization with
%  the augmented Lagrangian method of multipliers (ADMM). 
% 
% 
%         min  (1/2) ||M X-Y||^2_F + lambda ||V1||_{2,1} + i_+(V2)
%          X,Z              
%         subject to: sum(X) = ones(1,N)); V1 = X,  V2 = X;
%
%  Augmented Lagrangian (scaled version):
%
%       L(X,V1,V2,D1,D2) = (1/2) ||M X-Y||^2_F 
%                      + lambda ||V1||_{2,1} + mu/2||X-V1-D1||^2_F
%                      + i_+(V2) + mu/2||X-V2-D2||^2_F
%       
%  where D1 and D2 are the scaled Lagrange multipliers
%
%
%  ADMM:
%
%      do 
%        X  <-- arg min L(X,V1,V2,D1,D2)
%                    X,   s.t: sum(X) = ones(1,N));
%        V1  <-- arg min L(X,V1,V2,D1,D2)
%                     V1
%        V2  <-- arg min L(X,V1,V2,D1,D2)
%                     V2
%        D1  <-- D1 - (X-V1);
%        D2  <-- D2 - (X-V2);
%
%      while ~stop_rule
%  
% For details see
%
% M. D. Iordache, J. M. Bioucas-Dias and A. Plaza. "Collaborative Sparse 
% Regression for Hyperspectral Unmixing",  IEEE Transactions on Geoscience 
% and Remote Sensing, 2013  (accepted).
% http://dx.doi.org/10.1109/TGRS.2013.2240001 
%
% NOTE: this version  differs from that presented in the paper in two
% points: 
%     
%    a) the sum-to-one constraint is allowed
%
%    b) use just two spliting instead of the three proposed in the paper.
%       
% The version clsunsal_v1 implemts exactly the splitting scheme shown in
% the paper
%
% ------------------------------------------------------------------------
%  ===== Required inputs =============
%
%  M - [L(channels) x p(endmembers)] mixing matrix
%
%  y - matrix with  L(channels) x N(pixels).
%      each pixel is a linear mixture of p endmembers
%      signatures y = M*x + noise,
%
%      
%
%
%  ====================== Optional inputs =============================
%
%  'AL_ITERS' - Minimum number of augmented Lagrangian iterations
%               Default: 100;
%               
%  lambda - regularization parameter. 
%
%
%  'POSITIVITY'  = {'yes', 'no'}; Enforces the positivity constraint: 
%                   X >= 0
%                   Default 'no'
%
%  'ADDONE'  = {'yes', 'no'}; sum(X) = ones(1,N)
%              Default 'no'
% 
%   'TOL'    - tolerance for the primal and  dual residuals 
%              Default = 1e-4; 
%
%
%  'verbose'   = {'yes', 'no'}; 
%                 'no' - work silently
%                 'yes' - display warnings
%                  Default 'no'
%        
%  =========================== Outputs ==================================
%
% X  =  [pxN] estimated mixing matrix
%
%
% -------------------------------------------------------------------------
%
% Copyright (July, 2012):        José Bioucas-Dias (bioucas@lx.it.pt)
%
% CLSUNSAL is distributed under the terms of
% the GNU General Public License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ---------------------------------------------------------------------

%
%--------------------------------------------------------------
% test for number of required parametres
%--------------------------------------------------------------
if (nargin-length(varargin)) ~= 2
    error('Wrong number of required parameters');
end
% mixing matrixsize
[LM,p] = size(M);
% data set size
[L,N] = size(y);
if (LM ~= L)
    error('mixing matrix M and data set y are inconsistent');
end
% if (L<p)
%     error('Insufficient number of columns in y');
% end


%
%--------------------------------------------------------------
% Set the defaults for the optional parameters
%--------------------------------------------------------------
% maximum number of AL iteration
AL_iters = 1000;
% regularizatio parameter
lambda = 0.0;
% display only sunsal warnings
verbose = 'off';
% Positivity constraint
positivity = 'no';
% Sum-to-one constraint
addone = 'no';
% tolerance for the primal and dual residues
tol = 1e-4;
% initialization
x0 = 0;

%
%--------------------------------------------------------------
% Local variables
%--------------------------------------------------------------


%--------------------------------------------------------------
% Read the optional parameters
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'AL_ITERS'
                AL_iters = round(varargin{i+1});
                if (AL_iters <= 0 )
                       error('AL_iters must a positive integer');
                end
            case 'LAMBDA'
                lambda = varargin{i+1};
                if lambda < 0
                       error('lambda must be positive');
                end
            case 'POSITIVITY'
                positivity = varargin{i+1};
            case 'ADDONE'
                addone = varargin{i+1};
            case 'TOL'
                tol = varargin{i+1};
            case 'VERBOSE'
                verbose = varargin{i+1};
            case 'X0'
                x0 = varargin{i+1};
                if (size(x0,1) ~= p) || (size(x0,1) ~= N)
                    error('initial X is inconsistent with M or Y');
                end
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end


% compute mean norm
norm_y = sqrt(mean(mean(y.^2)));
% rescale M and Y and lambda
M = M/norm_y;
y = y/norm_y;
lambda = lambda/norm_y^2;

%
%---------------------------------------------
% just least squares
%---------------------------------------------
if sum(sum(lambda == 0)) &&  strcmp(positivity,'no') && strcmp(addone,'no')
    z = pinv(M)*y;
    % primal and dual residues
    res_p = 0;
    res_d = 0;
    return
end
%---------------------------------------------
% least squares constrained (sum(x) = 1)
%---------------------------------------------
SMALL = 1e-12;
B = ones(1,p);
a = ones(1,N);

if  strcmp(addone,'yes') && strcmp(positivity,'no') 
    F = M'*M;
    % test if F is invertible
    if rcond(F) > SMALL
        % compute the solution explicitly
        IF = inv(F);
        z = IF*M'*y-IF*B'*inv(B*IF*B')*(B*IF*M'*y-a);
        % primal and dual residues
        res_p = 0;
        res_d = 0;
        return
    end
end


%
%---------------------------------------------
%  Constants and initializations
%---------------------------------------------
mu_AL = 0.01;
mu = 10*mean(lambda(:)) + mu_AL;

%F = M'*M+mu*eye(p);
[UF,SF] = svd(M'*M);
sF = diag(SF);
IF = UF*diag(1./(sF+2*mu))*UF';
%IF = inv(F);
Aux = IF*B'*inv(B*IF*B');
x_aux = Aux*a;
IF1 = (IF-Aux*B*IF);


yy = M'*y;

%
%---------------------------------------------
%  Initializations
%---------------------------------------------

% no intial solution supplied
if x0 == 0
    x= IF*M'*y;
end

% auxiliary variables; (errors occur if x0 is provided)
v1 = x;
v2 = x;

% scaled Lagrange Multipliers
d1  = 0*v1;
d2  = 0*v2;

%
%---------------------------------------------
%  AL iterations - main body
%---------------------------------------------
tol1 = sqrt(N*p)*tol;
tol2 = sqrt(N*p)*tol;
i=1;
res_p = inf;
res_d = inf;
mu_changed = 0;
while (i <= AL_iters) && ((abs (res_p) > tol1) || (abs (res_d) > tol2)) 
    % save z to be used later
    if mod(i,10) == 1
        v10 = v1;
        v20 = v2;
    end
    % minimize with respect to v1
    v1 =  vector_soft_row(x-d1,lambda/mu);
    % minimize wrt v2
    v2 = x-d2;
    % teste for positivity
    if strcmp(positivity,'yes')
        v2 = max(0,v2);
    end
    % test for sum-to-one 
    if strcmp(addone,'yes')
       x = IF1*(yy + mu*(v1+d1)+ mu*(v2+d2))+x_aux;
    else
       x = IF*(yy + mu*(v1+d1)+ mu*(v2+d2));
    end
    
    % Lagrange multipliers update
    d1 = d1 -(x-v1);
    d2 = d2 -(x-v2);

    % update mu so to keep primal and dual residuals within a factor of 10
    if mod(i,10) == 1
        % primal residue
        res_p = sqrt(norm(x-v1,'fro')^2 + norm(x-v2,'fro')^2);
        % dual residue
        res_d = mu*norm(v1-v10+v2-v20,'fro');
        if  strcmp(verbose,'yes')
            fprintf(' i = %f, res_p = %f, res_d = %f\n',i,res_p,res_d)
        end
        % update mu
        if res_p > 10*res_d
            mu = mu*2;
            d1 = d1/2;
            d2 = d2/2;
            mu_changed = 1;
        elseif res_d > 10*res_p
            mu = mu/2;
            d1 = d1*2;
            d2 = d2*2;
            mu_changed = 1;
        end
        if  mu_changed
           % update IF and IF1

           IF = UF*diag(1./(sF+2*mu))*UF';
           Aux = IF*B'*inv(B*IF*B');
           x_aux = Aux*a;
           IF1 = (IF-Aux*B*IF);
           mu_changed = 0;
        end 
    end
    i=i+1;    
end
end

% ML tools
function [feats,model] = search_files_ML(files,library,y)
global cfg

% Tries to load the feats
feats_load = load(files(1,:));
feats = feats_load.feats;

% Tries to load a saved model
model_load = load(files(2,:));
model = model_load.model;

% Checks that the model's parameters fit the library
ind_X = length(feats) == size(model.X,2);
ind_y = isequal(model.ClassNames,unique(y));
ind_nb_tree = isequal(model.NumTrees,cfg.nb_tree);
ind_nb_entry = length(library.entry_names) == size(model.X,1);

if ~(ind_X && ind_y && ind_nb_tree && ind_nb_entry)
    error('Saved model does not fit the library')
end

end

function [feats,model] = train_ML(library_in,y,method,files)
global cfg

% Change lib into a ratios library
library = to_subs_ref(library_in,'all',[]);

% Training the model on the library with full features
if method == 1
    model_full = fitctree(library.att_data,y);
elseif method == 2
    model_full = TreeBagger(cfg.nb_tree,library.att_data,y);
end

% Select feats and apply them to the library
feats_ind = get_feats_Trees(model_full,method);
library_red = get_from_indexes(library,[],feats_ind);
feats = library_red.att_names;

% Training the model on the library with reduced features
if method == 1
    model = fitctree(library_red.att_data,y); 
elseif method == 2
    model = TreeBagger(cfg.nb_tree,library_red.att_data,y,'OOBPrediction','on','OOBPredictorImportance','on');
    mean(model.oobError)
end

% Saves the model
save(files(1,:),'feats');
save(files(2,:),'model');
end

function feats = get_feats_Trees(model,method)
%
%

all_ind = 1:size(model.X,2);

if method == 1
    feats = all_ind(model.Trees.predictorImportance > 0);
elseif method == 2
    feats = [];
    for i = 1:model.NumTrees
        feats = cat(2,feats, all_ind(model.Trees{i}.predictorImportance > 0));
        feats = unique(feats);
    end
end

end
%

function output = tmp_img_wl(img_modif_spec,wl)
    if wl < 1000
        xir = {'VNIR'};
    else
        xir = {'SWIR'};
    end

    inds = get_ind_from_gr_id(img_modif_spec,char(xir),wl);

    output = img_modif_spec.att_data(:,inds);
end

function res = tmp_img_multi(img,wl,diff,option,norm)
    % 
    % Parameters
    % ----------
    % img: matrix
    %       image to be modified
    % wl: 1xn array
    %       array of wavelengths
    % ref val: 1x(n-1) array
    %       ref value of absorbance
    % option: 1x(n-1) array
    %       options of comparison, 0: below ref val, 1: above ref val,
    %                              2: around 0, ref val is the tolerance
    %
    
    if nargin < 4
        option = sign(diff) == 1;
    end
    
    if isempty(option)
        option = sign(diff) == 1;
    end
    
    if nargin < 5
        norm = [];
    end
    
    if length(norm) == 1
        norm_uni = tmp_img_wl(img,wl(norm));
    elseif length(norm) == 2
        norm_sta = tmp_img_wl(img,wl(norm(1)));
        norm_end = tmp_img_wl(img,wl(norm(2)));
        a = (norm_end - norm_sta)./(wl(norm(2)) - wl(norm(1)));
    end
    
    res = true([size(img.att_data,1),1]);
    for i = 1:length(diff)
        % Normalizing by a line drawn between extreme wavelengths
        if length(norm) == 1
            img_wl_tmp_1 = tmp_img_wl(img,wl(i))./norm_uni;
            img_wl_tmp_2 = tmp_img_wl(img,wl(i+1))./norm_uni;
        elseif length(norm) == 2
            img_wl_tmp_1 = tmp_img_wl(img,wl(i))./(norm_sta + a.*(wl(i) - wl(1)));
            img_wl_tmp_2 = tmp_img_wl(img,wl(i+1))./(norm_sta + a.*(wl(i+1) - wl(1)));
        else
            img_wl_tmp_1 = tmp_img_wl(img,wl(i));
            img_wl_tmp_2 = tmp_img_wl(img,wl(i+1));
        end
        
        % Application of the options
        if option(i) == 0
            res = res & ((img_wl_tmp_1 - img_wl_tmp_2) < diff(i));
        elseif option(i) == 1
            res = res & ((img_wl_tmp_1 - img_wl_tmp_2) > diff(i));
        elseif option(i) == 2
            tmp = (img_wl_tmp_1 - img_wl_tmp_2);
            res = res & ((tmp > -diff(i)) & (tmp < diff(i)));
        else
            error('unknown option')
        end
    end
end

% Function for colors
function colors = color_classes(classes)
    % Association of minerals to a predefined color
    %
    %
    
    map_name_rgb = containers.Map(...
                {'Pyroxene',...
                 'Pyroxene-Serpentine',...
                 'Pyroxene-Serpentine-Smectite',...
                 'Pyroxene-Serpentine-Goethite',...
                 'Pyroxene-Smectite',...
                 'Pyroxene-Serpentine-Smectite-Goethite',...
                 'Pyroxene-Smectite-Goethite',...
                 'Pyroxene-Goethite',...
                 'Serpentine',...
                 'Serpentine-Smectite',...
                 'Serpentine-Goethite',...
                 'Smectite',...
                 'Serpentine-Smectite-Goethite',...
                 'Smectite-Goethite',...
                 'Goethite',...
                 'Goethite-Hematite',...
                 'Hematite',...
                 'Olivine',...
                 'Olivine-Serpentine',...
                 'Olivine-Serpentine-Smectite',...
                 'Olivine-Serpentine-Goethite',...
                 'Olivine-Smectite',...
                 'Olivine-Serpentine-Smectite-Goethite',...
                 'Olivine-Smectite-Goethite',...
                 'Olivine-Goethite',...
                 'Serpentine-Ni',...
                 'Serpentine-Smectite-Ni',...
                 'Serpentine-Goethite-Ni',...
                 'Smectite-Ni',...
                 'Serpentine-Smectite-Goethite-Ni',...
                 'Smectite-Goethite-Ni',...
                 'Kaolinite',...
                 'Kaolinite-Goethite',...
                 'Kaolinite-Ni',...
                 'Goethite-Ni',...
                 'Goethite-Hematite-Ni',...
                 'Hematite-Ni',...
                 'Olivine-Ni',...
                 'Olivine-Serpentine-Ni',...
                 'Olivine-Serpentine-Smectite-Ni',...
                 'Olivine-Serpentine-Goethite-Ni',...
                 'Olivine-Smectite-Ni',...
                 'Olivine-Serpentine-Smectite-Goethite-Ni',...
                 'Olivine-Smectite-Goethite-Ni',...
                 'Olivine-Goethite-Ni',...
                 'Quartz',...
                 'Quartz-Goethite',...
                 'LDH',...
                 'Ni',...
                 'Kerolite',...
                 'X'
                 },...
                {[150, 200, 200],...
                 [170, 200, 200],...
                 [170, 175, 175],...
                 [170, 162, 162],...
                 [170, 150, 150],...
                 [170, 137, 137],...
                 [170, 125, 125],...
                 [170, 100, 100],...
                 [  0, 200, 200],...
                 [  0, 175, 175],...
                 [  0, 162, 162],...
                 [  0, 150, 150],...
                 [  0, 137, 137],...
                 [  0, 125, 125],...
                 [  0, 100, 100],...
                 [  0,  75,  75],...
                 [  0,  50,  50],...
                 [ 85, 200, 200],...
                 [ 42, 200, 200],...
                 [ 42, 175, 175],...
                 [ 42, 162, 162],...
                 [ 42, 150, 150],...
                 [ 42, 137, 137],...
                 [ 42, 125, 125],...
                 [ 42, 100, 100],...
                 [  0, 250, 175],...
                 [  0, 225, 150],...
                 [  0, 212, 137],...
                 [  0, 200, 125],...
                 [  0, 187, 112],...
                 [  0, 175, 100],...
                 [  0,   0, 220],...
                 [  0,  45, 220],...
                 [ 90,  45, 220],...
                 [  0, 150,  75],...
                 [  0, 125,  50],...
                 [  0, 100,  25],...
                 [ 85, 250, 175],...
                 [ 42, 250, 175],...
                 [ 42, 225, 150],...
                 [ 42, 212, 137],...
                 [ 42, 200, 125],...
                 [ 42, 187, 112],...
                 [ 42, 175, 100],...
                 [ 42, 150,  75],...
                 [170,   0, 192],...
                 [  0,  29, 192],...
                 [ 25, 150, 150],...
                 [ 85, 255, 128],...
                 [  0, 200, 125],...
                 [  0,   0,   0]
                 });
    
    % Initialisation
    colors = NaN(length(classes),3);
    
    % Find the classes that are listed in the map
    ind_map_exist = ismember(classes,map_name_rgb.keys);
    
    % Creation of distinguishable colors to fill non referenced classes
    nb_dist_col = sum(~ind_map_exist);
    dis_col = distinguishable_colors(nb_dist_col);
    
    % Loop over the classes and filling the color matrix
    k = 1;
    for i = 1:length(classes)
        if ind_map_exist(i)
            colors(i,:) = hsv2rgb(map_name_rgb(classes{i})./255);
        else
            colors(i,:) = dis_col(k,:);
            k = k + 1;
        end
    end
end

function img_rgb = apply_color2map(mineral_map,colors)
    % Apply color to map - improve code
    a = mineral_map;
    b = a;
    c = a;
    
    indexes = unique(a);
    indexes(isnan(indexes)) = [];
    nb_col = length(indexes);
    for i = 1:nb_col
        idx = indexes(i);
        a(a==idx) = colors(i,1);
        b(b==idx) = colors(i,2);
        c(c==idx) = colors(i,3);
    end
    img_rgb(:,:,1) = a;
    img_rgb(:,:,2) = b;
    img_rgb(:,:,3) = c;
end
