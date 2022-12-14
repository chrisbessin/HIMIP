classdef Multi_Classifier
    properties
        ma_img
        prop_class
        error
        entry_names
        
        error_lib_rgb
        table_error_lib
        data_error_lib
        
        minerals_ind_rgb
        table_minerals_ind
        
        ROI_error_lib
        ROI_error_lib_rgb
        ROI_minerals_ind
        ROI_minerals_ind_rgb
        
        nb_cluster
    end
    
    methods
        function obj = Multi_Classifier(ma_img,prop_class,error,entry_names,nb_cluster)
            if nargin<1
                ma_img = [];
                prop_class = [];
                error = [];
                entry_names = [];
                nb_cluster = [];
            end
            obj.ma_img = ma_img;
            obj.prop_class = prop_class;
            obj.error = error;
            obj.entry_names = entry_names;
            obj.nb_cluster = nb_cluster;
        end
        
        function method = method(obj,method_name,img,library,att_names,spectrum_types,min_max_wl,min_max_xy)
            %{Calls the method based on its name
            
            %}
            global cfg
            method = obj;
            
            % saves the reconstructed RGB image
            img = create_img(img, cfg.red_default, cfg.green_default, cfg.blue_default, 2, []);
            
            if strcmp(method_name,'unmixing')
                method = method_1(method,img,library,att_names,spectrum_types,min_max_wl,min_max_xy);
            elseif strcmp(method_name,'dist')
                method = method_2(method,img,library,att_names,spectrum_types,min_max_wl,min_max_xy);
            elseif strcmp(method_name,'clustering')
                method = method_3(method,img,library,att_names,spectrum_types,min_max_wl,min_max_xy);
            elseif strcmp(method_name,'lib_tree')
                method = method_4(method,img,library,att_names,spectrum_types,min_max_wl,min_max_xy);
            elseif strcmp(method_name,'custom_tree')
                method = method_5(method,img,library,att_names,spectrum_types,min_max_wl,min_max_xy);
            else
                error('Cannot find the method name...')
            end
            
        end
        
        function method = method_1(obj,img,library,att_names,spectrum_types,min_max_wl,min_max_xy)
            % Method unmixing with library reduction
            %
            %
            global cfg
            method = obj;
            
            % Trims the spectra 
            att_names_C = att_names(strcmp(library.att_types(contains(att_names,library.att_groups)),'C'));
            for a = 1:length(att_names_C)
                min_val.(att_names_C{a}) = min_max_wl(2*a - 1);
                max_val.(att_names_C{a}) = min_max_wl(2*a);
                % Remove the unnecessary wavelengths
                % (Not applied to the img as it will be done in the
                % 'consis_multi_att' function)
                library = rm_min_max_wls(library,att_names_C{a},min_val.(att_names_C{a}),max_val.(att_names_C{a}));
            end
            
            % Crop the image
            img = img_modify(img,1,min_max_xy);
            [library,img] = consis_multi_att(library,img,att_names,1);
            
            % Creation of the classifier
            classif = Unit_Classifier(img);
            
            % Kmeans classification on all data
            clf_metric = 'sqeuclidean'; %'sqeuclidean' 'cityblock', 'cosine', 'correlation', 'hamming'
            nb_replicate = 1;
            max_iter = 5;
            classif = Classifier_UnSuprv_Kmeans(classif,clf_metric,cfg.nb_cluster,nb_replicate,max_iter,20);
            
            % Filling the missing values
            img.att_data = fillmissing(img.att_data,'linear',2);
            
            % Change spectra type
            classif.ma_img = apply_change_spec(classif.ma_img,att_names,spectrum_types,0);
            library = apply_change_spec(library,att_names,spectrum_types,0);
            
            % Keep endmembers through EIA_EIHA
            endmembers = EIA_NFINDR(classif.ma_img.att_data');
            
            % Classification of the endmembers with metric
            % Creation of the object
            nb_endmembers = size(endmembers,2);
            img_endmembers = img;
            img_endmembers.att_data = endmembers';
            img_endmembers.ind_img = linspace(1,nb_endmembers,nb_endmembers);
            classif_endmembers = Unit_Classifier(img_endmembers);
            
            % Classification
            metric = 'sqeuclidean';
            classif_endmembers = Classifier_Suprv_metric(classif_endmembers,metric,library);
            
            % Reduction of the library
            lib_entries_to_del = setdiff(library.entry_names,classif_endmembers.entry_names);
            library_endmembers = add_modif(library,1,lib_entries_to_del);
            
            % Unmixing
            unmix_method = 3;
            classif = Classifier_Suprv_Unmixing(classif,library_endmembers,unmix_method,1);
            [method.table_error_lib,method.error_lib_rgb] = build_table_img_res(classif,'cont');
            method.data_error_lib = classif.ma_img;
            
            % From subminerals to minerals
            classif = Classifier_Suprv_submin2min(classif,classif.prop_class);
            [method.table_minerals_ind,method.minerals_ind_rgb] = build_table_img_res(classif,'min');
            
            % Properties
            method.ma_img = classif.ma_img;
            method.prop_class = classif.prop_class;
            method.error = classif.error;
            method.entry_names = classif.entry_names;
            
            % ROI
            method = compute_ROI(method,1,1);
        end
        
        function method = method_2(obj,img,library,att_names,spectrum_types,min_max_wl,min_max_xy,params)
            % Classification based on the metric
            %
            %
            global cfg
            method = obj;
            
            if ~exist('params')
                params.nb_clust = 1000;
                params.nb_repl = 1;
            end
            
            % Trims the spectra 
            att_names_C = att_names(strcmp(library.att_types(contains(att_names,library.att_groups)),'C'));
            for a = 1:length(att_names_C)
                min_val.(att_names_C{a}) = min_max_wl(2*a - 1);
                max_val.(att_names_C{a}) = min_max_wl(2*a);
                % Remove the unnecessary wavelengths
                % (Not applied to the img as it will be done in the
                % 'consis_multi_att' function)
                library = rm_min_max_wls(library,att_names_C{a},min_val.(att_names_C{a}),max_val.(att_names_C{a}));
            end
            
            % Crop the image
            img = img_modify(img,1,min_max_xy);
            [library,img] = consis_multi_att(library,img,att_names,1);
            
            % Creation of the classifier
            classif = Unit_Classifier(img);
            
            % Change spectra type
            classif.ma_img = apply_change_spec(classif.ma_img,att_names,[1,1],0);
            library = apply_change_spec(library,att_names,[1,1],0);
            
            % Kmeans classification on all data
            nb_cluster = params.nb_clust;
            clf_metric = 'sqeuclidean'; %'sqeuclidean' 'cityblock', 'cosine', 'correlation', 'hamming'
            nb_replicate = params.nb_repl;
            max_iter = 5;
            classif = Classifier_UnSuprv_Kmeans(classif,clf_metric,cfg.nb_cluster,nb_replicate,max_iter,20);
            
            % Filling the missing values
            img.att_data = fillmissing(img.att_data,'linear',2);
            
            % Change spectra type
            classif.ma_img = apply_change_spec(classif.ma_img,att_names,spectrum_types,0);
            library = apply_change_spec(library,att_names,spectrum_types,0);
            
            % Kmeans classification on all data
            clf_metric = 'sqeuclidean';
            nb_cluster = 50;
            nb_replicate = 12;
            max_iter = 100;
            classif = Classifier_UnSuprv_Kmeans(classif,clf_metric,nb_cluster,nb_replicate,max_iter);

            % Classifier with metric
            metric = 'sqeuclidean';
            classif = Classifier_Suprv_metric(classif,metric,library);
            [method.table_error_lib,method.error_lib_rgb] = build_table_img_res(classif,'cont');
            method.data_error_lib = classif.ma_img;
            
            % From subminerals to minerals
            classif = Classifier_Suprv_submin2min(classif,classif.prop_class);
            [method.table_minerals_ind,method.minerals_ind_rgb] = build_table_img_res(classif,'min');
            
            % Properties
            method.ma_img = classif.ma_img;
            method.prop_class = classif.prop_class;
            method.error = classif.error;
            method.entry_names = classif.entry_names;
            
            % ROI
            method = compute_ROI(method,1,1);
        end
        
        function method = method_3(obj,img,library,att_names,spectrum_types,min_max_wl,min_max_xy)
            % Clustering without library
            %
            %
            global cfg
            method = obj;
            
            % Trims the spectra 
            att_names_C = att_names(strcmp(img.att_types(contains(att_names,img.att_groups)),'C'));
            for a = 1:length(att_names_C)
                min_val.(att_names_C{a}) = min_max_wl(2*a - 1);
                max_val.(att_names_C{a}) = min_max_wl(2*a);
                % Remove the unnecessary wavelengths
                img = rm_min_max_wls(img,att_names_C{a},min_val.(att_names_C{a}),max_val.(att_names_C{a}));
            end
            
            % Crop the image
            img = img_modify(img,1,min_max_xy);
            
            % Creation of the classifier
            classif = Unit_Classifier(img);
            
            % Kmeans classification on all data
            clf_metric = 'sqeuclidean'; %'sqeuclidean' 'cityblock', 'cosine', 'correlation', 'hamming'
            nb_replicate = 1;
            max_iter = 5;
            classif = Classifier_UnSuprv_Kmeans(classif,clf_metric,cfg.nb_cluster,nb_replicate,max_iter,20);
            
            % Filling the missing values
            img.att_data = fillmissing(img.att_data,'linear',2);
            
            % Change spectra type
            classif.ma_img = apply_change_spec(classif.ma_img,att_names,spectrum_types,0);
            
            % Hierarchical classification
            hierC_metric = 'correlation';
            cutoff = 0.03;
            classif = Classifier_UnSuprv_HierC(classif,hierC_metric,cutoff);
            
            % Kmeans classification
            clf_metric = 'sqeuclidean';
            nb_cluster = 20;
            nb_replicate = 12;
            max_iter = 100;
            classif = Classifier_UnSuprv_Kmeans(classif,clf_metric,nb_cluster,nb_replicate,max_iter);
            
            % Give names and errors
            classif.entry_names = cellstr(num2str([1:length(classif.prop_class)]'));
            classif.error = zeros(length(classif.prop_class),1);
            
            % Creation of the minerals table
            [method.table_minerals_ind,method.minerals_ind_rgb] = build_table_img_res(classif,'min');
            
            % Properties
            method.ma_img = classif.ma_img;
            method.prop_class = classif.prop_class;
            method.error = classif.error;
            method.entry_names = classif.entry_names;
            
            % ROI
            method = compute_ROI(method,1,1);
        end
        
        function method = method_4(obj,img,library,att_names,spectrum_types,min_max_wl,min_max_xy)
            % Application of ML by training on the library
            %
            %
            global cfg
            method = obj;
            
            % Check the spectrum_type
            spectrum_types = [3,3];
            
            % Trims the spectra 
            att_names_C = att_names(strcmp(...
                library.att_types(contains(att_names,...
                                           library.att_groups)), 'C'));
            for a = 1:length(att_names_C)
                min_val.(att_names_C{a}) = min_max_wl(2*a - 1);
                max_val.(att_names_C{a}) = min_max_wl(2*a);
                % Remove the unnecessary wavelengths
                % (Not applied to the img as it will be done in the
                % 'consis_multi_att' function)
                library = rm_min_max_wls(library,...
                                         att_names_C{a},...
                                         min_val.(att_names_C{a}),...
                                         max_val.(att_names_C{a}));
            end
            
            % Crop the image
            img = img_modify(img,1,min_max_xy);
            [library,img] = consis_multi_att(library,img,att_names,1);
            
            % Creation of the classifier
            classif = Unit_Classifier(img);
            
            % Kmeans classification on all data
            clf_metric = 'sqeuclidean'; %'sqeuclidean' 'cityblock', 'cosine', 'correlation', 'hamming'
            nb_replicate = 1;
            max_iter = 5;
            n_sub_group = 20;
            t_kmeans = tic;
            classif = Classifier_UnSuprv_Kmeans(classif,...
                                                clf_metric,...
                                                cfg.nb_cluster,...
                                                nb_replicate,...
                                                max_iter,...
                                                n_sub_group);
            if cfg.print_times
                toc(t_kmeans)
            end
            
            % Filling the missing values
            img.att_data = fillmissing(img.att_data,'linear',2);
            
            % Change spectra type
            classif.ma_img = apply_change_spec(classif.ma_img,att_names,spectrum_types,0);
            library = apply_change_spec(library,att_names,spectrum_types,0);
                        
            % Applying the trained the model
            t_ml = tic;
            classif = Classifier_Suprv_Tree(classif,library,2);
            if cfg.print_times
                toc(t_ml)
            end
            
            % Make the tables and images 
            [method.table_error_lib,method.error_lib_rgb] = build_table_img_res(classif,'cont');
            method.data_error_lib = classif.ma_img;
            
            % From subminerals to minerals
            classif = Classifier_Suprv_submin2min(classif,classif.prop_class);
            [method.table_minerals_ind,method.minerals_ind_rgb] = build_table_img_res(classif,'min');
            
            % Properties
            method.ma_img = classif.ma_img;
            method.prop_class = classif.prop_class;
            method.error = classif.error;
            method.entry_names = classif.entry_names;
            
            % ROI
            method = compute_ROI(method,1,1);
        end
        
        function method = method_5(obj,img,library,att_names,spectrum_types,min_max_wl,min_max_xy)
            % Application of a custom tree classifier
            %
            global cfg
            method = obj;
            
            % Check the spectrum_type
            spectrum_types = [3,3];
            
            % Trims the spectra 
            att_names_C = att_names(strcmp(img.att_types(contains(att_names,img.att_groups)),'C'));
            for a = 1:length(att_names_C)
                min_val.(att_names_C{a}) = min_max_wl(2*a - 1);
                max_val.(att_names_C{a}) = min_max_wl(2*a);
                % Remove the unnecessary wavelengths
                img = rm_min_max_wls(img,att_names_C{a},min_val.(att_names_C{a}),max_val.(att_names_C{a}));
            end
            
            % Crop the image
            img = img_modify(img,1,min_max_xy);
            
            % Creation of the classifier
            classif = Unit_Classifier(img);
            
            % Kmeans classification on all data
            clf_metric = 'sqeuclidean'; %'sqeuclidean' 'cityblock', 'cosine', 'correlation', 'hamming'
            nb_replicate = 1;
            max_iter = 5;
            classif = Classifier_UnSuprv_Kmeans(classif,clf_metric,cfg.nb_cluster,nb_replicate,max_iter,20);
            
            % Filling the missing values
            img.att_data = fillmissing(img.att_data,'linear',2);
            
            % Change spectra type
            classif.ma_img = apply_change_spec(classif.ma_img,att_names,spectrum_types,0);
            
            % Apply custom tree
            classif = Classifier_Suprv_Tree_custom(classif);
            [method.table_error_lib,method.error_lib_rgb] = build_table_img_res(classif,'cont');
            method.data_error_lib = classif.ma_img;
            
            % From subminerals to minerals
            classif = Classifier_Suprv_submin2min(classif,classif.prop_class);
            [method.table_minerals_ind,method.minerals_ind_rgb] = build_table_img_res(classif,'min');
            
            % Properties
            method.ma_img = classif.ma_img;
            method.prop_class = classif.prop_class;
            method.error = classif.error;
            method.entry_names = classif.entry_names;
            
            % ROI
            method = compute_ROI(method,0,1);
        end
        
        function obj = compute_ROI(obj,error_lib,minerals_ind)
            % Calculation of the ROI and definition of the images
            %
            % error_lib: boolean
            %       option to apply on the error
            % minerals_ind: boolean
            %       option to apply on the mineral proportion
            %
            global cfg
            
            if ~minerals_ind && ~error_lib
                error('No option selected.')
            end
            
            % Calculation of the ROI size
            size_window(1) = Inf;
            size_window(2) = 4*cfg.mass_ROI*10^6/(pi*cfg.diam_sample^2*cfg.average_density_samples); % size of the section in cm
            
            map = obj.ma_img.ind_img;
            
            % Size of the image in real coordinates
            x_img_min = min(obj.ma_img.axis.x(:));
            x_img_max = max(obj.ma_img.axis.x(:));
            y_img_min = min(obj.ma_img.axis.y(:));
            y_img_max = max(obj.ma_img.axis.y(:));
            l_img_x = x_img_max - x_img_min;
            l_img_y = y_img_max - y_img_min;
            
            % Size of the ROI in real coord
            % the size bounded by the size of the image
            l_roi_x = min(size_window(1),l_img_x);
            l_roi_y = min(size_window(2),l_img_y);

            % Size of the image in pixels
            [N_pxl_y, N_pxl_x] = size(map);
            res_x = l_img_x/N_pxl_x;
            res_y = l_img_y/N_pxl_y;
            
            % Size of the ROI in pixels
            l_pxl_x = round(l_roi_x/res_x,0);
            l_pxl_y = round(l_roi_y/res_y,0);
            
            % size of the moving window
            kernel = ones(l_pxl_y,l_pxl_x);
            max_nb = conv2(ones(N_pxl_y, N_pxl_x),kernel,'same');
            
            % Mineral -----------------------------------------------------
            if minerals_ind
                % Mineral map
                min_all = unique(map);
                min_all(isnan(min_all)) = [];
                nb_min = length(min_all);

                % Global proportion
                for el = 1:length(min_all)
                   mineral_count(el) = sum(sum(map == min_all(el)));
                end
                global_prop = mineral_count/N_pxl_y/N_pxl_x;
                
                %{
                % Dealing with low proportions
                % It is assumed that the proportions below 10% won't be
                % detected properly with ID2B, for low proportions we thus
                % artificially increase their values
                min_detect = 0.005;
                min_ID2B = 0.10;
                global_prop((global_prop < min_ID2B) & (global_prop > min_detect)) = min_ID2B;
                global_prop(global_prop < min_detect) = 0;    

                % The higher values are renormalized to obtain a sum of 100
                global_prop(global_prop > min_ID2B) = global_prop(global_prop > min_ID2B) * (1 - sum(global_prop(global_prop == min_ID2B))) / sum(global_prop(global_prop > min_ID2B));
                %}
                
                % Local proportions
                local_prop = zeros(N_pxl_y, N_pxl_x,nb_min);

                % local prop: maximum possible number per pixel
                for n = 1:nb_min
                    local_prop(:,:,n) = conv2(map == min_all(n), kernel, 'same')./max_nb;
                end
                
                % Chi2 test
                local_chi2 = zeros(N_pxl_y,N_pxl_x);
                for i = 1:N_pxl_y
                    for j = 1:N_pxl_x
                        local_chi2(i,j) = chi2_pro(global_prop', squeeze(local_prop(i,j,:)), l_pxl_x*l_pxl_y);
                    end
                end
                
                % Calculation of the ROI for prop
                [obj.ROI_minerals_ind, max_pval] = build_ROI_from_func(local_chi2, l_roi_x, l_roi_y, l_pxl_x, l_pxl_y, max_nb, obj.ma_img.coord);
            end
            
            % Error -------------------------------------------------------
            if error_lib
                % Average local error
                error_img_tmp = make_img_from_var(obj.ma_img,obj.error);
                error_img_tmp(isnan(error_img_tmp)) = max(error_img_tmp(:));  % Puts the NaN at the maximum value of the dataset to allow the convolution to run
                local_error = conv2(error_img_tmp,kernel,'same')./max_nb;
                
                % Calculation of the ROI for the error
                [obj.ROI_error_lib,max_err] = build_ROI_from_func(local_error,l_roi_x,l_roi_y,l_pxl_x,l_pxl_y,max_nb,obj.ma_img.coord);
            end
            
            % Creation of the ROI figure(s) -------------------------------
            obj = build_ROI_rgb(obj,error_lib,minerals_ind);
        end
        
        function obj = build_ROI_rgb(obj,error_lib,minerals_ind)
            % Creation of figure with the ROI 
            %
            % error_lib: boolean
            %       option to apply on the error
            % minerals_ind: boolean
            %       option to apply on the mineral proportion
            %
            %
            if ~error_lib && ~minerals_ind
                error('No option selected.')
            end
            
            fig = [];
            if error_lib
                fig = [fig, {'error_lib'}];
            end
            if minerals_ind
                fig = [fig, {'minerals_ind'}];
            end
            
            % export of the ROI images
            coord_x = obj.ma_img.axis.x;
            coord_y = obj.ma_img.axis.y;
            
            for f = 1:length(fig)
                if ~isempty(obj.(strcat('ROI_',fig{f})))
                    ROI_tmp = obj.(strcat('ROI_',fig{f}));
                else
                    ROI_tmp = [0,0,0,0];
                end
                
                % Figure creation
                fig_exp = figure('visible','off','Name','fig_exp');
                
                % Axes creation
                h1 = axes;
                cla;
                hold on;
                
                if coord_x(1) < coord_x(2)
                    h1.XDir = 'normal';
                    h1.XAxis.TickValues = coord_x;
                else
                    h1.XDir = 'reverse';
                    %coord_x = flip(coord_x);
                    h1.XAxis.TickValues = flip(coord_x);
                end
                h1.XTickMode = 'auto';
                if coord_y(1) < coord_y(2)
                    h1.YDir = 'reverse';
                    h1.YAxis.TickValues = coord_y;
                else
                    h1.YDir = 'normal';
                    %coord_y = flip(coord_y);
                    h1.YAxis.TickValues = flip(coord_y);
                end
                h1.YTickMode = 'auto';
                
                % Items of the axes
                img = obj.(strcat(fig{f}, '_rgb'));
                set(imagesc(img),'XData', coord_x,'YData',coord_y);
                %plot(ROI_tmp(1:2), ROI_tmp(3)*ones(1,2), 'w', 'LineWidth', 3);
                %plot(ROI_tmp(1:2), ROI_tmp(4)*ones(1,2), 'w', 'LineWidth', 3);
                
                h1.XLim = [min(obj.ma_img.axis.x), max(obj.ma_img.axis.x)];
                h1.YLim = [min(obj.ma_img.axis.y), max(obj.ma_img.axis.y)];
                h1.DataAspectRatio = [1,1,1];
                
                % Try to change the figure's size to fit the axis
                try
                    img_sz = size(fig_exp.Children.Children(3).CData);
                    pos = fig_exp.Position;
                    pos(4) = pos(3)/img_sz(2)*img_sz(1);
                    set(fig_exp,'Position',pos)
                catch
                end
                
                % Saving the file
                obj.(strcat('ROI_',fig{f},'_rgb')) = fig_exp;
            end
        end
        
        function export_key_param(obj,dir)
            % Exports the key images
            %
            %
            fig = [{'error_lib'}, {'minerals_ind'}];
            
            for f = 1:length(fig)
                % export the tables
                if ~isempty(obj.(strcat('table_',fig{f})))
                    writetable(obj.(strcat('table_',fig{f})),fullfile(dir,strcat(fig{f},'.csv')));
                end
                
                % export of the AOI values as txt
                if ~isempty(obj.(strcat('ROI_',fig{f})))
                    csvwrite(fullfile(dir,strcat('ROI_',fig{f},'.csv')),obj.(strcat('ROI_',fig{f})));
                end
                
                % save the minerals image with the ROI
                if ~isempty(obj.(strcat('ROI_',fig{f},'_rgb')))
                    img_tmp = copy(obj.(strcat('ROI_',fig{f},'_rgb')));
                    img_tmp.CurrentAxes.XTick = [];
                    img_tmp.CurrentAxes.YTick = [];
                    img_tmp.Children.Position = [0, 0, 1, 1];
                    img_tmp.Children.Units = 'pixels';
                    img_tmp.Position = img_tmp.Children.Position;
                    img_tmp.PaperPosition = [0, 0, 1, 1];
                    exportgraphics(img_tmp, fullfile(dir, strcat(fig{f}, '.png')))
                    
                    % Export Image legend
                    if strcmp(fig{f},'minerals_ind')
                        tab_exp = obj.(strcat('table_',fig{f}))(:,{'Minerals','Color'});
                        for i = 1:size(tab_exp,1)
                            tab_exp{i,'Color'}{1} = tab_exp{i,'Color'}{1}(41:47);
                        end
                    elseif strcmp(fig{f},'error_lib')
                        max_v = max(cellfun(@str2num,obj.(strcat('table_',fig{f})){:,'Error'}));
                        min_v = min(cellfun(@str2num,obj.(strcat('table_',fig{f})){:,'Error'}));
                        
                        list = min_v:((max_v-min_v)/10):max_v;
                        colors_classes = parula(length(list));
                        
                        % Color of the table
                        for c = 1:length(list)
                            color_column(c,1) = {rgb2hex(colors_classes(c,:))};
                        end
                        
                        tab_exp = table(list',color_column,'VariableNames',{'Error','Color'});
                    else
                        error('Wrong option...')
                    end
                    writetable(tab_exp,fullfile(dir,strcat('Colors_',fig{f},'.csv')));
                end
            end
            
            % save the RGB image
            if ~isempty(obj.ma_img.img_rgb)
                img_tmp = copy(obj.ma_img.img_rgb);
                img_tmp.CurrentAxes.XTick = [];
                img_tmp.CurrentAxes.YTick = [];
                img_tmp.Children.Position = [0, 0, 1, 1];
                img_tmp.Children.Units = 'pixels';
                img_tmp.Position = img_tmp.Children.Position;
                img_tmp.PaperPosition = [0, 0, 1, 1];
                exportgraphics(img_tmp, fullfile(dir,'HS_rgb_recons.png'))
                
                % Export Image scale 
                tab_exp_legend = table(img_tmp.Children.XLim',img_tmp.Children.YLim','VariableNames',{'X','Y'},'rowNames',{'Min','Max'});
                writetable(tab_exp_legend, fullfile(dir,strcat('Axis_XYlim_',fig{f},'.csv')),'WriteRowNames',true);
            end
        end
    end
end

function [ROI,max_val] = build_ROI_from_func(local_val,l_roi_x,l_roi_y,l_pxl_x,l_pxl_y,max_nb,coord)
    % 
    %
    %
    
    % Exclusion of the pixels too close to the border, by dividing
    % with a high number their value will not be a maximum
    exclu = max_nb;
    exclu(exclu < l_pxl_x*l_pxl_y) = -Inf;
    exclu(exclu >= l_pxl_x*l_pxl_y) = 0;
    local_val = local_val + exclu;
    
    % Find the maximum value
    [max_val,ind_xy] = max(local_val(:));
    [Id_x, Id_y] = ind2sub(size(local_val),ind_xy);
    
    % Build the window
    half_l_x = round(l_roi_x/2,0);
    half_l_y = round(l_roi_y/2,0);

    % Local prop
    coord_roi_xc = coord.x(Id_x,Id_y);
    coord_roi_yc = coord.y(Id_x,Id_y);
    x_roi_min = ceil(coord_roi_xc - half_l_x);
    x_roi_max = floor(coord_roi_xc + half_l_x);
    y_roi_min = ceil(coord_roi_yc - half_l_y);
    y_roi_max = floor(coord_roi_yc + half_l_y);

    %
    ROI = [x_roi_min,  x_roi_max, y_roi_min, y_roi_max];
end

function p_value = chi2_pro(prop_observed,prop_expected,N_tot)
% Returns the p-value of two categorical variables
%
% Parameters
% ----------
% prop_observed: array (1xn)
%       proportions of the n observed categories
% prop_expected: array (1xn)
%       proportions of the n expected categories
%
% Output
% ------
% p_value: double
%       p value
%

sel_no_zero = prop_expected>0;
prop_observed = prop_observed(sel_no_zero);
prop_expected = prop_expected(sel_no_zero);

N_observed = N_tot*prop_observed;
N_expected = N_tot*prop_expected;

chi2_rs = sum((N_observed - N_expected).^2./N_expected);

df = length(N_observed) - 1;

p_value = chi2cdf(chi2_rs, df, 'upper');

end

%{
% // ====================================================================
% // This file is part of the Endmember Induction Algorithms Toolbox for MATLAB 
% // Copyright (C) Grupo de Inteligencia Computacional, Universidad del 
% // País Vasco (UPV/EHU), Spain, released under the terms of the GNU 
% // General Public License.
% //
% // Endmember Induction Algorithms Toolbox is free software: you can redistribute 
% // it and/or modify it under the terms of the GNU General Public License 
% // as published by the Free Software Foundation, either version 3 of the 
% // License, or (at your option) any later version.
% //
% // Endmember Induction Algorithms Toolbox is distributed in the hope that it will
% // be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
% // of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% // General Public License for more details.
% //
% // You should have received a copy of the GNU General Public License
% // along with Endmember Induction Algorithms Toolbox. 
% // If not, see <http://www.gnu.org/licenses/>.
% // ====================================================================
%}

function [E,C] = EIA_NFINDR(data,p,maxit)
    % [E,C] = EIA_NFINDR(data,p,maxit)
    %
    % Manuel Grana <manuel.grana[AT]ehu.es>
    % Miguel Angel Veganzones <miguelangel.veganzones[AT]ehu.es>
    % Grupo de Inteligencia Computacional (GIC), Universidad del Pais Vasco /
    % Euskal Herriko Unibertsitatea (UPV/EHU)
    % http://www.ehu.es/computationalintelligence
    % 
    % Copyright (2011) Grupo de Inteligencia Computacional @ Universidad del Pais Vasco, Spain.
    % Copyright (2007) GRNPS group @ University of Extremadura, Spain. 
    %
    % N-FINDR endmembers induction algorithm.
    % ------------------------------------------------------------------------------
    % Input:   data      : column data matrix [nvariables x nsamples]
    %          p         : number of endmembers to be induced. If not provided it is calculated by HFC method with tol=10^(-5)
    %          maxit     : maximum number of iterations. Default = 3*p
    %
    % Output:  E         : set of induced endmembers [nvariables x p]
    %          C         : induced endmembers indexes vector [nsamples] with {0,1} values, where '1' indicates that the corresponding sample has been identified as an endmember.
    %
    % Bibliographical references:
    % [1] Winter, M. E., «N-FINDR: an algorithm for fast autonomous spectral end-member determination in hyperspectral data», presented at the Imaging Spectrometry V, Denver, CO, USA, 1999, vol. 3753, págs. 266-275.

    % Parameters
    if (nargin < 1)
        error('Insufficient parameters');
    end
    if (nargin < 2 || p <= 0)
        p = EIA_HFC(data,10^(-5));
    end
    if (nargin < 3 || maxit <= 0)
        maxit = 3*p;
    end

    % data size
    [nvariables,nsamples] = size(data);

    % Dimensionality reduction by PCA
    [~, zscores] = pca(data');
    data_pca = squeeze(zscores(:,1:p-1))';

    % Initialization
    E = zeros(nvariables,p);
    C = zeros(1,nsamples);
    IDX = zeros(1,p);
    TestMatrix = zeros(p);
    TestMatrix(1,:) = 1;
    for i = 1:p
        idx = floor(rand*nsamples) + 1;
        TestMatrix(2:p,i) = data_pca(:,idx);
        IDX(i) = idx;
    end
    actualVolume = abs(det(TestMatrix)); % instead of: volumeactual = abs(det(MatrixTest))/(factorial(p-1));
    it = 1;
    v1 = -1;
    v2 = actualVolume;

    % Algorithm
    while it<=maxit && v2>v1
        for k=1:p
            for i=1:nsamples
                actualSample = TestMatrix(2:p,k);
                TestMatrix(2:p,k) = data_pca(:,i);
                volume = abs(det(TestMatrix));  % instead of: volume = abs(det(MatrixTest))/(factorial(p-1));
                if volume > actualVolume
                    actualVolume = volume;
                    IDX(k) = i;
                else
                    TestMatrix(2:p,k) = actualSample;
                end
            end
        end
        it = it+1;
        v1 = v2;
        v2 = actualVolume;
    end
    for i = 1:p
        E(:,i) = data(:,IDX(i));
        C(IDX(i)) = 1;
    end
end