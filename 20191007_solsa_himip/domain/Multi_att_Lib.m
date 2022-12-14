classdef Multi_att_Lib < Multi_att
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        entry_names
        last_edited
        model_bag
        file_name
    end
    
    methods
        function obj = Multi_att_Lib(att_names,att_data,entry_names)
            % Creation of the object Multi_att_Lib
            %
            if nargin == 0
                att_names = [];
                att_data = [];
                entry_names = [{}];
            end
            obj = obj@Multi_att(att_names,att_data);
            obj.entry_names = entry_names;
        end
        
        function lib = load_csv(obj,file_path)
            %
            % file path : str
            %       path of the .csv file
            % data_path = 'D:\SOLSA\HIMIP\Hyper_Spectral_Libraries\test_library_fusion.csv';
            
            lib_table = readtable(file_path,'ReadVariableNames',true);
            
            colnames = lib_table.Properties.VariableNames;
            att_names = colnames(2:end);
            
            % Min names
            mineral_list = lib_table.Mineral_ID;
            
            % Att_data
            att_data = lib_table{:,2:end};
            
            % Check 
            if size(mineral_list,1) ~= size(att_data,1)
                error('Library has some entry in the wrong format...')
            end
            
            % Creation of the object
            lib = Multi_att_Lib(att_names,att_data,mineral_list);
            
            lib = build_prop_from_names(lib);
            
            % saving the file name
            file_spl = split(file_path,'\');
            file_spl = file_spl(end);
            file_spl = split(file_spl,'.');
            lib.file_name = file_spl{1};
        end
        
        function lib = load_usgs(obj,dir_path)
            % Load the library from the USGS Library
            %
            % Parameters
            % ----------
            % dir_path : string
            %        path of the library directory 
            %
            global cfg
            
            % Finding the wavelengths file and mineral directory
            home_dir = pwd;
            cd(dir_path);
            wls_file_path = fullfile(dir_path,getfield(dir('*Wavelengths*.txt'), 'name'));
            dir_spectra = fullfile(dir_path,getfield(dir('*Mineral*'), 'name'));
            cd(dir_spectra);
            spectra_files = dir('*');
            spectra_files = spectra_files(~[spectra_files.isdir]);
            spectra_files_names = {spectra_files.name}';
            cd(home_dir);
            
            % submineral list
            sepInd = strfind(spectra_files_names, '_');
            mineral_list = cellfun(@(x,y) x(y(2) + 1:end - 4), spectra_files_names, sepInd, 'un', 0);
            
            % wavelengths
            fileID = fopen(wls_file_path,'r');
            wls_load = textscan(fileID,'%f','HeaderLines',1);
            wavelengths = wls_load{1}*1000;     % conversion from µm to nm
            
            % Loads the spectra
            att_data = zeros(length(spectra_files),length(wavelengths));
            for i = 1:length(spectra_files)
                % spectra
                fileID = fopen(fullfile(dir_spectra,spectra_files(i).name),'r');
                spectrum_obj = textscan(fileID,'%f','HeaderLines',1);
                fclose(fileID);
                spectrum = spectrum_obj{1};
                spectrum(spectrum<0) = NaN;     % replace negative values by NaN
                att_data(i,:) = spectrum; 
            end
            
            % Attributes properties
            att_names = get_names_from_wvs(wavelengths);
            
            % Creation of the object
            lib = Multi_att_Lib(att_names,att_data,mineral_list);
            
            lib = build_prop_from_names(lib);
        end
        
        function lib = load_eco_lib(obj,dir_path)
            %{
            %}
            global cfg 
            lib = Multi_att_Lib();
            
            % Finding the spectrum files
            home_dir = pwd;
            cd(dir_path);
            spectra_files = dir('*spectrum.txt');
            spectra_files_names = {spectra_files.name}';
            
            % Loads the spectra
            for i = 1:length(spectra_files)
                % spectra and wavelength
                fileID = fopen(fullfile(dir_path,spectra_files_names{i}),'r');
                spectrum_obj = textscan(fileID,'%f %f','HeaderLines',21);
                fclose(fileID);
                wavelengths = spectrum_obj{1}*1000;
                spectrum = spectrum_obj{2}/100;
                if wavelengths(1) > wavelengths(2)
                    wavelengths = flip(wavelengths);
                    spectrum = flip(spectrum);
                end
                att_names = get_names_from_wvs(wavelengths);
                
                % Reload the file to get the name of the mineral 
                fileID = fopen(fullfile(dir_path,spectra_files_names{i}),'r');
                header = textscan(fileID,'%s');
                fclose(fileID);
                entry_name = [header{1}{2},'_',header{1}{7},'_',header{1}{12}];
                
                % Creation of the temporary library containing one spectrum
                lib_tmp = Multi_att_Lib(att_names,spectrum',entry_name);
                
                % Add to the export library
                lib = lib + lib_tmp;
            end
            cd(home_dir);
            
            lib = build_prop_from_names(lib);
        end
        
        function lib = get_from_indexes(obj,ind_lin,ind_col)
            %{
            
            %}
            lib = get_from_indexes@Multi_att(obj,ind_lin,ind_col);
            
            % Case of empty indices variable
            if isempty(ind_lin)
                ind_lin = true(1,size(obj.att_data,1));
            end
            
            lib.entry_names = lib.entry_names(ind_lin);
        end
        
        function lib = plus(obj,lib_add)
            %{
            %}
            
            lib = plus@Multi_att(obj,lib_add);
            
            lib.entry_names = obj.entry_names;
            if iscellstr(lib_add.entry_names)
                nb_lib = length(lib_add.entry_names);
            else
                nb_lib = 1;
            end
            for i = 1:nb_lib
                % check is cellstr or str
                if iscellstr(lib_add.entry_names)
                    add_name_tmp = lib_add.entry_names{i};
                else
                    add_name_tmp = lib_add.entry_names;
                end
                % check if 'sample' already exists
                sep_add_n_tmp = strfind(add_name_tmp, '_');
                min_tmp = add_name_tmp(1:sep_add_n_tmp(1) - 1);
                test_ret_min = add_name_tmp(sep_add_n_tmp(1) + 1:end);
                cont = contains(lib.entry_names,test_ret_min);
                if max(cont)
                    names_dup = lib.entry_names(cont);
                    sep_add_n_tmp = strfind(names_dup, '_');
                    ids = cellfun(@(x,y) str2double(x(y(end) + 1:end)), names_dup, sep_add_n_tmp, 'un', 0);
                    lib.entry_names = cellstr([lib.entry_names; [min_tmp,'_',test_ret_min,'_',num2str(max(cell2mat(ids)) + 1)]]);
                else
                    ids = add_name_tmp(sep_add_n_tmp(end) + 1:end);
                    if all(ismember(ids,'0123456789'))
                        % if new entry already has an ID
                        lib.entry_names = cellstr([lib.entry_names; add_name_tmp]);
                    else
                        % if new entry doesn't have an ID
                        lib.entry_names = cellstr([lib.entry_names; [add_name_tmp,'_1']]);
                    end
                end
            end
        end
        
        function lib = add_modif(obj,delete,entry_names,att_add)
            %{ 
            Addition or modification of an entry of the database
            
            Parameters 
            ----------
            del : boolean
                   option to delete the entry
            entry_names: list of string
                    names of entries to be deleted / added
            att_add : :obj:
                    multi_att object to be added
            %}
            
            % Test if the entry name and ID are already in the library
            comp_entry = ismember(obj.entry_names,entry_names);
            test_exist = max(comp_entry);
            
            if ~isempty(entry_names)
                if test_exist
                    if delete
                        disp('Entry to be deleted is found. Deleting...')
                    else
                        error('Entry name already in the Library. Aborting...')
                    end
                else
                    if delete
                        error('Entry to be deleted is not found. Aborting...')
                    else
                        disp('Additing the entry...')
                    end
                end

                if delete
                    % parameters
                    obj.att_data(comp_entry,:) = [];

                    % names
                    obj.entry_names(comp_entry) = [];
                else
                    % parameters
                    obj = obj + att_add;
                    % ! In the case of an empty Library the first element will
                    % define the Library parameters !! 
                end
            end
            
            % Creation of the object
            lib = obj;
        end
        
        function lib_table = build_lib_table(obj)
            % Builds the table of the Library
            %
            % Output
            % ------
            % lib_table : table
            %        table containing all the info of the library
            %
            
            % Create colname
            colnames = [{'Mineral_ID'}, obj.att_names ];
            
            % Size of the data
            [nb_line, nb_col] = size(obj.att_data);
            
            % Merge all data
            lib_table = cell2table(cell(nb_line,nb_col+1),'VariableNames',colnames);
            lib_table.Mineral_ID = obj.entry_names;
            if nb_line > 0
                lib_table{:,2:end} = num2cell(obj.att_data);
            end
        end
        
        function lib = to_subs_ref(obj,option,param)
            %
            %
            %
            ma = to_subs_ref@Multi_att(obj,option,param);
            lib = copy(ma,'lib');
            prop = setdiff(properties(Multi_att_Lib),properties(Multi_att));
            for i = 1:length(prop)
                lib.(prop{i}) = obj.(prop{i});
            end
        end
        
        function model = reduce_lib(obj,split_perc,model_option,multi_class)
            %{
            
             Parameters
             ----------
             split_perc: float, optionnal
                   percentage of the dataset to be used for the training
                   default 70%
             model_option: int, optionnal
                   model to be trained (1: decision tree (default), 2:
                   random_forest)
            multi_class: boolean, optionnal
                    option to compute the regression on multi classes
            
             Returns
             -------
             model: :obj:
                   model fitted on the training dataset
            %}
            
            if nargin <=1
                split_perc = 70;
            end
            if nargin <= 2
                model_option = 1;
            end
            if nargin <=3
                multi_class = 1;
            end
            
            names = obj.entry_names;
            
            % Extract the names of the minerals in the entry_names
            % attributes
            names_sep = cellfun(@(x) split( x , ["-", "_"] ), names, 'un',0);  % separate the entry_names regarding - and _
            names_bool_low = cellfun(@(x) isstrprop(x,'lower'),names_sep,'un',0);  % identification of the lower case letters
            % the char kept have first a capital letter and following only lower cases 
            names_fin = cellfun(@(x1,y1) x1(cell2mat(cellfun(@(x1,y1) y1(1)==0 & sum(y1)>=1 & sum(y1)==length(y1)-1 ...
                                                                ,x1,y1,'un',0))),names_sep,names_bool_low,'un',0);
            
            % All the possible minerals
            all_names = names_fin{1};
            for i = 2:length(names_fin)
                all_names = [all_names; names_fin{i}(~ismember(names_fin{i},all_names))];
            end
            
            % Indicators of minerals
            ind_min = zeros(length(names_fin),length(all_names));
            for i = 1:length(names_fin)
                ind_min(i,:) = contains(all_names',names_fin{i});
            end
            
            % Application of the spectra change
            att_names = [{'VNIR'}, {'SWIR'}];
            spectrum_types = [7,7];
            obj = apply_change_spec(obj,att_names,spectrum_types,0);
            
            % Definition of the data and targets matrix
            X = obj.att_data;
            Y = ind_min;
            
            % Random split
            N_samples = size(X,1);
            k = randperm(N_samples);
            split_val = round(N_samples*split_perc/100);
            
            % Application of the split
            X_train = X(k(1:split_val),:);
            X_test = X(k(split_val+1:end),:);
            Y_train = Y(k(1:split_val),:);
            Y_test = Y(k(split_val+1:end),:);
            
            if model_option == 1
                % Apply a decision tree
                for k = 1:length(all_names)
                    model = fitctree(X_train,Y_train(:,k));
                    view(model,'mode','graph')
                end
            elseif model_option == 2
                % Training
                nTrees = 100;
                
                F1_score = zeros(1,length(all_names));
                % Apply a random forest
                for k = 1:length(all_names)
                    model = TreeBagger(nTrees,X_train,Y_train(:,k), 'Method', 'classification'); 
                    Y_predict = str2double(model.predict(X_test));  % Predictions is a char though. We want it to be a number.
                    
                    % Metric
                    confMat = confusionmat(Y_predict,Y_test(:,k));
                    
                    % Precision and recall
                    recall = zeros(1:size(confMat,1));
                    precision = zeros(1:size(confMat,1));
                    for i = 1:size(confMat,1)
                        recall(i) = confMat(i,i)/sum(confMat(i,:));
                        precision(i) = confMat(i,i)/sum(confMat(:,i));
                    end
                    recall(isnan(recall)) = [];
                    Recall = sum(recall)/size(confMat,1);
                    Precision = sum(precision)/size(confMat,1);
                    
                    % F1-score
                    F1_score(k) = 2*Recall*Precision/(Precision + Recall);
                end
                [all_names,cellstr(num2str(F1_score'))]
            end
        end
            
        function lib = train_model(obj,model_option)
            %{
            
            Parameters
            ----------
            model_option: int, optionnal
                    model to be trained (1: decision tree, 2:
                    random_forest (default))
            %}
            
            lib = obj;
            
            if nargin == 1
                model_option = 2;
            end
                                    
            % Initialisations
            lib_spectra = lib.att_data;
            
            % Get the minerals (labels)
            tmp_entry_names = split(lib.entry_names{1},'_');
            y = cellstr(tmp_entry_names{1});
            for i = 2:length(lib.entry_names)
                tmp_entry_names = split(lib.entry_names{i},'_');
                y = [y; tmp_entry_names{1}];
            end
            
            % Training the model
            if model_option == 1
                lib.model = fitctree(lib_spectra,y); %,'MaxDepth',10
            elseif model_option == 2
                lib.model = TreeBagger(100,lib_spectra,y);
            end
        end
        
        function export(obj, dir_path, file_name, matlab_file)
            % Export of the Library
            %
            % Parameters
            % ----------
            % dir_path : string
            %        directory to export the file
            % file_name : string
            %        name of the exported file
            % matlab_file : boolean
            %        option to export as a matlab file (default .csv)
            %
            lib_table_exp = build_lib_table(obj);
            
            if nargin <4
                matlab_file = 0;
            end
            
            if matlab_file
                save(fullfile(dir_path,strcat(file_name,'.mat')),'lib_table_exp');
            else
                writetable(lib_table_exp, fullfile(dir_path,strcat(file_name,'.csv')));
            end
        end
        
        function show(obj)
            % Print the table
            %
            lib_table_exp = build_lib_table(obj)
        end
    end
end

        
function att_names = get_names_from_wvs(wavelengths)
    global cfg

    att_names = [];
    if min(wavelengths) < cfg.lim_vnir_swir
        att_names = strcat({'VNIR_C_'},num2str(round(wavelengths(wavelengths < cfg.lim_vnir_swir),1)))';
    end
    if max(wavelengths) > cfg.lim_vnir_swir
        att_names = [att_names, strcat({'SWIR_C_'},num2str(round(wavelengths(wavelengths >= cfg.lim_vnir_swir),1)))'];
    end
end
