classdef Multi_att
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        att_names
        att_data
        att_groups
        att_types
        att_ids
    end
    % properties
    % ----------
    % att_names: vector of str, 1xm
    %       names of the attributes
    % att_data: matrix nxm
    %       values of the attributes
    %
    
    methods
        function obj = Multi_att(att_names,att_data)
            if nargin == 0
                att_names = [];
                att_data = [];
            end
            obj.att_names = att_names;
            obj.att_data = att_data;
        end
        
        function multi_att = build_prop_from_names(obj)
            % Computes the properties (att_groups, att_types, att_ids) 
            % from the names (att_names) of a multi_att object
            % 
            
            sepInd = strfind(obj.att_names, '_');
            
            % Groups
            [obj.att_groups,groups_ind] = unique(cellfun(@(x,y) x(1:y(1)-1), obj.att_names, sepInd, 'un', 0),'stable');   % option stable to keep the same order as in the file
            
            % Types
            att_types_tmp = cellfun(@(x,y) x(y(1)+1:y(2)-1), obj.att_names, sepInd, 'un', 0);
            obj.att_types = att_types_tmp(groups_ind);
            
            % Definition of the numerical values of the continuous values
            for g = 1:length(obj.att_groups)
                tag = obj.att_groups{g};
                names_tmp = obj.att_names(startsWith(obj.att_names,tag))';
                if strcmp(obj.att_types(g),'C')
                    obj.att_ids.(obj.att_groups{g}) = str2double(strrep(erase(names_tmp,strcat(tag,'_C_')),'_','.')).';
                elseif strcmp(obj.att_types(g),'N')
                    obj.att_ids.(obj.att_groups{g}) = erase(names_tmp,strcat(tag,'_N_')).';
                end
            end
            
            % Output
            multi_att = obj;
        end

        function multi_att = build_names_from_prop(obj)
            % Computes the names (att_names) from the properties
            % (att_groups, att_types and att_ids)
            %
            % Output
            % ------
            % multi_att: :obj:
            %       copy of the multi_att object 
            %
            
            names = [];
            for g = 1:length(obj.att_groups)
                type = obj.att_types{g};
                group = obj.att_groups{g};
                try
                    names = [names, strcat(group,'_',type,'_',strrep(strrep(cellstr(num2str(obj.att_ids.(group)')),'.','_'),' ',''))'];
                catch
                    names = [names, strcat(group,'_',type,'_',cellstr(obj.att_ids.(group)')')];
                end
            end
            
            obj.att_names = names;
            multi_att = obj;
        end

        function multi_att = copy(obj,type)
            % Copy of a multi_att obj
            % 
            % Parameters
            % ----------
            % type: string
            %       either 'ma', 'img', or 'lib' depending on the output
            %       subobject to be created
            %
            % Output
            % ------
            % multi_att: :obj:
            %       copy of the multi_att object 
            %
            
            % Creation of the object
            if strcmp(type,'ma')
                multi_att = Multi_att();
            elseif strcmp(type,'img')
                multi_att = Multi_att_Img();
            elseif strcmp(type,'lib')
                multi_att = Multi_att_Lib();
            end
            
            % Copy of the properties
            multi_att.att_names = obj.att_names;
            multi_att.att_data = obj.att_data;
            multi_att.att_groups = obj.att_groups;
            multi_att.att_types = obj.att_types;
            multi_att.att_ids = obj.att_ids;
        end
        
        function ind = get_ind_from_gr_id(obj,group_name,id)
            % Finds the indice of a colum based on the group and the ID
            % For continuous variables the closest value is selected
            %
            % Parameters
            % ----------
            % group_name: string
            %       name of the group
            % id: string
            %       id to be found
            % 
            % Ouput
            % -----
            % ind: int
            %       indice matching the group_name and id
            %
            
            mem = ismember(obj.att_groups,group_name);
            
            if max(mem)
                if strcmp(obj.att_types(mem),'N')
                    % Numerical variables
                    ind = find(strcmp(obj.att_names, strcat(group_name,'_N_',id)));
                elseif strcmp(obj.att_types(mem),'C')
                    % Continuous variables
                    % Index of the first element of the group
                    sepInd = strfind(obj.att_names, '_');
                    ind_gr = contains(obj.att_groups,group_name);
                    [~,ind_st] = unique(cellfun(@(x,y) x(1:y(1)-1), obj.att_names, sepInd, 'un', 0),'stable');
                    % We searched for the closest value and use its indice
                    ind = ind_st(ind_gr) - 1 + index_close_value_reg_array(obj.att_ids.(group_name), id);
                end
            else
                error('Group name not in the object...')
            end
            
        end
        
        function multi_att = get_from_indexes(obj,ind_lin,ind_col)
            % Extract part of a multi att object based on line and column
            % indices
            %
            % Parameters
            % ----------
            % ind_lin: array of int
            %       indices of the lines to be extracted
            % ind_col: array of int
            %       indices of the columns to be extracted
            %
            % Output
            % ------
            % multi_att: :obj:
            %       multi_att object containing only the selected indices
            %
            
            multi_att = obj;
            
            % Case of empty indices variable
            if isempty(ind_lin)
                ind_lin = true(1,size(obj.att_data,1));
            end
            
            if isempty(ind_col)
                ind_col = true(1,size(obj.att_data,2));
            end
            
            % Application on the multi_att parameters
            multi_att.att_data = obj.att_data(ind_lin,ind_col);
            multi_att.att_names = obj.att_names(ind_col);
            multi_att = build_prop_from_names(multi_att);
        end
        
        function multi_att = get_from_group(obj,group_names)
            % Extract the data of a property group 
            %
            % Parameters
            % ----------
            % group_names: cell string array
            %        name of the groups to be kept
            %
            % Output
            % ------
            % multi_att: :obj:
            %       multi_att object containing only the selected groups
            %
           
            multi_att = obj;
            
            % Check if group_names contain valid names
            valid_name_ind = ismember(obj.att_groups,group_names);
            valid_gr_names = obj.att_groups(valid_name_ind);
            
            % Check if there are groups to be selected
            if ~isempty(valid_gr_names)
                % Init
                sel_col_to_keep = false(1,length(obj.att_names));
                multi_att.att_ids = [];
                
                % Loop over the selected groups
                for g = 1:length(valid_gr_names)
                    sel_col_to_keep = sel_col_to_keep | contains(obj.att_names,valid_gr_names{g});
                    multi_att.att_ids.(valid_gr_names{g}) = obj.att_ids.(valid_gr_names{g});
                end

                % Application of the selection to parameters
                multi_att.att_data = obj.att_data(:,sel_col_to_keep);
                multi_att.att_names = obj.att_names(sel_col_to_keep);
                multi_att.att_groups = valid_gr_names;
                multi_att.att_types = obj.att_types(valid_name_ind);
            elseif ~isempty(group_names)
                error('No valid group names in the input...')
            end
        end
        
        function multi_att = add_col(obj,name_col,data_col,ind_col)
            % Adds a column in a multi_att object
            %
            % Parameters
            % ----------
            % name_col: string
            %       name of the column to add
            % data_col: array of double (nx1)
            %       data, the size has to be consistent with the object
            % ind_col: integer
            %       indice of the column to insert the new column
            %
            %
            % Output
            % ------
            % multi_att: :obj:
            %       multi_att object with added column
            %
            
            insert = @(a, x, n)cat(2,  x(:,1:n-1), a, x(:,n:end));
            
            obj.att_names = insert(name_col,obj.att_names,ind_col);
            obj = build_prop_from_names(obj);
            
            obj.att_data = insert(data_col,obj.att_data,ind_col);
            
            multi_att = obj;
        end
        
        function multi_att = add_group(obj,group_name,ids,data,cont,ind)
            % Adds a group to a multi_att
            %
            % Parameters
            % ----------
            % group_name: tring
            %       name of the new group
            % ids: array of string
            %       IDs of the columns
            %       if the object is considered continuous the values must 
            %       be in increasing order and of the same size as the IDs
            % data: matrix nxm
            %       data to be added, they must match  the object size and
            %       the id size
            % cont: boolean
            %       option to consider the IDs as continuous
            % ind: integer
            %       place to insert the new group
            %
            % Output
            % ------
            % multi_att: :obj:
            %       multi_att object with the added group
            %
            
            if nargin <= 4
                cont = 0;
                ind = length(obj.att_groups) + 1;
            end
            
            
            % If the ma contains the group, it is remove first
            if contains(group_name,obj.att_groups)
                % Remove the different columns
                sel_att = ~contains(obj.att_names,group_name);
                obj = get_from_indexes(obj,[],sel_att);
            end
            
            insert = @(a, x, n)cat(2,  x(:,1:n-1), a, x(:,n:end));
            
            obj.att_groups = insert(group_name,obj.att_groups,ind);
            
            if cont
                obj.att_types = insert({'C'},obj.att_types,ind);
                if isa(ids,'double')
                    obj.att_ids.(group_name) = ids;
                elseif isa(ids,'cellstr')
                    obj.att_ids.(group_name) = str2num(strrep(ids,'_','.'));
                end
            else
                obj.att_types = insert({'N'},obj.att_types,ind);
                obj.att_ids.(group_name) = ids;
            end
            obj = build_names_from_prop(obj);
            
            % build att_data
            [~,ind_data] = find(contains(obj.att_names,group_name),1) ;
            obj.att_data = insert(data,obj.att_data,ind_data);
            
            multi_att = obj;
        end
        
        function [multi_att, ind_nan_lin] = remove_NaN(obj)
            % Removes the columns and line that contain only NaNs
            %
            % Output
            % ------
            % multi_att: :obj:
            %       object without lines or columns of NaN
            % ind_nan_lin: array of boolean, optional
            %       lines deleted, used for the function of the subclass
            %       Multi_att_img
            %
            
            % Selection of the NaN columns and lines
            ind_nan_data = isnan(obj.att_data);
            ind_nan_lin = min(ind_nan_data,[],2)==1;
            ind_nan_col = min(ind_nan_data,[],1)==1;
            
            % Extraction of the matrix without NaN lines and columns
            multi_att = get_from_indexes(obj,~ind_nan_lin,~ind_nan_col);
        end
        
        function multi_att = to_subs_ref(obj,option,param)
            %
            %
            %
            
            if strcmp(option,'all')
                groups = obj.att_groups;
                att_names_tmp = [];
                for i = 1:length(groups)
                    len_group = length(obj.att_ids.(groups{i}));
                    len_group_sub = len_group*(len_group - 1)/2;
                    
                    % Init output var
                    ids1.(groups{i}) = [];
                    ids2.(groups{i}) = [];
                    names_grp = cell(1,len_group_sub);
                    names_grp(:) = {''};
                    
                    % Building all pairs indices and their names
                    id = 0;
                    for j = 1:len_group - 1
                        for k = j + 1:len_group
                            id = id + 1;
                            names_grp(id) = cellstr(strcat(groups{i},'_C_',num2str(round(obj.att_ids.(groups{i})(j))),'_',num2str(round(obj.att_ids.(groups{i})(k)))));
                            ids1.(groups{i}) = [ids1.(groups{i}), j];
                            ids2.(groups{i}) = [ids2.(groups{i}), k];
                        end
                    end
                    att_names_tmp = [att_names_tmp, names_grp];
                end
                
            elseif strcmp(option,'feats')
                % groups
                sep_und = cellfun(@(x) split(x,'_'),param,'un',0);
                group_names = cellfun(@(x) x{1},sep_und,'un',0);
                groups_un = unique(group_names);
                
                % reordering the group names according to the ma
                groups = [];
                for g = 1:length(obj.att_groups)
                    if contains(obj.att_groups(g),groups_un)
                        groups = [groups, obj.att_groups(g)];
                    end
                end
                
                % Building all pairs indices
                for i = 1:length(groups)
                    % Init
                    sep_und_gr = sep_und(strcmp(group_names,groups{i}));
                    ids1.(groups{i}) = str2double(cellfun(@(x) x{3} ,sep_und_gr,'un',0));
                    ids2.(groups{i}) = str2double(cellfun(@(x) x{4} ,sep_und_gr,'un',0));
                    ma_group = obj.get_from_group(groups{i});
                    
                    % Loop over unique wavelengths to find their indices
                    uni_all = unique([ids1.(groups{i}),ids2.(groups{i})]);
                    for j = 1:length(uni_all)
                        tmp_ids = get_ind_from_gr_id(ma_group,groups{i},uni_all(j));
                        ids1.(groups{i})(ids1.(groups{i}) == uni_all(j)) = tmp_ids;
                        ids2.(groups{i})(ids2.(groups{i}) == uni_all(j)) = tmp_ids;
                    end
                end
                
                % Names
                att_names_tmp = param;
            end
            
            % Creation of the data from the pairs of wavelengths
            att_data_tmp = [];
            for i = 1:length(groups)
                ma_group = obj.get_from_group(groups{i});
                
                % Init
                data_grp = NaN(size(obj.att_data,1),length(ids1.(groups{i})));
                
                % Creation of all substraction of wavelengths
                for j = 1:length(ids1.(groups{i}))
                    data_grp(:,j) = ma_group.att_data(:,ids1.(groups{i})(j)) - ma_group.att_data(:,ids2.(groups{i})(j));
                end
                att_data_tmp = [att_data_tmp, data_grp];
            end
            multi_att = Multi_att(att_names_tmp,att_data_tmp);
            multi_att = build_prop_from_names(multi_att);
        end
        
        function multi_att = rm_min_max_wls(obj,prop_name,min_val,max_val)
            % This function removes the values of a continuous outside of a
            % given range
            %
            % Parameters
            % ----------
            % obj: :obj:
            %       Multi_att object or subobject
            % prop_name: string
            %       Name of the property
            % min_val: double
            %       minimum value of the continuous property
            % max_val: double
            %       maximum values of the continuous property
            %
            %
            % Output
            % ------
            % multi_att: :obj:
            %       multi_att object without the wls outside of the bounds
            %
            
            mem_prop = ismember(obj.att_groups,prop_name);
            
            % Check if the parameters are consistent
            if max(mem_prop)==0
                error('The property was not found in the object, aborting...')
            elseif ~strcmp(obj.att_types,'C')
                error('This is not declared as a continuous variable, cannot proceed');
            end
            
            % Removes the values outside of the min max bounds
            min_val(isnan(min_val)) = -Inf;
            max_val(isnan(max_val)) = Inf;
            sel = obj.att_ids.(prop_name) <= max_val & obj.att_ids.(prop_name) >= min_val;
            
            sel_data = true(1,size(obj.att_data,2));
            sel_data(contains(obj.att_names,prop_name)) = sel;
            
            % Updates the parameters
            multi_att = get_from_indexes(obj,[],sel_data);
        end
        
        function multi_att = apply_change_spec(obj,prop_names,sp_trans,mean_opt)
            % Application of spectra transformation to the object
            %
            % Parameters
            % ----------
            % prop_names: array of string
            %       names of the groups to apply the transformation
            % sp_trans: array of integers
            %       index of the transformation
            % mean_opt: boolean
            %       option to compute the mean of the curves
            %
            % Output
            % ------
            % multi_att: :obj:
            %       multi_att object with the transformed group
            %
            
            multi_att = obj;
            
            % Option mean
            if mean_opt
                multi_att.att_data = mean(multi_att.att_data,1);
            end
            
            % Remove NaN
            multi_att = remove_NaN(multi_att);
            
            % Loop over the properties
            for i = 1:length(prop_names)
                mem = ismember(multi_att.att_groups,prop_names{i});
                if max(mem) == 1 && strcmp(multi_att.att_types(mem), 'C')
                    % Extraction the group
                    mul_att_group = get_from_group(multi_att,prop_names{i});
                    
                    % Creation of the 'change_spectra' object
                    grp_att_data_trans = Change_spectra_type(multi_att.att_ids.(prop_names{i}), mul_att_group.att_data);
                    
                    % Modification of the variable
                    grp_ind = contains(multi_att.att_names,prop_names{i});
                    multi_att.att_data(:,grp_ind) = grp_att_data_trans.apply(sp_trans(i) - 1);
                end
            end
        end
        
        function [ref_ma,sec_ma] = consis_multi_att(obj1,obj2,att_names,ref)
            % Make two objects Multi_att consistent
            % Either the common parameters of both objects are taken, 
            % or one object is the reference and its parameters are kept
            % and applied to the other object
            %
            % Parameters
            % ----------
            % obj1: :obj:
            %        first multi_att object
            % obj2: :obj:
            %        second multi_att object
            % att_names: cell str
            %       names of the groups to be kept
            % ref: boolean
            %        option to consider obj1 as the reference, we keep all
            %        its fields even if they don't exist in obj2 (set to
            %        NaN)
            %
            % Output
            % ------
            % ref_ma: :obj:
            %       multi_att reference object or first object rendered
            %       consistent with the second one
            % sec_ma: :obj:
            %       multi_att secondary object rendered
            %       consistent with the first one
            %
            
            % Remove NaN
            %obj1 = remove_NaN(obj1);
            %obj2 = remove_NaN(obj2);
            
            % Only keep the groups specified
            if ~isempty(att_names)
                obj1 = get_from_group(obj1,att_names);
                obj2 = get_from_group(obj2,att_names);
            end
            
            % Continuous values made consistent
            [ref_ma,sec_ma] = consis_def_contin(obj1,obj2,ref);
            
            if ref
                % The reference is kept as is
                ref_ma = obj1;
                
                % Fills the non existing columns in the sec_ma with NaN
                data = NaN(size(obj2.att_data,1),1);
                for i = 1:length(ref_ma.att_names)
                    name = ref_ma.att_names{i};
                    if ~contains(name,sec_ma.att_names)
                        sec_ma = add_col(sec_ma,name,data,find(contains(ref_ma.att_names,name),1));
                    end
                end
                
                % Remove the different columns
                sel_att_sec = ismember(sec_ma.att_names,ref_ma.att_names);
                sec_ma = get_from_indexes(sec_ma,[],sel_att_sec);
            else
                % Numerical variables
                % Selection of columns with the same name
                sel_att_ref = ismember(ref_ma.att_names,sec_ma.att_names);
                sel_att_sec = ismember(sec_ma.att_names,ref_ma.att_names);
                
                % Keeps columns that are the same in both objects
                ref_ma = get_from_indexes(ref_ma,[],sel_att_ref);
                sec_ma = get_from_indexes(sec_ma,[],sel_att_sec);
            end
            
            % Filling the missing values
            ref_ma.att_data = fillmissing(ref_ma.att_data,'linear',2);
            sec_ma.att_data = fillmissing(sec_ma.att_data,'linear',2);
        end
        
        function [ref_ma,sec_ma] = consis_def_contin(obj1,obj2,ref)
            % Makes the continuous variables consistent between two multi
            % att
            % It takes the variable with the highest resolution as a
            % reference, unless specified otherwise in the parameter
            %
            % Parameters
            % ----------
            % ref: boolean
            %       option to select the first object as the reference
            %
            % Output
            % ------
            % ref_ma: :obj:
            %       multi_att reference object or first object rendered
            %       consistent with the second one
            % sec_ma: :obj:
            %       multi_att secondary object rendered
            %       consistent with the first one
            %
            
            
            group1_C = obj1.att_groups(strcmp(obj1.att_types,'C'));
            group2_C = obj2.att_groups(strcmp(obj2.att_types,'C'));
            common_groups  = group1_C(ismember(group1_C,group2_C));
            
            for g = 1:length(common_groups)
                grp = common_groups{g};
                c_var1 = obj1.att_ids.(grp);
                c_var2 = obj2.att_ids.(grp);
                
                % Removing intervals not shared
                if (min(c_var1) ~= min(c_var2)) || (max(c_var1) ~= max(c_var2))
                    min_val = max(min(c_var1),min(c_var2));
                    max_val = min(max(c_var1),max(c_var2));
                    obj1 = rm_min_max_wls(obj1,grp,min_val,max_val);
                    obj2 = rm_min_max_wls(obj2,grp,min_val,max_val);
                end
                
                % Update the continuous values
                c_var1 = obj1.att_ids.(grp);
                c_var2 = obj2.att_ids.(grp);
                
                % Check if the variables are not consistent
                if ~isequal(c_var1,c_var2)
                    if ref 
                        % If the option is activated the first is
                        % considered as the reference
                        ref_ind = 1;
                    else
                        % Check resolutions
                        res1 = (max(c_var1) - min(c_var1))/(length(c_var1) - 1);
                        res2 = (max(c_var2) - min(c_var2))/(length(c_var2) - 1);
                        
                        % The lowest resolution is considered as the reference
                        if res1 < res2
                            ref_ind = 1;
                        else
                            ref_ind = 2;
                        end
                    end
                    
                    % Interpolation of the second continuous values on the
                    % points of the reference
                    if ref_ind == 1
                        % Case where the first variable is the reference
                        data_interp = NaN(size(obj2.att_data,1),length(obj1.att_ids.(grp)));
                        
                        % Skip data containing only NaNs
                        ind_tot = (1:size(obj2.att_data,1))';
                        ind_no_nan = ind_tot(sum(isnan(obj2.att_data),2)~=size(obj2.att_data,2));
                        for i = 1:length(ind_no_nan)
                            data_interp(ind_no_nan(i),:) = interp1(obj2.att_ids.(grp),...
                                                           obj2.att_data(ind_no_nan(i),contains(obj2.att_names,grp)),...
                                                           obj1.att_ids.(grp),...
                                                           'spline');
                        end
                        obj2 = add_group(obj2,grp,obj1.att_ids.(grp),data_interp,1,g);
                    elseif ref_ind == 2
                        % Case where the second variable is the reference
                        data_interp = NaN(size(obj2.att_data,1),length(obj1.att_ids.(grp)));
                        
                        % Skip data containing only NaNs
                        ind_tot = (1:size(obj1.att_data,1))';
                        ind_no_nan = ind_tot(sum(isnan(obj1.att_data),2)~=size(obj1.att_data,2));
                        for i = 1:length(ind_no_nan)
                            data_interp(ind_no_nan(i),:) = interp1(obj1.att_ids.(grp),...
                                                                   obj1.att_data(ind_no_nan(i),contains(obj1.att_names,grp)),...
                                                                   obj2.att_ids.(grp),...
                                                                   'spline');
                        end
                        obj1 = add_group(obj1,grp,obj2.att_ids.(grp),data_interp,1,g);
                    end
                end
            end
            ref_ma = obj1;
            sec_ma = obj2;
        end
        
        function multi_att = plus(obj,ma_add)
            % Addition of two multi att objects
            % The first one is consisdered as the reference unless it's
            % empty
            
            if isempty(obj.att_data)
                multi_att = ma_add;
            else
                [ref_cons,ma_cons] = consis_multi_att(obj,ma_add,[],1);
                multi_att = ref_cons;                
                multi_att.att_data = cat(1,ref_cons.att_data,ma_cons.att_data);
            end
        end
    end
end

function index = index_close_value_reg_array(array_val, value)
    % Gives the index of a given wavelength in a file
    %
    % Parameters
    % ----------
    % array_val : array of double
    %       regular values
    % value : double
    %       value to be found in array_val
    %
    % Output
    % ------
    % index : integer
    %       index of the values in the array_val
    %
    
    % Finding the indexes
    dist = abs(array_val - value);
    step = (max(array_val) - min(array_val))/(length(array_val) - 1); % average step
    if min(dist) < step*200
        [~,index] = min(dist);
    else
        error('Could not find a element close enough')
    end

end
