classdef Multi_att_Img < Multi_att
    
    properties
        coord
        axis
        surf
        sel
        ind_img
        img_rgb
    end
    
    methods
        function obj = Multi_att_Img(att_names,att_data,coord,surf,sel,ind_img)
            % Creation of the object Multi_att_Img
            %
            if nargin == 0
                att_names = [];
                att_data = [];
                coord = [];
                surf = [];
                sel = [];
                ind_img = [];
            end
            
            obj = obj@Multi_att(att_names,att_data);
            
            obj.coord = coord;
            obj.surf = surf;
            obj.sel = sel;
            obj.ind_img = ind_img;
            
            if nargin > 0
                obj = update_axis_val(obj);
            end
        end
        
        function multi_att_img = update_classes(obj,new_classes)
            % Updates the indices of the image with new ones
            % the attribute att_data must be updated beforehand 
            % 
            % Parameters
            % ----------
            % obj: :obj:
            %       Multi_att_img
            % new_classes: array of int
            %       value of the new classes, the vector must have the same
            %       length as the old one 
            %
            % Output
            % ------
            % multi_att_img: :obj:
            %       multi_att_img containing the new_classes
            %
            
            multi_att_img = obj;
            
            old_classes = unique(multi_att_img.ind_img);
            old_classes(isnan(old_classes)) = [];
            uni_new = unique(new_classes);
            uni_new(isnan(uni_new)) = [];
            
            % Check the consistency between the new variable and the object
            % 
            if (length(new_classes) ~= length(old_classes)) || (length(uni_new) ~= size(multi_att_img.att_data,1))
                error('The variable new_classes is not consistent with the multi_att_img object.');
            end
            
            input_img = multi_att_img.ind_img;
            for i = 1:length(old_classes)
                input_img(input_img==i) = -double(new_classes(i)); 
            end
            multi_att_img.ind_img = -input_img;
        end
        
        function multi_att = remove_NaN(obj)
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
            
            [multi_att, ind_nan_lin] = remove_NaN@Multi_att(obj);
            
            % Update the ind_img
            %{
            new_classes = [];
            ite = 1;
            for i = 1:length(ind_nan_lin)
                if ind_nan_lin(i)
                    new_classes = [new_classes, NaN];
                else
                    new_classes = [new_classes, ite];
                    ite = ite+1;
                end
            end
            multi_att = update_classes(multi_att,new_classes);
            %}
        end
        
        function multi_att = get_from_indexes(obj,ind_lin,ind_col,skip_update)
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
            
            if nargin <= 3 
                skip_update = 0;
            end
            
            multi_att = get_from_indexes@Multi_att(obj,ind_lin,ind_col);
            
            % Update the ind_img if the ind_lin is not empty or is not
            % containing only ones (all row selected)
            if ~isempty(ind_lin) && ~(size(obj.att_data,1) == length(ind_lin) && min(ind_lin) == 1) && ~skip_update
                if size(obj.att_data,1) == length(ind_lin)
                    ind_nan_lin = ~ind_lin;
                else
                    ind_nan_lin = ~ismember(1:size(obj.att_data,1),ind_lin);
                end
                
                % Creation of the vector of indices of replacement
                % containing Nans where the row is deleted, and a updated
                % indice otherwise
                new_classes = [];
                ite = 1;
                for i = 1:length(ind_nan_lin)
                    if ind_nan_lin(i)
                        new_classes = [new_classes, NaN];
                    else
                        new_classes = [new_classes, ite];
                        ite = ite + 1;
                    end
                end
                
                % Update of the change 
                multi_att = update_classes(multi_att,new_classes);
            end
        end
        
        function img = load_envi_files(obj, dir_path)
            % Loads envi or rgb files into Multi_att objects
            %
            % Parameters
            % ----------
            % dir_path: string
            %       path of the folder that contains the envi files
            %
            
            % Get the type of data
            pathparts = strsplit(dir_path,'\');
            dir = pathparts{end};
            
            % Load the data from the file
            load_data = Load_Data_File(fullfile(pathparts{1:end-1}),[],lower(dir));
            names.(dir) = load_data.data_names;
            data.(dir) = load_data.data_loaded;
            att_names = [names.(dir)];
            
            % Reverse the axis for the SWIR - due to a difference in the
            % acquisition 
            if strcmp(dir,'SWIR')
                data.('SWIR') = flip(data.('SWIR'),2);
            end
            
            % Creation of fictitious coordinate variables
            [N_line,N_samp, ~] = size(data.(dir));
            coord_def.x = repmat(linspace(1,N_samp,N_samp), [N_line,1]);
            coord_def.y = repmat(linspace(1,N_line,N_line)', [1,N_samp]);
            coord_def.z = zeros([N_line,N_samp]);
            
            % Binning - Re-Sampling
            x_samp = 1;
            y_samp = 1;
            
            % Application of the binning
            coord_def.x = coord_def.x(1:y_samp:end, 1:x_samp:end);
            coord_def.y = coord_def.y(1:y_samp:end, 1:x_samp:end);
            coord_def.z = coord_def.z(1:y_samp:end, 1:x_samp:end);
            data.(dir) = data.(dir)(1:y_samp:end, 1:x_samp:end, :);
                        
            % Definition of the selection
            sel_def = ones(size(coord_def.x));
            
            % Application of the selection to the Coordinates
            coord_def.x(~sel_def) = NaN;
            coord_def.y(~sel_def) = NaN;
            coord_def.z(~sel_def) = NaN;
            
            % Fusion of the data
            nb_ind = length(coord_def.x(:));
            [N_line,N_samp, N_C] = size(data.(dir));
            att_data = reshape(data.(dir),[N_samp*N_line,N_C]);
            
            % All indices of the image
            ind = 1:nb_ind;
            
            % Properties of the object
            ind_img_def = reshape(ind,size(coord_def.x));   % flip + transpose : to mimic the linear indexing of Matlab
            
            % TEMP: constant surface element, use the real surface obtained
            % by proper calculation
            id_x_mid = round(size(coord_def.x,2)/2,0);
            id_y_mid = round(size(coord_def.x,1)/2,0);
            surf_el = (coord_def.x(id_y_mid,id_x_mid+1) - coord_def.x(id_y_mid,id_x_mid))*(coord_def.y(id_y_mid+1,id_x_mid) - coord_def.y(id_y_mid,id_x_mid)); % temp, use the calculated surface
            surf_def = surf_el*ones(size(coord_def.x));
            
            % Definition of the object
            img = Multi_att_Img(att_names,att_data,coord_def,surf_def,sel_def,ind_img_def);
            img = build_prop_from_names(img);
        end
        
        function img = load_fused_data(obj, dir_path_acquisition,dir_path_fusion)
            % Loads fused data into Multi_att objects
            %
            % Parameters
            % ----------
            % dir_path: string
            %       path of the folder that contains the fused data
            %
            global cfg 
            
            % Binning - Re-Sampling
            x_samp = cfg.corsen_x;
            y_samp = cfg.corsen_y;
            
            % Getting the files for the fusion
            [indexes, coord_def, att_groups] = ... 
                load_fusion_files(dir_path_fusion,x_samp,y_samp);
            
            % size of the image
            [sz_x,sz_y] = size(coord_def.x);
            
            % Get the selection of the data 
            [sel_def, min_max_idx_xy] = get_selection(coord_def);
            
            % Application of the selection to remove indexes outside
            for d = 1:length(att_groups)
                indexes.(att_groups{d}).index_x(~sel_def) = 0;
                indexes.(att_groups{d}).index_y(~sel_def) = 0;
            end
            
            % Application of the selection to the Coordinates
            coord_def.x(~sel_def) = NaN;
            coord_def.y(~sel_def) = NaN;
            coord_def.z(~sel_def) = NaN;
            
            % Linear filling of the missiing values
            coord_def.x = fillmissing(coord_def.x,'linear',2);
            coord_def.y = fillmissing(coord_def.y,'linear',2);
            coord_def.z = fillmissing(coord_def.z,'linear',2);
            
            % All indices of the image
            nb_ind = sz_x*sz_y;
            ind = 1:nb_ind;
            
            % Load data and merge them into the data file
            att_data = [];
            att_names = [];
            for d = 1:length(att_groups)
                [att_names_tmp, att_data_tmp] = ...
                    create_att_fus(dir_path_acquisition,dir_path_fusion,...
                                   att_groups{d}, indexes.(att_groups{d}),...
                                   sel_def, x_samp, y_samp,ind);
                att_names = [att_names, att_names_tmp];
                att_data = [att_data, att_data_tmp];
            end
            
            surf_el = nanmedian(abs(diff(coord_def.x(~isnan(coord_def.x)))))...
                *nanmedian(abs(diff(coord_def.x(~isnan(coord_def.y))))); % temp, use the calculated surface
            surf_def = surf_el*ones([sz_x,sz_y]);
            
            % Definition of the object
            img = Multi_att_Img(att_names,att_data,coord_def,surf_def,sel_def,reshape(ind,[sz_x,sz_y]));
            img = build_prop_from_names(img);
            
            % Application of the square selection (sides are removed)
            img = img_modify(img,1,min_max_idx_xy);
        end
        
        function img = update_axis_val(obj)
            % Updates the value of the axis attributes based on the
            % coordinates
            %
            %
            
            % Takes the middle row and column of the coordinates
            coord_mid = size(obj.coord.x);
            coord_mid = round(coord_mid/4,0);
            raw_axis_x = obj.coord.x([coord_mid(1), 2*coord_mid(1), 3*coord_mid(1)],:);
            raw_axis_y = obj.coord.y(:,[coord_mid(2), 2*coord_mid(2), 3*coord_mid(2)]);
            raw_axis_x = nanmean(raw_axis_x, 1);
            raw_axis_y = nanmean(raw_axis_y, 2);
            
            % Keeps the min and max points
            min_x = nanmin(raw_axis_x);
            min_y = nanmin(raw_axis_y);
            max_x = nanmax(raw_axis_x);
            max_y = nanmax(raw_axis_y);
            
            raw_axis_x(raw_axis_x ~= min_x & raw_axis_x ~= max_x) = NaN;
            raw_axis_y(raw_axis_y ~= min_y & raw_axis_y ~= max_y) = NaN;
            
            % Creates a line between the min and max points
            % Keeps the order increasing
            interp_x = fillmissing(raw_axis_x,'linear');
            obj.axis.x = interp_x;
            
            interp_y = fillmissing(raw_axis_y,'linear');
            obj.axis.y = interp_y;
                        
            img = obj;
        end
        
        function rgb_img = rgb_image_builder(obj, R_WL, G_WL, B_WL, norm_opt, output_file)
            % Shows a hyperspectral image with three RGB channels
            % If necessary saves it
            %
            % Parameters
            % ----------
            % R_WL : wavelength attributed to Red
            % G_WL : wavelength attributed to Green
            % B_WL : wavelength attributed to Blue
            % norm_opt : int
            %           0 : no normalisation
            %           1 : linear normalisation
            %           2 : histeq normalisation
            % output_file : name of the output file to write
            %
            % Output
            % ------
            % rgb_img : rgb image
            %
            global cfg
            
            WLs = [cellstr(num2str(R_WL)), cellstr(num2str(G_WL)), cellstr(num2str(B_WL))];
            rgb_indexes = ones(1,3);
            
            for w = 1:3
                if contains('rgb',WLs{w})
                    ind_gr = contains(obj.att_groups,'RGB');
                    if max(ind_gr)
                        if strcmp(WLs{w},'r')
                            rgb_indexes(w) = get_ind_from_gr_id(obj,'RGB','R');
                        elseif strcmp(WLs{w},'g')
                            rgb_indexes(w) = get_ind_from_gr_id(obj,'RGB','G');
                        elseif strcmp(WLs{w},'b')
                            rgb_indexes(w) = get_ind_from_gr_id(obj,'RGB','B');
                        end
                    else
                        fprinf('Error: File does not contain RGB data.')
                    end
                elseif str2double(WLs{w}) < cfg.lim_vnir_swir
                    ind_gr = contains(obj.att_groups,'VNIR');
                    if max(ind_gr)
                        rgb_indexes(w) = get_ind_from_gr_id(obj,'VNIR',str2double(WLs{w}));
                    end
                elseif str2double(WLs{w}) > cfg.lim_vnir_swir
                    ind_gr = contains(obj.att_groups,'SWIR');
                    if max(ind_gr)
                        rgb_indexes(w) = get_ind_from_gr_id(obj,'SWIR',str2double(WLs{w}));
                    end
                else
                    error('wrong input value to build the RGB image')
                end
            end
            % Display the false-color image of the sample
            
            % RGB image building
            rgb_img = reshape(obj.att_data(:,rgb_indexes),[size(obj.ind_img), 3]);
            
            % Normalisation if required
            if exist('norm_opt')
                if norm_opt == 1
                    n = size(rgb_img);
                    minX = repmat(min(rgb_img), [n(1), 1]);
                    maxX = repmat(max(rgb_img), [n(1), 1]);
                    rgb_img = (rgb_img - minX)./(maxX - minX);
                elseif norm_opt == 2
                    rgb_img = histeq(rgb_img);
                end
            end
            
            % Set the 'value' to 1 for all
            hsv_img = rgb2hsv(rgb_img);
            rgb_img = hsv2rgb(hsv_img);
            
            % Saving it as a .png if required
            if (exist('output_file') ~= 0)
                if ~isempty(output_file)
                    imwrite(rgb_img, output_file);
                end
            end
        end
        
        function img = create_img(obj, R_WL, G_WL, B_WL, norm_opt, output_file)
            %{Creation of the attributes img_rgb
            %}
            
            img = obj;
            
            % Figure creation
            fig_img = figure('visible','off','Name','fig_img');
            coord_x = img.axis.x;
            coord_y = img.axis.y;
            
            % Axes creation
            h1 = axes;
            cla;
            hold on;

            if coord_x(1) < coord_x(2)
                h1.XDir = 'normal';
                h1.XAxis.TickValues = coord_x;
            else
                h1.XDir = 'reverse';
                h1.XAxis.TickValues = flip(coord_x);
            end
            h1.XTickMode = 'auto';
            if coord_y(1) < coord_y(2)
                h1.YDir = 'reverse';
                h1.YAxis.TickValues = coord_y;
            else
                h1.YDir = 'normal';
                h1.YAxis.TickValues = flip(coord_y);
            end
            h1.YTickMode = 'auto';

            % Items of the axes
            rgb_img = rgb_image_builder(img, R_WL, G_WL, B_WL, norm_opt, output_file);
            set(imagesc(rgb_img),'XData', coord_x,'YData',coord_y);
            
            h1.XLim = [min(coord_x), max(coord_x)];
            h1.YLim = [min(coord_y), max(coord_y)];
            h1.DataAspectRatio = [1,1,1];
            
            % Change the figure's size to fit the axis
            img_sz = size(fig_img.Children.Children.CData);
            pos = fig_img.Position;
            pos(4) = pos(3)/img_sz(2)*img_sz(1);
            set(fig_img,'Position',pos)
            
            % Saving the file
            img.img_rgb = fig_img;            
        end
        
        function data_img = reshape_att_data(obj)
            % Reshape the att_data back into the image frame
            %
            %
            
            data_img = reshape(obj.ind_img,[length(obj.ind_img(:)),1]);
            simg = length(data_img);
            data_img(:,2) = 1:simg;
            data_img = sortrows(data_img,1);
            
            [sx,sy] = size(obj.att_data);
            data_img(1:sx,3:2 + sy) = obj.att_data;
            data_img(sx+1:end,3:2 + sy) = NaN(simg-sx,sy);
            
            data_img = sortrows(data_img,2);
            data_img = data_img(:,3:end);
        end
        
        function img = img_modify(obj,corsen_rate,min_max_idx_xy)
            % Edition of Hyper Spectral Data
            % It allows to modify the spectrum type and to downscale the
            % image.
            %
            % Parameters
            % ----------
            % corsen_rate: int
            %       corsening rate
            % min_max_idx_xy: array of 4 doubles
            %       min and max indexes of X and Y
            %
                        
            % Corsening the image
            if ~isempty(min_max_idx_xy)
                xMin = min_max_idx_xy(1);
                xMax = min_max_idx_xy(2);
                yMin = min_max_idx_xy(3);
                yMax = min_max_idx_xy(4);
            else
                size_hs_data = size(obj.ind_img);
                xMin = 1;
                xMax = size_hs_data(2);
                yMin = 1;
                yMax = size_hs_data(1);
            end
            xRate = corsen_rate;
            yRate = corsen_rate;
            
            % Application of Crop_downscale - TODO : test if necessary
            obj.ind_img = obj.ind_img(yMin:yRate:yMax, xMin:xRate:xMax);
            obj.surf = obj.surf(yMin:yRate:yMax, xMin:xRate:xMax);
            obj.sel = obj.sel(yMin:yRate:yMax, xMin:xRate:xMax);
            
            % Removes the data
            ind_keep = reshape(obj.ind_img,[size(obj.ind_img,1)*size(obj.ind_img,2),1]);
            obj.att_data = obj.att_data(ind_keep,:);
            
            % Update the indexes 
            obj.ind_img = reshape(1:length(obj.ind_img(:)),size(obj.ind_img));
            
            % Update the coordinates
            obj.coord.x = obj.coord.x(yMin:yRate:yMax, xMin:xRate:xMax);
            obj.coord.y = obj.coord.y(yMin:yRate:yMax, xMin:xRate:xMax);
            obj.coord.z = obj.coord.z(yMin:yRate:yMax, xMin:xRate:xMax);
            
            % Update the axis values
            img = update_axis_val(obj);
        end
        
        function image_sc = make_img_from_var(obj,var)
            %
            %
            % Parameters
            % ----------
            % var: array of double
            %       values to be 
            %       the size of the vector should be the same as the number
            %       of index in the obj
            %
            
            % Apply attribute to map
            image_sc = obj.ind_img;
            indexes = unique(image_sc);
            indexes(isnan(indexes)) = [];
            
            % Loop over all indexes that are replaced
            if length(indexes) == length(image_sc(:))
                image_sc = reshape(var,size(obj.ind_img));
            else
                nb_col = length(indexes);
                for i = 1:nb_col
                    idx = indexes(i);
                    image_sc(image_sc==idx) = var(i);
                end
            end
        end
        
        function img = to_subs_ref(obj,option,param)
            %
            %
            %
            ma = to_subs_ref@Multi_att(obj,option,param);
            img = copy(ma,'img');
            prop = setdiff(properties(Multi_att_Img),properties(Multi_att));
            for i = 1:length(prop)
                img.(prop{i}) = obj.(prop{i});
            end
        end
    end
end

function [idx_sta,idx_end] = find_idx_longest_interval(vec)
%
% Returns the indexes of the longest interval of a boolean array
%

if size(vec,2)==1
    vec = vec';
end

vec = diff([0,vec]);

% Finds the starts and ends of segments
idx_sta_all = strfind(vec, 1);
idx_end_all = strfind(vec, -1);

% Add values for segments touching the border
% Empty values
if isempty(idx_sta_all)
    idx_sta_all = 1;
end
if isempty(idx_end_all)
    idx_end_all = length(vec);
end
% Inconsistencies
if idx_sta_all(1) > idx_end_all(1)
    idx_sta_all = [1, idx_sta_all];
end
if idx_sta_all(end) > idx_end_all(end)
    idx_end_all = [idx_end_all, length(vec)];
end

% Finds the longest segment
len = idx_end_all - idx_sta_all;
[~,id] = max(len);

idx_sta = idx_sta_all(id);
idx_end = idx_end_all(id);
end

function [sel_def, min_max_idx_xy] = get_selection(coord_def)

    % size of the image
    [sz_x,sz_y] = size(coord_def.x);
    
    % Selection of the pixels with data 
    % based on the coordinates (!! definition can change, to remove when all image have their filter)
    val_x_unk = [-25, 25];
    val_y_unk = [0, Inf];
    sel_def_init = (coord_def.x > val_x_unk(1) &... 
                    coord_def.x < val_x_unk(2)) &...
                   (coord_def.y > val_y_unk(1) &...
                    coord_def.y < val_y_unk(2));

    % Square selection : removes the useless sides of the image
    opt_sq = 0;
    if opt_sq
        se = strel('disk',10);
        sel_def_init_closed = imclose(sel_def_init,se);
        x_middle = sel_def_init_closed(:, round(size(sel_def_init_closed,2)/2));
        y_middle = sel_def_init_closed(round(size(sel_def_init_closed,1)/2), :);
        [val_y(1), val_y(2)] = find_idx_longest_interval(x_middle);
        [val_x(1), val_x(2)] = find_idx_longest_interval(y_middle);
    else
        [indy,indx] = find(sel_def_init == 1);

        val_x = [min(indx),max(indx)];
        val_y = [min(indy),max(indy)];
    end
    min_max_idx_xy = [val_x, val_y];

    % Definition of the selection
    try
        tmp = load(fullfile(dir_path_fusion,'shadow_filter.mat'),...
                    'shadow_filter');
        sel_square_def = rot90(tmp.shadow_filter,-1);
        sel_square_def = sel_square_def(1:y_samp:end,1:x_samp:end);
    catch
        sel_square_def = zeros([sz_x,sz_y]);
        sel_square_def(val_y(1):val_y(2),val_x(1):val_x(2)) = 1;
        sel_square_def = logical(sel_square_def);
    end

    sel_def = sel_def_init & sel_square_def;
end

function [att_names, att_data] = create_att_fus(dir_path_acquisition, dir_path_fusion, group, indexes, sel_def, x_samp, y_samp, ind)
    load_data = Load_Data_File(dir_path_acquisition, dir_path_fusion, group);
    att_names = load_data.data_names;
    att_data_grp = load_data.data_loaded;
    
    if strcmp(group,'rgb')
        att_data = att_data_grp(1:y_samp:end, 1:x_samp:end,:);
        sel_def_3 = repmat(sel_def,[1, 1, 3]);
        att_data(~sel_def_3) = NaN;

        % RGB case : the indexes are not used, it is assumed that the
        % reconstructed data exists in the folder
        size_rgb = size(att_data);
        if size_rgb(3) > 3
            att_data = att_data(:, :, end - 2:end);
        end
        att_data(:,1:3) = reshape(att_data,...
                                    [size_rgb(1)*size_rgb(2),3]);
    else
        % Indices in the selection
        len_img = length(ind);
        att_len = size(att_data_grp,3);
        data_grp_sz = size(att_data_grp);
        
        % Get the indexes of the spectra in the fusion image
        ind_x_resh = reshape(indexes.index_x, 1, len_img);
        ind_y_resh = reshape(indexes.index_y, 1, len_img);
        
        % Making sure the indexes are within the allowed bounds
        ind_x_resh(ind_x_resh > data_grp_sz(2)) = data_grp_sz(2);
        ind_x_resh(ind_x_resh < 1) = 1;
        ind_y_resh(ind_y_resh > data_grp_sz(1)) = data_grp_sz(1);
        ind_y_resh(ind_y_resh < 1) = 1;
        
        % Get the indexes
        glob_ind_len = sub2ind(data_grp_sz,...
                               reshape(repmat(ind_y_resh,att_len,1), 1, att_len*len_img),... % y index
                               reshape(repmat(ind_x_resh,att_len,1), 1, att_len*len_img),... % x index 
                               repmat(1:att_len,1,len_img)); % z index 
        ind_fin = reshape(glob_ind_len,att_len,len_img)';
        att_data = att_data_grp(ind_fin);
        
        % Set the Nan where the mask is applied
        ind_sel = ind(~sel_def);
        att_data(ind_sel,:) = NaN(length(ind_sel),att_len);
    end
end

