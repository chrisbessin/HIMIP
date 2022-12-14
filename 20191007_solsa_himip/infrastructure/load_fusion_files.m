function [indexes,coord_def,att_groups] = load_fusion_files(dir_path_fusion,x_samp,y_samp,check)
    %
    %
    % Parameters
    % ----------
    % dir_path_fusion: string
    %       path of the fusion directory
    % x_samp: int
    %       resampling along X
    % y_samp: int
    %       resampling along Y
    % check: boolean
    %       check if the fusion file load is possible
    %
    
    if nargin < 4
        check = false;
    end
    
    % Load the coordinates
    load(fullfile(dir_path_fusion,'XYZ_Data'),'X','Y','Z')
    
    % Rotation of the coordinates
    X = rot90(X,-1);
    Y = rot90(Y,-1);
    Z = rot90(Z,-1);
    
    % Loading the indexes
    dirs = [{'VNIR'}, {'SWIR'}, {'XRF'}];
    for d = 1 : length(dirs)
        name = lower(dirs(d));
        file_tmp = fullfile(dir_path_fusion,'S2P',dirs(d),'index_x+index_y.mat');
        if isfile(file_tmp{1})
            indexes.(name{1}) = load(file_tmp{1});
            indexes.(name{1}).index_x = rot90(indexes.(name{1}).index_x,-1);
            indexes.(name{1}).index_y = rot90(indexes.(name{1}).index_y,-1);
        end
    end
    
    % Test if indexes exist in the folder 
    if exist('indexes','var')
        att_groups = fieldnames(indexes);
    else
        error('Not a proper directory, choose a directory with fused data')
    end
    
    if ~check
        % Coordinates
        coord_def.x = X(1:y_samp:end, 1:x_samp:end);
        coord_def.y = Y(1:y_samp:end, 1:x_samp:end);
        coord_def.z = Z(1:y_samp:end, 1:x_samp:end);

        % Indexes of the properties
        for d = 1:length(att_groups)
            tmp_mat_x = getfield(indexes,att_groups{d},'index_x');
            tmp_mat_y = getfield(indexes,att_groups{d},'index_y');
            indexes.(att_groups{d}).index_x = tmp_mat_x(1:y_samp:end, 1:x_samp:end);
            indexes.(att_groups{d}).index_y = tmp_mat_y(1:y_samp:end, 1:x_samp:end);
        end
    end
end

