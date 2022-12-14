function run_HyperSpectral_local(varargin) %CoreScan_Experiment_uuid,Fusion_Experiment_uuid,CoreScan_Config_ID,HyperSpectral_Config_ID,DB_URI,CoreScan_Table_Name,Fusion_Table_Name,HyperSpectral_Experiment_uuid,check)
%{
7. GUI runs ID2A HyperSpectral: 
    Check_Inputs_HyperSpectral(
			CoreScan_Experiment_uuid: string, 
			Fusion_Experiment_uuid: string, 
			CoreScan_Config_ID: int, 
			HyperSpectral_Config_ID: int, 
			DB_URI: string,
			CoreScan_Table_Name: string /*optional, default: standard name in the Solsa Sample DB Schema*/,
			Fusion_Table_Name: string /*optional, default: standard name in the Solsa Sample DB Schema*/,
			HyperSpectral_Experiment_uuid: string,
            check: boolean /*optional: checks if the run is possible*/)
    
	Response:
			Check_Config: string /*PASS, FAIL*/
			Check_Config_Info: string
			Check_File: string /*PASS, FAIL*/
			Check_File_Info: string
			Check_Data: string /*PASS, FAIL*/
			Check_Data_Info: string
			Check_Global: string /*PASS, FAIL*/
			Estimated_Time: int /*seconds*/
			Results: string /*PASS, FAIL*/
			Results_info: string
    
	Comments:
			Info fields are strings describing the error, format to be defined

    cmd:
/d/SOLSA/HIMIP/20191007_solsa_himip/interface/run_HyperSpectral/for_redistribution_files_only/run_HyperSpectral.exe CoreScan_Experiment_uuid 25d60a86-2be1-4948-b1ca-5a0c29a7b2a7 Fusion_Experiment_uuid 7e1052ab-036b-4818-b814-a30e2bad3bbc CoreScan_Config_ID 1 HyperSpectral_Config_ID 23 DB_URI https://solsa.crystallography.net/db/test_samples/ HyperSpectral_Experiment_uuid 09d15fe4-f8a8-4100-98f1-0bfbb9fb34c1 check 0
/d/SOLSA/HIMIP/20191007_solsa_himip/interface/run_HyperSpectral/for_redistribution_files_only/run_HyperSpectral.exe
acquisition_path D:\SOLSA\Quick_interface\data\results\ER-NC00-0072
config_path D:\SOLSA\HIMIP\Hyperspectral_GUI\Config.m check 0
%}


warning ('off','all');

% Init values of the arguments

v.acquisition_path = 'D:\SOLSA\Quick_interface\data\results\ER-NC00-0072_1';
v.config_path = 'D:\SOLSA\HIMIP\Hyperspectral_GUI\Config.m';
v.check = '0';

[~,SolsaID,~] = fileparts(v.acquisition_path);
SolsaID = 'ER-NC00-0072';

% Check if the number of argument is even (keyword + arg)
known_vars = fieldnames(v);
if mod(nargin,2) ~= 0
    error('Expecting an even number of arguments');
end

% Replace init values by the arg if it exists
for idx = 1 : 2 : nargin - 1
    if ~ismember(varargin{idx}, known_vars)
        error('Argument "%s" is not a variable name I know', varargin{idx});
    end
    val = varargin{idx + 1};
    v.(varargin{idx}) = val;
end

% Changes the format of the check if it was entered with command line 

if ischar(v.check)
    v.check = str2num(v.check);
end


try
    % Get the config
    fileID = fopen(v.config_path,'r');
    config_txt = fscanf(fileID,'%c');
    eval(config_txt)
    
    % Check Config
    names_expected_cfg = [
        {'dir_proj'               };...
        {'dir_path_library'       };...
        {'path_library'           };...
        {'dir_data_fus'           };...
        {'dir_data'               };...
        {'file_data'              };...
        {'darkref_tag'            };...
        {'whiteref_tag'           };...
        {'red_default'            };...
        {'green_default'          };...
        {'blue_default'           };...
        {'min_wl_vnir'            };...
        {'max_wl_vnir'            };...
        {'min_wl_swir'            };...
        {'max_wl_swir'            };...
        {'lim_vnir_swir'          };...
        {'average_density_samples'};...
        {'diam_sample'            };...
        {'mass_ROI'               };...
        {'spectrum_types'         };...
        {'corsen_x'               };...
        {'corsen_y'               };...
        {'att_names'              };...
        {'nb_cluster'             };...
        {'process_sample'         };...
        {'nb_tree'                };...
        {'print_times'            };...
        {'elements'               };...
        {'elements_masse_mol'     };...
        {'elts_mass_dict'         };...
        {'res_x_swir'             };...
        {'res_y_swir'             }];
    
    ind_check_config = all(ismember(names_expected_cfg, fieldnames(cfg)));
    
    % config
    if ind_check_config
        print_check_config = 'PASS\nConfig files are found and complete\n';
    else
        print_check_config = 'FAIL\nConfig files are either not found or incomplete\n';
    end
catch
    print_check_config = 'FAIL\nCould not load the config from the DB.\n';
end

% Check File
% Corescan
try
    % Get global info from the DB
    % Get experiment corescan
    path_corescan_exp = v.acquisition_path;
    
    % Checking the files
    files_to_check_corescan = [
        {[path_corescan_exp,'\SWIR\',cfg.darkref_tag,SolsaID,'_SWIR*.raw']},...
        {[path_corescan_exp,'\SWIR\',cfg.darkref_tag,SolsaID,'_SWIR*.hdr']},...
        {[path_corescan_exp,'\SWIR\',cfg.whiteref_tag,SolsaID,'_SWIR*.raw']},...
        {[path_corescan_exp,'\SWIR\',cfg.whiteref_tag,SolsaID,'_SWIR*.hdr']},...
        {[path_corescan_exp,'\SWIR\',SolsaID,'_SWIR*.raw']},...
        {[path_corescan_exp,'\SWIR\',SolsaID,'_SWIR*.hdr']},...
        {[path_corescan_exp,'\VNIR\',cfg.darkref_tag,SolsaID,'_FX10*.raw']},...
        {[path_corescan_exp,'\VNIR\',cfg.darkref_tag,SolsaID,'_FX10*.hdr']},...
        {[path_corescan_exp,'\VNIR\',cfg.whiteref_tag,SolsaID,'_FX10*.raw']},...
        {[path_corescan_exp,'\VNIR\',cfg.whiteref_tag,SolsaID,'_FX10*.hdr']},...
        {[path_corescan_exp,'\VNIR\',SolsaID,'_FX10*.raw']},...
        {[path_corescan_exp,'\VNIR\',SolsaID,'_FX10*.hdr']}...
        ];
    file_names_corescan = [
        {'darkref_SWIR_file.raw \n'},... 
        {'darkref_SWIR_file.hdr \n'},...
        {'whiteref_SWIR_file.raw \n'},... 
        {'whiteref_SWIR_file.hdr \n'},...
        {'SWIR_file.raw \n'},... 
        {'SWIR_file.hdr \n'},...
        {'darkref_VNIR_file.raw \n'},... 
        {'darkref_VNIR_file.hdr \n'},...
        {'whiteref_VNIR_file.raw \n'},... 
        {'whiteref_VNIR_file.hdr \n'},...
        {'VNIR_file.raw \n'},... 
        {'VNIR_file.hdr \n'}...
        ];
    
    ind_exist_corescan_list = true(1,length(files_to_check_corescan));
    for f = 1:length(files_to_check_corescan)
        ind_exist_corescan_list(f) = ~isempty(dir(files_to_check_corescan{f}));
    end
    
    if min(ind_exist_corescan_list)
        ind_exist_corescan = true;
        fail_check_corescan = '';
    else
        fail_check_corescan = ['Some corescan files are missing: \n', file_names_corescan{~ind_exist_corescan_list}];
        ind_exist_corescan = false;
    end
catch
    fail_check_corescan = 'Files not found. Could not load the corescan experiment from the DB. ';
    ind_exist_corescan = false;
end

% Library 
if ~exist(cfg.path_library, 'file')
    fprintf(['Could not find the Library at: ', cfg.path_library, '\n'])
end

% Fusion
try
    % Get experiment fusion
    path_fusion = [v.acquisition_path,'\XYZ_DFusion'];
    
    % Checking the files
    files_to_check_fusion = [
        {[path_fusion,'\XYZ_Data.mat']},...
        {[path_fusion,'\shadow_filter.mat']},...
        {[path_fusion,'\S2P\SWIR\index_x+index_y.mat']},...
        {[path_fusion,'\S2P\VNIR\index_x+index_y.mat']}];
    file_names_fusion = [
        {'XYZ_Data.mat '},...
        {'shadow_filter.mat '},...
        {'SWIR\index_x+index_y.mat '},...
        {'VNIR\index_x+index_y.mat '}];

    test_exist_fusion = true(1,length(files_to_check_fusion));
    for f = 1:length(files_to_check_fusion)
        test_exist_fusion(f) = ~isempty(dir(files_to_check_fusion{f}));
    end
    
    if min(test_exist_fusion)
        test_exist_fusion = true;
        fail_check_fusion = '';
    else
        fail_check_fusion = ['Some fusion files are missing: ',file_names_fusion{~test_exist_fusion},'\n'];
        test_exist_fusion = false;
    end
catch
    fail_check_fusion = 'Files not found. Could not load the fusion experiment from the DB. \n';
    test_exist_fusion = false;
end

% Check files merge
if ind_exist_corescan && test_exist_fusion
    print_check_files = 'PASS\nAll corescan and fusion files are found\n';
else
    print_check_files = ['FAIL\n',fail_check_corescan,fail_check_fusion,'\n'];
end

% Check Load Data
% SWIR VNIR
try
    if ind_exist_corescan
        fprintf('Check load Data\n')
        img = Multi_att_Img();
        img = load_envi_files(img, [path_corescan_exp,'\SWIR']);
        img = load_envi_files(img, [path_corescan_exp,'\VNIR']);
        ind_data_swir_vnir = true;
    else
        error('fail')
    end
catch
    ind_data_swir_vnir = false;
end

% Fusion
try
    if ind_exist_corescan
        x_samp = 8;
        y_samp = 8;
        fus_data = load_fusion_files(path_fusion,x_samp,y_samp,true);
        ind_data_fusion = isfield(fus_data,'vnir') && isfield(fus_data,'swir');
    else
        error('fail')
    end
catch
    ind_data_fusion = false;
end

% Output
if ind_data_swir_vnir && ind_data_fusion
    print_check_data_load = 'PASS\nAll data can be loaded.\n';
else
    data_ind = [{'SWIR & VNIR '}, {'Fusion files '}];
    print_check_data_load = ['FAIL\nSome data cannot be loaded: ',data_ind{[~ind_data_swir_vnir, ~ind_data_fusion]},'\n'];
end

% Check Global
ind_global = min(ind_exist_corescan) && min(test_exist_fusion) && ind_check_config && ind_data_swir_vnir && ind_data_fusion;

if ind_global
    print_check_glob = 'PASS\n';
else
    print_check_glob = 'FAIL\n';
end

% Estimated time
try
    sx = size(fus_data.swir.index_x,1);
    time_seconds = 7.24E-07*sx^2 + 0.007346188*sx +  3.326796916;
catch
    time_seconds = '##';
end

% Print outputs
if v.check 
    fprintf([num2str(time_seconds),' seconds\n'])
    fprintf(print_check_glob)
    fprintf(print_check_files)
    fprintf(print_check_config)
    fprintf(print_check_data_load)
end


% Run the process
if ~v.check 
    if ind_global
        try            
            % Init
            min_max_wl = [cfg.min_wl_vnir,cfg.max_wl_vnir,cfg.min_wl_swir,cfg.max_wl_swir];

            % Load the data 
            fprintf('Load Data\n')
            img = Multi_att_Img();
            img = load_fused_data(img,path_corescan_exp,path_fusion);

            % Load the library
            fprintf('Load Library\n')
            library = Multi_att_Lib();
            library = load_csv(library,cfg.path_library);
            
            % Application of the method 
            fprintf('Apply Classification\n')
            multi_class = Multi_Classifier();
            multi_class = method(...
                multi_class,cfg.process_sample,...
                img,...
                library,...
                cfg.att_names,...
                cfg.spectrum_types,...
                min_max_wl,...
                []...
            );
            
            % Exporting images and tables
            export_path = path_corescan_exp;
            export_key_param(multi_class,export_path)
            
            % Check exported files
            files_to_check_exp = [
                {'HS_rgb_recons.png'},...
                {'ROI_error_lib.csv'},...
                {'ROI_minerals_ind.csv'},...
                {'error_lib.csv'},...
                {'minerals_ind.csv'},...
                {'error_lib.png'},...
                {'minerals_ind.png'},...
             ];
            
            exported_files = '[';
            for f = 1:length(files_to_check_exp)
                if ~isempty(dir(fullfile(export_path,files_to_check_exp{f})))
                    exported_files = [exported_files, files_to_check_exp{f}, ', '];
                end
            end
            exported_files = [exported_files(1:end-2),']'];
            
            ind_process = 1;
        catch
            ind_process = 0;
        end
        
        % Output message
        if ind_process
            fprintf('PASS\n')
            fprintf('Processed and updated successfully\n')
        else
            fprintf('FAIL\n')
            fprintf('Cannot process: ')
            if ~ind_process
                fprintf('error in the hyperspectral processing ')
            end
            fprintf('\n')
        end
    else
        fprintf('FAIL\nCheck failed, abort...\n')
    end
end
end
