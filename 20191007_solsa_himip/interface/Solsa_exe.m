function Solsa_exe(varargin) %CoreScan_Experiment_uuid,Fusion_Experiment_uuid,CoreScan_Config_ID,HyperSpectral_Config_ID,DB_URI,CoreScan_Table_Name,Fusion_Table_Name,HyperSpectral_Experiment_uuid,check)
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

%}


warning ('off','all');

% Curl exe, if in conda
test_curl = evalc('!curl');
if strlength(test_curl) >= 53
    if strcmp(test_curl(1:53),"'curl' n'est pas reconnu en tant que commande interne")
        conda_exe = 'C:\Users\christophe.bessin\AppData\Local\Continuum\anaconda3\Scripts\activate.bat C:\Users\christophe.bessin\AppData\Local\Continuum\anaconda3';
    else
        conda_exe = '';
    end
else
    conda_exe = '';
end

if ~isempty(conda_exe)
    conda_exe = [conda_exe, ' & '];
end

% Init values of the arguments
v.CoreScan_Experiment_uuid = '25d60a86-2be1-4948-b1ca-5a0c29a7b2a7';
v.Fusion_Experiment_uuid = '7e1052ab-036b-4818-b814-a30e2bad3bbc';
v.CoreScan_Config_ID = '1';
v.HyperSpectral_Config_ID = '26';
v.DB_URI = 'https://solsa.crystallography.net/db/test_samples/';
v.CoreScan_Table_Name = 'hhh';
v.Fusion_Table_Name = 'lkj';
v.HyperSpectral_Experiment_uuid = '09d15fe4-f8a8-4100-98f1-0bfbb9fb34c1';
v.check = '0';

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
    config_json_entry = evalc(['!',conda_exe,'curl -sS ',v.DB_URI,'hyperspectral_config/',v.HyperSpectral_Config_ID,'?format=json']);
    
    config_struct_entry = jsondecode(config_json_entry);
    config_json = config_struct_entry.data.attributes.config_json;
    
    global cfg
    cfg = jsondecode(config_json);
    
    % Check Config
    names_expected_cfg =   [{'dir_proj'               };...
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
    
    ind_check_config = all(ismember(names_expected_cfg,fieldnames(cfg)));
    
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
    get_json = evalc(['!',conda_exe,'curl -sS ',v.DB_URI,'experiment_corescan/',v.CoreScan_Experiment_uuid,'?format=json']);
    json_corescan = jsondecode(get_json);
    corescan_id = json_corescan.data.attributes.id;
    SolsaID = json_corescan.data.attributes.sample_id.attributes.SolsaID;
    uuid_sample = json_corescan.data.attributes.sample_id.attributes.uuid;
    path_corescan_exp = json_corescan.data.attributes.corescan_dir_path;
    
    % Checking the files
    files_to_check_corescan = [{[path_corescan_exp,'\SWIR\',cfg.darkref_tag,SolsaID,'_SWIR*.raw']},{[path_corescan_exp,'\SWIR\',cfg.darkref_tag,SolsaID,'_SWIR*.hdr']},...
                               {[path_corescan_exp,'\SWIR\',cfg.whiteref_tag,SolsaID,'_SWIR*.raw']},{[path_corescan_exp,'\SWIR\',cfg.whiteref_tag,SolsaID,'_SWIR*.hdr']},...
                               {[path_corescan_exp,'\SWIR\',SolsaID,'_SWIR*.raw']},{[path_corescan_exp,'\SWIR\',SolsaID,'_SWIR*.hdr']},...
                               {[path_corescan_exp,'\VNIR\',cfg.darkref_tag,SolsaID,'_FX10*.raw']},{[path_corescan_exp,'\VNIR\',cfg.darkref_tag,SolsaID,'_FX10*.hdr']},...
                               {[path_corescan_exp,'\VNIR\',cfg.whiteref_tag,SolsaID,'_FX10*.raw']},{[path_corescan_exp,'\VNIR\',cfg.whiteref_tag,SolsaID,'_FX10*.hdr']},...
                               {[path_corescan_exp,'\VNIR\',SolsaID,'_FX10*.raw']},{[path_corescan_exp,'\VNIR\',SolsaID,'_FX10*.hdr']}];
    file_names_corescan = [{'darkref_SWIR_file.raw \n'}, {'darkref_SWIR_file.hdr \n'},...
                            {'whiteref_SWIR_file.raw \n'}, {'whiteref_SWIR_file.hdr \n'},...
                            {'SWIR_file.raw \n'}, {'SWIR_file.hdr \n'},...
                            {'darkref_VNIR_file.raw \n'}, {'darkref_VNIR_file.hdr \n'},...
                            {'whiteref_VNIR_file.raw \n'}, {'whiteref_VNIR_file.hdr \n'},...
                            {'VNIR_file.raw \n'}, {'VNIR_file.hdr \n'}];
    
    ind_exist_corescan_list = true(1,length(files_to_check_corescan));
    for f = 1:length(files_to_check_corescan)
        ind_exist_corescan_list(f) = ~isempty(dir(files_to_check_corescan{f}));
    end
    
    if min(ind_exist_corescan_list)
        ind_exist_corescan = true;
        fail_check_corescan = '';
    else
        fail_check_corescan = ['Some corescan files are missing: \n',file_names_corescan{~ind_exist_corescan_list}];
        ind_exist_corescan = false;
    end
catch
    fail_check_corescan = 'Files not found. Could not load the corescan experiment from the DB. ';
    ind_exist_corescan = false;
end

% Fusion
try
    % Get experiment fusion
    get_json = evalc(['!',conda_exe,'curl -sS ',v.DB_URI,'experiment_fusion/',v.Fusion_Experiment_uuid,'?format=json']);
    json_fusion = jsondecode(get_json);
    path_fusion = json_fusion.data.attributes.data_fusion_dir_path;
    
    % Checking the files
    files_to_check_fusion = [{[path_fusion,'\XYZ_Data.mat']},...
                             {[path_fusion,'\S2P\SWIR\index_x+index_y.mat']},...
                             {[path_fusion,'\S2P\VNIR\index_x+index_y.mat']}];
    file_names_fusion = [{'XYZ_Data.mat '},{'Indexes_SWIR.mat '},{'Indexes_VNIR.mat '}];

    test_exist_fusion = true(1,length(files_to_check_fusion));
    for f = 1:length(files_to_check_fusion)
        test_exist_fusion(f) = ~isempty(dir(files_to_check_fusion{f}));
    end
    
    if min(test_exist_fusion)
        test_exist_fusion = true;
        fail_check_fusion = '';
    else
        fail_check_fusion = ['Some fusion files are missing: ',file_names_fusion{~test_exist_fusion}];
        test_exist_fusion = false;
    end
catch
    fail_check_fusion = 'Files not found. Could not load the fusion experiment from the DB. ';
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
        if isfield(fus_data,'vnir') && isfield(fus_data,'swir')
            ind_data_fusion = true;
        else
            ind_data_fusion = false;
        end
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
    data_ind = [{'SWIR & VNIR '},{'Fusion files '}];
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
            img = Multi_att_Img();
            img = load_fused_data(img,path_corescan_exp,path_fusion);

            % Load the library
            library = Multi_att_Lib();
            library = load_csv(library,cfg.path_library);
            
            % Application of the method 
            multi_class = Multi_Classifier();
            multi_class = method(multi_class,cfg.process_sample,img,library,cfg.att_names,cfg.spectrum_types,min_max_wl,[]);
            
            % Exporting images and tables
            export_path = path_corescan_exp;
            export_key_param(multi_class,export_path)
            
            % Check exported files
            files_to_check_exp = [{'HS_rgb_recons.png'},...
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

        try
            % EXPORT the results
            % Sample struct
            sample_id.attributes.SolsaID = SolsaID;
            sample_id.attributes.uuid = uuid_sample;
            sample_id.type = "sample";

            % Make the table
            mine = [];
            err = [];
            col = [];
            prop = [];
            for i = 1: length(multi_class.table_minerals_ind.Minerals)
                mine = [mine, multi_class.table_minerals_ind.Minerals{i},', '];
                err = [err, multi_class.table_minerals_ind.Error{i},', '];
                col = [col, multi_class.table_minerals_ind.Color{i},', '];
                prop = [prop, multi_class.table_minerals_ind.Prop{i},', '];
            end
            mine = ['[', mine(1:end-2), ']'];
            err = ['[', err(1:end-2), ']'];
            col = ['[', col(1:end-2), ']'];
            prop = ['[', prop(1:end-2), ']'];

            % Experiment details
            experiment_details_json.results_global_table.Minerals = mine;
            experiment_details_json.results_global_table.Error = err;
            experiment_details_json.results_global_table.Color = col;
            experiment_details_json.results_global_table.Proportion = prop;
            experiment_details_json.method = cfg.process_sample;
            experiment_details_json.path_export_path = replace(export_path,'\','\\');
            experiment_details_json.exported_files = exported_files;

            % Structure results Hyperspectral
            results_struct_hs.data.attributes.('acquisition_status') = 'done';
            results_struct_hs.data.attributes.('config_id') = v.HyperSpectral_Config_ID;
            results_struct_hs.data.attributes.('experiment_details_json') = experiment_details_json;
            results_struct_hs.data.attributes.('corescan_id') = corescan_id;
            results_struct_hs.data.attributes.('measurement_datetime') = char(datetime(datetime,'Format','yyyy-MM-dd HH:mm:ss'));
            results_struct_hs.data.attributes.('special_details') = 'test';
            results_struct_hs.data.attributes.('status_message') = 'tt';
            results_struct_hs.data.attributes.('uuid') = char(java.util.UUID.randomUUID);
            results_struct_hs.data.type = 'experiment_hyperspectral';

            % Structure results ROI
            results_struct_roi.data.attributes.('uuid') = char(java.util.UUID.randomUUID);
            results_struct_roi.data.attributes.('region_top') = multi_class.ROI_minerals_ind(3);
            results_struct_roi.data.attributes.('region_bottom') = multi_class.ROI_minerals_ind(4);
            results_struct_roi.data.attributes.('sample_id') = sample_id;
            results_struct_roi.data.attributes.('special_details') = '';
            results_struct_roi.data.type = 'region_of_interest';

            % Create json file and correct it 
            json_output_exp_hs = jsonencode(results_struct_hs);
            json_output_exp_hs = strrep(json_output_exp_hs,'"data":{','"data":[{');
            json_output_exp_hs = strrep(json_output_exp_hs,'}}','}]}');
            
            % Export tmp file
            path_file_hs = fullfile(cfg.dir_data,'tmp_hs.json');
            filejson = fopen(path_file_hs,'w');
            fprintf(filejson, json_output_exp_hs);
            fclose(filejson);
            
            % PUT hs results
            ret_put = evalc(['!',conda_exe,'curl -sS --header "Content-Type: application/json" -X PUT -d @', path_file_hs,' ',v.DB_URI,'experiment_hyperspectral']);
            if contains(ret_put,'401 Unauthorized')
                error('No credentials')
            end
            
            % Create json file and correct it 
            json_output_roi = jsonencode(results_struct_roi);
            json_output_roi = strrep(json_output_roi,'"data":{','"data":[{');
            json_output_roi = strrep(json_output_roi,'}}','}]}');
            
            % Export tmp file
            path_file_roi = fullfile(cfg.dir_data,'tmp_roi.json');
            filejson = fopen(path_file_roi,'w');
            fprintf(filejson, json_output_roi);
            fclose(filejson);
            
            % PUT ROI results
            ret_put = evalc(['!',conda_exe,'curl -sS --header "Content-Type: application/json" -X PUT -d @', path_file_roi,' ',v.DB_URI,'region_of_interest']);
            if contains(ret_put,'401 Unauthorized')
                error('No credentials')
            end

            ind_DB_update = 1;
        catch
            ind_DB_update = 0;
        end

        % Output message
        if ind_process && ind_DB_update
            fprintf('PASS\n')
            fprintf('Processed and updated successfully\n')
        else
            fprintf('FAIL\n')
            fprintf('Cannot process: ')
            if ~ind_process
                fprintf('error in the hyperspectral processing ')
            elseif ~ind_DB_update
                fprintf('error when updating the database')
            end
            fprintf('\n')
        end
    else
        fprintf('FAIL\nCheck failed, abort...\n')
    end
end
end
