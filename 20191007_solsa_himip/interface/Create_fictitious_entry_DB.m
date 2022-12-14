conda_exe = 'C:\Users\christophe.bessin\AppData\Local\Continuum\anaconda3\Scripts\activate.bat C:\Users\christophe.bessin\AppData\Local\Continuum\anaconda3';

config_path = 'D:\SOLSA\HIMIP\20191007_solsa_himip\utils\Config.m';
run(config_path)

config_struct.data.attributes.('uuid') = char(java.util.UUID.randomUUID);
config_struct.data.attributes.('label') = 'corsen3';
config_struct.data.attributes.('config_json') = cfg;
config_struct.data.attributes.('special_details') = '';
config_struct.data.type = 'hyperspectral_config';

% Create json file and correct it 
json_config = jsonencode(config_struct);
json_config = strrep(json_config,'"data":{','"data":[{');
json_config = strrep(json_config,'}}','}]}');

% Export tmp file
path_file_hs = fullfile(cfg.dir_data,'tmp_hs.json');
filejson = fopen(path_file_hs,'w');
fprintf(filejson, json_config);
fclose(filejson);

% PUT hs results
evalc(['!',conda_exe,' & curl -sS --header "Content-Type: application/json" -X PUT -d @', path_file_hs,' https://solsa.crystallography.net/db/test_samples/hyperspectral_config'])


% Sample struct
sample_id.attributes.SolsaID = 'ER-NC00-0056';
sample_id.attributes.uuid = 'a4495996-288c-11e9-882a-e77ce34a9608';
sample_id.type = "sample";


% Structure experiment_corescan // corescan_config.uuid
exp_corescan_struct.data.attributes.('sample_id') = sample_id;
exp_corescan_struct.data.attributes.('uuid') = char(java.util.UUID.randomUUID);
exp_corescan_struct.data.attributes.('corescan_dir_path') = 'D:\\SOLSA\\DATA\\111_fused\\Ech_ref\\ER-NC00-0056';
exp_corescan_struct.data.attributes.('acquisition_status') = 'done';
exp_corescan_struct.data.attributes.('experiment_details_json') = '';
exp_corescan_struct.data.type = 'experiment_corescan';

% Create json file and correct it 
json_output_exp_cs = jsonencode(exp_corescan_struct);
json_output_exp_cs = strrep(json_output_exp_cs,'"data":{','"data":[{');
json_output_exp_cs = strrep(json_output_exp_cs,'}}','}]}');

% Export tmp file
path_file_hs = fullfile(cfg.dir_data,'tmp_hs.json');
filejson = fopen(path_file_hs,'w');
fprintf(filejson, json_output_exp_cs);
fclose(filejson);

% PUT hs results
evalc(['!',conda_exe,' & curl -sS --header "Content-Type: application/json" -X PUT -d @', path_file_hs,' ', DB_URL,'/',v.DB_URI,'experiment_corescan'])



% Structure experiment_fusion // fusion_config.uuid
exp_fusion_struct.data.attributes.('uuid') = char(java.util.UUID.randomUUID);
exp_fusion_struct.data.attributes.('corescan_id') = 3;
exp_fusion_struct.data.attributes.('data_fusion_dir_path') = 'D:\\SOLSA\\DATA\\111_fused\\Ech_ref\\ER-NC00-0056\\XYZ_DFusion';
exp_fusion_struct.data.attributes.('acquisition_status') = 'done';
exp_fusion_struct.data.attributes.('experiment_details_json') = '';
exp_fusion_struct.data.type = 'experiment_fusion';

% Create json file and correct it 
json_output_exp_fus = jsonencode(exp_fusion_struct);
json_output_exp_fus = strrep(json_output_exp_fus,'"data":{','"data":[{');
json_output_exp_fus = strrep(json_output_exp_fus,'}}','}]}');

% Export tmp file
path_file_hs = fullfile(cfg.dir_data,'tmp_hs.json');
filejson = fopen(path_file_hs,'w');
fprintf(filejson, json_output_exp_fus);
fclose(filejson);

% PUT hs results
evalc(['!',conda_exe,' & curl -sS --header "Content-Type: application/json" -X PUT -d @', path_file_hs,' ', DB_URL,'/',v.DB_URI,'experiment_fusion'])




%get_json = evalc(['!',conda_exe,' & curl -sS ',DB_URL,'/',v.DB_URI,'experiment_corescan/',v.CoreScan_Experiment_uuid,'?format=json']);