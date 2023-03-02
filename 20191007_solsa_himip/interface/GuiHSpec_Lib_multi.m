function varargout = GuiHSpec_Lib_multi(varargin)
% HSPEC_LIB_MAIN_DISPLAY MATLAB code for GuiHSpec_Lib_multi.fig
%      HSPEC_LIB_MAIN_DISPLAY, by itself, creates a new HSPEC_LIB_MAIN_DISPLAY or raises the existing
%      singleton*.
%
%      H = HSPEC_LIB_MAIN_DISPLAY returns the handle to a new HSPEC_LIB_MAIN_DISPLAY or the handle to
%      the existing singleton*.
%
%      HSPEC_LIB_MAIN_DISPLAY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HSPEC_LIB_MAIN_DISPLAY.M with the given input arguments.
%
%      HSPEC_LIB_MAIN_DISPLAY('Property','Value',...) creates a new HSPEC_LIB_MAIN_DISPLAY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GuiHSpec_Lib_multi_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GuiHSpec_Lib_multi_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GuiHSpec_Lib_multi

% Last Modified by GUIDE v2.5 27-Jul-2020 15:46:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GuiHSpec_Lib_multi_OpeningFcn, ...
                   'gui_OutputFcn',  @GuiHSpec_Lib_multi_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before GuiHSpec_Lib_multi is made visible.
function GuiHSpec_Lib_multi_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GuiHSpec_Lib_multi (see VARARGIN)

% Choose default command line output for GuiHSpec_Lib_multi
handles.output = hObject;

% Initialisations
try
    % Run from Matlab 
    %run('D:\SOLSA\HIMIP\20191007_solsa_himip\utils\Config.m');
    run(fullfile(pwd, 'HS_Analysis', 'HIMIP', '20191007_solsa_himip', 'utils','Config.m'));
    
    % Add paths
    folders = [{'infrastructure'},{'domain'},{'application'}, {'interface'},{'utils'}];
    for i = 1:length(folders)
        %addpath(['D:\SOLSA\HIMIP\20191007_solsa_himip\',folders{i}])
        addpath([pwd, folders{i}])
    end
catch
    % Run from the GUI
    path = fullfile(pwd,'Config.m');
    fileID = fopen(path,'r');
    config_txt = fscanf(fileID,'%c');
    eval(config_txt)
end

% Objects initialisations
% Library
set(handles.file_library, 'String', cfg.path_library);
set(handles.file_ext_library, 'String',cfg.path_library);

% Data
set(handles.data_file, 'String', fullfile(cfg.dir_data,cfg.file_data));

% RGB init
handles.red = cfg.red_default;
set(handles.Red_value, 'String', handles.red);
handles.green = cfg.green_default;
set(handles.Green_value, 'String', handles.green);
handles.blue = cfg.blue_default;
set(handles.Blue_value, 'String', handles.blue);

% dir init
handles.dir_data = cfg.dir_data_fus;
handles.path_library = cfg.path_library;
handles.path_ext_library = cfg.path_library;

% Wavelengths
handles.VNIR.min_wl_val = cfg.min_wl_vnir;
handles.VNIR.max_wl_val = cfg.max_wl_vnir;
handles.SWIR.min_wl_val = cfg.min_wl_swir;
handles.SWIR.max_wl_val = cfg.max_wl_swir;
handles.lim_vnir_swir = cfg.lim_vnir_swir;

set(handles.min_wl_vnir,'String',cfg.min_wl_vnir)
set(handles.max_wl_vnir,'String',cfg.max_wl_vnir)
set(handles.min_wl_swir,'String',cfg.min_wl_swir)
set(handles.max_wl_swir,'String',cfg.max_wl_swir)

handles.hs_lib = [];
handles.ind_sel_table_lib = [];
handles.ind_sel_table_ext_lib = [];
handles.ind_sel_table_results = [];
handles.ind_sel_table_lib_check = [];
handles.mulatt_sel_pts_img = Multi_att_Lib();
handles.sel_pts_img = [];

%
handles.button_name = [{'reflectance'} {'log10'} {'Stav_Golay'} {'hull_removal'} {'hull_curve'} {'noise_hull_rm'} {'hull_noise_rm'} {'diff_1'} {'diff_2'} {'bollinger'}];

handles.img_axes4 = 'minerals_ind';
handles.ROI = [];
handles.process_sample = cfg.process_sample;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = GuiHSpec_Lib_multi_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% =========================================================================
% ============================== Library ==================================

function browse_library_Callback(hObject, eventdata, handles)
global cfg

% Gets the library from file
[library_file, dir_path, oCancel] = uigetfile('*.*', 'Select a file',handles.path_library);
if ~oCancel
    disp('User selects Cancel')
    return
end
handles.path_library = fullfile(dir_path,library_file);
t_lib = tic;

% Loading the library
handles.hs_lib = load_library(dir_path,library_file);

% Display in the table
update_lib_table(handles,'lib')

% Update handles structure
guidata(hObject, handles);

if cfg.print_times
    fprintf('Lib ')
    toc(t_lib)
end

function file_library_Callback(hObject, eventdata, handles)

% Default library file
file_library = strtrim(get(hObject, 'String'));
if isempty(file_library)
    file_library = handles.path_library;
end
set(hObject, 'String', file_library);

% Update handles structure
guidata(hObject, handles);

function file_library_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file_library (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function table_lib_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to table_lib (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

handles.ind_sel_table_lib = eventdata.Indices(:,1);

% Update handles structure
guidata(hObject, handles);

update_curves_spectra_viz(handles)

% ============================ Edit Library ===============================

function create_library_Callback(hObject, eventdata, handles)
% hObject    handle to create_library (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Gets the library from file
[library_file, dir_path, oCancel] = uiputfile('*.csv','Select a file',handles.path_library);
if ~oCancel
    disp('User selects Cancel')
    return
end

% Loading the library
handles.path_library = fullfile(dir_path,library_file);

% Creating the library
hs_lib = Multi_att_Lib();

handles.hs_lib = hs_lib;

% Update handles structure
guidata(hObject, handles);

update_lib_table(handles,'lib')

function delete_selected_entries_Callback(hObject, eventdata, handles)
% hObject    handle to delete_selected_entries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Deletion of the selected entries
entry_names = handles.hs_lib.entry_names(handles.ind_sel_table_lib);
handles.hs_lib = add_modif(handles.hs_lib,1,entry_names);

update_lib_table(handles,'lib')

% Update handles structure
guidata(hObject, handles);

function save_library_Callback(hObject, eventdata, handles)
% hObject    handle to save_library (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

csv_table = build_lib_table(handles.hs_lib);

writetable(csv_table,handles.path_library)

function edit_name_Callback(hObject, eventdata, handles)
% hObject    handle to edit_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tmp = split(handles.hs_lib.entry_names(handles.ind_sel_table_lib),'_');
old_name_suff = [tmp{2},'_',tmp{3}];

handles.hs_lib.entry_names(handles.ind_sel_table_lib) = cellstr([get(handles.new_name_entry,'String'),'_',old_name_suff]);

update_lib_table(handles,'lib')

% Update handles structure
guidata(hObject, handles);

function new_name_entry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to new_name_entry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% =========================================================================
% ========================= External Library ==============================

function browse_ext_library_Callback(hObject, eventdata, handles)

% Gets the library from file
[library_file, dir_path, oCancel] = uigetfile('*.*', 'Select a file',handles.path_ext_library);
if ~oCancel
    disp('User selects Cancel')
    return
end
handles.path_ext_library = fullfile(dir_path,library_file);

% Loading the library
handles.hs_ext_lib = load_library(dir_path,library_file);

% Display in the table
update_lib_table(handles,'ext_lib')

% Update handles structure
guidata(hObject, handles);

function file_ext_library_Callback(hObject, eventdata, handles)

% Default library file
file_ext_library = strtrim(get(hObject, 'String'));
if isempty(file_ext_library)
    file_ext_library = handles.path_library;
end
set(hObject, 'String', file_ext_library);

% Update handles structure
guidata(hObject, handles);

function file_ext_library_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file_ext_library (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function table_ext_lib_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to table_ext_lib (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

handles.ind_sel_table_ext_lib = eventdata.Indices(:,1);

% Update handles structure
guidata(hObject, handles);

update_curves_spectra_viz(handles)

% ================== Add Entry from External Library ======================

function add_spec_to_lib_from_ext_Callback(hObject, eventdata, handles)
% hObject    handle to add_spec_to_lib_from_ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


idx_entry = handles.ind_sel_table_ext_lib;

if ~isempty(idx_entry)
    handles.entry = get_from_indexes(handles.hs_ext_lib,idx_entry,[]);
    
    % Add the modification to the Library
    handles.hs_lib = add_modif(handles.hs_lib,0,handles.entry.entry_names,handles.entry);
    
    update_lib_table(handles,'lib')
end

% Update handles structure
guidata(hObject, handles);

% =========================================================================

function hs_lib = load_library(dir_path,library_file)
% Load a library
% 
%
try
    path_library = fullfile(dir_path,library_file);

    hs_lib = Multi_att_Lib();
    extension_file = library_file(end-2:end);
    if strcmp(extension_file, 'csv')
        % csv case
        hs_lib = load_csv(hs_lib,path_library);
    else
        % other cases
        % Trying to read it as a usgs library
        dir_ref = pwd;
        try
            hs_lib = load_usgs(hs_lib,dir_path);
        catch
            cd(dir_ref)
            try
                hs_lib = load_eco_lib(hs_lib,dir_path);
            catch
                cd(dir_ref)
                error('Could not load the library. File selected is not a valid library, or its folder does not contain a USGS Library...')
            end
        end
    end
catch
    ErrorMessage = lasterr;
    errordlg(ErrorMessage,'Error');
end

% =========================================================================
% ================ Edit Library from Hyper-Spectral Image =================

% ============================= Load file =================================

function browse_data_Callback(hObject, eventdata, handles)
% hObject    handle to browse_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global cfg 
try
    % Gets the spectra from file
    dir_data = uigetdir(handles.dir_data);

    % Update of the dir_data
    handles.dir_data = dir_data;

    browse_t = tic;
    
    % Load the data 
    img = Multi_att_Img();
    handles.img = load_fused_data(img,dir_data,fullfile(dir_data,'XYZ_DFusion'));

    % Empty the selected points table
    handles.sel_pts_img = [];
    handles.mulatt_sel_pts_img = Multi_att_Lib();
    handles.sel_spec_for_edit = [];

    % Set the sample name
    sp_tmp = split(dir_data,'\');
    handles.sample_name = sp_tmp{end};
    set(handles.data_file, 'String', handles.sample_name);

    % Update handles structure
    guidata(hObject, handles);

    % Display the image
    update_curves_axes1(hObject, handles)
    
    if cfg.print_times
        fprintf('Load data ')
        toc(browse_t)
    end
catch
    ErrorMessage = lasterr;
    errordlg(ErrorMessage,'Error');
end

function data_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ========================= Table management ==============================
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

% Get the coordinate of cursor when mouse clicked
coordinate = get(gca, 'CurrentPoint');
coor_real = round(coordinate(1, 1:2));

% From coordinates to indices
% !! warning : the coordinate X is displayed in reverse (- imshow cannot
% display otherwise), the sign is reversed manually to recover the correct
% value
min_mat = abs(handles.img.coord.x - coor_real(1)) + abs(handles.img.coord.y - coor_real(2)); 
[~, index] = min(min_mat(:));

%[corr(1), corr(2)] = ind2sub(size(handles.img.coord.x),index);

% Set coordinates in the table
size_sel_pts_img = size(handles.sel_pts_img);
handles.sel_pts_img(size_sel_pts_img(1) + 1,:) = [coor_real, rand, index];

% Get the data from the image
multi_att_tmp = get_from_indexes(handles.img,handles.sel_pts_img(:,4),[],1);

handles.mulatt_sel_pts_img = copy(handles.img,'lib');
handles.mulatt_sel_pts_img.att_data = multi_att_tmp.att_data;
handles.mulatt_sel_pts_img.entry_names = strrep(cellstr(strcat(num2str(handles.sel_pts_img(:,1)),'_',num2str(handles.sel_pts_img(:,2)))),' ','');

% Spectrum to be added
if handles.compute_spectra_average.Value
    handles.mulatt_add = handles.mulatt_sel_pts_img;
    handles.mulatt_add.att_data = mean(handles.mulatt_sel_pts_img.att_data,1);
else
    handles.mulatt_add = [];
end

% Update handles
guidata(hObject, handles)

% Display spectrum
update_curves_axes1(hObject, handles)
update_curves_spectra_viz(handles)

function delete_all_edit_spectra_Callback(hObject, eventdata, handles)
% hObject    handle to delete_all_edit_spectra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.mulatt_sel_pts_img = Multi_att_Lib();
handles.sel_pts_img = [];

% Update handles
guidata(hObject, handles)

update_curves_axes1(hObject, handles)
update_curves_spectra_viz(handles)

function delete_items_Callback(hObject, eventdata, handles)
% hObject    handle to delete_items (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Deletion in the multi att
name = strcat(num2str(handles.sel_pts_img(handles.sel_spec_for_edit,1)),'_',num2str(handles.sel_pts_img(handles.sel_spec_for_edit,2)));
handles.mulatt_sel_pts_img = add_modif(handles.mulatt_sel_pts_img,1,name);

% Deletion in the table
handles.sel_pts_img(handles.sel_spec_for_edit,:) = [];

% Update handles
guidata(hObject, handles)

update_curves_axes1(hObject, handles)
update_curves_spectra_viz(handles)

function table_pts_for_edit_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to table_pts_for_edit (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

handles.sel_spec_for_edit = eventdata.Indices(:,1);

% Update handles structure
guidata(hObject, handles);


% ======================= Min-Max Wavelengths =============================
function min_wl_vnir_Callback(hObject, eventdata, handles)
% hObject    handle to min_wl_vnir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

min_wl_vnir_input = str2num(get(handles.min_wl_vnir,'String'));

if min_wl_vnir_input < handles.VNIR.max_wl_val
    handles.VNIR.min_wl_val = min_wl_vnir_input;
else
    set(handles.min_wl_vnir,'String',handles.VNIR.min_wl_val)
end

% Update handles
guidata(hObject, handles)

update_curves_spectra_viz(handles)

function max_wl_vnir_Callback(hObject, eventdata, handles)
% hObject    handle to max_wl_vnir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

max_wl_vnir_input = str2num(get(handles.max_wl_vnir,'String'));

if (max_wl_vnir_input > handles.VNIR.min_wl_val) && (max_wl_vnir_input < handles.lim_vnir_swir)
    handles.VNIR.max_wl_val = max_wl_vnir_input;
else
    set(handles.max_wl_vnir,'String',handles.VNIR.max_wl_val)
end

% Update handles
guidata(hObject, handles)

update_curves_spectra_viz(handles)

function min_wl_swir_Callback(hObject, eventdata, handles)
% hObject    handle to min_wl_swir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

min_wl_swir_input = str2num(get(handles.min_wl_swir,'String'));

if (min_wl_swir_input < handles.SWIR.max_wl_val) && (min_wl_swir_input > handles.lim_vnir_swir)
    handles.SWIR.min_wl_val = min_wl_swir_input;
else
    set(handles.min_wl_swir,'String',handles.SWIR.min_wl_val)
end

% Update handles
guidata(hObject, handles)

update_curves_spectra_viz(handles)

function max_wl_swir_Callback(hObject, eventdata, handles)
% hObject    handle to max_wl_swir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

max_wl_swir_input = str2num(get(handles.max_wl_swir,'String'));

if max_wl_swir_input > handles.SWIR.min_wl_val
    handles.SWIR.max_wl_val = max_wl_swir_input;
else
    set(handles.max_wl_swir,'String',handles.SWIR.max_wl_val)
end

% Update handles
guidata(hObject, handles)

update_curves_spectra_viz(handles)

% =========================================================================
% ==================== Add Entry from HS Image ============================
function add_spec_to_Lib_Callback(hObject, eventdata, handles)
% hObject    handle to add_spec_to_Lib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    if size(handles.table_pts_for_edit.Data,1) == 1
        idx_entry = 1;
    else
        idx_entry = handles.sel_spec_for_edit;
    end

    if ~isempty(idx_entry)
        handles.entry = get_from_indexes(handles.mulatt_sel_pts_img,idx_entry,[]);
        entry_name = [get(handles.mineral_new_entry,'String'),'_', handles.sample_name];

        % Renaming the entries
        handles.entry.entry_names = repmat({entry_name},[length(handles.entry.entry_names),1]);

        % Add the modification to the Library
        handles.hs_lib = add_modif(handles.hs_lib,0,entry_name,handles.entry);

        update_lib_table(handles,'lib')
    end

    % Update handles structure
    guidata(hObject, handles);
catch
    ErrorMessage = lasterr;
    errordlg(ErrorMessage,'Error');
end


% =========================================================================
% ====================== Figures management ===============================

% ======================== Buttons display spectra ========================

function compute_spectra_average_Callback(hObject, eventdata, handles)
% hObject    handle to compute_spectra_average (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of compute_spectra_average

update_curves_spectra_viz(handles)

function spectrum_type_vnir_Callback(hObject, eventdata, handles)
% hObject    handle to spectrum_type_vnir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.spect_trans_vnir = str2num(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

update_curves_spectra_viz(handles)

function spectrum_type_vnir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spectrum_type_vnir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

init_val = [1,7];

set(hObject,'String',replace(num2str(init_val),'  ',','));

handles.spect_trans_vnir = init_val;

% Update handles structure
guidata(hObject, handles);

function spectrum_type_swir_Callback(hObject, eventdata, handles)
% hObject    handle to spectrum_type_swir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.spect_trans_swir = str2num(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

update_curves_spectra_viz(handles)

function spectrum_type_swir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spectrum_type_swir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

init_val = [1,7];

set(hObject,'String',replace(num2str(init_val),'  ',','));

handles.spect_trans_swir = init_val;

% Update handles structure
guidata(hObject, handles);

% ========================== Display figures ==============================
function spectra_viz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spectra_viz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate spectra_viz
xlabel('Wavelength (nm)')
title('Spectra')
ylim([0 1])

% Set the following of the motion of the mouse
set(gcf, 'WindowButtonMotionFcn', @mouseMove);

function mouseMove(object, eventdata, handles)

chil_f = get(gcf,'Children');
tag_axe = 'spectra_viz';

for i = 1:length(chil_f)
    chil_tmp = chil_f(i).Children;
    for j = 1:length(chil_tmp)
        axe_test = chil_f(i).Children(j);
        if strcmp(axe_test.Tag,tag_axe)
            axe = axe_test;
        end
    end
end

try
    C = get (axe, 'CurrentPoint');
    Xlim = get(axe, 'XLim');
    Ylim = get(axe, 'YLim');

    % Create an object text and save it in the axes in UserData
    user_data = get(axe, 'UserData');
    if ~isfield(user_data, 'text_handle') || ~ishandle(user_data.text_handle)
        user_data.text_handle = text(nan, nan,'','Parent', axe, 'FontSize',8,'Color','r');
        set(axe, 'UserData', user_data); % save to axes
    end

    pos_x = Xlim(2) - (Xlim(2)-Xlim(1))/9;

    % Edit the text object
    obj_text = get(axe,'UserData');

    % Displays only if above the axes
    if Xlim(1) < C(1,1) && Xlim(2) > C(1,1) && ...
       Ylim(1) < C(1,2) && Ylim(2) > C(1,2)
        set(obj_text.text_handle,'Position', [pos_x, .1], 'String', ['x= ',num2str(round(C(1,1),1)),newline,'y= ',num2str(round(C(1,2),3))]);
    else
        set(obj_text.text_handle,'Position', [pos_x, .1], 'String', '');
    end
catch
end

function update_curves_spectra_viz(handles)

% Init
axes(handles.spectra_viz);
hold on
cla

% Display the spectra selected in the library
if ~isempty(handles.ind_sel_table_lib)
    draw_curves_spectra_viz(get_from_indexes(handles.hs_lib,handles.ind_sel_table_lib,[]),0,handles)
end

% Display the spectra selected in the external library
if ~isempty(handles.ind_sel_table_ext_lib)
    draw_curves_spectra_viz(get_from_indexes(handles.hs_ext_lib,handles.ind_sel_table_ext_lib,[]),0,handles)       
end

% Display the spectra of the points digitalized on the image 
if ~isempty(handles.mulatt_sel_pts_img.att_data)
    draw_curves_spectra_viz(handles.mulatt_sel_pts_img,1,handles)
end

% Display the spectra of the library check
if ~isempty(handles.ind_sel_table_lib_check)
    draw_curves_spectra_viz(get_from_indexes(handles.find_mineral.data_error_lib,handles.ind_sel_table_lib_check,[]),0,handles)
end

% Display the spectra of the results
if ~isempty(handles.ind_sel_table_results)
    draw_curves_spectra_viz(get_from_indexes(handles.find_mineral.ma_img,handles.ind_sel_table_results,[]),0,handles)
end

handles.spectra_viz.XLim = [handles.VNIR.min_wl_val, handles.SWIR.max_wl_val];

function draw_curves_spectra_viz(multi_att,mean_calc,handles)
% Update the curves visualisation when a parameter is changed
%
% Parameters
% ---------
% mul_att: :obj:
%       object multiattribute, its spectra will be displayed
% sp_types: array of int
%       types of the spectra transformation to apply
% mean_calc: boolean
%       option to display the average of the 
% handles
%

% Initialisation 
prop_names = [{'VNIR'},{'SWIR'}];
sp_trans.VNIR = handles.spect_trans_vnir;
sp_trans.SWIR = handles.spect_trans_swir;

% Parameters of the additional curves
curve_style = [{'-'} {'-'} {'-.'} {'-'} {'--'} {'--'} {'-.'} {'-'} {'-'} {'-.'} {'--'}];
colors      = [{'k'} {'y'} {'b'}  {'r'} {'r'}  {'k'}  {'k'}  {'b'} {'b'} {'k'}  {'r'} ];

mean_opt = mean_calc && get(handles.compute_spectra_average,'Value');

for i = 1:length(prop_names)
    % Remove NaNs
    multi_att = remove_NaN(multi_att);
    
    % Remove the parts of the spectra outside of the bounds
    mem = ismember(multi_att.att_groups,prop_names{i});
    if max(mem)==1
        multi_att = rm_min_max_wls(multi_att,prop_names{i},...
                                 handles.(prop_names{i}).min_wl_val,...
                                 handles.(prop_names{i}).max_wl_val);
    end
    
    % Loop over the different spectra transformation options
    % Skip if all the wavelengths have been removed in  the previous step
    mem = ismember(multi_att.att_groups,prop_names{i});
    if max(mem)==1
        for j = 1:length(sp_trans.(prop_names{i}))
            s = sp_trans.(prop_names{i})(j);
            multi_att_mod = apply_change_spec(multi_att,{prop_names{i}},s,0);
            ma_group = get_from_group(multi_att_mod,prop_names{i});
            if mean_opt
                multi_att_mean = apply_change_spec(multi_att,{prop_names{i}},s,1);
                ma_group_mean = get_from_group(multi_att_mean,prop_names{i});
                plot(multi_att.att_ids.(prop_names{i}),ma_group.att_data,'k-.')
                plot(multi_att_mean.att_ids.(prop_names{i}),ma_group_mean.att_data,'r-')
            else
                plot(multi_att.att_ids.(prop_names{i}),ma_group.att_data,[colors{s} curve_style{s}])
            end
        end
    end
end

function update_curves_axes1(hObject, handles)

handles = guidata(hObject);

% Display
handles.norm_opt = 2;
rgb_img = rgb_image_builder(handles.img, handles.red, handles.green, handles.blue, handles.norm_opt);

% Finds the middle indexes of the image to avoid edge effect
coord_x = handles.img.axis.x;
coord_y = handles.img.axis.y;

axes(handles.axes1); 
hold on
cla

h1 = handles.axes1;
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

% When click on the figure, jump to click function
set(imagesc(rgb_img), 'ButtonDownFcn', {@axes1_ButtonDownFcn, handles.axes1}, 'XData', coord_x, 'YData', coord_y)

if ~isempty(handles.sel_pts_img)
    scatter(handles.sel_pts_img(:,1),handles.sel_pts_img(:,2),[],handles.sel_pts_img(:,3), '.')
    table_coord_edit = table(handles.sel_pts_img(:,1:2));
    set(handles.table_pts_for_edit,'Data',table_coord_edit{:,:},'ColumnName',[{'Ang. Samp.'},{'Line'}])
else
    table_coord_edit = table([]);
    set(handles.table_pts_for_edit,'Data',table_coord_edit{:,:},'ColumnName',[{'Ang. Samp.'},{'Line'}])
end

% Setting the aspect ratio
h1.XLim = [min(coord_x), max(coord_x)];
h1.YLim = [min(coord_y), max(coord_y)];
h1.DataAspectRatio = [1,1,1];

function update_lib_table(handles,lib)
hs_lib = handles.(strcat('hs_',lib));

% Display in the table
table_min_id = table(hs_lib.entry_names);
set(handles.(strcat('table_',lib)), 'Data', table_min_id{:,:},...
                                    'ColumnName',[{'Mineral_ID'}],...
                                    'ColumnWidth',{200})

function update_image_axes4(handles)
% Set the image axes4
%

% Possible options
fig = [{'error_lib'}, {'minerals_ind'}];

if contains(handles.img_axes4,fig) % Check if availaible option
    if ~isempty(handles.find_mineral.(strcat('ROI_',handles.img_axes4,'_rgb'))) % Check if the attribute exists
        % Copy of the axes in the attribute
        h1 = handles.find_mineral.(strcat('ROI_',handles.img_axes4,'_rgb')).Children;
        
        % Copy of the properties to the GUI axes
        handles.axes4.XDir = h1.XDir;
        handles.axes4.XAxis.TickValues = h1.XAxis.TickValues;
        handles.axes4.XTickMode = h1.XTickMode;
        handles.axes4.XLim = h1.XLim;

        handles.axes4.YDir = h1.YDir;
        handles.axes4.YAxis.TickValues = h1.YAxis.TickValues;
        handles.axes4.YTickMode = h1.YTickMode;
        handles.axes4.YLim = h1.YLim;

        handles.axes4.DataAspectRatio = h1.DataAspectRatio;
        
        copyobj(h1.Children,handles.axes4)
    else
        % Display empty axes
        axes(handles.axes4)
        cla
    end
end


% ============================ Color RGB ==================================

function Red_value_Callback(hObject, eventdata, handles)
% hObject    handle to Red_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Red_value as text
%        str2double(get(hObject,'String')) returns contents of Red_value as a double

% Default red value
red = get(hObject, 'String');

red_t = str2double(red);
if ~isnan(red_t)
    red = red_t;
end

if isempty(red)
    red = handles.red;
end
set(hObject, 'String', red);
handles.red = red;

% Update handles structure
guidata(hObject, handles);

function Red_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Red_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Green_value_Callback(hObject, eventdata, handles)
% hObject    handle to Green_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Green_value as text
%        str2double(get(hObject,'String')) returns contents of Green_value as a double

% Default green value
green = get(hObject, 'String');

green_t = str2double(green);
if ~isnan(green_t)
    green = green_t;
end

if isempty(green)
    green = handles.green;
end
set(hObject, 'String', green);
handles.green = green;

% Update handles structure
guidata(hObject, handles);

function Green_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Green_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Blue_value_Callback(hObject, eventdata, handles)
% hObject    handle to Blue_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Blue_value as text
%        str2double(get(hObject,'String')) returns contents of Blue_value as a double

% Default green value
blue = get(hObject, 'String');

blue_t = str2double(blue);
if ~isnan(blue_t)
    blue = blue_t;
end

if isempty(blue)
    blue = handles.blue;
end
set(hObject, 'String', blue);
handles.blue = blue;

% Update handles structure
guidata(hObject, handles);

function Blue_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Blue_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function update_rgb_values_Callback(hObject, eventdata, handles)
% hObject    handle to update_rgb_values (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_curves_axes1(hObject,handles)


% =========================================================================
% ====================== Apply Mineral Mapping ============================

function table_error_lib_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to table_error_lib (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

handles.ind_sel_table_lib_check = eventdata.Indices(:,1);

% Update handles structure
guidata(hObject, handles);

update_curves_spectra_viz(handles)

function table_minerals_ind_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to table_minerals_ind (see GCBO)

handles.ind_sel_table_results = eventdata.Indices(:,1);

% Update handles structure
guidata(hObject, handles);

update_curves_spectra_viz(handles)

function display_results_error_Callback(hObject, eventdata, handles)
% hObject    handle to display_results_error (see GCBO)

handles.img_axes4 = 'error_lib';
update_image_axes4(handles)

% Update handles structure
guidata(hObject, handles);

function display_results_minerals_Callback(hObject, eventdata, handles)
% hObject    handle to display_results_minerals (see GCBO)

handles.img_axes4 = 'minerals_ind';
update_image_axes4(handles)

% Update handles structure
guidata(hObject, handles);

function display_ratio_wls_Callback(hObject, eventdata, handles)
% hObject    handle to display_ratio_wls (see GCBO)

handles.img_axes4 = 'ratio_wls';
update_image_axes4(handles)

% Update handles structure
guidata(hObject, handles);

function process_dist_Callback(hObject, eventdata, handles)
handles.process_sample = 'dist';

% Update handles structure
guidata(hObject, handles);

function process_unmixing_Callback(hObject, eventdata, handles)
handles.process_sample = 'unmixing';

% Update handles structure
guidata(hObject, handles);

function process_clustering_Callback(hObject, eventdata, handles)
handles.process_sample = 'clustering';

% Update handles structure
guidata(hObject, handles);

function process_custom_DTree_Callback(hObject, eventdata, handles)
handles.process_sample = 'custom_tree';

% Update handles structure
guidata(hObject, handles);

function process_library_DTree_Callback(hObject, eventdata, handles)
handles.process_sample = 'lib_tree';

% Update handles structure
guidata(hObject, handles);


% ---- 
function test_process_sample_Callback(hObject, eventdata, handles)
% hObject    handle to test_process_sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global cfg
proc_t = tic;
try
    spectrum_types = [6,6];
    att_names = [{'VNIR'}, {'SWIR'}];

    min_max_wl = [];
    for a = 1:length(att_names)
        min_max_wl = [min_max_wl handles.(att_names{a}).min_wl_val];
        min_max_wl = [min_max_wl handles.(att_names{a}).max_wl_val];
    end
    
    % Application of the method 
    multi_class = Multi_Classifier();

    multi_class = method(multi_class,handles.process_sample,handles.img,handles.hs_lib,att_names,spectrum_types,min_max_wl,[]);

    % Creation of mineral table
    col_minerals_ind = multi_class.table_minerals_ind.Properties.VariableNames;
    set(handles.table_minerals_ind,'Data', multi_class.table_minerals_ind{:,:},'ColumnName',col_minerals_ind,'ColumnWidth',{120,35,50,50});

    if contains(['dist','unmixing','custom_tree','lib_tree'],handles.process_sample)
        % Creation and display of error table
        col_error_lib = multi_class.table_error_lib.Properties.VariableNames;
        set(handles.table_error_lib,'Data', multi_class.table_error_lib{:,:},'ColumnName',col_error_lib, 'ColumnWidth',{120,35,50,50});
    else
        % Empty the error
        empty = table();
        set(handles.table_error_lib,'Data',empty{:,:})
    end

    % Exporting results
    export_key_param(multi_class,handles.dir_data)
    
    if cfg.print_times
        fprintf('Process ')
        toc(proc_t)
    end
    handles.find_mineral = multi_class;

    % Display Image
    update_image_axes4(handles)

    % Update handles structure
    guidata(hObject, handles);
catch
    ErrorMessage = lasterr;
    errordlg(ErrorMessage,'Error');
end
