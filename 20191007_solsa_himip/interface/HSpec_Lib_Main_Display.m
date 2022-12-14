function varargout = HSpec_Lib_Main_Display(varargin)
% HSPEC_LIB_MAIN_DISPLAY MATLAB code for HSpec_Lib_Main_Display.fig
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
%      applied to the GUI before HSpec_Lib_Main_Display_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HSpec_Lib_Main_Display_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HSpec_Lib_Main_Display

% Last Modified by GUIDE v2.5 12-Feb-2020 17:53:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HSpec_Lib_Main_Display_OpeningFcn, ...
                   'gui_OutputFcn',  @HSpec_Lib_Main_Display_OutputFcn, ...
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


% --- Executes just before HSpec_Lib_Main_Display is made visible.
function HSpec_Lib_Main_Display_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HSpec_Lib_Main_Display (see VARARGIN)

% Choose default command line output for GuiHSpec_Lib_multi
handles.output = hObject;

% Initialisations
run('D:\SOLSA\HIMIP\20191007_solsa_himip\utils\Config.m');

% Objects initialisations
% Library
set(handles.file_library, 'String', cfg.path_library);
set(handles.file_ext_library, 'String', cfg.path_library);

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

handles.dir_data = cfg.dir_data_fus;
handles.path_library = cfg.path_library;
handles.path_ext_library = cfg.path_library;

% Wavelengths
handles.min_wl_val = cfg.min_wl_swir;
handles.max_wl_val = cfg.max_wl_swir;
handles.lim_vnir_swir = cfg.lim_vnir_swir;

set(handles.min_wl,'String',cfg.min_wl_swir)
set(handles.max_wl,'String',cfg.max_wl_swir)

handles.ind_sel_table_lib = [];
handles.ind_sel_table_ext_lib = [];
handles.ind_sel_table_results = [];
handles.ind_sel_table_lib_check = [];
handles.mulatt_sel_pts_img = Multi_att_Lib();
handles.sel_pts_img = [];

% 
handles.button_name = [{'reflectance'},...
                       {'log10'},...
                       {'Stav_Golay'},...
                       {'hull_removal'},...
                       {'hull_curve'},...
                       {'noise_hull_rm'},...
                       {'hull_noise_rm'},...
                       {'diff_1'},...
                       {'diff_2'},...
                       {'bollinger'},...
                       {'norm_max_1'}];

handles.img_axes4 = 'minerals_ind';
handles.ROI = [];

compute_mass_weight(hObject,handles)

% UIWAIT makes HSpec_Lib_Main_Display wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HSpec_Lib_Main_Display_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% =========================================================================
% ============================== Library ==================================

function browse_library_Callback(hObject, eventdata, handles)

% Gets the library from file
[library_file, dir_path, oCancel] = uigetfile('*.*', 'Select a file',handles.path_library);
if ~oCancel
    disp('User selects Cancel')
    return
end
handles.path_library = fullfile(dir_path,library_file);

% Loading the library
handles.hs_lib = load_library(dir_path,library_file);

% Display in the table
update_lib_table(handles,'lib')

% Update handles structure
guidata(hObject, handles);

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

path_library = fullfile(dir_path,library_file);

hs_lib = Multi_att_Lib();
extension_file = library_file(end-2:end);
if extension_file == 'csv'
    % csv case
    hs_lib = load_csv(hs_lib,path_library);
else
    % other cases
    % Trying to read it as a usgs library
    try
        hs_lib = load_usgs(hs_lib,dir_path);
    catch
        error('Could not load the library. File selected is not a valid library, or its folder does not contain a USGS Library...')
    end
end

% =========================================================================
% ================ Edit Library from Hyper-Spectral Image =================

% ============================= Load file =================================

function browse_data_Callback(hObject, eventdata, handles)
% hObject    handle to browse_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Gets the spectra from file
dir_data = uigetdir(handles.dir_data);

% Update of the dir_data
handles.dir_data = dir_data;

% Load the data 
img = Multi_att_Img();
handles.img = load_envi_files(img,dir_data);

% Init coord min max
handles.xmin_val = min(handles.img.coord.x(:));
handles.xmax_val = max(handles.img.coord.x(:));
handles.ymin_val = min(handles.img.coord.y(:));
handles.ymax_val = max(handles.img.coord.y(:));
set(handles.xmin,'String',handles.xmin_val);
set(handles.xmax,'String',handles.xmax_val);
set(handles.ymin,'String',handles.ymin_val);
set(handles.ymax,'String',handles.ymax_val);

% Empty the selected points table
handles.sel_pts_img = [];
handles.mulatt_sel_pts_img = Multi_att_Lib();
handles.sel_spec_for_edit = [];

% Get the type of data
pathparts = strsplit(dir_data,'\');
handles.data_type = pathparts{end};

% Update handles structure
guidata(hObject, handles);

% Display the image
update_curves_axes1(hObject, handles)

function data_file_Callback(hObject, eventdata, handles)
% hObject    handle to data_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Default data file
data_path = strtrim(get(hObject, 'String'));

if isempty(data_path)
    data_path = fullfile(handles.dir_data,handles.file_data);
end
set(hObject, 'String', data_path);

% Update handles structure
guidata(hObject, handles);

function data_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ========================= Image Extension ===============================

function xmin_Callback(hObject, eventdata, handles)
% hObject    handle to ee (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ee as text
%        str2double(get(hObject,'String')) returns contents of ee as a double

val = str2double(get(hObject, 'String'));
if (val >= handles.xmin_val) && (val < handles.xmax_val)
   handles.xmin_val = val;
end

set(handles.xmin,'String',handles.xmin_val)

% Update handles structure
guidata(hObject, handles);

modify_hs_data(hObject,handles)
update_curves_axes1(hObject, handles)

function xmax_Callback(hObject, eventdata, handles)
% hObject    handle to fef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fef as text
%        str2double(get(hObject,'String')) returns contents of fef as a double

val = str2double(get(hObject, 'String'));
if (val > handles.xmin_val) && (val <= handles.xmax_val)
   handles.xmax_val = val;
end

set(handles.xmax,'String',handles.xmax_val)

% Update handles structure
guidata(hObject, handles);

modify_hs_data(hObject,handles)
update_curves_axes1(hObject, handles)

function ymin_Callback(hObject, eventdata, handles)
% hObject    handle to ff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ff as text
%        str2double(get(hObject,'String')) returns contents of ff as a double

val = str2double(get(hObject, 'String'));
if (val >= handles.ymin_val) && (val < handles.ymax_val)
   handles.ymin_val = val;
end

set(handles.ymin,'String',handles.ymin_val)

% Update handles structure
guidata(hObject, handles);

modify_hs_data(hObject,handles)
update_curves_axes1(hObject, handles)

function ymax_Callback(hObject, eventdata, handles)
% hObject    handle to f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f as text
%        str2double(get(hObject,'String')) returns contents of f as a double

val = str2double(get(hObject, 'String'));
if (val > handles.ymin_val) && (val <= handles.ymax_val)
   handles.ymax_val = val;
end

set(handles.ymax,'String',handles.ymax_val)

% Update handles structure
guidata(hObject, handles);

modify_hs_data(hObject,handles)
update_curves_axes1(hObject, handles)

function modify_hs_data(hObject,handles)
% Applies modification to the hyperspectral data
%

tmp = abs(handles.img.axis.x - handles.xmin_val);
xmin_idx = find(tmp == min(tmp),1);
tmp = abs(handles.img.axis.x - handles.xmax_val);
xmax_idx = find(tmp == min(tmp),1);

tmp = abs(handles.img.axis.y - handles.ymin_val);
ymin_idx = find(tmp == min(tmp),1);
tmp = abs(handles.img.axis.y - handles.ymax_val);
ymax_idx = find(tmp == min(tmp),1);

min_max_idx_xy = [xmin_idx, xmax_idx, ymin_idx, ymax_idx];

handles.img = img_modify(handles.img,1,min_max_idx_xy);

% Update handles structure
guidata(hObject, handles);


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
%coor_real(1) = - coor_real(1);
min_mat = abs(handles.img.coord.x - coor_real(1)) + abs(handles.img.coord.y - coor_real(2)); 
[~, index] = min(min_mat(:));

%[corr(1), corr(2)] = ind2sub(size(handles.img.coord.x),index)

% Set coordinates in the table
size_sel_pts_img = size(handles.sel_pts_img);
handles.sel_pts_img(size_sel_pts_img(1) + 1,:) = [coor_real, rand, index];

% Get the data from the image
multi_att_tmp = get_from_indexes(handles.img,handles.sel_pts_img(:,4),[]);

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
function min_wl_Callback(hObject, eventdata, handles)
% hObject    handle to min_wl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

min_wl_input = str2num(get(handles.min_wl,'String'));

if min_wl_input < handles.max_wl_val
    handles.min_wl_val = min_wl_input;
else
    set(handles.min_wl,'String',handles.min_wl_val)
end

% Update handles
guidata(hObject, handles)

update_curves_spectra_viz(handles)

function max_wl_Callback(hObject, eventdata, handles)
% hObject    handle to max_wl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

max_wl_input = str2num(get(handles.max_wl,'String'));

if (max_wl_input > handles.min_wl_val)
    handles.max_wl_val = max_wl_input;
else
    set(handles.max_wl,'String',handles.max_wl_val)
end

% Update handles
guidata(hObject, handles)

update_curves_spectra_viz(handles)

% =========================================================================
% =================== Mass proportion of Mineral ==========================

function st_Fe_Callback(hObject, eventdata, handles)
% hObject    handle to st_Fe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of st_Fe as text
%        str2double(get(hObject,'String')) returns contents of st_Fe as a double

compute_mass_weight(hObject,handles)

function st_Mg_Callback(hObject, eventdata, handles)
% hObject    handle to st_Mg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of st_Mg as text
%        str2double(get(hObject,'String')) returns contents of st_Mg as a double

compute_mass_weight(hObject,handles)

function st_Ni_Callback(hObject, eventdata, handles)
% hObject    handle to st_Ni (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of st_Ni as text
%        str2double(get(hObject,'String')) returns contents of st_Ni as a double

compute_mass_weight(hObject,handles)

function st_Co_Callback(hObject, eventdata, handles)
% hObject    handle to st_Co (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of st_Co as text
%        str2double(get(hObject,'String')) returns contents of st_Co as a double

compute_mass_weight(hObject,handles)

function st_Cr_Callback(hObject, eventdata, handles)
% hObject    handle to st_Cr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of st_Cr as text
%        str2double(get(hObject,'String')) returns contents of st_Cr as a double

compute_mass_weight(hObject,handles)

function st_Mn_Callback(hObject, eventdata, handles)
% hObject    handle to st_Mn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of st_Mn as text
%        str2double(get(hObject,'String')) returns contents of st_Mn as a double

compute_mass_weight(hObject,handles)

function st_Na_Callback(hObject, eventdata, handles)
% hObject    handle to st_Na (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of st_Na as text
%        str2double(get(hObject,'String')) returns contents of st_Na as a double

compute_mass_weight(hObject,handles)

function st_Ca_Callback(hObject, eventdata, handles)
% hObject    handle to st_Ca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of st_Ca as text
%        str2double(get(hObject,'String')) returns contents of st_Ca as a double

compute_mass_weight(hObject,handles)

function st_Al_Callback(hObject, eventdata, handles)
% hObject    handle to st_Al (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of st_Al as text
%        str2double(get(hObject,'String')) returns contents of st_Al as a double

compute_mass_weight(hObject,handles)

function st_Si_Callback(hObject, eventdata, handles)
% hObject    handle to st_Si (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of st_Si as text
%        str2double(get(hObject,'String')) returns contents of st_Si as a double

compute_mass_weight(hObject,handles)

function st_O_Callback(hObject, eventdata, handles)
% hObject    handle to st_O (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of st_O as text
%        str2double(get(hObject,'String')) returns contents of st_O as a double

compute_mass_weight(hObject,handles)

function st_H_Callback(hObject, eventdata, handles)
% hObject    handle to st_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of st_H as text
%        str2double(get(hObject,'String')) returns contents of st_H as a double

compute_mass_weight(hObject,handles)

function st_C_Callback(hObject, eventdata, handles)
% hObject    handle to st_C (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of st_C as text
%        str2double(get(hObject,'String')) returns contents of st_C as a double

compute_mass_weight(hObject,handles)

function st_S_Callback(hObject, eventdata, handles)
% hObject    handle to st_S (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of st_S as text
%        str2double(get(hObject,'String')) returns contents of st_S as a double

compute_mass_weight(hObject,handles)

function compute_mass_weight(hObject,handles)
%
%
global cfg

% Getting the Stoechiometric numbers
Nb_stoe = zeros(1, length(cfg.elements));
pwt_x_Nb_stoe = zeros(1, length(cfg.elements));
for e = 1:length(cfg.elements)
   Nb_stoe(e) = str2num(get(handles.(strcat('st_',cfg.elements{e})),'String'));
   pwt_x_Nb_stoe(e) = Nb_stoe(e) * cfg.elts_mass_dict(cfg.elements{e});
end

% Total molar mass
sum_mass_mol = sum(pwt_x_Nb_stoe);

% Computing the mass percentage and setting it in the interface
pwt_min_norm = zeros(1, length(cfg.elements));
for e = 1:length(cfg.elements)
   pwt_min_norm(e) = pwt_x_Nb_stoe(e)/sum_mass_mol*100;
   set(handles.(strcat('pwt_',cfg.elements{e})),'String',pwt_min_norm(e));
end
handles.pwt_min_norm_dict = containers.Map(cfg.elements,pwt_min_norm);

% Update handles
guidata(hObject,handles)


% ==================== Add Entry from HS Image ============================
function add_spec_to_Lib_Callback(hObject, eventdata, handles)
% hObject    handle to add_spec_to_Lib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

idx_entry = handles.sel_spec_for_edit;

if ~isempty(idx_entry)
    handles.entry = get_from_indexes(handles.mulatt_sel_pts_img,idx_entry,[]);
    entry_name = get(handles.mineral_new_entry,'String');
    
    % Add the modification to the Library
    handles.hs_lib = add_modif(handles.hs_lib,0,entry_name,handles.entry);
    
    update_lib_table(handles,'lib')
end

% Update handles structure
guidata(hObject, handles);


% =========================================================================
% ====================== Figures management ===============================

% ======================== Buttons display spectra ========================
function reflectance_Callback(hObject, eventdata, handles)
% hObject    handle to reflectance (see GCBO)

% Hint: get(hObject,'Value') returns toggle state of reflectance
update_curves_spectra_viz(handles)

function log10_Callback(hObject, eventdata, handles)
% hObject    handle to log10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of log10
update_curves_spectra_viz(handles)

function Stav_Golay_Callback(hObject, eventdata, handles)
% hObject    handle to Stav_Golay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Stav_Golay
update_curves_spectra_viz(handles)

function hull_removal_Callback(hObject, eventdata, handles)
% hObject    handle to hull_removal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hull_removal
update_curves_spectra_viz(handles)

function hull_curve_Callback(hObject, eventdata, handles)
% hObject    handle to hull_curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hull_curve
update_curves_spectra_viz(handles)

function noise_hull_rm_Callback(hObject, eventdata, handles)
% hObject    handle to noise_hull_rm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noise_hull_rm
update_curves_spectra_viz(handles)

function hull_noise_rm_Callback(hObject, eventdata, handles)
% hObject    handle to hull_noise_rm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hull_noise_rm
update_curves_spectra_viz(handles)

function diff_1_Callback(hObject, eventdata, handles)
% hObject    handle to diff_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of diff_1
update_curves_spectra_viz(handles)

function diff_2_Callback(hObject, eventdata, handles)
% hObject    handle to diff_2 (see GCBO)
% Hint: get(hObject,'Value') returns toggle state of diff_2

update_curves_spectra_viz(handles)

function bollinger_Callback(hObject, eventdata, handles)
% hObject    handle to bollinger (see GCBO)
% Hint: get(hObject,'Value') returns toggle state of bollinger

update_curves_spectra_viz(handles)

function compute_spectra_average_Callback(hObject, eventdata, handles)
% hObject    handle to compute_spectra_average (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of compute_spectra_average

update_curves_spectra_viz(handles)


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

function update_curves_spectra_viz(handles)

% Init
axes(handles.spectra_viz);
hold on
cla

% Building the transformation vector
handles.spect_trans = [];
for s = 1 : length(handles.button_name)
    if getfield(handles, handles.button_name{s}, 'Value')
        handles.spect_trans = [handles.spect_trans, s];
    end
end

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
    draw_curves_spectra_viz(get_from_indexes(handles.ma_clusters,handles.ind_sel_table_lib_check,[]),0,handles)
end

% Display the spectra of the results
if ~isempty(handles.ind_sel_table_results)
    draw_curves_spectra_viz(get_from_indexes(handles.find_mineral.data_error_lib,handles.ind_sel_table_results,[]),0,handles)
end

handles.spectra_viz.XLim = [handles.min_wl_val, handles.max_wl_val];

function draw_curves_spectra_viz(multi_att,mean_calc,handles)
% Update the curves visualisation when a parameter is updated
%
% Paramters
% ---------
% wavelengths: range of double (1xn)
%       wavelengths of the spectra to be displayed
% spectra: matrix of double (mxn)
%       spectra to display
% handles
%

% Initialisation 
prop_names = multi_att.att_groups(contains(multi_att.att_types,'C'));
sp_trans = handles.spect_trans;

% Parameters of the additional curves
curve_style = [{'-'} {'-'} {'-.'} {'-'} {'--'} {'--'} {'-.'} {'-'} {'-'} {'-.'}];
colors = [{'k'} {'y'} {'b'} {'r'} {'r'} {'k'} {'k'} {'b'} {'b'} {'k'}];

mean_opt = mean_calc && get(handles.compute_spectra_average,'Value');

for i = 1:length(prop_names)
    % Remove the parts of the spectra outside of the bounds
    multi_att = rm_min_max_wls(multi_att,prop_names{i},...
                               handles.min_wl_val,...
                               handles.max_wl_val);
    
    % Loop over the different spectra transformation options
    % Skip if all the wavelengths have been removed in  the previous step
    mem = ismember(multi_att.att_groups,prop_names{i});
    if max(mem)==1
        for j = 1:length(sp_trans)
            s = sp_trans(j);
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
h1.DataAspectRatio = [1,1,1];
h1.XLim = [min(coord_x), max(coord_x)];
h1.YLim = [min(coord_y), max(coord_y)];

function update_lib_table(handles,lib)
hs_lib = handles.(strcat('hs_',lib));
% Display in the table
table_min_id = table(hs_lib.entry_names);
set(handles.(strcat('table_',lib)), 'Data', table_min_id{:,:},...
                           'ColumnName',[{'Mineral_ID'}],...
                           'ColumnWidth',{200})

function update_image_axes4(handles)

% Finds the middle indexes of the image to avoid edge effect
coord_x = handles.find_mineral.ma_img.axis.x;
coord_y = handles.find_mineral.ma_img.axis.y;

fig = [{'error_lib'}, {'minerals_ind'}];

axes(handles.axes4);
hold on
cla

h1 = handles.axes4;
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

if contains(handles.img_axes4,fig)
    set(imagesc(handles.find_mineral.(strcat(handles.img_axes4,'_rgb'))),'XData', coord_x,'YData',coord_y);
    ROI = handles.find_mineral.(strcat('ROI_',handles.img_axes4));
    if ~isempty(ROI)
        plot(ROI(1:2),ROI(3)*ones(1,2),'k','LineWidth',3);
        plot(ROI(1:2),ROI(4)*ones(1,2),'k','LineWidth',3);
    end
end

% Setting the aspect ratio
h1.DataAspectRatio = [1,1,1];
h1.XLim = [min(coord_x), max(coord_x)];
h1.YLim = [min(coord_y), max(coord_y)];

export_key_param(handles.find_mineral,'D:\')

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
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

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

function test_process_sample_Callback(hObject, eventdata, handles)
% hObject    handle to test_process_sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spectrum_type = [6];
att_names = [{handles.data_type}];

min_max_wl = [handles.min_wl_val handles.max_wl_val];

% Application of the method 
find_mineral = Multi_Classifier();
find_mineral = method_2(find_mineral,handles.img,handles.hs_lib,att_names,spectrum_type,min_max_wl,[]);

% ROI
find_mineral = compute_ROI(find_mineral);

% Saving the spectra of the clusters
handles.ma_clusters = find_mineral.data_error_lib;

export_key_param(find_mineral,handles.dir_data)

% Display Tables
col_minerals_ind = find_mineral.table_minerals_ind.Properties.VariableNames;
col_error_lib = find_mineral.table_error_lib.Properties.VariableNames;
set(handles.table_minerals_ind,'Data', find_mineral.table_minerals_ind{:,:},'ColumnName',col_minerals_ind,'ColumnWidth',{120,35,50,50});
set(handles.table_error_lib,'Data', find_mineral.table_error_lib{:,:},'ColumnName',col_error_lib, 'ColumnWidth',{120,35,50,50});

handles.find_mineral = find_mineral;

% Display Image
update_image_axes4(handles)

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in process_unmixing.
function process_unmixing_Callback(hObject, eventdata, handles)
% hObject    handle to process_unmixing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spectrum_type = [6];
att_names = [{handles.data_type}];

min_max_wl = [handles.min_wl_val handles.max_wl_val];

% Application of the method 
find_mineral = Multi_Classifier();
find_mineral = method_3(find_mineral,handles.img,[],att_names,spectrum_type,min_max_wl,[]);

% Display Tables
col_minerals_ind = find_mineral.table_minerals_ind.Properties.VariableNames;
set(handles.table_minerals_ind,'Data', find_mineral.table_minerals_ind{:,:},'ColumnName',col_minerals_ind,'ColumnWidth',{120,35,50,50});
set(handles.table_error_lib,'Data',[]);
handles.find_mineral = find_mineral;

% Display Image
update_image_axes4(handles)

% Update handles structure
guidata(hObject, handles);


