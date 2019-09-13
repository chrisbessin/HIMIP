
function varargout = GuiComputeReflectance(varargin)
% GuiComputeReflectance MATLAB code for GuiComputeReflectance.fig
%      GuiComputeReflectance, by itself, creates a new GuiComputeReflectance or raises the existing
%      singleton*.
%
%      H = GuiComputeReflectance returns the handle to a new GuiComputeReflectance or the handle to
%      the existing singleton*.
%
%      GuiComputeReflectance('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GuiComputeReflectance.M with the given input arguments.
%
%      GuiComputeReflectance('Property','Value',...) creates a new GuiComputeReflectance or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GuiComputeReflectance_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GuiComputeReflectance_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GuiComputeReflectance
% Last Modified by GUIDE v2.5 27-Dec-2018 15:48:32

% Author: Thanh Bui (thanh.bui@erametgroup.com)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GuiComputeReflectance_OpeningFcn, ...
                   'gui_OutputFcn',  @GuiComputeReflectance_OutputFcn, ...
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


% --- Executes just before GuiComputeReflectance is made visible.
function GuiComputeReflectance_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GuiComputeReflectance (see VARARGIN)

% Choose default command line output for GuiComputeReflectance
handles.output = hObject;
% Add paths
if (~isdeployed)
    addpath('D:\Matlab\GUI_SOLSA\Utilities')
end


% Display solsa logo
if isdeployed
    img = imread('solsa_logo1.jpg');
else
    img = imread('images/solsa_logo1.jpg');
end
axes(handles.axes1); imshow(img)


handles.dataPath = ' ';
handles.reflectanceFullPath = ' ';
handles.samplePath = ' ';
handles.whiteRefPath = ' ';
handles.darkRefPath = ' ';
handles.reflectancePath = ' ';
handles.currentPath = 'C:\';


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GuiComputeReflectance wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GuiComputeReflectance_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%

% --- Executes on button press in dataPb.
function dataPb_Callback(hObject, eventdata, handles)
% hObject    handle to dataPb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dataPath = uigetdir(handles.currentPath);
if (isequal(dataPath, 0))
    disp('User selects Cancel')
    return
end
handles.dataPath = dataPath;
set(handles.dataText, 'String', dataPath)

handles.currentPath = dataPath;
% Update handles
guidata(hObject, handles);


function dataText_Callback(hObject, eventdata, handles)
% hObject    handle to dataText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataText as text
%        str2double(get(hObject,'String')) returns contents of dataText as a double


% --- Executes during object creation, after setting all properties.
function dataText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in reflectanceFullPb.
function reflectanceFullPb_Callback(hObject, eventdata, handles)
% hObject    handle to reflectanceFullPb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, path, oCancel] = uiputfile('.raw', 'Specify a file', handles.currentPath);
if ~oCancel
    disp('User selects Cancel')
    return
end

set(handles.reflectanceFullText, 'String', fullfile(path, filename))
handles.reflectanceFullPath = fullfile(path, filename);

handles.currentPath = path;
% Update handles
guidata(hObject, handles);

function reflectanceFullText_Callback(hObject, eventdata, handles)
% hObject    handle to reflectanceFullText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reflectanceFullText as text
%        str2double(get(hObject,'String')) returns contents of reflectanceFullText as a double


% --- Executes during object creation, after setting all properties.
function reflectanceFullText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reflectanceFullText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% -------------------- For separate sample paths --------------------------
% --- Executes on button press in samplePb.
function samplePb_Callback(hObject, eventdata, handles)
% hObject    handle to samplePb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, path, oCancel] = uigetfile('.raw', 'Select a file', handles.currentPath);
if ~oCancel
    return
end


set(handles.sampleText, 'String', fullfile(path, filename))
handles.samplePath = fullfile(path, filename);

handles.currentPath = path;
% Update handles
guidata(hObject, handles);


function sampleText_Callback(hObject, eventdata, handles)
% hObject    handle to sampleText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sampleText as text
%        str2double(get(hObject,'String')) returns contents of sampleText as a double


% --- Executes during object creation, after setting all properties.
function sampleText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sampleText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in whitePb.
function whitePb_Callback(hObject, eventdata, handles)
% hObject    handle to whitePb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, path, oCancel] = uigetfile('.raw', 'Select a file', handles.currentPath);
if ~oCancel
    disp('User selects Cancel')
    return
end

set(handles.whiteText, 'String', fullfile(path, filename))
handles.whiteRefPath = fullfile(path, filename);

handles.currentPath = path;
% Update handles
guidata(hObject, handles);


function whiteText_Callback(hObject, eventdata, handles)
% hObject    handle to whiteText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whiteText as text
%        str2double(get(hObject,'String')) returns contents of whiteText as a double


% --- Executes during object creation, after setting all properties.
function whiteText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whiteText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in darkPb.
function darkPb_Callback(hObject, eventdata, handles)
% hObject    handle to darkPb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, path, oCancel] = uigetfile('.raw', 'Select a file', handles.currentPath);
if ~oCancel
    disp('User selects Cancel')
    return
end

set(handles.darkText, 'String', fullfile(path, filename))
handles.darkRefPath = fullfile(path, filename);

handles.currentPath = path;
% Update handles
guidata(hObject, handles);

function darkText_Callback(hObject, eventdata, handles)
% hObject    handle to darkText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of darkText as text
%        str2double(get(hObject,'String')) returns contents of darkText as a double


% --- Executes during object creation, after setting all properties.
function darkText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to darkText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in reflectancePb.
function reflectancePb_Callback(hObject, eventdata, handles)
% hObject    handle to reflectancePb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, path, oCancel] = uiputfile('.raw', 'Specify a file', handles.currentPath);
if ~oCancel
    disp('User selects Cancel')
    return
end
set(handles.reflectanceText, 'String', fullfile(path, filename))
handles.reflectancePath = fullfile(path, filename);

handles.currentPath = path;
% Update handles
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function reflectancePb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reflectancePb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function reflectanceText_Callback(hObject, eventdata, handles)
% hObject    handle to darkText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of darkText as text
%        str2double(get(hObject,'String')) returns contents of darkText as a double


% --- Executes during object creation, after setting all properties.
function reflectanceText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to darkText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in computePb.
function computePb_Callback(hObject, eventdata, handles)
% hObject    handle to computePb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.statusText, 'String', 'Busy ...')
pause(1)
short_range = 1;
hullR = get(handles.hullRemovalCb, 'Value');
if (handles.dataPath ~= ' ') % Using datapath
    fprintf('Automatic selection of files \n')
    fprintf('Computing reflectance ..., please wait ...\n')
    [sampleFile, whiterefFile, darkrefFile] = SearchDataFiles (handles.dataPath);
    reflFile = handles.reflectanceFullPath;
    [refl, wavelength, rgb_img, rgb_img_he] = compute_reflectance(sampleFile, whiterefFile, darkrefFile, short_range, hullR, reflFile);
    fprintf('... Computation finished \n')
elseif (handles.samplePath ~= ' ') % Using sample, dark and white reference paths
    fprintf('Computing reflectance ..., please wait ...\n')
    sampleFile = handles.samplePath;
    whiterefFile = handles.whiteRefPath;
    darkrefFile = handles.darkRefPath;
    reflFile = handles.reflectancePath;
    [refl, wavelength, rgb_img, rgb_img_he] = compute_reflectance(sampleFile, whiterefFile, darkrefFile, short_range, hullR, reflFile);
    fprintf('... Computation finished \n')
else
    fprintf('Please specify correct paths \n')
end
set(handles.statusText, 'String', 'Free');



% --- Executes on button press in resetPb.
function resetPb_Callback(hObject, eventdata, handles)
% hObject    handle to resetPb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dataPath = ' ';
handles.reflectanceFullPath = ' ';
handles.samplePath = ' ';
handles.whiteRefPath = ' ';
handles.darkRefPath = ' ';
handles.reflectancePath = ' ';

set(handles.dataText, 'String', ' ')
set(handles.reflectanceFullText, 'String', ' ')
set(handles.sampleText, 'String', ' ')
set(handles.whiteText, 'String', ' ')
set(handles.darkText, 'String', ' ')
set(handles.reflectanceText, 'String', ' ')

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in hullRemovalCb.
function hullRemovalCb_Callback(hObject, eventdata, handles)
% hObject    handle to hullRemovalCb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hullRemovalCb


%% 3) Continuum removal computation

% --- Executes on button press in reflCRPb.
function reflCRPb_Callback(hObject, eventdata, handles)
% hObject    handle to reflCRPb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileName, filePath, oCheck] = uigetfile('*.raw', 'Select a reflectance file', handles.currentPath);
if ~oCheck
    return
end
set(handles.reflCREdit, 'String', fullfile(filePath, fileName))
handles.reflForCR_file = fullfile(filePath, fileName);

handles.currentPath = filePath;
% Update handles
guidata(hObject, handles)

function reflCREdit_Callback(hObject, eventdata, handles)
% hObject    handle to reflCREdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reflCREdit as text
%        str2double(get(hObject,'String')) returns contents of reflCREdit as a double


% --- Executes during object creation, after setting all properties.
function reflCREdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reflCREdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in hullRemovalPb.
function hullRemovalPb_Callback(hObject, eventdata, handles)
% hObject    handle to hullRemovalPb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fileName, filePath, oCheck] = uiputfile('*.raw', 'Specify a file', handles.currentPath);
if ~oCheck % User pressed Cancel
    return
end
set(handles.hullRemovalEdit, 'String', fullfile(filePath, fileName))
handles.hullRemoval_file = fullfile(filePath, fileName);

currentPath = filePath;
% Update handles
guidata(hObject, handles)

function hullRemovalEdit_Callback(hObject, eventdata, handles)
% hObject    handle to hullRemovalEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hullRemovalEdit as text
%        str2double(get(hObject,'String')) returns contents of hullRemovalEdit as a double


% --- Executes during object creation, after setting all properties.
function hullRemovalEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hullRemovalEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in hullRemovalCompPb.
function hullRemovalCompPb_Callback(hObject, eventdata, handles)
% hObject    handle to hullRemovalCompPb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    reflFile = handles.reflForCR_file;
    reflCRFile = handles.hullRemoval_file;
catch
    return
end
set(handles.statusText, 'String', 'Busy ...')
pause(1)
[refl_cr, wavelength, rgb_img_cr] = compute_hull_removal(reflFile, reflCRFile);
set(handles.statusText, 'String', 'Free')
