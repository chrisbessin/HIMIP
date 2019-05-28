function varargout = GuiRegistration(varargin)
% GUIREGISTRATION MATLAB code for GuiRegistration.fig
%      GUIREGISTRATION, by itself, creates a new GUIREGISTRATION or raises the existing
%      singleton*.
%
%      H = GUIREGISTRATION returns the handle to a new GUIREGISTRATION or the handle to
%      the existing singleton*.
%
%      GUIREGISTRATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIREGISTRATION.M with the given input arguments.
%
%      GUIREGISTRATION('Property','Value',...) creates a new GUIREGISTRATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GuiRegistration_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GuiRegistration_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GuiRegistration

% Last Modified by GUIDE v2.5 22-Mar-2019 14:50:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GuiRegistration_OpeningFcn, ...
                   'gui_OutputFcn',  @GuiRegistration_OutputFcn, ...
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


% --- Executes just before GuiRegistration is made visible.
function GuiRegistration_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GuiRegistration (see VARARGIN)

% Choose default command line output for GuiRegistration
handles.output = hObject;

% Add paths
if(~isdeployed)
    addpath('D:\Matlab\GUI_SOLSA\Utilities')
    addpath('D:\Matlab\GUI_SOLSA\myprivate')
end

% Display solsa logo
if isdeployed
    img = imread('solsa_logo1.jpg');
else
    img = imread('images\solsa_logo1.jpg');
end
axes(handles.axes1); imshow(img)

% Root path
handles.currentPath = 'C:\';
handles.hyperData = 0;
handles.loadedCP = 0;

% Turn off axes
set(handles.axesFix, 'XColor','none','YColor','none');
set(handles.axesMoving, 'XColor','none','YColor','none');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GuiRegistration wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GuiRegistration_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in fixOpenPb.
function fixOpenPb_Callback(hObject, eventdata, handles)
% hObject    handle to fixOpenPb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileName, filePath, notCancel] = uigetfile({'*.raw; *.png; *.jpg'}, 'Select a file', handles.currentPath);
if (~notCancel)
    disp('Please select corresponding data');
    return
end
handles.currentPath = filePath;
dataFile = fullfile(filePath, fileName);
if strcmp(fileName(end-2:end), 'raw')
    fprintf('Loading hyperspectral data \n');
    [fixedData, fixedInfo, fixedImg] = access_spectra_data(dataFile);
    handles.fixedData = fixedData;
    handles.fixedInfo = fixedInfo;
    assignin('base', 'fixedInfo', fixedInfo);
    handles.hyperData = 1;
    
else
    fprintf('Loading rgb image \n');
    fixedImg = imread(fullfile(filePath, fileName));
    
end
set(handles.fixedDataText, 'String', fileName)
axes(handles.axesFix); imshow(fixedImg); axis on
handles.fixedImg = fixedImg;


% Update handles
guidata(hObject, handles)


% --- Executes on button press in movingOpenPb.
function movingOpenPb_Callback(hObject, eventdata, handles)
% hObject    handle to movingOpenPb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileName, filePath, notCancel] = uigetfile({'*.raw; *.png; *.jpg'}, 'Select a file', handles.currentPath);
if (~notCancel)
    disp('Please select corresponding data');
    return
end
handles.currentPath = filePath;
dataFile = fullfile(filePath, fileName);
if strcmp(fileName(end-2:end), 'raw')
    fprintf('Loading hyperspectral data \n');
    [movingData, movingInfo, movingImg] = access_spectra_data(dataFile);
    handles.movingData = movingData;
    handles.movingInfo = movingInfo;

    assignin('base', 'movingData', movingData)
    assignin('base', 'movingInfo', movingInfo)
    handles.hyperData = 1;
else
    fprintf('Loading rgb image \n');
    movingImg = imread(fullfile(filePath, fileName));   
end
set(handles.movingDataText, 'String', fileName)
axes(handles.axesMoving); imshow(movingImg); axis on
handles.movingImg = movingImg;


% Update handles
guidata(hObject, handles)

% --- Executes on button press in executePb.
function executePb_Callback(hObject, eventdata, handles)
% hObject    handle to executePb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(get(handles.existingCPRb, 'Value'))
    if (handles.loadedCP == 0)
        fprintf('Using existing control points \n');
        if isdeployed
            load ControlPoints.mat
        else
            load Data/ControlPoints.mat
        end
        controlPoints.fixedPoints = fixedPoints;
        controlPoints.movingPoints = movingPoints;
        handles.controlPoints = controlPoints; 
    else
        fprintf('Using loaded control points \n')
    end
end

try
    controlPoints = handles.controlPoints;
    fixedImg = handles.fixedImg;
    movingImg = handles.movingImg;
catch
    return
end


fprintf('Computing the transformation \n')
movingPoints = controlPoints.movingPoints;
fixedPoints = controlPoints.fixedPoints;
% Compute the transform
if length(movingPoints) < 4
    msgbox('At least 4 control points must be selected!', 'Error', 'error')
    return
end

% Select a transformation type
contents = cellstr(get(handles.transTypePM, 'String'));
transType = contents{get(handles.transTypePM, 'Value')}

tform = fitgeotrans(movingPoints, fixedPoints, transType);

% Constrain the transformed image to the same number of rows, columns and spatial limits as the fixed image. 
Rfixed = imref2d(size(fixedImg));
registeredImg = imwarp(movingImg,tform,'FillValues', max(movingImg(:)), 'OutputView',Rfixed);
handles.registeredImg = registeredImg;
if (handles.hyperData)
    registeredData = imwarp(handles.movingData, tform, 'FillValues', max(handles.movingData(:)), 'OutputView', Rfixed);
    handles.registeredData = registeredData;
end


if(get(handles.existingCPRb, 'Value'))
    figure(1);
else
    figure(2);
end
subplot 121, imshow(registeredImg), axis on, title('Registered Image')
subplot 122, imshow(fixedImg), axis on, title('Fixed Image')


% Update handles
guidata(hObject, handles)

% --- Executes on button press in existingCPRb.
function existingCPRb_Callback(hObject, eventdata, handles)
% hObject    handle to existingCPRb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of existingCPRb

if(get(hObject, 'Value'))
    fprintf('Using existing control points \n');
    if isdeployed
        load ControlPoints.mat
        set(handles.cpEdit, 'String', 'ControlPoints')
    else
        load Data/ControlPoints.mat
        set(handles.cpEdit, 'String', 'Data/ControlPoints')
    end
    controlPoints.fixedPoints = fixedPoints;
    controlPoints.movingPoints = movingPoints;
    handles.controlPoints = controlPoints;   
else
    return
end

% Update handles
guidata(hObject, handles)

% --- Executes on button press in selectCPRb.
function selectCPRb_Callback(hObject, eventdata, handles)
% hObject    handle to selectCPRb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of selectCPRb

try 
    movingImg = handles.movingImg;
    fixedImg = handles.fixedImg;
catch
    return
end

if (get(hObject, 'Value'))
    fprintf('Select control points from images\n');
    % Select control points
    [movingPoints, fixedPoints] = Guicpselect(movingImg, fixedImg, 'Wait',true);
    if length(movingPoints) < 4
        msgbox('At least 4 control points must be selected!', 'Error', 'error')
    end
    fprintf('Control point selection is finished \n')
    controlPoints.fixedPoints = fixedPoints;
    controlPoints.movingPoints = movingPoints;
    handles.controlPoints = controlPoints; 
           
else
    return
end

% Update handles
guidata(hObject, handles)


% --- Executes on button press in closePb.
function closePb_Callback(hObject, eventdata, handles)
% hObject    handle to closePb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get(handles.output, 'Tag') is the 'Tag' of the GUI
Figures = findobj('Type','Figure','-not','Tag',get(handles.output,'Tag'));
close(Figures)


% --- Executes on button press in savepb.
function savePb_Callback(hObject, eventdata, handles)
% hObject    handle to savepb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.hyperData) % Hyperspectral data
    [fileName, filePath, oCancel] = uiputfile('.raw', 'Specify a file', handles.currentPath);
    if ~oCancel
        disp('User selects Cancel')
        return
    end

    handles.currentPath = filePath; % Update current path
    registeredDataPath = fullfile(filePath, fileName);
    registeredHdrPath = strcat(registeredDataPath(1:end-4), '.hdr');


    info = handles.movingInfo;
    info.samples = size(handles.registeredData,2);
    precision = 'single';
    info.data_type = precision2datatype(precision);
    if (~get(handles.fusedDataCb, 'Value'))
        data = single(handles.registeredData);  
    else  % Save fused data
        data = cat(3, handles.registeredData, handles.fixedData);
        if(mean(handles.movingInfo.Wavelength) <= mean(handles.fixedInfo.Wavelength))
            info.Wavelength = [handles.movingInfo.Wavelength; handles.fixedInfo.Wavelength];
        else
            info.Wavelength = [handles.fixedInfo.Wavelength; handles.movingInfo.Wavelength];
        end
        info.bands = handles.movingInfo.bands + handles.fixedInfo.bands;
    end
    % Write file
    multibandwrite(data, registeredDataPath, 'bil', 'precision', precision);
    write_envihdr(info, registeredHdrPath);  
else
    [fileName, filePath, oCancel] = uiputfile({'*.jpg'; '*.png'}, 'Specify a file', handles.currentPath);
    if ~oCancel
        disp('User selects Cancel')
        return
    end

    handles.currentPath = filePath; % Update current path
    registeredDataPath = fullfile(filePath, fileName);
    imwrite(handles.registeredImg, registeredDataPath)
end
    
    
msgbox('Registered data have been saved successfully.')

% Update handles structures
guidata(hObject, handles)


% --- Executes on button press in loadCPPb.
function loadCPPb_Callback(hObject, eventdata, handles)
% hObject    handle to loadCPPb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.existingCPRb, 'Value')
    [fileName, filePath, oCancel] = uigetfile('.mat', 'Select control points', './');
    if(~oCancel)
        disp ('User selects Cancel')
        return
    end
    cpPath = fullfile(filePath, fileName);
    load(cpPath)
    set(handles.cpEdit, 'String', cpPath)
    
    controlPoints.fixedPoints = fixedPoints;
    controlPoints.movingPoints = movingPoints;
    handles.controlPoints = controlPoints; 
    handles.loadedCP = 1;
   
else
    return
end

% Update structures
guidata(hObject, handles)




function cpEdit_Callback(hObject, eventdata, handles)
% hObject    handle to cpEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cpEdit as text
%        str2double(get(hObject,'String')) returns contents of cpEdit as a double


% --- Executes during object creation, after setting all properties.
function cpEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cpEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveCpPb.
function saveCpPb_Callback(hObject, eventdata, handles)
% hObject    handle to saveCpPb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (get(handles.selectCPRb, 'Value'))

    [fileName, filePath, oCancel] = uiputfile('.mat', 'Specify a file', './');
    if ~oCancel
        disp ('User selects Cancel')
        return
    end
    try
        cpPath = fullfile(filePath, fileName);
        movingPoints = handles.controlPoints.movingPoints;
        fixedPoints = handles.controlPoints.fixedPoints;
        save(cpPath, 'movingPoints', 'fixedPoints');
    catch
        return
    end
else
    return
end


% --- Executes on button press in fusedDataCb.
function fusedDataCb_Callback(hObject, eventdata, handles)
% hObject    handle to fusedDataCb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fusedDataCb


% --- Executes on selection change in transTypePM.
function transTypePM_Callback(hObject, eventdata, handles)
% hObject    handle to transTypePM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns transTypePM contents as cell array
%        contents{get(hObject,'Value')} returns selected item from transTypePM
contents = cellstr(get(hObject, 'String'));
transType = contents{get(hObject, 'Value')}

% --- Executes during object creation, after setting all properties.
function transTypePM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to transTypePM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
