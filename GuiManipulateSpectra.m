function varargout = GuiManipulateSpectra(varargin)
% GUIMANIPULATESPECTRA MATLAB code for GuiManipulateSpectra.fig
%      GUIMANIPULATESPECTRA, by itself, creates a new GUIMANIPULATESPECTRA or raises the existing
%      singleton*.
%
%      H = GUIMANIPULATESPECTRA returns the handle to a new GUIMANIPULATESPECTRA or the handle to
%      the existing singleton*.
%
%      GUIMANIPULATESPECTRA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIMANIPULATESPECTRA.M with the given input arguments.
%
%      GUIMANIPULATESPECTRA('Property','Value',...) creates a new GUIMANIPULATESPECTRA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GuiManipulateSpectra_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GuiManipulateSpectra_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GuiManipulateSpectra

% Last Modified by GUIDE v2.5 21-May-2019 09:29:32

% Author: Thanh Bui (thanh.bui@erametgroup.com)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GuiManipulateSpectra_OpeningFcn, ...
                   'gui_OutputFcn',  @GuiManipulateSpectra_OutputFcn, ...
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


% --- Executes just before GuiManipulateSpectra is made visible.
function GuiManipulateSpectra_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GuiManipulateSpectra (see VARARGIN)

% Choose default command line output for GuiManipulateSpectra
handles.output = hObject;
% Add paths
if(~isdeployed)
    addpath('D:\Matlab\GUI_SOLSA\Utilities')
end

% Display solsa logo
if isdeployed
    img = imread('solsa_logo1.jpg');
else
    img = imread('images/solsa_logo1.jpg');
end
axes(handles.axes2); imshow(img)

% Store data to handles for exchanging among different functions
handles.data = 0;
handles.info = struct('Wavelength', ' ');
handles.rgb_img_org = 0;
handles.rgb_img = 0;
handles.rgb_img_disp = 0;
handles.spectraFig = ' ';
handles.key = ' ';
handles.currentPath = 'C:\';

% For average data
handles.xStep = 1;
handles.yStep = 1;

handles.spectra = struct('pos', ' ', 'spectrum', 0, 'wavelength', 0);
handles.fileName = ' ';

set(handles.axes1, 'XColor','none','YColor','none');


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GuiManipulateSpectra wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GuiManipulateSpectra_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% =======================================================================
% 1) Open a hyperspectral image and display it
% --- Executes on button press in openPb.
function openPb_Callback(hObject, eventdata, handles)
% hObject    handle to openPb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.statusText, 'String', 'Loading data ...')

[filename, path, oCheck] = uigetfile('.raw', 'Select a file', handles.currentPath);
if ~oCheck
    disp('User selected Cancel')
    set(handles.statusText, 'String', ' ')
    return
end
handles.currentPath = path;
dataFile = fullfile(path, filename);
set(handles.fileNameST, 'String', filename);

[data, info, rgb_img_org] = access_spectra_data(dataFile);
handles.data = data;
handles.info = info;
handles.rgb_img_org = rgb_img_org;
handles.rgb_img = rgb_img_org;
handles.fileName = filename;

% Band selection
if mean(info.Wavelength > 1100)
    handles.r_wl = 2000; 
    handles.g_wl = 2200; 
    handles.b_wl = 2350;
else
    handles.r_wl = 700; 
    handles.g_wl = 600; 
    handles.b_wl = 500;
end
set(handles.rSlider, 'Value', (handles.r_wl-min(info.Wavelength))/(max(info.Wavelength)-min(info.Wavelength)));
set(handles.rEdit, 'String', handles.r_wl);
set(handles.gSlider, 'Value', (handles.g_wl-min(info.Wavelength))/(max(info.Wavelength)-min(info.Wavelength)));
set(handles.gEdit, 'String', handles.g_wl);
set(handles.bSlider, 'Value', (handles.b_wl-min(info.Wavelength))/(max(info.Wavelength)-min(info.Wavelength)));
set(handles.bEdit, 'String', handles.b_wl);

% img slider
set(handles.imgSlider, 'Value', 0.5);

%handles.spectra = struct();
handles.spectra = struct('pos', ' ', 'spectrum', 0, 'wavelength', 0);

%axes_position = get(handles.axes1, 'Position')
%set(hObject, 'Units', 'pixels');
%[imHeight, imWidth, ~] = size(rgb_img);
%position = get(hObject, 'Position')
%set(hObject, 'Position', [position(1:2) imWidth + 100 imHeight + 100]);

% Display image
axes(handles.axes1);
cla(handles.axes1, 'reset')
handles.rgb_img_disp = rgb_img_org;
hImg = imshow(handles.rgb_img_disp); axis on;


% set(handles.axes1, ...
%     'Visible', 'on', ...
%     'Units', 'pixels', ...
%     'Position', [0 0 imWidth imHeight]);


%get(img, 'CurrentCharacter')
% When click on the figure, jump to click function
set(hImg, 'ButtonDownFcn', {@axes1_ButtonDownFcn, hObject})

% Finish loading the data
set(handles.statusText, 'String', ' ')

%Update handles structure
guidata(hObject, handles)


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
coordinate = get(gca, 'CurrentPoint');
data = handles.data;
info = handles.info;
xStep = handles.xStep; yStep = handles.yStep;

dataSize = size(data);


% Get the coordinate of cursor when mouse clicked
corr = round(coordinate(1, 1:2));
% Extract spectrum and compute its continuum removal
wavelength = info.Wavelength;
if strcmp(wavelength, ' ')
    return
end

startX = corr(1) - xStep;
if startX < 1
    startX = 1;
end
endX = corr(1) + xStep;
if endX > dataSize(2)
    endX = dataSize(2);
end

startY = corr(2) - yStep;
if startY < 1
    startY = 1;
end
endY = corr(2) + yStep;
if endY > dataSize(1)
    endY = dataSize(1);
end

temp = data(startY:endY, startX:endX, :);
size(temp)
refl = squeeze(mean(mean(temp,1),2));


%refl = squeeze(data(corr(2), corr(1), :));

% Reduce noise
order = get(handles.orderPm, 'Value');
frameLen = round(get(handles.framelenSlider, 'Value')) + 1;
if frameLen <= order
    frameLen = order + mod(order,2) + 1;
elseif ~mod(frameLen, 2)
    frameLen = frameLen + 1;
end
try
    spectrum = sgolayfilt(refl, order, frameLen);
    [spectrum_cr, ~, ~] = ContinuumRemovalSeg (wavelength, spectrum, 10);
catch
    msgbox('Please load data first')
    return
end

% Store wavelength and spectra
spectra.pos = sprintf('Pos: %d, %d', corr(1), corr(2));
spectra.spectrum = spectrum;
spectra.wavelength = wavelength;

%ctrlPressed = strcmp(handles.key, 'control');
ctrlPressed = 1;
% Store to handles
if (ctrlPressed)
    fprintf('Control is pressed\n')
    handles.spectra = [handles.spectra,spectra];
    %resetPb_Callback(handles.resetPb, eventdata, handles);

else      
    handles.spectra = struct('pos', ' ', 'spectrum', 0, 'wavelength', 0);
    handles.spectra(2) = spectra;
    
end

% Display variable to workspace
assignin('base', 'spectra', handles.spectra)
% Display on the spectra list, start from 2 to avoid the first empty data
set(handles.specListLb, 'String', {handles.spectra(2:end).pos}, 'Value', 1)


% Plot a point when left button clicked
hold on, scatter(corr(1), corr(2), '.');

% Display spectrum and its continuum removal
handles.spectraFig = figure(1);
if (ctrlPressed)
    hold on, subplot 211, plot(wavelength, spectrum, 'DisplayName', sprintf('Pos: %d, %d', corr(1), corr(2))),
else
    hold off, subplot 211, plot(wavelength, spectrum, 'DisplayName', sprintf('Pos: %d, %d', corr(1), corr(2))),
end

ylabel('Reflectance'), 
legend('show')
if (ctrlPressed)
    hold on, subplot 212, plot(wavelength, spectrum_cr, 'DisplayName', sprintf('Pos: %d, %d', corr(1), corr(2))), 
else
    hold off, subplot 212, plot(wavelength, spectrum_cr, 'DisplayName', sprintf('Pos: %d, %d', corr(1), corr(2))), 
end

ylabel('Hull removal'), 
xlabel('Wavelength (nm)')
legend('show')
% Display the cursor location
set(handles.spectraFig, 'WindowButtonMotionFcn', @(obj, event)cursorLocation(obj, event, 'BottomLeft', ' (%.1f, %.2f)', 'k'))

% Update handles
guidata(hObject, handles)


% --- Executes on button press in resetPb.
function resetPb_Callback(hObject, eventdata, handles)
% hObject    handle to resetPb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Display image

set(handles.framelenSlider, 'Value', 3);
set(handles.orderPm, 'Value', 2);
axes(handles.axes1)
cla(handles.axes1, 'reset')
handles.rgb_img_disp = handles.rgb_img_org;
img = imshow(handles.rgb_img_disp); axis on

% reset struct of spectra
handles.spectra = struct('pos', ' ', 'spectrum', 0, 'wavelength', 0);

% Reset spectral list
set(handles.specListLb, 'String', ' ', 'Value', 1)
set(handles.imgSlider, 'Value', 0.5);

% When click on the figure, jump to click function
set(img, 'ButtonDownFcn', {@axes1_ButtonDownFcn, hObject})


% Reset the band selection
wavelength = handles.info.Wavelength;
if mean(wavelength > 1100)
    handles.r_wl = 2000; 
    handles.g_wl = 2200; 
    handles.b_wl = 2350;
else
    handles.r_wl = 700; 
    handles.g_wl = 600; 
    handles.b_wl = 500;
end
set(handles.rSlider, 'Value', (handles.r_wl-min(wavelength))/(max(wavelength)-min(wavelength)));
set(handles.rEdit, 'String', handles.r_wl);
set(handles.gSlider, 'Value', (handles.g_wl-min(wavelength))/(max(wavelength)-min(wavelength)));
set(handles.gEdit, 'String', handles.g_wl);
set(handles.bSlider, 'Value', (handles.b_wl-min(wavelength))/(max(wavelength)-min(wavelength)));
set(handles.bEdit, 'String', handles.b_wl);

set(handles.winSizeX, 'String', 1);
set(handles.winSizeY, 'String', 1);
handles.xStep = 1;
handles.yStep = 1;


% Close all figures
Figures = findobj('Type','Figure','-not','Tag',get(handles.output,'Tag'));
close(Figures)

% Update handles
guidata(hObject, handles)


% --- Executes on button press in savePb.
function savePb_Callback(hObject, eventdata, handles)
% hObject    handle to savePb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
spectra = handles.spectra;
if length(spectra) < 2
    msgbox(sprintf('Please select spectra'), 'Error', 'error')
    return
end
imgFile = handles.fileName;
saveAll = get(handles.saveAllCb, 'Value');
% Plus 1 because the checklist start from spectra 2, refer line 192
spectraIndices = get(handles.specListLb, 'Value') + 1; 
[filename, path] = uiputfile('.txt', 'Specify a file', handles.currentPath);
if ~filename  % filename is empty
    return
end
handles.currentPath = path;
filePath = fullfile(path, filename);
write_spectra(filePath, imgFile, spectra, spectraIndices, saveAll);

% Update handles
guidata(hObject, handles)




% --- Executes on button press in saveAllCb.
function saveAllCb_Callback(hObject, eventdata, handles)
% hObject    handle to saveAllCb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of saveAllCb


% --- Executes on button press in openSpectraPb.
function openSpectraPb_Callback(hObject, eventdata, handles)
% hObject    handle to openSpectraPb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, path] = uigetfile('.txt', 'Select a file', handles.currentPath);
if ~filename  % filename is empty, do nothing
    return
end
handles.currentPath = path;
filePath = fullfile(path, filename);
singleSpec = get(handles.singleSpecCB, 'Value');
[data, legendTxt] = read_spectra(filePath);
if mean(data(:,1) < 10)
    data(:,1) = 1000*data(:,1);
end
if singleSpec
    legendTxt = filename(1:end-4);
end
overlay = get(handles.overlayCb, 'Value');
% Overlay or not
if overlay
    f = figure(1); % To superimpose several spectra on the same figure
else
    f = figure(2);
end
for i = 2:size(data,2)
    if singleSpec
        hold on, subplot 211, plot(data(:,1), data(:,i), 'DisplayName', legendTxt)  
    else
        hold on, subplot 211, plot(data(:,1), data(:,i), 'DisplayName', legendTxt{i-1})  
    end
    ylabel('Reflectance'), 
    l = legend('show');
    set(l, 'Interpreter','none')
    
    [spectrum_cr, ~, ~] = ContinuumRemovalSeg (data(:,1), data(:,i), 10);
    if singleSpec
        hold on, subplot 212, plot(data(:,1), spectrum_cr, 'DisplayName', legendTxt), 
    else
        hold on, subplot 212, plot(data(:,1), spectrum_cr, 'DisplayName', legendTxt{i-1}), 
    end
    ylabel('Hull removal'), 
    xlabel('Wavelength (nm)')
end
lcr = legend('show');
set(lcr,'Interpreter', 'none')

set(f, 'WindowButtonMotionFcn', @(obj, event)cursorLocation(obj, event, 'BottomLeft', ' (%.1f, %.2f)', 'k'))

% Updata handles
guidata(hObject, handles)

% --- Executes on slider movement.
function framelenSlider_Callback(hObject, eventdata, handles)
% hObject    handle to framelenSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

if ~ishandle(handles.spectraFig)
    return
end

frameLen = round(get(hObject,'Value')) + 1;
order = get(handles.orderPm, 'Value');
spectra = handles.spectra;
refl = spectra(length(spectra)).spectrum;
wavelength = spectra(length(spectra)).wavelength;

if frameLen <= order
    frameLen = order + mod(order,2) + 1;
elseif ~mod(frameLen, 2)
    frameLen = frameLen + 1;
end
spectrum = sgolayfilt(refl, order, frameLen);

[spectrum_cr, ~, ~] = ContinuumRemovalSeg (wavelength, spectrum, 10);

handles.spectraFig = figure(1); 
hold on, subplot 211, plot(wavelength, spectrum, 'DisplayName', sprintf('O:%d, FL:%d', order, frameLen)),
ylabel('Reflectance'), 
hold on, subplot 212, plot(wavelength, spectrum_cr, 'DisplayName', sprintf('O:%d, FL:%d', order, frameLen)), 
ylabel('Hull removal'), 
xlabel('Wavelength (nm)')
% Display the cursor location
set(handles.spectraFig, 'WindowButtonMotionFcn', @(obj, event)cursorLocation(obj, event, 'BottomLeft', ' (%.1f, %.2f)', 'k'))


% --- Executes during object creation, after setting all properties.
function framelenSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to framelenSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in orderPm.
function orderPm_Callback(hObject, eventdata, handles)
% hObject    handle to orderPm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns orderPm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from orderPm



% --- Executes during object creation, after setting all properties.
function orderPm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to orderPm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function imgSlider_Callback(hObject, eventdata, handles)
% hObject    handle to imgSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Display image
curVal = get(hObject, 'Value');
if curVal <= 0.5
    curVal = curVal*2;
else
    curVal = curVal*4;
end
axes(handles.axes1)
%cla(handles.axes1, 'reset')
handles.rgb_img_disp = curVal*handles.rgb_img;
img = imshow(handles.rgb_img_disp); axis on;
% When click on the figure, jump to click function
set(img, 'ButtonDownFcn', {@axes1_ButtonDownFcn, hObject})

% Update handles
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function imgSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imgSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% ====================================================
% 4) Change hyperspectral bands to obtain rgb image


% --- Executes on slider movement.
function rSlider_Callback(hObject, eventdata, handles)
% hObject    handle to rSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
try
    wavelength = handles.info.Wavelength;
catch
    return
end

r_wl = min(wavelength) + (max(wavelength) - min(wavelength))*get(hObject, 'Value');
set(handles.rEdit, 'String', r_wl);
handles.r_wl = r_wl;

% Update handles structure
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function rSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function rEdit_Callback(hObject, eventdata, handles)
% hObject    handle to rEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rEdit as text
%        str2double(get(hObject,'String')) returns contents of rEdit as a double

try
    wavelength = handles.info.Wavelength;
catch
    return
end
r_wl = str2double(get(hObject, 'String'));
if (r_wl >= min(wavelength)) && (r_wl <= max(wavelength))
    handles.r_wl = r_wl;
    r_wl_norm = (r_wl - min(wavelength))/(max(wavelength) - min(wavelength));
    set(handles.rSlider, 'Value', r_wl_norm);
else
    msgbox('Input value out of the wavelength range')
    return
end

% Update handles
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function rEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on slider movement.
function gSlider_Callback(hObject, eventdata, handles)
% hObject    handle to gSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

try
    wavelength = handles.info.Wavelength;
catch
    return
end

g_wl = min(wavelength) + (max(wavelength) - min(wavelength))*get(hObject, 'Value');
set(handles.gEdit, 'String', g_wl);
handles.g_wl = g_wl;

% Update handles
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function gSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function gEdit_Callback(hObject, eventdata, handles)
% hObject    handle to gEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gEdit as text
%        str2double(get(hObject,'String')) returns contents of gEdit as a double

try
    wavelength = handles.info.Wavelength;
catch
    return
end
g_wl = str2double(get(hObject, 'String'));
if (g_wl >= min(wavelength)) && (g_wl <= max(wavelength))
    handles.g_wl = g_wl;
    g_wl_norm = (g_wl - min(wavelength))/(max(wavelength) - min(wavelength));
    set(handles.gSlider, 'Value', g_wl_norm);
else
    msgbox('Input value out of the wavelength range')
    return
end

% Update handles
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function gEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function bSlider_Callback(hObject, eventdata, handles)
% hObject    handle to bSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

try
    wavelength = handles.info.Wavelength;
catch
    return
end

b_wl = min(wavelength) + (max(wavelength) - min(wavelength))*get(hObject, 'Value');
set(handles.bEdit, 'String', b_wl);
handles.b_wl = b_wl;

% Update handles
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function bSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function bEdit_Callback(hObject, eventdata, handles)
% hObject    handle to bEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bEdit as text
%        str2double(get(hObject,'String')) returns contents of bEdit as a double

try
    wavelength = handles.info.Wavelength;
catch
    return
end
b_wl = str2double(get(hObject, 'String'));
if (b_wl >= min(wavelength)) && (b_wl <= max(wavelength))
    handles.b_wl = b_wl;
    b_wl_norm = (b_wl - min(wavelength))/(max(wavelength) - min(wavelength));
    set(handles.bSlider, 'Value', b_wl_norm);
else
    msgbox('Input value out of the wavelength range')
    return
end

% Update handles
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function bEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in updatePb.
function updatePb_Callback(hObject, eventdata, handles)
% hObject    handle to updatePb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    data = handles.data;
    wavelength = handles.info.Wavelength;
    r_wl = handles.r_wl;
    g_wl = handles.g_wl;
    b_wl = handles.b_wl;
    rgb_img = hyperspectraldata2rgbimg(data, wavelength, r_wl, g_wl, b_wl);
catch
    msgbox('Load data first!')
    return
end


handles.rgb_img = rgb_img;

set(handles.imgSlider, 'Value', 0.5);

axes(handles.axes1)
handles.rgb_img_disp = rgb_img;
img = imshow(handles.rgb_img_disp); axis on

% When click on the figure, jump to click function
set(img, 'ButtonDownFcn', {@axes1_ButtonDownFcn, hObject})

% Update handles
guidata(hObject, handles)


% --- Executes on selection change in specListLb.
function specListLb_Callback(hObject, eventdata, handles)
% hObject    handle to specListLb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns specListLb contents as cell array
%        contents{get(hObject,'Value')} returns selected item from specListLb

get(hObject, 'Value')

% --- Executes during object creation, after setting all properties.
function specListLb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to specListLb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in overlayCb.
function overlayCb_Callback(hObject, eventdata, handles)
% hObject    handle to overlayCb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of overlayCb


% --- Executes on button press in saveImgPb.
function saveImgPb_Callback(hObject, eventdata, handles)
% hObject    handle to saveImgPb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileName, pathName, oCancel] = uiputfile({'*.png'; '*.jpg'}, 'Save as an image', handles.currentPath);
if ~oCancel
    return
end
handles.currentPath = pathName;
imgSave = handles.rgb_img_disp;
imwrite(imgSave, fullfile(pathName, fileName));

% Update handles
guidata(hObject, handles)


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
eventdata.Key;
handles.key = eventdata.Key;

% Update handles
guidata(hObject, handles)


% --- Executes on key release with focus on figure1 or any of its controls.
function figure1_WindowKeyReleaseFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was released, in lower case
%	Character: character interpretation of the key(s) that was released
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) released
% handles    structure with handles and user data (see GUIDATA)
if(strcmp(eventdata.Key, 'control'))
    handles.key = ' ';
end

% Update handles
guidata(hObject, handles)


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if(strcmp(eventdata.Key, 'control'))
    handles.key = ' ';
end



function winSizeX_Callback(hObject, eventdata, handles)
% hObject    handle to winSizeX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of winSizeX as text
%        str2double(get(hObject,'String')) returns contents of winSizeX as a double
xStep = str2double(get(hObject, 'String'));
if(isnan(xStep) || xStep < 0)
    msgbox('Please enter a positive integer', 'Error')
    return
end
handles.xStep = round((xStep-1)/2);

% Update handle structure
guidata(hObject, handles)



% --- Executes during object creation, after setting all properties.
function winSizeX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to winSizeX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function winSizeY_Callback(hObject, eventdata, handles)
% hObject    handle to winSizeY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of winSizeY as text
%        str2double(get(hObject,'String')) returns contents of winSizeY as a double

yStep = str2double(get(hObject, 'String'));
if (isnan(yStep) || yStep < 0)
    msgbox('Please enter a positive integer', 'Error')
    return
end
handles.yStep = round((yStep-1)/2);

% Update handle structure
guidata(hObject, handles) 

% --- Executes during object creation, after setting all properties.
function winSizeY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to winSizeY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function fileNameST_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in clearPB, clear spectral list
function clearPB_Callback(hObject, eventdata, handles)
% hObject    handle to clearPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes1)
cla(handles.axes1, 'reset')
img = imshow(handles.rgb_img_disp); axis on

% reset struct of spectra
handles.spectra = struct('pos', ' ', 'spectrum', 0, 'wavelength', 0);

% Reset spectral list
set(handles.specListLb, 'String', ' ', 'Value', 1)

% When click on the figure, jump to click function
set(img, 'ButtonDownFcn', {@axes1_ButtonDownFcn, hObject})

% Close all figures
Figures = findobj('Type','Figure','-not','Tag',get(handles.output,'Tag'));
close(Figures)

% Update handles
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function clearPB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clearPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in singleSpecCB.
function singleSpecCB_Callback(hObject, eventdata, handles)
% hObject    handle to singleSpecCB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of singleSpecCB


% --- Executes on button press in generateSLPb.
function generateSLPb_Callback(hObject, eventdata, handles)
% hObject    handle to generateSLPb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spectraPath = uigetdir(handles.currentPath)
if (isequal(spectraPath, 0))
    disp('User selects Cancel')
    return
end
handles.currentPath = spectraPath;
[fileName, filePath, oCheck] = uiputfile('*.mat', 'Specify a file', handles.currentPath);
if ~oCheck
    disp('User selects Cancel')
    return
end
handles.currentPath = filePath;

spectralLibFile = fullfile(filePath, fileName)

generate_spectral_library(spectraPath, spectralLibFile)

% Update handle structures
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function generateSLPb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to generateSLPb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
