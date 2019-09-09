function varargout = GuiHimip(varargin)
% GUIHIMIP MATLAB code for GuiHimip.fig
%      GUIHIMIP, by itself, creates a new GUIHIMIP or raises the existing
%      singleton*.
%
%      H = GUIHIMIP returns the handle to a new GUIHIMIP or the handle to
%      the existing singleton*.
%
%      GUIHIMIP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIHIMIP.M with the given input arguments.
%
%      GUIHIMIP('Property','Value',...) creates a new GUIHIMIP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GuiHimip_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GuiHimip_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GuiHimip

% Last Modified by GUIDE v2.5 14-Jan-2019 09:35:01

% Author: Thanh Bui (thanh.bui@erametgroup.com)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GuiHimip_OpeningFcn, ...
                   'gui_OutputFcn',  @GuiHimip_OutputFcn, ...
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


% --- Executes just before GuiHimip is made visible.
function GuiHimip_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GuiHimip (see VARARGIN)

% Choose default command line output for GuiHimip
handles.output = hObject;
if (~isdeployed)
    addpath('D:\Matlab\GUI_SOLSA');
    addpath('D:\Matlab\GUI_SOLSA\Unmixing');
    addpath('D:\Matlab\GUI_SOLSA\Utilities');
end
% Display solsa logo
if isdeployed
    img = imread('solsa_logo1.jpg');
else
    img = imread('images/solsa_logo1.jpg');
end
axes(handles.axes1); imshow(img)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GuiHimip wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GuiHimip_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function reflComp_Callback(hObject, eventdata, handles)
% hObject    handle to reflComp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call reflectance computation
GuiComputeReflectance

% --------------------------------------------------------------------
function processing_Callback(hObject, eventdata, handles)
% hObject    handle to processing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function helptag_Callback(hObject, eventdata, handles)
% hObject    handle to helptag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uManual_Callback(hObject, eventdata, handles)
% hObject    handle to uManual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isdeployed
    winopen('UserManual.pdf')
else
    winopen('Data/UserManual.pdf')
end

% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GuiAbout

% --------------------------------------------------------------------
function closeTag_Callback(hObject, eventdata, handles)
% hObject    handle to closeTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all


% --------------------------------------------------------------------
function sManipulation_Callback(hObject, eventdata, handles)
% hObject    handle to sManipulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call manipulation gui
GuiManipulateSpectra


% --------------------------------------------------------------------
function registration_Callback(hObject, eventdata, handles)
% hObject    handle to registration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call registration GUI
GuiRegistration


% --------------------------------------------------------------------
function hUnmixing_Callback(hObject, eventdata, handles)
% hObject    handle to hUnmixing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call hyperspectral unmixing
GuiUnmixing


% --------------------------------------------------------------------
function reflComputation_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to reflComputation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call reflectance computation
GuiComputeReflectance


% --------------------------------------------------------------------
function spectraMan_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to spectraMan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call manipulation gui
GuiManipulateSpectra


% --------------------------------------------------------------------
function hUnmixing_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to hUnmixing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Call hyperspectral unmixing
GuiUnmixing


% --------------------------------------------------------------------
function closeTag_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to closeTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get(handles.output, 'Tag') is the 'Tag' of the GUI
Figures = findobj('Type','Figure','-not','Tag',get(handles.output,'Tag'));
close(Figures);


% --------------------------------------------------------------------
function specRegistration_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to specRegistration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call registration GUI
GuiRegistration
