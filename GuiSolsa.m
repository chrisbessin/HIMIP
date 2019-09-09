function varargout = GuiSolsa(varargin)
% GUISOLSA MATLAB code for GuiSolsa.fig
%      GUISOLSA, by itself, creates a new GUISOLSA or raises the existing
%      singleton*.
%
%      H = GUISOLSA returns the handle to a new GUISOLSA or the handle to
%      the existing singleton*.
%
%      GUISOLSA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUISOLSA.M with the given input arguments.
%
%      GUISOLSA('Property','Value',...) creates a new GUISOLSA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GuiSolsa_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GuiSolsa_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GuiSolsa

% Last Modified by GUIDE v2.5 28-Dec-2018 15:37:21

% Author: Thanh Bui (thanh.bui@erametgroup.com)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GuiSolsa_OpeningFcn, ...
                   'gui_OutputFcn',  @GuiSolsa_OutputFcn, ...
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


% --- Executes just before GuiSolsa is made visible.
function GuiSolsa_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GuiSolsa (see VARARGIN)

% Choose default command line output for GuiSolsa
handles.output = hObject;
addpath('D:\Matlab\GUI_SOLSA');
addpath('D:\Matlab\GUI_SOLSA\Unmixing');
addpath('D:\Matlab\GUI_SOLSA\Utils');

% Display solsa logo
img = imread('images/solsa_logo1.jpg');
axes(handles.axes1); imshow(img)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GuiSolsa wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GuiSolsa_OutputFcn(hObject, eventdata, handles) 
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


% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
