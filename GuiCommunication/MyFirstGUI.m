function varargout = MyFirstGUI(varargin)
% MYFIRSTGUI MATLAB code for MyFirstGUI.fig
%      MYFIRSTGUI, by itself, creates a new MYFIRSTGUI or raises the existing
%      singleton*.
%
%      H = MYFIRSTGUI returns the handle to a new MYFIRSTGUI or the handle to
%      the existing singleton*.
%
%      MYFIRSTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MYFIRSTGUI.M with the given input arguments.
%
%      MYFIRSTGUI('Property','Value',...) creates a new MYFIRSTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MyFirstGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MyFirstGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MyFirstGUI

% Last Modified by GUIDE v2.5 03-Dec-2016 10:46:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MyFirstGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MyFirstGUI_OutputFcn, ...
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


% --- Executes just before MyFirstGUI is made visible.
function MyFirstGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MyFirstGUI (see VARARGIN)

handles.parameters = SecondGUIStuff; % Create the parameters object which will be shared by reference.
handles.parameterListener = addlistener( handles.parameters, 'StuffChanged', @(src,obj) updateGUI( handles.MyFirstGUI ) ); % Listen for the StuffChanged event and trigger the updateGUI callback in response
guidata( hObject, handles ) % Update handles structure so we can call updateGUI

updateGUI( handles.MyFirstGUI ) % Update the GUI initial state from the default parameter object - or if the parameter object were given non-default values this will alwas ensure UI is in sync at start

% Choose default command line output for MyFirstGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MyFirstGUI wait for user response (see UIRESUME)
% uiwait(handles.MyFirstGUI);


% --- Outputs from this function are returned to the command line.
function varargout = MyFirstGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pbLaunch.
function pbLaunch_Callback(hObject, eventdata, handles)
% hObject    handle to pbLaunch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MySecondGUI( handles.parameters ) % Launch the 2nd GUI passing in a reference copy of the parameters


function updateGUI( hGUI )

handles = guidata( hGUI ); % Get the handles from the main GUI handle - never pass handles structure directly into a callback like this
set( handles.textParameter1, 'String', num2str( handles.parameters.parameter1 ) )
set( handles.textParameter2, 'String', handles.parameters.parameter2 )
