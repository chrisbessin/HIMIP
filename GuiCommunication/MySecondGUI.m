function varargout = MySecondGUI(varargin)
% MYSECONDGUI MATLAB code for MySecondGUI.fig
%      MYSECONDGUI, by itself, creates a new MYSECONDGUI or raises the existing
%      singleton*.
%
%      H = MYSECONDGUI returns the handle to a new MYSECONDGUI or the handle to
%      the existing singleton*.
%
%      MYSECONDGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MYSECONDGUI.M with the given input arguments.
%
%      MYSECONDGUI('Property','Value',...) creates a new MYSECONDGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MySecondGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MySecondGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MySecondGUI

% Last Modified by GUIDE v2.5 03-Dec-2016 10:49:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MySecondGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MySecondGUI_OutputFcn, ...
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


% --- Executes just before MySecondGUI is made visible.
function MySecondGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MySecondGUI (see VARARGIN)

% A bit of validation - optional, but I like to do it, for documentation
% also
assert( numel( varargin ) == 1 && isa( varargin{1}, 'SecondGUIStuff' ),...
    'MySecondGUI:InvalidInput', 'MySecondGUI must be called with a SecondGUIStuff object as its argument' )

handles.parameters = varargin{1}; % Set the reference copy of the parameters object on handles.

% Initialise the GUI from the current state of the parameters object
set( handles.editParameter1, 'String', handles.parameters.parameter1 );
set( handles.sliderParameter2, 'Value', handles.parameters.parameter2 );
set( handles.textParameter2, 'String', num2str( handles.parameters.parameter2 ) );

% Choose default command line output for MySecondGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MySecondGUI wait for user response (see UIRESUME)
% uiwait(handles.MySecondGUI);


% --- Outputs from this function are returned to the command line.
function varargout = MySecondGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editParameter1_Callback(hObject, eventdata, handles)
% hObject    handle to editParameter1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editParameter1 as text
%        str2double(get(hObject,'String')) returns contents of editParameter1 as a double


% --- Executes during object creation, after setting all properties.
function editParameter1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editParameter1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderParameter2_Callback(hObject, eventdata, handles)
% hObject    handle to sliderParameter2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set( handles.textParameter2, 'String', num2str( get( hObject, 'Value' ) ) );

% --- Executes during object creation, after setting all properties.
function sliderParameter2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderParameter2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pbApply.
function pbApply_Callback(hObject, eventdata, handles)
% hObject    handle to pbApply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Simply assign the paremeters from the UI to the stored parameters object
handles.parameters.parameter1 = get( handles.editParameter1, 'String' );
handles.parameters.parameter2 = get( handles.sliderParameter2, 'Value' );
handles.parameters.notify( 'StuffChanged' ); % Notify anything listening that something has changed so that they can react appropriately.
