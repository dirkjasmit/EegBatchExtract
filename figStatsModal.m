function varargout = figStatsModal(varargin)
% FIGSTATSMODAL MATLAB code for figStatsModal.fig
%      FIGSTATSMODAL, by itself, creates a new FIGSTATSMODAL or raises the existing
%      singleton*.
%
%      H = FIGSTATSMODAL returns the handle to a new FIGSTATSMODAL or the handle to
%      the existing singleton*.
%
%      FIGSTATSMODAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIGSTATSMODAL.M with the given input arguments.
%
%      FIGSTATSMODAL('Property','Value',...) creates a new FIGSTATSMODAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before figStatsModal_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to figStatsModal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help figStatsModal

% Last Modified by GUIDE v2.5 28-Mar-2025 10:46:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @figStatsModal_OpeningFcn, ...
                   'gui_OutputFcn',  @figStatsModal_OutputFcn, ...
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


% --- Executes just before figStatsModal is made visible.
function figStatsModal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to figStatsModal (see VARARGIN)

% Choose default command line output for figStatsModal
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes figStatsModal wait for user response (see UIRESUME)
% uiwait(handles.figureAnalysis);


% --- Outputs from this function are returned to the command line.
function varargout = figStatsModal_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonOK.
function pushbuttonOK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function editVar1_Callback(hObject, eventdata, handles)
% hObject    handle to editVar1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editVar1 as text
%        str2double(get(hObject,'String')) returns contents of editVar1 as a double


% --- Executes during object creation, after setting all properties.
function editVar1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editVar1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editVar2_Callback(hObject, eventdata, handles)
% hObject    handle to editVar2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editVar2 as text
%        str2double(get(hObject,'String')) returns contents of editVar2 as a double


% --- Executes during object creation, after setting all properties.
function editVar2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editVar2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuAnalysi.
function popupmenuAnalysi_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuAnalysi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuAnalysi contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuAnalysi


% --- Executes during object creation, after setting all properties.
function popupmenuAnalysi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuAnalysi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figureAnalysis.
function figureAnalysis_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figureAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    data = guidata(hObject);
    data.Var1 = data.editVar1.String;
    data.Var2 = data.editVar2.String;
    data.Analysis = data.popupmenuAnalysis.Value;
    guidata(hObject, data);
catch
    warning('Failed to save the data')
end

% Hint: delete(hObject) closes the figure
delete(hObject);
