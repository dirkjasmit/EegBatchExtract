
function varargout = EegBatchReview(varargin)
% EEGBATCHREVIEW MATLAB code for EegBatchReview.fig
%      EEGBATCHREVIEW, by itself, creates a new EEGBATCHREVIEW or raises the existing
%      singleton*.
%
%      H = EEGBATCHREVIEW returns the handle to a new EEGBATCHREVIEW or the handle to
%      the existing singleton*.
%
%      EEGBATCHREVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EEGBATCHREVIEW.M with the given input arguments.
%
%      EEGBATCHREVIEW('Property','Value',...) creates a new EEGBATCHREVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EegBatchReview_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EegBatchReview_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EegBatchReview

% Last Modified by GUIDE v2.5 06-May-2025 20:12:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EegBatchReview_OpeningFcn, ...
                   'gui_OutputFcn',  @EegBatchReview_OutputFcn, ...
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

end

% --- Executes just before EegBatchReview is made visible.
function EegBatchReview_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EegBatchReview (see VARARGIN)

% Choose default command line output for EegBatchReview
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EegBatchReview wait for user response (see UIRESUME)
% uiwait(handles.MainGUI);

end

% --- Outputs from this function are returned to the command line.
function varargout = EegBatchReview_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = guidata(hObject);

FilterSpec = {'*.set', 'EEGLAB'};
fid = fopen('.EegBatchReview_FilePath.ini','r');
defaultpath = ".";
if fid>0
    try
        defaultpath = fgetl(fid);
        fclose(fid);
    catch
    end
end 
[filenames, pathname, ~] = uigetfile(FilterSpec, ...
    'Select an EEGLAB file', defaultpath, ...
    'multiselect', 'on');

% Check if the user selected files or canceled
if isequal(filenames, 0)
    disp('User canceled file selection.');
    return;
end

% Ensure filenames is a cell array (single selection returns a char)
if ischar(filenames)
    filenames = {filenames};
end

% save the path when not cancelled
fid = fopen('.EegBatchReview_FilePath.ini','w');
if fid>0
    fprintf(fid,'%s',pathname);
    fclose(fid);
end

% save in struct and wait for next buttonpress
data.filenames = filenames;
data.pathname = pathname;

data.editSelectPath.String = pathname;
data.editSelectFileNum.String = sprintf('%d files', length(filenames));

guidata(hObject,data)

end

% end of function ---------------------------------------------------------




function editSelectPath_Callback(hObject, eventdata, handles)
% hObject    handle to editSelectPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSelectPath as text
%        str2double(get(hObject,'String')) returns contents of editSelectPath as a double

end

% --- Executes during object creation, after setting all properties.
function editSelectPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSelectPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function editSelectFileNum_Callback(hObject, eventdata, handles)
% hObject    handle to editSelectFileNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSelectFileNum as text
%        str2double(get(hObject,'String')) returns contents of editSelectFileNum as a double

end

% --- Executes during object creation, after setting all properties.
function editSelectFileNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSelectFileNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = guidata(hObject);

% Store batch state
data.batch.index = 1;
data.batch.total = length(data.filenames);
data.batch.hObject = hObject;  % for later use in callbacks

guidata(hObject, data);

% Start batch process
startReviewForNextFile(hObject);

end
