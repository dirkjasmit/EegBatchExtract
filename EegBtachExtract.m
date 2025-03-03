function varargout = EegBtachExtract(varargin)
% EEGBTACHEXTRACT MATLAB code for EegBtachExtract.fig
%      EEGBTACHEXTRACT, by itself, creates a new EEGBTACHEXTRACT or raises the existing
%      singleton*.
%
%      H = EEGBTACHEXTRACT returns the handle to a new EEGBTACHEXTRACT or the handle to
%      the existing singleton*.
%
%      EEGBTACHEXTRACT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EEGBTACHEXTRACT.M with the given input arguments.
%
%      EEGBTACHEXTRACT('Property','Value',...) creates a new EEGBTACHEXTRACT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EegBtachExtract_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EegBtachExtract_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EegBtachExtract

% Last Modified by GUIDE v2.5 02-Mar-2025 11:18:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EegBtachExtract_OpeningFcn, ...
                   'gui_OutputFcn',  @EegBtachExtract_OutputFcn, ...
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


% --- Executes just before EegBtachExtract is made visible.
function EegBtachExtract_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EegBtachExtract (see VARARGIN)

% Choose default command line output for EegBtachExtract
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EegBtachExtract wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% initialise data
data = guidata(hObject);
data.EEG = eeg_emptyset();

% set the values of the uicontrols
try
    strlist = readtable(sprintf('%s.ini',get(hObject,'name')),'delimiter','\t','filetype','text');
    SetUIControlData(hObject, strlist);
catch
    warning('Initialization file not found. Will be created on close.')
end

% read in the chanlocs
try
    fn = data.edit2.String;
    if strlength(fn)>3
        fn = fn(1:strfind(fn,'(')-2);
        EEG = pop_loadset([data.edit1.String '/' fn]);
    end
    data.chanlocs = EEG.chanlocs;
catch
    warning('could not open chanlocs file')
end

% get the most recently used filenames
try
fid = fopen('.EegBatchExtract_SelectedFiles.ini','r');
if fid>0
    try
        pathname = [];
        filenames = {};
        
        pathname = fgetl(fid);
        if pathname==-1
            pathname = [];
        else
            filenames{1} = fgetl(fid);
            while filenames{end} ~= -1
                filenames{end+1} = fgetl(fid);
            end
            if filenames{end} == -1
                filenames = filenames(1:end-1);
            end
        end
    catch
        warning('Init file for file names not present (or failed reading)')
    end
    data.pathname = pathname;
    data.filenames = filenames;
    
    data.edit4.String = pathname;
    data.edit5.String = sprintf("%d files",length(filenames));
end 
catch
    warning('No such file')
end

guidata(hObject, data)



% --- Outputs from this function are returned to the command line.
function varargout = EegBtachExtract_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonSelect.
function pushbuttonSelect_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = guidata(hObject);

FilterSpec = {'*.set', 'EEGLAB'};
fid = fopen('.EegBatchExtract_FilePath.ini','r');
defaultpath = ".";
if fid>0
    try
        defaultpath = fgetl(fid);
        fclose(fid);
    catch
    end
end 
[filenames, pathname, filterindex] = uigetfile(FilterSpec, ...
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
fid = fopen('.EegBatchExtract_FilePath.ini','w');
if fid>0
    fprintf(fid,'%s',pathname);
    fclose(fid);
end

% save in struct and wait for next buttonpress
data.filenames = filenames;
data.pathname = pathname;

data.edit4.String = pathname;
data.edit5.String = sprintf('%d files', length(filenames));

guidata(hObject,data)

% end of function ---------------------------------------------------------


% --- Executes on button press in pushbuttonChanlocs.
function pushbuttonChanlocs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonChanlocs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = guidata(hObject);

FilterSpec = {'*.set', 'EEGLAB'};
fid = fopen('.EegBatchExtract_DefaultPath.ini','r');
DefaultPath = ".";
if fid>0
    try
        DefaultPath = fgetl(fid);
        fclose(fid);
    catch
    end
end 
[FileName,PathName,FilterIndex] = uigetfile(FilterSpec, ...
    'Select an EEGLAB file with channel locations', DefaultPath);

% Check if the user selected files or canceled
if isequal(FileName, 0)
    disp('User canceled file selection.');
    return;
end

% save the returned pathname
fid = fopen('.EegBatchExtract_DefaultPath.ini','w');
if fid>0
    fprintf(fid,'%s',PathName);
    fclose(fid);
end 

% put filename in editbox
data.edit1.String = PathName;
data.edit2.String = sprintf('%s (%d channels)', FileName, 0);
pause(.1)

EEG = pop_loadset([PathName '/' FileName]);
data.chanlocs = EEG.chanlocs;
data.edit2.String = sprintf('%s (%d channels)', FileName, length(data.chanlocs));
pause(.05)

guidata(hObject, data);



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonBatch.
function pushbuttonBatch_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonBatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = guidata(hObject);

% predefined
freqnames = {'alpha','beta','theta','delta','alphalo','alphahi','betalo','betahi'};
freqs = {[8 13], [13 30], [4 8], [1 3], [7 9.5], [9.5 13], [13 21], [21 30]};


% what analysis was selected?
selectedButton = data.uibuttongroup.SelectedObject;
var = selectedButton.String;

if ismember(var, {'DFA','Power'}) && data.popupmenuFreq.Value==1
    errordlg('No frequency selected!')
    return
end

pars = struct;
pars.freqname = freqnames{data.popupmenuFreq.Value-1};
pars.freq = freqs{data.popupmenuFreq.Value-1};
pars.pathname = data.pathname;
pars.var = var;
pars.chanlocs = data.chanlocs;

% DFA, IAF, and Power have the same output format.
if ismember(var, {'DFA','IAF','Power'})
    RunSingleParameterExtraction(data.filenames, pars)

elseif ismember(var, {'Coherence'})
    RunCoherence(data.filenames, pars)
    
elseif ismember(var, {'PSD'})
    RunPSD(data.filenames, pars)
    
end


function RunPSD(filenames, pars)
% function for analysis PSD

% ask for frequency bounds
dims = [1 35];  % Textbox dimensions
definput = {'1.0', '45.0'};  % Default values
answer = inputdlg({'Enter lower bound:', 'Enter upper bound:'}, ...
    'input lower and upper bound for PSD', dims, definput);
% Convert the cell array to numbers
if ~isempty(answer) % Check if user didn't cancel
    lo = round(str2double(answer{1}));
    hi = round(str2double(answer{2}));
else
    return
end

% initialize output. Save the data into a glocal variable for cohoerence.
global AllPSDfreqs
AllPSDfreqs = lo:hi;
global AllPSD
AllPSD = nan(length(lo:hi), length(pars.chanlocs));

% loop thru the files selected
for f=1:length(filenames)
    EEG = pop_loadset([pars.pathname '/' filenames{f}]);
    % remove channels not in chanlocs
    remove = find(~ismember({EEG.chanlocs.labels}, {pars.chanlocs.labels}));
    if ~isempty(remove)
        EEG = pop_select(EEG, 'nochannel', remove);
    end

    % impute the rest
    EEG = pop_interp(EEG, pars.chanlocs, 'spherical');

    % this takes a bit of time
    [P, fs] = pfft(EEG.data(:,:)', EEG.srate, hanning(EEG.srate*4), .5);

    cnt = 0;
    for frq=lo:(hi-1)
        cnt = cnt+1;
        ndx = fs>=frq & fs<(frq+1);
        AllPSD(frq,:,f) = mean(P(ndx,:),1);
    end
end

disp('Data stored in global AllPSD (freqs in AllPSDfreqs)')



function RunCoherence(filenames, pars)
% function for analysis coherence

prefix = ['COH_' pars.var];

% initialize table with no values in the table (rowcount = 0). One
% string for Id and the rest are doubles
% T = table('Size', [0, length(pars.chanlocs)+1], ...
%     'VariableTypes', ["string" repelem("double",length(pars.chanlocs))], ...
%     'VariableNames', [{'Id'} cellfun(@(x)sprintf('%s_%s',prefix,x),{pars.chanlocs.labels},'uni',0)]);

% initialize output. Save the data into a glocal variable for cohoerence.
global MeanCOH
MeanCOH = nan(length(pars.chanlocs));

% loop thru the files selected
for f=1:length(filenames)
    EEG = pop_loadset([pars.pathname '/' filenames{f}]);
    % remove channels not in chanlocs
    remove = find(~ismember({EEG.chanlocs.labels}, {pars.chanlocs.labels}));
    if ~isempty(remove)
        EEG = pop_select(EEG, 'nochannel', remove);
    end

    % impute the rest
    EEG = pop_interp(EEG, pars.chanlocs, 'spherical');

    % this takes a bit of time
    [COH, fs] = fastcoherence(EEG.data(:,:)', 'srate', EEG.srate, 'window', hanning(EEG.srate*4));

    ndx = fs>=pars.freq(1) & fs<pars.freq(2);
    
    MeanCOH(:,:,f) = mean(COH(:,:, ndx),3);
end



function RunSingleParameterExtraction(filenames, pars)
% function for analysis of EEG extracting a single parameter per
% channel

% initialize based on the analysis selected (var)
switch pars.var
    case 'IAF'
        prefix = pars.var;
        % if var is IAF, ask for lower and upper values for alpha as that
        % may be different from the dropdown menu
        dims = [1 35];  % Textbox dimensions
        definput = {'7.0', '13.0'};  % Default values
        answer = inputdlg({'Enter lower bound:', 'Enter upper bound:'}, ...
            'input frequency band for IAF estimation', dims, definput);
        % Convert the cell array to numbers
        if ~isempty(answer) % Check if user didn't cancel
            lo = str2double(answer{1});
            hi = str2double(answer{2});
        else
            return
        end

    case {'Power','DFA'}
        prefix = sprintf('%s_$s', pars.var, pars.freqname);
        % determine frequency for Power and DFA
end

% initialize table with no values in the table (rowcount = 0). One
% string for Id and the rest are doubles
T = table('Size', [0, length(pars.chanlocs)+1], ...
    'VariableTypes', ["string" repelem("double",length(pars.chanlocs))], ...
    'VariableNames', [{'Id'} cellfun(@(x)sprintf('%s_%s',prefix,x),{pars.chanlocs.labels},'uni',0)]);

% loop thru the files selected
for f=1:length(filenames)
    EEG = pop_loadset([pars.pathname '/' filenames{f}]);
    % remove channels not in chanlocs
    remove = find(~ismember({EEG.chanlocs.labels}, {pars.chanlocs.labels}));
    if ~isempty(remove)
        EEG = pop_select(EEG, 'nochannel', remove);
    end


    % impute the rest
    EEG = pop_interp(EEG, pars.chanlocs, 'spherical');

    % do calculations
    switch pars.var
        case 'DFA'
            val = dfa(abs(hilbert(filter_fir(EEG.data, EEG.srate, pars.freq(1), pars.freq(2), 3.0, true)')), EEG.srate);

        case 'IAF'
            [P,fs] = pfft(EEG.data(:,:)', EEG.srate, hanning(EEG.srate*4), .5);
            ndx = fs>=lo & fs<hi;
            val = sum(P(ndx,:).*(fs(ndx)'),1) ./ sum(P(ndx,:),1);

        case 'Power'
            [P,fs] = pfft(EEG.data(:,:)', EEG.srate, hanning(EEG.srate*4), .5);
            ndx = fs>=pars.freq(1) & pars.fs<freq(2);
            val = mean(P(ndx,:),1);
    end

    T.Id(f) = pars.filenames{f};
    T(f,2:end) = array2table(val);
end

% save the table
[outfilename, outpathname] = uiputfile('*.txt', 'Save As Tab-delimted output');
if outfilename ~= 0
    writetable(T,[outpathname '/' outfilename], 'Delimiter', '\t');
    disp('file written')
else
    warning('No file written')
end

% --------  end of function -----------------------




function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuFreq.
function popupmenuFreq_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuFreq contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuFreq


% --- Executes during object creation, after setting all properties.
function popupmenuFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

data = guidata(hObject);

% save all the settings
strlist = GetUIControlData(hObject);
writetable(strlist,sprintf('%s.ini',get(hObject,'name')),'delimiter','\t','filetype','text');

if isfield(data, 'filenames')
    fid = fopen('.EegBatchExtract_SelectedFiles.ini','w');
    if fid>0
        fprintf(fid,'%s\n',data.pathname);
        for f=1:length(data.filenames)
            fprintf(fid,'%s\n',data.filenames{f});
        end
        fclose(fid);
    end
end

delete(hObject);


% pass a handle to the gui, and it will extract all the UIControl object
% data (slider, checkbox value and edit strings).
function strlist = GetUIControlData(hObject)

strlist = struct;

ch = get(hObject,'ch');
count = 0;
for c=1:length(ch)
    if strcmpi(get(ch(c),'Type'),'uicontrol')
        skip=false;
        switch get(ch(c),'Style')
            case 'edit'
                tmp1 = sprintf('%s',get(ch(c),'tag'));
                tmp2 = sprintf('%s',get(ch(c),'string'));
            case {'checkbox','slider'}
                tmp1 = sprintf('%s',get(ch(c),'tag'));
                tmp2 = sprintf('%.4f',get(ch(c),'value'));
            case {'popupmenu'}
                tmp1 = sprintf('%s',get(ch(c),'tag'));
                tmp2 = sprintf('%d',get(ch(c),'value'));
            otherwise
                skip=true;
        end
        if ~skip
            count = count + 1;
            strlist(count).key = tmp1;
            strlist(count).val = tmp2;
        end
    end
end             

strlist = struct2table(strlist);

% 
function SetUIControlData(hObject, strlist)

ch = get(hObject,'ch');
for c=1:length(ch)
    if strcmpi(get(ch(c),'Type'),'uicontrol')
        switch get(ch(c),'Style')
            
            case 'edit'
                ndx = find(strcmpi(strlist.key, get(ch(c),'tag')));
                if length(ndx)==1
                    if isnumeric(strlist.val)
                        set(ch(c),'string', sprintf('%.4f', strlist.val(ndx)));
                    else
                        set(ch(c),'string', sprintf('%s', strlist.val{ndx}));
                    end
                end
                pause(0.005)
                
            case {'checkbox','popupmenu'}
                ndx = find(strcmpi(strlist.key, get(ch(c),'tag')));
                if length(ndx)==1
                    if isnumeric(strlist.val)
                        set(ch(c),'value', strlist.val(ndx));
                    else
                        set(ch(c),'value', str2num(strlist.val{ndx}));
                    end
                end
                pause(0.005);
                
            case {'slider'}
                ndx = find(strcmpi(strlist.key, get(ch(c),'tag')));
                if length(ndx)==1
                    if isnumeric(strlist.val)
                        set(ch(c),'value', strlist.val(ndx));
                    else
                        set(ch(c),'value', double(strlist.val{ndx}));
                    end
                end
                ch(c).Callback(ch(c),[])
                pause(0.005);
        end
    end
end             


% --- Executes on button press in checkboxMeta.
function checkboxMeta_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxMeta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxMeta
