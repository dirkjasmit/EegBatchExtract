function varargout = EegBatchExtract(varargin)
% EEGBATCHEXTRACT MATLAB code for EegBatchExtract.fig
%      EEGBATCHEXTRACT, by itself, creates a new EEGBATCHEXTRACT or raises the existing
%      singleton*.
%
%      H = EEGBATCHEXTRACT returns the handle to a new EEGBATCHEXTRACT or the handle to
%      the existing singleton*.
%
%      EEGBATCHEXTRACT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EEGBATCHEXTRACT.M with the given input arguments.
%
%      EEGBATCHEXTRACT('Property','Value',...) creates a new EEGBATCHEXTRACT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EegBatchExtract_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EegBatchExtract_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EegBatchExtract

% Last Modified by GUIDE v2.5 13-Mar-2025 17:02:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EegBatchExtract_OpeningFcn, ...
                   'gui_OutputFcn',  @EegBatchExtract_OutputFcn, ...
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


% --- Executes just before EegBatchExtract is made visible.
function EegBatchExtract_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EegBatchExtract (see VARARGIN)

% Choose default command line output for EegBatchExtract
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EegBatchExtract wait for user response (see UIRESUME)
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
    fn = data.editImputeFilename.String;
    if strlength(fn)>3
        fn = fn(1:strfind(fn,'(')-2);
        EEG = pop_loadset([data.editImputePath.String '/' fn]);
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
    
    data.editSelectPath.String = pathname;
    data.editSelectFileNum.String = sprintf("%d files",length(filenames));
end 
catch
    warning('No such file')
end

guidata(hObject, data)



% --- Outputs from this function are returned to the command line.
function varargout = EegBatchExtract_OutputFcn(hObject, eventdata, handles) 
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

data.editSelectPath.String = pathname;
data.editSelectFileNum.String = sprintf('%d files', length(filenames));

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
data.editImputePath.String = PathName;
data.editImputeFilename.String = sprintf('%s (%d channels)', FileName, 0);
pause(.1)

EEG = pop_loadset([PathName '/' FileName]);
data.chanlocs = EEG.chanlocs;
data.editImputeFilename.String = sprintf('%s (%d channels)', FileName, length(data.chanlocs));
pause(.05)

guidata(hObject, data);



function editImputePath_Callback(hObject, eventdata, handles)
% hObject    handle to editImputePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editImputePath as text
%        str2double(get(hObject,'String')) returns contents of editImputePath as a double


% --- Executes during object creation, after setting all properties.
function editImputePath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editImputePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editImputeFilename_Callback(hObject, eventdata, handles)
% hObject    handle to editImputeFilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editImputeFilename as text
%        str2double(get(hObject,'String')) returns contents of editImputeFilename as a double


% --- Executes during object creation, after setting all properties.
function editImputeFilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editImputeFilename (see GCBO)
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
freqnames = {'alpha','beta','theta','delta','alphalo','alphahi','betalo','betahi','SENSEalpha','SENSEtheta'};
freqs = {[8 13], [13 30], [4 8], [1 3], [7 9.5], [9.5 13], [13 21], [21 30], [8.5 12.0], [3.0 8.5]};


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
pars.dB = data.checkboxdB.Value;

% DFA, IAF, and Power have the same output format.
if ismember(var, {'DFA','IAF','Power'})
    RunSingleParameterExtraction(data.filenames, pars);

elseif ismember(var, {'Coherence'})
    RunCoherence(data.filenames, pars);
    
elseif ismember(var, {'PSD'})
    RunPSD(data.filenames, pars);
    
elseif ismember(var, {'FAA'})
    RunFAA(data.filenames, pars);
    
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
            if pars.dB
                AllPSD(frq,:,f) = 10*log10(mean(P(ndx,:),1));
            else
                AllPSD(frq,:,f) = mean(P(ndx,:),1);
            end
        end
    end

    disp('Data stored in global AllPSD (freqs in AllPSDfreqs)')
    
    % open a figure for the PSD
    chanlocs = pars.chanlocs;

    fig = figure('Name', 'plot Power Spectrum', 'Tag', 'figPSD', 'NumberTitle', 'off', ...
                 'Position', [100, 100, 700, 500]); % Adjust size as needed
    % Create a frame (uipanel), checkbox, and create axes inside the panel
    frame = uipanel('Parent', fig, 'Title', 'Power spectrum', 'Tag', 'framePSD', ...
                    'Position', [0.005, 0.1, .98, 0.9]); % [x, y, width, height]
    ax = axes('Parent', frame, 'Position', [0.07, 0.1, .91, 0.88], 'Tag', 'axesPSD'); % Adjust within panel
    checkbox = uicontrol('Style', 'checkbox', 'Parent', fig, 'Tag', 'cbSummary', ...
                         'String', 'Summarize into regions', 'Fontsize', 11, ...
                         'Units', 'normalized', ...
                         'Position', [0.05, 0.03, 0.4, 0.05], ...
                         'Callback', @(src, event) checkbox_callback(src, ax, AllPSDfreqs, AllPSD, chanlocs, 10, pars.dB));
                     
    checkbox_callback(checkbox, ax, AllPSDfreqs, AllPSD, chanlocs, 10, pars.dB);




% Callback function for checkbox
function checkbox_callback(hObject, ax, fs, P, chanlocs, fontsize, dB)

    % hardcoded limits
    ndx = fs>=1 & fs<45;  
    
    % plot either all channels or a summary
    if get(hObject,'Value')==0
        x = fs(ndx);
        if dB
            y = 10*log10(P(ndx,:));
        else 
            y = P(ndx,:);
        end
        newx = min(fs(ndx)):.25:max(fs(ndx));
        res = arrayfun(@(c)spline(x, y(:,c), newx)', 1:size(y,2), 'UniformOutput', false);
        newy = cell2mat(res);
        plot(ax, newx, newy);
        
    else
        numlabels = {'theta','radius','X','Y','Z','sph_theta','sph_phi','sph_radius'};
        tab = struct2table(chanlocs);
        % repair: convert cell array of double to double with missings
        for lab=1:length(numlabels)
            if ismember(numlabels{lab}, tab.Properties.VariableNames) && ~isnumeric(tab.(numlabels{lab}))
                values = cellfun(@(x)ifthen(isempty(x), nan, double(x)), tab.(numlabels{lab}));
                tab.(lab) = values;
            end
        end
        
        relX = tab.X ./ sqrt(tab.X.^2+tab.Y.^2+tab.Z.^2);
        relY = tab.Y ./ sqrt(tab.X.^2+tab.Y.^2+tab.Z.^2);
        % relZ = tab.Z ./ sqrt(tab.X.^2+tab.Y.^2+tab.Z.^2);
        
        ant    = relX>=-1E-5;
        post   = ~ant;
        medial = abs(relY)<.41;
        left   = relY>=.41;
        right  = relY<=.41;
        
        regP = nan(size(P,1),6);
        regP(:,1) = mean(P(:,ant & left),2);
        regP(:,2) = mean(P(:,ant & medial),2);
        regP(:,3) = mean(P(:,ant & right),2);
        regP(:,4) = mean(P(:,post & left),2);
        regP(:,5) = mean(P(:,post & medial),2);
        regP(:,6) = mean(P(:,post & right),2);

        x = fs(ndx);
        if dB
            y = 10*log10(regP(ndx,:));
        else 
            y = regP(ndx,:);
        end
        newx = min(fs(ndx)):.25:max(fs(ndx));
        res = arrayfun(@(c)spline(x, y(:,c), newx)', 1:size(y,2), 'UniformOutput', false);
        newy = cell2mat(res);
        plot(ax, newx, newy);
        legend('ant left','ant medial','ant right','post left','post medial','post right')
 
    end
    set(gca, 'fontsize', fontsize+2)
    xlabel('frequency (Hz)')
    ylabel('Power({\mu}V^2/Hz)')   


function RunFAA(filenames, pars)
    % function for analysis coherence

    % initialize table with no values in the table (rowcount = 0). One
    % string for Id and the rest are doubles
    % T = table('Size', [0, length(pars.chanlocs)+1], ...
    %     'VariableTypes', ["string" repelem("double",length(pars.chanlocs))], ...
    %     'VariableNames', [{'Id'} cellfun(@(x)sprintf('%s_%s',prefix,x),{pars.chanlocs.labels},'uni',0)]);

    % initialize output. Save the data into a glocal variable for cohoerence.
    global FAA
    FAA = nan(length(filenames),1);

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
        ndxF3 = ismember({pars.chanlocs.labels}, {'F3'});
        ndxF4 = ismember({pars.chanlocs.labels}, {'F3'});
        if ~sum(ndxF3 | ndxF4)
            % determine not by name but by approximate location... [sph_theta
            % sph_phi]
            ideal_F3 = [43.5 29.3];
            ideal_F4 = [-43.5 29.3];
            delta = nan(1,length(pars.chanlocs));
            for ch=1:length(pars.chanlocs)
                deltaF3(ch) = angularDifference(ideal_F3(1), ideal_F3(2), pars.chanlocs(ch).sph_theta, pars.chanlocs(ch).sph_phi);
                deltaF4(ch) = angularDifference(ideal_F4(1), ideal_F4(2), pars.chanlocs(ch).sph_theta, pars.chanlocs(ch).sph_phi);
            end
            ndxF3 = minindex(deltaF3);
            ndxF4 = minindex(deltaF4);
        end        

        [P, fs] = pfft(EEG.data([ndxF3 ndxF4],:)', EEG.srate, hanning(EEG.srate*4), .5);
        ndx = fs>=pars.freq(1) & fs<pars.freq(2);



        tmp = log(mean(P(ndx,:),1));
         % index 1 is F3, index 2 is F4. This is the correct way so that higher
         % scores mean greater right brain alpha, greter left brain activity,
         % and therefore apporach related behaviors (positive affect).
        FAA(f) = tmp(2) - tmp(1);
    end

    disp('FAA data stored in global variable FAA')


function angle_diff = angularDifference(theta1, phi1, theta2, phi2)
    % Convert degrees to radians
    theta1 = deg2rad(theta1);
    phi1 = deg2rad(phi1);
    theta2 = deg2rad(theta2);
    phi2 = deg2rad(phi2);

    % Convert to Cartesian unit vectors
    A = [cos(phi1) * cos(theta1), cos(phi1) * sin(theta1), sin(phi1)];
    B = [cos(phi2) * cos(theta2), cos(phi2) * sin(theta2), sin(phi2)];

    % Compute the angular difference using the dot product
    cos_alpha = dot(A, B); % No need to divide by magnitude since they are unit vectors
    cos_alpha = max(min(cos_alpha, 1), -1); % Ensure numerical stability

    % Compute the angle in degrees
    angle_diff = rad2deg(acos(cos_alpha));


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

    disp('Coherence data stored in global variable MeanCOH')


function RunSingleParameterExtraction(filenames, pars)
% function for analysis of EEG extracting a single parameter per
% channel

    global T

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
            prefix = sprintf('%s_%s', pars.var, pars.freqname);
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
                ndx = fs>=pars.freq(1) & fs<pars.freq(2);
                val = mean(P(ndx,:),1);
                if pars.dB
                    val = 10*log10(val);
                end
        end

        T.Id(f) = filenames{f};
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

fprintf('Data Table also stored in global variable T\n')

% plotting the data in T
toPlot = nanmean(table2array(T(:,2:end)));
figure; 
cmin = min(toPlot(:));
cmax = max(toPlot(:));
topoplot((toPlot), pars.chanlocs, 'maplimits', [cmin*.98 cmax*.98])
colorbar


% --------  end of function -----------------------




function editSelectPath_Callback(hObject, eventdata, handles)
% hObject    handle to editSelectPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSelectPath as text
%        str2double(get(hObject,'String')) returns contents of editSelectPath as a double


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



function editSelectFileNum_Callback(hObject, eventdata, handles)
% hObject    handle to editSelectFileNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSelectFileNum as text
%        str2double(get(hObject,'String')) returns contents of editSelectFileNum as a double


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
                try
                    tmp1 = sprintf('%s',get(ch(c),'tag'));
                    tmp2 = sprintf('%s',get(ch(c),'string'));
                catch
                end
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


% --- Executes on button press in checkboxdB.
function checkboxdB_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxdB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxdB
