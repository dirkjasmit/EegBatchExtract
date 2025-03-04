function [COH, fs] = fastcoherence(data,varargin)

% function COH = fastcoherence(Data,Rate,Window,Overlap)
%
% Calculates the coherence using cross spectral density and FFT
% input:
% - Data:   matrix of signals in columns
%
% NOTE: call fastcoherence(data1,data2,...) for cross coherence between two
% datasets
%
% varargin:
% - srate/rate:  sampling rate (default 1)
% - window:  scalar of fft window size or the taper (e.g. hanning(256))
%            (default: size of first dimension of data)
% - verbose: (true/false) track progress (default: false)
% - nyqvist: (true/false) return values up to nyqvist (default: false,
%            return the full spectrum)
%
% output:
% - COH: normalized cross-coherenceP:
% - fs: frequencies

if isnumeric(varargin{1})
    % second dataset
    data2 = varargin{1};
    varlist = varargin(2:end);
else
    data2 = data;
    varlist = varargin(:);
end

if size(data,3)==1 && size(data,2)>size(data,1)*2
    warning('Data are continuous signals in rows. Transposing...')
    data = data';
end
if size(data2,3)==1 && size(data2,2)>size(data2,1)*2
    warning('Data are continuous signals in rows. Transposing...')
    data2 = data2';
end


% read parameters
srate = 1;
window = size(data,1);
firstcall = true;
verbose = false;
nyqvist = false;
for v=1:2:length(varlist)
    if strcmpi(varlist{v},'srate') ||strcmpi(varlist{v},'rate')
        srate = varlist{v+1};
    elseif strcmpi(varlist{v},'window')
        window = varlist{v+1};
    elseif strcmpi(varlist{v},'firstcall')
        firstcall = varlist{v+1};
    elseif strcmpi(varlist{v},'verbose')
        verbose = varlist{v+1};
    elseif strcmpi(varlist{v},'nyqvist')
        nyqvist = varlist{v+1};
    end
end

% check parameters
if isscalar(window)
    window = hanning(window);
end

len = length(window);



% If data are continuous, they will be organized in columns. Chop into
% epochs of 'Len' size. reshape data into epochs of size Len. First
% transpose data (signals in columns, chop up, permute back to (samples x
% channels x epochs)
if size(data,3)==1
    nepochs = floor(size(data,1)./len);
    nchans = size(data,2);
    nchans2 = size(data2,2);
    data = permute(reshape(data(1:len*nepochs,:)',nchans,len,[]),[2,1,3]);
    data2 = permute(reshape(data2(1:len*nepochs,:)',nchans2,len,[]),[2,1,3]);
end

nchans = size(data,2);
nchans2 = size(data2,2);

% you can easily determine the frequencies from the parameters...
fs = linspace(0,srate./2,len/2+1);

% check if the data are too big. If so, call this function back with a
% subset of the data. Critical number is 40, but only check for the FIRST
% call.
crit = 40;
if firstcall && ~isnumeric(varargin{1}) && nchans>crit
    % too many channels in the first call: recusrsive call!
    COH = nan(nchans,nchans,length(fs));
    for c1=1:40:size(data,2)
        for c2=1:40:size(data,2)
            tmpvar = [varargin  {'firstcall' false}];
            c1max = ifthen(c1+39>nchans,nchans,c1+39);
            c2max = ifthen(c2+39>nchans2,nchans2,c2+39);
            tmp = fastcoherence(data(:,c1:c1max,:),data(:,c2:c2max,:),tmpvar{:});
            try
                COH(c1:c1max,c2:c2max,:) = tmp(:,:,1:size(COH,3));
            catch E
                error(E.msg)
            end
        end
    end
    return
end

% start workinmg on data and data2 coherence.         
nepochs = size(data,3);
nepochs2 = size(data2,3);
if nepochs~=nepochs2
    error('data1 and data2 are not the same length')
end

% start looping thru windows and get fft's and cross
if verbose
    textprogressbar('progress: ')
    try
        textprogressbar(0);
    catch
        textprogressbar('progress: ')
    end
end
nchans = size(data,2);
nchans2 = size(data2,2);
psd1 = nan(len,nchans*nchans2,nepochs);
psd2 = nan(len,nchans2*nchans,nepochs);
cross = nan(len,nchans*nchans2,nepochs);
for e=1:size(data,3)
    if verbose
        textprogressbar(100*e/size(data,3))
    end
    F1 = fft(data(:,:,e).*repmat(window,1,nchans));
    F2 = fft(data2(:,:,e).*repmat(window,1,nchans2));
    psd1(:,:,e) = repmat(F1.*conj(F1),1,nchans2)./len;
    psd2(:,:,e) = repelem(F2.*conj(F2),1,nchans)./len;
    cross(:,:,e) = repmat(F1,1,nchans2).*repelem(conj(F2),1,nchans) ./ len;
end
if verbose
    textprogressbar('done')
end

% calculate coherence
tmp = abs((sum(cross,3).^2)./(sum(psd1,3).*sum(psd2,3)));
COH = permute(reshape(tmp,[],nchans,nchans2),[2,3,1]);
if nyqvist
    COH = COH(:,:,1:end/2+1);
end
if size(data,2)==1 && size(data2,2)==1
    COH = squeeze(COH);
end
