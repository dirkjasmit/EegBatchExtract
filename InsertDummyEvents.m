function [EEGOut] = InsertDummyEvents(EEG, sec, event)
% function [EEGOut] = InsertDummyEvents(EEG, EpochLen, Rate, EventType);
% Inserts events at every EpochLen seconds. Uses Rate (default 250 Hz) to
% calculate the correct sample to insert events at.
% EventType is the event type field value for all events inserted

EEGOut = EEG;

if nargin<3
    event = 'dummy';
end
if nargin < 2
    EpochLen = 1;
end

% check if not already epoched
if EEG.trials>1
    error('Only for EEG structs with 2D data')
end

epochs = floor(EEG.pnts / (sec * EEG.srate));

ne=length(EEG.event);
nu=length(EEG.urevent);
disp(sprintf('Appending %d events.', epochs));
for e=1:epochs
    EEGOut.urevent(nu+e).type    = event;
    EEGOut.urevent(nu+e).latency = 1+((e-1)*(sec*EEG.srate));
    EEGOut.event(ne+e).type     = EEGOut.urevent(nu+e).type;
    EEGOut.event(ne+e).latency  = EEGOut.urevent(nu+e).latency;
    EEGOut.event(ne+e).duration = 1;
    EEGOut.event(ne+e).urevent  = nu+e;
end

EEGOut = eeg_checkset(EEGOut);

    