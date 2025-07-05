function startReviewForNextFile(hObject)
    data = guidata(hObject);

    if data.batch.index > data.batch.total
        msgbox('Batch review complete!', 'Done');
        system(sprintf('open %s', data.editSelectPath.String))
        return;
    end

    % Load the EEG dataset
    filename = data.filenames{data.batch.index};
    fullpath = fullfile(data.pathname, filename);
    EEG = pop_loadset('filename', filename, 'filepath', data.pathname);

    % Save EEG to guidata for rejection handling
    data.batch.origtrials = EEG.trials;
    if (EEG.trials==1)
        EEG = InsertDummyEvents(EEG, 2);
        EEG = pop_epoch(EEG, {'dummy'}, [0 2]);
    end
    data.batch.currentEEG = EEG;
    guidata(hObject, data);

    % Show EEG for manual rejection
    pop_eegplot(EEG, 1, 1, 0, [], 'title', 'Review EEG', ...
        'winlength', 8, ...
                'command', ['h = findobj(''Tag'', ''MainGUI'');' ...
                'data = guidata(h);' ...
                'data.batch.currentEEG = eeg_eegrej(data.batch.currentEEG, TMPREJ);' ...
                'guidata(h, data);' ...
                'batchReviewCallback(h);']);
end