function batchReviewCallback(figHandle)
    % Retrieve data from original GUI
    data = guidata(figHandle);

    % Get the cleaned EEG
    EEG = data.batch.currentEEG;
    if data.batch.origtrials == 1
        EEG.data = EEG.data(:,:);
        EEG.pnts = EEG.trials * EEG.pnts;
        EEG.xmax = (EEG.pnts-1)/EEG.srate;
        EEG.trials = 1;
        EEG = eeg_checkset(EEG);
    end

    % Create output filename
    originalFile = data.filenames{data.batch.index};
    [~, name, ext] = fileparts(originalFile);
    outputFile = ['Cln' name ext];

    % Save cleaned dataset
    pop_saveset(EEG, 'filename', outputFile, 'filepath', data.pathname, 'savemode', 'onefile');
    disp(['Saved cleaned file: ', outputFile]);

    % Move to next file
    data.batch.index = data.batch.index + 1;
    guidata(data.batch.hObject, data);

    % Continue batch
    startReviewForNextFile(data.batch.hObject);
end
