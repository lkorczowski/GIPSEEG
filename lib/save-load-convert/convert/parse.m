function [OUT_EEG, OUT_STM] = parse( EEG, STM )
    % OUT_EEG = chans x samples
    % OUT_STM = numEvents x 2
    %       col 1, index to EEG data
    %       col 2, 1 TA / 0 NT

    STM_TA  = find(STM.data(:,2)==33285);
    STM_TA_SAMP = zeros(length(STM_TA), 1);
    for t = 1:length(STM_TA_SAMP)
        [ ~, data_idx ] = min(abs(EEG.data(:,1)-STM.data(STM_TA(t),1)));
        STM_TA_SAMP(t)  = data_idx;
    end
    
    STM_NT  = find(STM.data(:,2)==33286);
    STM_NT_SAMP = zeros(length(STM_NT), 1);
    for t = 1:length(STM_NT_SAMP)
        [ ~, data_idx ] = min(abs(EEG.data(:,1)-STM.data(STM_NT(t),1)));
        STM_NT_SAMP(t)  = data_idx;
    end
    
    OUT_EEG = EEG.data(:, 2:17)';
    OUT_STM = unique([STM_TA_SAMP; STM_NT_SAMP]);
    OUT_STM = [OUT_STM zeros(length(OUT_STM),1)];
    for t = 1:length(STM_TA_SAMP)
       idx = find(STM_TA_SAMP(t) == OUT_STM(:,1));
       OUT_STM(idx,2) = 1;
    end
end

