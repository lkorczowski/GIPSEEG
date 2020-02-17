function [ ] = writeStim( filename, stim )
% this function generate a txt file 
% it writes in filename the content of stim (1:=TA ;else:=NA)


    FID = fopen( filename, 'wt' );
        % Manage different stimulation
    for i = 1:length(stim)-1
        if (stim(i,end) == 1)
            % TA (TARGET)
            fprintf(FID, '%d %d TA\r\n', stim(i,1),stim(i,2));
        else
            % NT (NON-TARGET)
            fprintf(FID, '%d %d NT\r\n', stim(i,1),stim(i,2));
        end
    end
    
    if (stim(end,end) == 1)
        % TA (TARGET)
        fprintf(FID, '%d %d TA\r\n', stim(end,1),stim(end,2));
    else
        % NT (NON-TARGET)
        fprintf(FID, '%d %d NT\r\n', stim(end,1),stim(end,2));
    end
    fclose(FID);
end

