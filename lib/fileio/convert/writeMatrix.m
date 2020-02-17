function [ ] = writeMatrix( filename, signal, format )

    if (nargin < 3)
        format = '%15.6f ';
    end
    [numChans, numSamples] = size(signal);
    charToPrint = [repmat(format,1,numSamples) '\r\n'];
    lastCharToPrint = [repmat(format,1,numSamples)];
    
    FID = fopen( filename, 'w' );
    for chan = 1:numChans-1
        fprintf(FID, charToPrint, signal(chan,:));
    end
    fprintf(FID, lastCharToPrint, signal(numChans,:));
    fclose(FID);
    
end

