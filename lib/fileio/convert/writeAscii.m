function [ ] = writeAscii( filename, text )

    FID = fopen( filename, 'w' );
    for i = 1:length(text)-1
        fprintf(FID, '%s\r\n', text{i});
    end
    fprintf(FID, '%s', text{length(text)});
    fclose(FID);
end

