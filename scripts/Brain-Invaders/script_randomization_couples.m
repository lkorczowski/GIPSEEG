%calcul combinaisons
clear all
couples=22;

Permu={'12a','1a2','21a','2a1','a12','a21'};

Permu=repmat(Permu,[1, ceil(couples/length(Permu))]);

Permu=Permu(randperm(length(Permu)))

% output filename
outfile = 'myTextFile';
filename = [outfile,'.txt'];
% initialize/open the file
fid = fopen(filename, 'w');

% write each cell to the text file
ncols= length(Permu');
for row=1:ncols
            fprintf(fid, 'couple#%d :',row);
        fprintf(fid, '%s \n', Permu{row});
    fprintf(fid, '\n');
end

% close file when done
fclose(fid);