function compt=Flash2Serie(F, Fs)
%%
if nargin<2
    Fs=128;%sample rate, default 128Hz
end
Indice=find(F);
length(Indice);
DiffInd=diff([Indice;length(F)]);

compt=zeros(size(DiffInd));
for gInd=1:length(DiffInd)
    ADD=0;
    
    % first case: if delay between 2 stimulations is higher than 10
    % secondes it is a break
    if DiffInd(gInd)>(Fs*8-1) %new break (+1000)
        ADD=1000;
    elseif DiffInd(gInd)>(Fs*3-1) %new level (+100)
        % second case: if delay between 2 stimulations is higher than 3
        % secondes it is a new level
        ADD=100;
    elseif DiffInd(gInd)>(Fs*1.1-1) %new trial (+10)
        % third case: if delay between 2 stimulations is higher than 1.25
        % secondes it is a life level
        ADD=10;
    elseif DiffInd(gInd)>0 %same serie (+1)
            ADD=1;
    end
    
        compt(gInd)=ADD;

    
end

%semilogy(compt,'*')