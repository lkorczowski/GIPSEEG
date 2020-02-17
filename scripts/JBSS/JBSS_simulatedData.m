clear all;
close all;
% number of iterations
nbIter = 100;

freqs  = [3,7,11,17,23];
ampls  = [0.5,1,1.5,2,2.5];
phases = [0,0,0,0,0;
          0,pi/6,pi/4,pi/3,pi/2];
nbChannels = 5;
conds = [1,2;2,5;5,10;10,50;50,100];
SNRs = [100,10,5,2,1,0.5,0.1];
weightsInter = [0.5,0.6,0.7,0.8,0.9];

folder = '\\filesrv4\home$\bouchafl\.windows\Bureau\simDatNew\';

for it=6:nbIter
    disp(['iteration: ' int2str(it)]);
    disp('one cond');
    results1Cond = JBSS_sim(1,freqs,ampls,phases,nbChannels,conds,SNRs,weightsInter);
    filename = [folder 'results1Cond' int2str(it) '.mat'];
    save(filename,'results1Cond','-v7.3');
    clear results1Cond;
    disp('two cond');
    results2Cond = JBSS_sim(2,freqs,ampls,phases,nbChannels,conds,SNRs,weightsInter);
    filename = [folder 'results2Cond' int2str(it) '.mat'];
    save(filename,'results2Cond','-v7.3');
    clear results2Cond;
    disp('ERP');
    resultsERP   = JBSS_sim(3,freqs,ampls,phases,nbChannels,conds,SNRs,weightsInter);
    filename = [folder 'resultsERP' int2str(it) '.mat'];
    save(filename,'resultsERP','-v7.3');
    clear resultsERP;
end