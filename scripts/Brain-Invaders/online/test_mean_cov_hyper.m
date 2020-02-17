clear all
REP='D:\gipsa-lab-extensions\scenarios\share\MDM\P300\brain-invaders-mdm-adaptive\config\';
C1csv='mdm-generic-C1';
C0csv='mdm-generic-C0';
REP4p='D:\gipsa-lab-extensions\scenarios\share\MDM\P300\brain-invaders-mdm-adaptive-4players\config\';

C1=csvread([REP C1csv '.csv']);
C0=csvread([REP C0csv '.csv']);

distance_riemann(C0,C1)