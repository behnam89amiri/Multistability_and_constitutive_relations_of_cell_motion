% convert High res simulation data to low res simulation data 
clear all
close all
clc
%%%% only simulation data sets
nameOfData='1_ctrl_30s_Sim';
% nameOfData='2_lat_30s_Sim';
% nameOfData='3_blebb_30s_Sim';

Folder='StandardDataSets';
ResultsFolder='LowResDataSets';
currentFolde=pwd;
if strcmp(strrep(nameOfData,'Sim',''),nameOfData)==0
    Simulation=1;
else
    Simulation=0;
end
cd(Folder);
load(nameOfData);
cd(currentFolde)
if Simulation==1 %find the low temporal resolution data (10min time frame) from complete simulation data sets
    TimeRes=30; % every ten minutes
    AllID=Data(:,1);AllT=Data(:,2);AllX=Data(:,3);AllFN=Data(:,4);
    AllCells=unique(AllID);
    IDlowtmp=[]; Tlowtmp=[]; Xlowtmp=[]; FNlowtmp=[];
    for ii=1:numel(AllCells)
        thisCell = find(AllID==AllCells(ii));
        thisCellLowTemp = thisCell(1:TimeRes:end);
        thisID=AllID(thisCellLowTemp);thisT=AllT(thisCellLowTemp);
        thisX=AllX(thisCellLowTemp);thisFN=AllFN(thisCellLowTemp);
        IDlowtmp=[IDlowtmp;thisID]; Tlowtmp=[Tlowtmp;thisT-thisT(1)];
        Xlowtmp=[Xlowtmp;thisX]; FNlowtmp=[FNlowtmp;thisFN];
    end
    DATA=[IDlowtmp Tlowtmp/60 Xlowtmp FNlowtmp];
end

ID=DATA(:,1); T=DATA(:,2); X=DATA(:,3); FN=DATA(:,4);
%clean the data 
deleted= union(find(ID==0),find(X==0 & FN==0));
DATA(deleted,:)=[];FN(deleted)=[];X(deleted)=[];T(deleted)=[];ID(deleted)=[];
cd(ResultsFolder)

newNameOfData=strrep(nameOfData,'30s','10min');
save(newNameOfData,'DATA')