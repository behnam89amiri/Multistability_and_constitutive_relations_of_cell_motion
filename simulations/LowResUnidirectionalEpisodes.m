clc
clear all
close all
% % %simulation
nameOfData='1_ctrl_10min_Sim';
% nameOfData='2_lat_10min_Sim';
% nameOfData='3_blebb_10min_Sim';
% % %experiment
% nameOfData='4_ctrl_10min';
% nameOfData='5_lat_10min';
% nameOfData='6_blebb_10min';
% nameOfData='7_untreated_10min';

currFolder=pwd;
Folder='LowResDataSets';
cd(Folder)
load(nameOfData)
cd(currFolder)
ID=DATA(:,1); T=DATA(:,2); X=DATA(:,3); FN=DATA(:,4);
cells=unique(ID);
sm=1; %smoothing window for velocity calculation
minDelXEpisode=0;minDelTEpisode=9; %episode of movement that are shorter than this period or displacement will be ignored
% change nan elements of FN vector of each cell to the mean value of FN
for ii=1:numel(cells)
    thisCell = find(ID==cells(ii));
    FN2(thisCell)=nanmean(FN(thisCell))*ones(numel(thisCell),1) ;
end
FN=FN2;DATA(:,4)=FN2;
FN(isnan(FN)==1)=0;
AllV=zeros(size(FN));EpisodePersistence=[];EpisodeV=[];EpisodeB=[];i=1;
%% find unidirectional moving episodes
for ii=1:numel(cells)
    thisCell = find(ID==cells(ii));
    thisX=X(thisCell);thisXsm=smooth(thisX,sm);thisT=T(thisCell);
    thisV=gradient(thisXsm,thisT);
    AllV(thisCell)=thisV;
    %%% cell-averaged properties
    cellFN(i)=nanmean(FN(thisCell));
    celltrajectoryLength(i)=sum(abs(diff(thisXsm)));
    cellPersPath1(i)=abs(thisXsm(end)-thisXsm(1))/celltrajectoryLength(i);
    cellMaxDisplacement(i)=max(abs(thisXsm-thisXsm(1)));
    cellPersPath2(i)=cellMaxDisplacement(i)/celltrajectoryLength(i);
    cellInstV(i)=celltrajectoryLength(i)/(10*numel(thisT));
    cellAvgV(i)=abs(thisXsm(end)-thisXsm(1))/(thisT(end)-thisT(1));
    cellMeasuredTime(i)=10*numel(thisX);
    cellID(i)=cells(ii);
    %%% persistence time
    Vsign2=sign(thisV)';
    thisChangeDirPoints=[1  find(diff(Vsign2)~=0)+1  numel(thisT)];
    episodePersistence=[];episodeV=[];episodeB=[];
    for j=1:numel(thisChangeDirPoints)-1
        episode=thisChangeDirPoints(j):thisChangeDirPoints(j+1);
        DeltaXepisode=abs(thisXsm(episode(end))-thisXsm(episode(1)));
        DeltaTepisode=thisT(episode(end))-thisT(episode(1));
        if DeltaXepisode>minDelXEpisode &&  DeltaTepisode>minDelTEpisode
            episodePersistence=[episodePersistence; DeltaTepisode];
            episodeV=[episodeV; DeltaXepisode/DeltaTepisode];
            episodeB=[episodeB;cellFN(i)];
        end
    end
    EpisodePersistence=[EpisodePersistence;episodePersistence];EpisodeV=[EpisodeV;episodeV];EpisodeB=[EpisodeB;episodeB];
    cellPersTime(i)=mean(episodePersistence);% cell-averagedPersTime
    i=i+1;
end
AllV=AllV/60;cellAvgV=cellAvgV/60;cellInstV=cellInstV/60;EpisodeV=EpisodeV/60;
cd(Folder)
save([nameOfData 'PersistenceProps'],'DATA','AllV','cellFN','celltrajectoryLength','cellPersPath1','cellMaxDisplacement',...
    'cellPersPath2','cellInstV','cellAvgV','EpisodePersistence','EpisodeV','EpisodeB','cellPersTime',...
    'cellMeasuredTime','cellID')
