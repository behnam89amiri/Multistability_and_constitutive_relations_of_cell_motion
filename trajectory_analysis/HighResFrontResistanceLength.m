%% This subroutine gets the time point of front collpases and the avg velocity before that.
clear all
clc
close all
warning off
currentfolder=pwd;

% % % % % Experiments
DataName='1_ctrl_30s';TimeRes=30;
% DataName='2_lat_30s';TimeRes=30;
% DataName='3_blebb_30s';TimeRes=30;

% % % %Simulation
DataName='1_ctrl_30s_Sim';TimeRes=20;
% DataName='2_lat_30s_Sim';TimeRes=20;
% DataName='3_blebb_30s_Sim';TimeRes=20;

cd('StandardDataSets') 
load([DataName '.mat']);
load([DataName 'StateEpisodes'])
cd(currentfolder)

CellIndx=Data(:,1);
Time=Data(:,2);
XC=Data(:,3);
XF=Data(:,5);
XB=Data(:,6);
AllCPs=Data(:,9);
Cells=unique(CellIndx);
Lnuc=Data(:,7)-Data(:,8);
if strcmp(strrep(DataName,'Sim',''),DataName)==0 %if simulation
    sme=1; smc=1;smV=1;
else
    sme=10; smc=10; smV=1; 
end
minDelXEpisode=0.2; minDelTEpisode=60;%second
LdotTresh=-0;             VelWindow=120;dilationWind=1;
sm=15;sm2=25;
DeltaLFRet=[];EpisodeVFRet=[];EpisodePersistenceFRet=[];PHASEVFRet=[];EpisodeCellID=[];
iCol=1;
for icell=1:length(Cells)
    ThisCelltmp=find(CellIndx==Cells(icell));
    if numel(ThisCelltmp)>1
        ThisCell=ThisCelltmp;
        XCthisCell=smooth(XC(ThisCell),sm);XFthisCell=smooth(XF(ThisCell),sm2);XBthisCell=smooth(XB(ThisCell),sm2);
        
        MovDirthisCell=zeros(numel(ThisCell),1);
        ThisCellEpisodes=find(CellNum==Cells(icell));
        ChangePoints=[1 ;find(AllCPs(ThisCell)); numel(ThisCell)];
        for iii=1:numel(ChangePoints)-1
            if velindx(ThisCellEpisodes(iii))>0.002
        MovDirthisCell(ChangePoints(iii):ChangePoints(iii+1))=1;
            elseif velindx(ThisCellEpisodes(iii))<-0.002
        MovDirthisCell(ChangePoints(iii):ChangePoints(iii+1))=-1;
            end
        end   
        MovDirthisCellExtended=abs(MovDirthisCell);
        MovDirthisCellExtended(dilationWind:end)= MovDirthisCellExtended(dilationWind:end)+  MovDirthisCellExtended(1:end-dilationWind+1);
        MovDirthisCellExtended(1:end-dilationWind+1)= MovDirthisCellExtended(dilationWind:end)+  MovDirthisCellExtended(1:end-dilationWind+1);
        MovDirthisCellExtended(MovDirthisCellExtended~=0)=1;
        Xcent=smooth(XCthisCell-XCthisCell(1),30);
        Xbent=XBthisCell-XBthisCell(1);Xfent=XFthisCell-XFthisCell(1);
        time=(Time(ThisCell)-Time(ThisCell(1)));
        vC=gradient(Xcent,time);    vB=gradient(Xbent,time); vF=gradient(Xfent,time);
        
        LFthisCell=XFthisCell-XCthisCell;LBthisCell=XCthisCell-XBthisCell; LBdot=gradient(LBthisCell);LFdot=gradient(LFthisCell);
        LthisCellsm=smooth(XF(ThisCell),10)-smooth(XB(ThisCell),10);
        %% persistence time of edges:
        XCthisCellsm=smooth(XC(ThisCell)+0.000000*(rand(numel(time),1)-0.5),smc);
        XFthisCellsm=smooth(XF(ThisCell)+0.000000*(rand(numel(time),1)-0.5),sme);
        XBthisCellsm=smooth(XB(ThisCell)+0.000000*(rand(numel(time),1)-0.5),sme);
        LFthisCellsm=max(XFthisCellsm-XCthisCellsm,0); LBthisCellsm=max(XCthisCellsm-XBthisCellsm,0);
        VCsm=smooth(gradient(XCthisCellsm,time),smV);VFsm=smooth(gradient(XFthisCellsm,time),smV);VBsm=smooth(gradient(XBthisCellsm,time),smV);
        VCsign=sign(VCsm);VFsign=sign(VFsm);VBsign=sign(VBsm);
        thisChangeDirPointsC=[1;  find(diff(VCsign)~=0)+1;  numel(time)];
        thisChangeDirPointsF=[1;  find(diff(VFsign)~=0)+1;  numel(time)];
        thisChangeDirPointsB=[1;  find(diff(VBsign)~=0)+1;  numel(time)];  
        deltaLFRet=[];episodeVFRet=[];episodePersistenceFRet=[];PhaseVFRet=[];
        LastEpisodeVF=nan;
        for j=1:numel(thisChangeDirPointsF)-1
            episodeF=thisChangeDirPointsF(j):thisChangeDirPointsF(j+1);
            VwindowBefore=abs(XCthisCellsm(episodeF(1))-XCthisCellsm(max(1,episodeF(1)-VelWindow)))/(time(episodeF(1))-time(max(1,episodeF(1)-VelWindow)));
            DeltaXepisodeF=abs(XCthisCellsm(episodeF(end))-XCthisCellsm(episodeF(1)));%delta x of the center
            DeltaTepisodeF=time(episodeF(end))-time(episodeF(1));
            LdotepisodeF=(LFthisCell(episodeF(end))-LFthisCell(episodeF(1)))/DeltaTepisodeF ;
            DeltaXFepisodeF=abs(XFthisCellsm(episodeF(end))-XFthisCellsm(episodeF(1)));%delta x of the f edge
            DeltaLFepisodeF=abs(LFthisCellsm(episodeF(end))-LFthisCellsm(episodeF(1)));%delta L of the f edge
            RelVsignF=sign((XCthisCellsm(episodeF(end))-XCthisCellsm(episodeF(1)))*(XF(ThisCell(episodeF(end)))-XF(ThisCell(episodeF(1)))));
            if numel(find(MovDirthisCell(episodeF)==1))>numel(episodeF)*0.75 %if half of the episode is in the moving up state(f is real front)
                %LdotepisodeF+ means expanding cell,   RelVsignF+ means same direction movement with smooth center
                if DeltaXepisodeF>minDelXEpisode &&  DeltaTepisodeF>minDelTEpisode
                    if RelVsignF<0 && LdotepisodeF<LdotTresh %front collapse
                        deltaLFRet=[deltaLFRet;DeltaLFepisodeF];
                        episodePersistenceFRet=[episodePersistenceFRet;DeltaTepisodeF];
                        episodeVFRet=[episodeVFRet; LastEpisodeVF];
                        PhaseVFRet=[PhaseVFRet ; VwindowBefore];
                        EpisodeCellID=[EpisodeCellID ; Cells(icell)];
                    end
                end
            end
            LastEpisodeVF=DeltaXFepisodeF/DeltaTepisodeF;
        end
        LastEpisodeVB=nan;
        for j=1:numel(thisChangeDirPointsB)-1
            episodeB=thisChangeDirPointsB(j):thisChangeDirPointsB(j+1);
            VwindowBefore=abs(XCthisCellsm(episodeB(1))-XCthisCellsm(max(1,episodeB(1)-VelWindow)))/(time(episodeB(1))-time(max(1,episodeB(1)-VelWindow)));
            DeltaXepisodeB=abs(XCthisCellsm(episodeB(end))-XCthisCellsm(episodeB(1)));
            DeltaTepisodeB=time(episodeB(end))-time(episodeB(1));
            LdotepisodeB=(LBthisCell(episodeB(end))-LBthisCell(episodeB(1)))/DeltaTepisodeB;
            DeltaXBepisodeB=abs(XBthisCellsm(episodeB(end))-XBthisCellsm(episodeB(1)));%delta x of the b edge
            DeltaLBepisodeB=abs(LBthisCellsm(episodeB(end))-LBthisCellsm(episodeB(1)));%delta L of the b edge
            RelVsignB=sign((XCthisCellsm(episodeB(end))-XCthisCellsm(episodeB(1)))*(XB(ThisCell(episodeB(end)))-XB(ThisCell(episodeB(1)))));
            if numel(find(MovDirthisCell(episodeB)==-1))>numel(episodeB)*0.75 %if half of the episode is in the moving down state(b is real front)
                if DeltaXepisodeB>minDelXEpisode &&  DeltaTepisodeB>minDelTEpisode
                    if  RelVsignB<0 && LdotepisodeB<LdotTresh %front collapse
                        deltaLFRet=[deltaLFRet;DeltaLBepisodeB];
                        episodePersistenceFRet=[episodePersistenceFRet;DeltaTepisodeB];
                        episodeVFRet=[episodeVFRet; LastEpisodeVB];
                        PhaseVFRet=[PhaseVFRet ; VwindowBefore];
                        EpisodeCellID=[EpisodeCellID ; Cells(icell)];
                    end
                end
            end
            LastEpisodeVB=DeltaXBepisodeB/DeltaTepisodeB;
        end       
        DeltaLFRet=[DeltaLFRet;deltaLFRet];EpisodeVFRet=[EpisodeVFRet;episodeVFRet]; EpisodePersistenceFRet=[EpisodePersistenceFRet;episodePersistenceFRet];
        PHASEVFRet=[PHASEVFRet;PhaseVFRet];    
    else  
    end  
end
cd('StandardDataSets')
save([ DataName 'FrontResLength.mat'],'EpisodePersistenceFRet','DeltaLFRet','PHASEVFRet','EpisodeVFRet','EpisodeCellID')
cd(currentfolder)
