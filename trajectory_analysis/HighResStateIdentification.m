%% This subroutine identifies the states between change points
% [CellIndx Time Xc FN Xf Xb Xnf Xnb stateChange vrf vrb vcr]
clear all
% close all
clc
warning off
currentfolder=pwd;

% % %Experiments
DataName='1_ctrl_30s';TimeRes=30;
% DataName='2_lat_30s';TimeRes=30;
% DataName='3_blebb_30s';TimeRes=30;
% % %Simulation
DataName='1_ctrl_30s_Sim';TimeRes=20;
% DataName='2_lat_30s_Sim';TimeRes=20;
% DataName='3_blebb_30s_Sim';TimeRes=20;


cd('StandardDataSets')
load([DataName '.mat']);
cd(currentfolder)

CellIndx=Data(:,1);
Time=Data(:,2);
XC=Data(:,3);
FNc=Data(:,4);
XF=Data(:,5);
XB=Data(:,6);
AllCPs=Data(:,9);
L=XF-XB;
Lf=XF-XC; Lb=XC-XB;
Cells=unique(CellIndx);
Lnuc=Data(:,7)-Data(:,8);

Fs=3600/TimeRes;
trshExcitationDuration=4; sm3=3; sm2=3;
ChangeofPhase=[0; find(AllCPs); numel(Time)];%all detected Phases even S to O
minEpisode=0;minEpisode=20;
osclTR=5; MovTr=0.002;
AvgvelWindow=240;
ChangeofCell=[find(diff(CellIndx))];%all detected Phases
ChanngePoints=union(ChangeofPhase,ChangeofCell);
iphase=1;WLs=[];Vbefors=[];
VbeforBackmoves=[];  AlldurationOfBackmoves=[];
for  Iphase=1:numel(ChanngePoints)-1
    Thisphase=[ChanngePoints(Iphase)+1:ChanngePoints(Iphase+1)];
    if numel(Thisphase)>minEpisode
        Phases(iphase,:)=[ChanngePoints(Iphase)+1 , ChanngePoints(Iphase+1)];
        thisL=L(Thisphase);thisXF=XF(Thisphase);thisXB=XB(Thisphase);thisXC=XC(Thisphase);thisFN=FNc(Thisphase);
        time=(Time(Thisphase)-Time(Thisphase(1))); thisLf=Lf(Thisphase);thisLb=Lb(Thisphase);
        avgL=smooth(thisL,300);
        highFiltL=smooth(thisL,20);
        BandPassFiltL=highFiltL-avgL;
        avgX=smooth(thisXC,300);
        highFiltX=smooth(thisXC,20);
        BandPassFiltX=highFiltX-avgX;  
        %remove a part at begining and end because of smoothing error
        OCLindx(iphase)=mean(abs((BandPassFiltL(1+round(minEpisode/2):end-round(minEpisode/2)))))+1*mean(abs((BandPassFiltX(1+round(minEpisode/2):end-round(minEpisode/2)))));
        %%% based on bandpass filter
        %         y=bpfilt(smooth(thisL,5),0.02,0.1,2,1);%signal, freq1 (1/min), freq2 (1/min), sampling frequency (1/min), isplot
        % OCLindx(iphase)=mean(abs(gradient(y)));
        % OCLindx(iphase)=mean(abs((y)));
        % OCLindx(iphase)=std((gradient(y)));  
        velindx(iphase)=(thisXC(end)-thisXC(1))/time(end);
        FNphase(iphase)=nanmean(thisFN);
        TotalTime=time(end)/3600;
        PhaseLength(iphase)=TotalTime;
        %% Duration of back excitations
        thisXCsm=smooth(thisXC,sm2);
        thisXFsm=smooth(thisXF,sm2);
        thisXBsm=smooth(thisXB,sm2);
        MovDir=sign(velindx(iphase));%for all phases
        if abs(velindx(iphase))<0.002
            MovDir=0;
        end
        thisLsm=smooth(thisXF-thisXB,20);  thisfLsm=thisXF-thisXC;  thisbLsm=thisXC-thisXB;
        deltaLf=[];
        %only real moving phases
        if MovDir>0
            durationOfBackmoves=[];
            tmpVb=(smooth(gradient(thisXBsm),sm3)-eps*ones(size(thisXFsm)));
            tmpVbsign=sign(tmpVb);
            tmpBchange=find(diff(tmpVbsign))+1;
            durationOfBackmoves=diff([1; tmpBchange;numel(thisXBsm)]);
            durationOfBackmoves(tmpVbsign([1; tmpBchange])>0)=[]; %remove forward movingback episodes
            ShortEpisode=find(durationOfBackmoves<trshExcitationDuration);
            durationOfBackmoves(ShortEpisode)=[];%remove short episodes
            %saving individual back excitation data(duration and velocity before)
            IndsVavdg=[max(1,tmpBchange-AvgvelWindow) tmpBchange-1]; %indices of a window before episodes between sign changes
            vbeforBackmoves=(thisXFsm(IndsVavdg(:,2))-thisXFsm(IndsVavdg(:,1)))./(time(IndsVavdg(:,2))-time(IndsVavdg(:,1)));
            vbeforBackmoves=[nan;vbeforBackmoves];
            vbeforBackmoves(tmpVbsign([1; tmpBchange])>0)=[];
            vbeforBackmoves(ShortEpisode)=[];%remove short episodes
            VbeforBackmoves=[VbeforBackmoves;vbeforBackmoves];
            AlldurationOfBackmoves=[AlldurationOfBackmoves;durationOfBackmoves];
            
            %%%front
            tmpVf=(smooth(gradient(thisXFsm),sm3)+eps*ones(size(thisXFsm)));
            tmpVfsign=sign(tmpVf);
            tmpFchange=find(diff(tmpVfsign))+1;
            if isempty(tmpFchange)~=1
                if tmpVfsign(tmpFchange(1))<0
                    tmpFchangetmp=tmpFchange;
                    if floor(numel(tmpFchangetmp)/2)~=(numel(tmpFchangetmp)/2);% if the number of change points is odd
                        tmpFchangetmp(end)=[];
                    end
                    Inds=reshape(tmpFchangetmp,[2,numel(tmpFchangetmp)/2]); %indices of start offront collapse(row1) and end of front collapse (row2)
                elseif tmpVfsign(tmpFchange(1))>0
                    tmpFchangetmp=tmpBchange(2:end);
                    if floor(numel(tmpFchangetmp)/2)~=(numel(tmpFchangetmp)/2);% if the number of change points is odd
                        tmpFchangetmp(end)=[];
                    end
                    Inds=reshape(tmpFchangetmp,[2,numel(tmpFchangetmp)/2]);
                end
                deltaLf=thisfLsm(Inds(1,:))-thisfLsm(Inds(2,:));
            end
            % front persistence time   
        elseif MovDir<0
            durationOfBackmoves=[];
            tmpVb=(smooth(gradient(thisXFsm),sm3)+eps*ones(size(thisXFsm)));
            tmpVbsign=sign(tmpVb);
            tmpBchange=find(diff(tmpVbsign))+1;
            durationOfBackmoves=diff([1; tmpBchange;numel(thisXBsm)]);
            durationOfBackmoves(tmpVbsign([1; tmpBchange])<0)=[]; %remove forward movingback episodes
            ShortEpisode=find(durationOfBackmoves<trshExcitationDuration);
            durationOfBackmoves(ShortEpisode)=[];%remove short episodes
            %saving individual back excitation data(duration and velocity before)
            IndsVavdg=[max(1,tmpBchange-AvgvelWindow) tmpBchange-1]; %indices of a window before episodes between sign changes
            vbeforBackmoves=(thisXBsm(IndsVavdg(:,2))-thisXBsm(IndsVavdg(:,1)))./(time(IndsVavdg(:,2))-time(IndsVavdg(:,1)));
            vbeforBackmoves=[nan;vbeforBackmoves];
            vbeforBackmoves(tmpVbsign([1; tmpBchange])<0)=[];
            vbeforBackmoves(ShortEpisode)=[];%remove short episodes
            VbeforBackmoves=[VbeforBackmoves;vbeforBackmoves];
            AlldurationOfBackmoves=[AlldurationOfBackmoves;durationOfBackmoves];  
            %%%front
            tmpVf=(smooth(gradient(thisXBsm),sm3)+eps*ones(size(thisXBsm)));
            tmpVfsign=sign(tmpVf);
            tmpFchange=find(diff(tmpVfsign))+1;
            if isempty(tmpFchange)~=1
                if tmpVfsign(tmpFchange(1))>0
                    tmpFchangetmp=tmpFchange;
                    if floor(numel(tmpFchangetmp)/2)~=(numel(tmpFchangetmp)/2);% if the number of change points is odd
                        tmpFchangetmp(end)=[];
                    end
                    Inds=reshape(tmpFchangetmp,[2,numel(tmpFchangetmp)/2]); %indices of start offront collapse(row1) and end of front collapse (row2)
                elseif tmpVfsign(tmpFchange(1))<0
                    tmpFchangetmp=tmpBchange(2:end);
                    if floor(numel(tmpFchangetmp)/2)~=(numel(tmpFchangetmp)/2);% if the number of change points is odd
                        tmpFchangetmp(end)=[];
                    end
                    Inds=reshape(tmpFchangetmp,[2,numel(tmpFchangetmp)/2]);
                end
                deltaLf=thisbLsm(Inds(1,:))-thisbLsm(Inds(2,:));
            end
        else
            durationOfBackmoves=nan;
        end
        AvgBackExcitationDuration(iphase)=nanmean(durationOfBackmoves);
        Freq(iphase)=numel(durationOfBackmoves)/TotalTime;
        DeltaLf=deltaLf;DeltaLf(DeltaLf<0)=nan;
        AvgdeltaLf(iphase)=nanmean(DeltaLf);
        minpeakDist=20; minpeakPromL=7;
        [pksL,locsL,wL,pL]=findpeaks(thisLsm,'MinPeakProminence',minpeakPromL,'MinPeakDistance',minpeakDist);
        FreqOfPeaksL(iphase)=(numel(pksL))/TotalTime;
        AvgPeaksWidthL(iphase)=mean(wL)*30/3600;
        CellNum(iphase)=CellIndx(Thisphase(1));
        iphase=iphase+1;
    else
        pause
    end
end
%%%%%find states based on ocl index and velocity
MovIndx=(abs(velindx));
for  iphase=1:numel(MovIndx)
    if MovIndx(iphase)<MovTr && OCLindx(iphase)<=osclTR
        State(iphase)=1;stateName='SS';
    elseif MovIndx(iphase)<MovTr && OCLindx(iphase)>osclTR
        State(iphase)=2;stateName='SO';
    elseif MovIndx(iphase)>=MovTr && OCLindx(iphase)<=osclTR
        if velindx(iphase)>0
            State(iphase)=3;stateName='MS';
        else
            State(iphase)=-3;stateName='MS';
        end
    elseif MovIndx(iphase)>=MovTr && OCLindx(iphase)>osclTR
        if velindx(iphase)>0
            State(iphase)=4;stateName='MO';
        else
            State(iphase)=-4;stateName='MO';
        end
    end
end
cd('StandardDataSets')
save([DataName 'StateEpisodes'],'OCLindx','velindx','FNphase','PhaseLength','AvgBackExcitationDuration',...
    'Freq','FreqOfPeaksL','CellNum','State','MovIndx');
cd(currentfolder);