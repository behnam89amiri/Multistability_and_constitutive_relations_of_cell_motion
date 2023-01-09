%% This subroutine finds the state change points
% input is Data file (folder: StandardDataSets) with the format:
% [CellIndx Time Xc FN Xf Xb Xnf Xnb stateChange vrf vrb vcr]

clear all
clc
close all
warning off
currentfolder=pwd;

%Experiment
% DataName='DataSetsContManual2'; TimeRes=30;
% DataName='DataSetsLatManual2'; TimeRes=30;
DataName='DataSetsBlebManual2'; TimeRes=30;

%Simulation
DataName='1_ctrl_30s_Sim';TimeRes=20;
% DataName='2_lat_30s_Sim';TimeRes=20;
% DataName='3_blebb_30s_Sim';TimeRes=20;


MethodofCP=1; %Method of finding state change points; 2 is automatic and 1 is manual
cd('StandardDataSets')
load([DataName '.mat']);
cd(currentfolder)
CellIndx=Data(:,1);Time=Data(:,2);
XC=Data(:,3);XF=Data(:,5);XB=Data(:,6);
Cells=unique(CellIndx);
Lth=0.98;sm=1;Nperm=1000;filterOut=round(600/TimeRes);
for icell=1:length(Cells)
    ThisCelltmp=find(CellIndx==Cells(icell));
    if numel(ThisCelltmp)>1
        SamplePoints=[1:filterOut:numel(ThisCelltmp)];
        ThisCell=ThisCelltmp(SamplePoints);
        XCthisCell=XC(ThisCell);XFthisCell=XF(ThisCell);XBthisCell=XB(ThisCell);
        XCsmthisCell=smooth(XCthisCell,sm);
        Xcent=XCsmthisCell-XCsmthisCell(1);
        Xbent=XBthisCell-XBthisCell(1);Xfent=XFthisCell-XFthisCell(1);
        time=(Time(ThisCell)-Time(ThisCell(1)));
        vC=gradient(Xcent,time);    vB=gradient(Xbent,time); vF=gradient(Xfent,time);
        initpoint=1;endpoint=length(vC);
        CPs=[initpoint endpoint];
        maxintwithnoPC=0; % intervals that has been checked and contain no more PC
        endloop=0;repeatfind=0;
        while endloop==0
            intstart=CPs(maxintwithnoPC+1);
            intend=CPs(maxintwithnoPC+2);
            vseg=vC(intstart:intend);
            meanvseg=mean(vseg);
            S=cumsum(vseg-meanvseg);
            Sdiff0=max(S)-min(S);
            Sdiff=[];
            for i=1:Nperm
                vsegperm=vseg(randperm(length(vseg)));
                Sperm=cumsum(vsegperm-meanvseg);
                Sdiff(i)=max(Sperm)-min(Sperm);
            end
            Lconf=length(find(Sdiff<Sdiff0))/Nperm;
            if Lconf>Lth
                [Scp,indScp]=max(abs(S(2:end-1)));
                IndScp=indScp+1+intstart-1; %index of the max S in the original vC vector(indScp is the index in vseg)
                if ismember(IndScp,CPs)~=1
                    CPsnotsorted=[CPs IndScp];
                    CPs=sort(CPsnotsorted);
                end
            else
                maxintwithnoPC=maxintwithnoPC+1;
            end
            if maxintwithnoPC==length(CPs)-1 %if the index of maximum interval with no CP is equal to all current interval.
                endloop=1;
            end
        end
        if numel(CPs)>2
            shortloop=0;j=2;
            while shortloop==0
                if CPs(j)-CPs(j-1)<15*2/filterOut
                    CPs(j)=[];
                else
                    j=j+1;
                end
                if j==numel(CPs)
                    shortloop=1;
                end
            end
        end
        if MethodofCP==1 %manual CP
            tmpFig=figure; plot(XFthisCell);hold on;plot(XBthisCell);plot(XCthisCell);plot([CPs;CPs],[min(XBthisCell)*ones(size(CPs)) ;max(XFthisCell)*ones(size(CPs))],'r')
            axis tight;axf=gcf;axf.Position=[20 100 1000*time(end)/(3600*10) max([XBthisCell;XFthisCell])-min([XBthisCell;XFthisCell])];
            [CPstmp,~]=getpts;
            CPs=sort([initpoint round(CPstmp)' endpoint]);
            close
        end
        allCPs=zeros(size(vC));
        allCPs(CPs)=1;
        tmpCP=repmat(allCPs',[filterOut,1]);
        tmpCP2=tmpCP(:);
        allCPs=tmpCP2(1:numel(ThisCelltmp));
        allCPS=zeros(size(allCPs));
        allCPS(find(diff(allCPs)>0)+floor(filterOut/2))=1;
        allCPS=allCPS(1:numel(ThisCelltmp));
    else
         allCPs=0;
    end 
    AllCPs(ThisCelltmp)=allCPS;
end
Data(:,9)=AllCPs;
cd('StandardDataSets')
NewDataName=DataName;
save(NewDataName,'Data','sm','Lth')
% save([ DataName '.mat'],'Data','ExperimentName','Names','Vtresh','sm','Lth')
cd(currentfolder)
%% merge all  consecutive similar phases and % exclude very short different phases in the middle of a long phase
% if an episode is shorter than 1h merge it with previous episode
currentfolder=pwd;
NewDataName='BlebExp5'; TimeRes=30;
clearvars -except NewDataName currentfolder TimeRes
cd('StandardDataSets')
load(NewDataName)
cd(currentfolder)
CellIndx=Data(:,1);Time=Data(:,2);
XC=Data(:,3);XF=Data(:,5);XB=Data(:,6);
AllCPs=Data(:,9);
Cells=unique(CellIndx);
short=60*60/TimeRes;% 30 min is the shortest acceptable state %%%%%%%%%%%%%%
%If the episode is short merge it with the previous episode
for icell=1:length(Cells)
    ThisCelltmp=find(CellIndx==Cells(icell));
    tmpX=XC(ThisCelltmp);
    tmpT=Time(ThisCelltmp);
    ThisPhaseChanges=[0 find(AllCPs(ThisCelltmp))' numel(ThisCelltmp)];
    ThisLengthOfPhases=diff(ThisPhaseChanges);
    tmp=find(ThisLengthOfPhases(2:end)<short); %from the second episode to the last
    AllCPs(ThisCelltmp(ThisPhaseChanges(tmp+1)))=0;
end
Data(:,9)=AllCPs;
cd('StandardDataSets')
save([NewDataName '6'],'Data')
cd(currentfolder)