clear all
clc
% close all
warning off
currentfolder=pwd;

% % %Experiments
DataName='1_ctrl_30s'; TimeRes=30;
% % %Simulation
% DataName='1_ctrl_30s_Sim';TimeRes=20;
cd('StandardDataSets')
load([DataName '.mat']);
cd(currentfolder)
numberforBackProtrude=1;smSlip=4;Simulation=0;velsm=3;
sm=20;minlengthOfPhase=60;window=minlengthOfPhase+10; Locs=[]; treshhold2=0.1; sm2=3;Ylim=[-0.06 0.05]; 
if strcmp(strrep(DataName,'Sim',''),DataName)==0 
Simulation=1;numberforBackProtrude=0.8;smSlip=5;velsm=1;
Data=Data(1:300000,:);
end
cellINDX=Data(:,1);
Time=Data(:,2);
XC=Data(:,3);
FNc=Data(:,4);
XF=Data(:,5);
XB=Data(:,6);
Xfnuc=Data(:,7);
Xbnuc=Data(:,8); 
AllCPs=Data(:,9);
L=XF-XB;
Lftmp=abs(XF-Xfnuc); Lbtmp=abs(XB-Xbnuc);
Cells=unique(cellINDX);
if Simulation==1
Vrf=Data(:,10);Vrb=Data(:,11);Vcr=Data(:,12);
end
CellDir=sign(XF-XB);
Vf=zeros(size(XF));Vb=zeros(size(XF));MoveDir=zeros(size(XF));
for i=1:numel(Cells)
    tmp=find(cellINDX==Cells(i));
    Vf(tmp)=gradient(smooth(XF(tmp),velsm),Time(tmp));
    Vb(tmp)=gradient(smooth(XB(tmp),velsm),Time(tmp));
    tmp2=[0 ;find(AllCPs(tmp)) ;numel(tmp)];
    for j=1:numel(tmp2)-1
        thisPhase=tmp(tmp2(j)+1:tmp2(j+1));
        avgV=(XC(thisPhase(end))-XC(thisPhase(1)))/(Time(thisPhase(end))-Time(thisPhase(1)));
        Dir=0;
        if avgV>0.002
            Dir=1;
        elseif avgV<-0.002
             Dir=-1;
        end
    MoveDir(tmp(tmp2(j)+1:tmp2(j+1)))=Dir;
    end
end
ChangeofCell=[0; find(diff(cellINDX)); numel(Time)];
ChanngePoints=ChangeofCell;
%% STATISTICS OVER THE CELL RETURNS
PropAroundReturn1=[]; PropAroundReturn2=[];PropAroundReturn3=[];PropAroundReturn4=[]; PropAroundReturn5=[];PropAroundReturn6=[];
 AllLags=[];
for iphase=1:numel(ChanngePoints)-1
    Thisphase=[ChanngePoints(iphase)+1:ChanngePoints(iphase+1)];
    CellDirPhase=sign(CellDir(Thisphase(2)));
    %%%% RETURNS
    Vcsm= gradient(smooth(XC(Thisphase),sm),Time(Thisphase))+0.000001;
    Vfsm= gradient(smooth(XF(Thisphase),sm),Time(Thisphase))+0.000001;
    Vbsm= gradient(smooth(XB(Thisphase),sm),Time(Thisphase))+0.000001;
    returns=find(diff(sign(Vcsm)));
    tmpp=diff([0;returns;numel(Vcsm)]); tmpp2=find(tmpp<minlengthOfPhase);
    tmpp3=union(tmpp2,tmpp2-1);tmpp3(tmpp3<1)=[];tmpp3(tmpp3==numel(returns)+1)=[];
    Returns=returns;Returns(tmpp3)=[];
    %find exact return place!
     returnsb=find(diff(sign(Vbsm)));
     returnsf=find(diff(sign(Vfsm)));
    ttmp=diff(sign(Vcsm));DirChange=ttmp(Returns);
    propAroundReturn1=nan(numel(Returns),2*window+1);
    propAroundReturn2=nan(numel(Returns),2*window+1);
    propAroundReturn3=nan(numel(Returns),2*window+1);
    propAroundReturn4=nan(numel(Returns),2*window+1);
    propAroundReturn5=nan(numel(Returns),2*window+1);
    propAroundReturn6=nan(numel(Returns),2*window+1);
    for i=1:numel(Returns)
        if DirChange(i)>0 & CellDirPhase<0
            frontV=Vf;backV=Vb; frontL=Lftmp;backL=Lbtmp;
            tmmp=Returns(i)-returnsb;tmmp(tmmp<0)=[];
            Lag=min(tmmp); 
        elseif DirChange(i)<0 & CellDirPhase<0
            frontV=Vb;backV=Vf; frontL=Lbtmp;backL=Lftmp;
            tmmp=Returns(i)-returnsf;tmmp(tmmp<0)=[];
            Lag=min(tmmp);  
        elseif DirChange(i)>0 & CellDirPhase>0
            frontV=Vb;backV=Vf; frontL=Lbtmp;backL=Lftmp;
            tmmp=Returns(i)-returnsf;tmmp(tmmp<0)=[];
            Lag=min(tmmp);  
            if Simulation==1
              frontVr=Vrb;backVr=Vrf;  
            end
        elseif DirChange(i)<0 & CellDirPhase>0
            frontV=Vf;backV=Vb; frontL=Lftmp;backL=Lbtmp;
            tmmp=Returns(i)-returnsb;tmmp(tmmp<0)=[];
            Lag=min(tmmp);
            if Simulation==1
              frontVr=smooth(Vrf,3);backVr=smooth(Vrb,3);  
            end
        end
        %front velocity
        whichProp1=-DirChange(i)*frontV; % so that the velocity before return should be positive
        tmpAft=whichProp1(Thisphase(Returns(i):min(numel(Thisphase),Returns(i)+window)));
        tmpAft2=nan(1,window+1);tmpAft2(1,1:numel(tmpAft))=tmpAft;
        tmpBef=whichProp1(Thisphase(max(1,Returns(i)-window):max(1,Returns(i)-1)));
        tmpBef2=nan(1,window);tmpBef2(1,numel(tmpBef2)-numel(tmpBef)+1:numel(tmpBef2))=tmpBef;
        propAroundReturn1(i,:)=[tmpBef2,tmpAft2];
        [pks,locs]=findpeaks(-propAroundReturn1(i,:),'MinPeakHeight',treshhold2); % - sign is to find the slipp retractions of front
        Locs=[Locs,locs];
        %rear velocity
        whichProp2=-DirChange(i)*backV;
        SlipDir=1*CellDirPhase; %+1 for back and -1 for front
        tmpAft=whichProp2(Thisphase(Returns(i):min(numel(Thisphase),Returns(i)+window)));
        tmpAft2=nan(1,window+1);tmpAft2(1,1:numel(tmpAft))=tmpAft;
        tmpBef=whichProp2(Thisphase(max(1,Returns(i)-window):max(1,Returns(i)-1)));
        tmpBef2=nan(1,window);tmpBef2(1,numel(tmpBef2)-numel(tmpBef)+1:numel(tmpBef2))=tmpBef;
        propAroundReturn2(i,:)=[tmpBef2,tmpAft2];
        %front length
        whichProp3=frontL;
        tmpAft=whichProp3(Thisphase(Returns(i):min(numel(Thisphase),Returns(i)+window)));
        tmpAft2=nan(1,window+1);tmpAft2(1,1:numel(tmpAft))=tmpAft;
        tmpBef=whichProp3(Thisphase(max(1,Returns(i)-window):max(1,Returns(i)-1)));
        tmpBef2=nan(1,window);tmpBef2(1,numel(tmpBef2)-numel(tmpBef)+1:numel(tmpBef2))=tmpBef;
        propAroundReturn3(i,:)=[tmpBef2,tmpAft2];
        %rear length
        whichProp4=backL;
        tmpAft=whichProp4(Thisphase(Returns(i):min(numel(Thisphase),Returns(i)+window)));
        tmpAft2=nan(1,window+1);tmpAft2(1,1:numel(tmpAft))=tmpAft;
        tmpBef=whichProp4(Thisphase(max(1,Returns(i)-window):max(1,Returns(i)-1)));
        tmpBef2=nan(1,window);tmpBef2(1,numel(tmpBef2)-numel(tmpBef)+1:numel(tmpBef2))=tmpBef;
        propAroundReturn4(i,:)=[tmpBef2,tmpAft2];
    
        if Simulation==1
        whichProp5=frontVr;
        tmpAft=whichProp5(Thisphase(Returns(i):min(numel(Thisphase),Returns(i)+window)));
        tmpAft2=nan(1,window+1);tmpAft2(1,1:numel(tmpAft))=tmpAft;
        tmpBef=whichProp5(Thisphase(max(1,Returns(i)-window):max(1,Returns(i)-1)));
        tmpBef2=nan(1,window);tmpBef2(1,numel(tmpBef2)-numel(tmpBef)+1:numel(tmpBef2))=tmpBef;
        propAroundReturn5(i,:)=[tmpBef2,tmpAft2];
        %rear length
        whichProp6=backVr;
        tmpAft=whichProp6(Thisphase(Returns(i):min(numel(Thisphase),Returns(i)+window)));
        tmpAft2=nan(1,window+1);tmpAft2(1,1:numel(tmpAft))=tmpAft;
        tmpBef=whichProp6(Thisphase(max(1,Returns(i)-window):max(1,Returns(i)-1)));
        tmpBef2=nan(1,window);tmpBef2(1,numel(tmpBef2)-numel(tmpBef)+1:numel(tmpBef2))=tmpBef;
        propAroundReturn6(i,:)=[tmpBef2,tmpAft2];
        end       
        AllLags=[AllLags;Lag];       
    end
    PropAroundReturn1=[PropAroundReturn1;propAroundReturn1];
    PropAroundReturn2=[PropAroundReturn2;propAroundReturn2];
    PropAroundReturn3=[PropAroundReturn3;propAroundReturn3];
    PropAroundReturn4=[PropAroundReturn4;propAroundReturn4];
            if Simulation==1
       PropAroundReturn5=[PropAroundReturn5;propAroundReturn5];
     PropAroundReturn6=[PropAroundReturn6;propAroundReturn6];
            end
end
figure;hold on;TimeRes=0.5;Tshift=2.5;
for ij=1:size(PropAroundReturn2,1)
SMPropAroundReturn2(ij,:)=smooth(PropAroundReturn2(ij,:),sm2);
end
for ij=1:size(PropAroundReturn1,1)
SMPropAroundReturn1(ij,:)=smooth(PropAroundReturn1(ij,:),sm2);
end
for i=1:size(PropAroundReturn1,1)
    plot([-window:window]*TimeRes+Tshift,SMPropAroundReturn1(i,:),'LineWidth',0.5,'color',[0.8 0.8 0.95])
    if i<=size(PropAroundReturn2,1)
        plot([-window:window]*TimeRes+Tshift,SMPropAroundReturn2(i,:),'LineWidth',0.5,'color',[0.95 0.8 0.8])
    end
end
pb=plot([-window:window]*TimeRes+Tshift,nanmean(SMPropAroundReturn2),'color',[0.85 0.325 0.1],'Linewidth',3.5);
ax=gca;ax.FontSize=17;
ax.LabelFontSizeMultiplier=1.4;ax.FontName='Times New Roman';
xlabel('t - t_{rev} (min)');ylabel(['v (' char(181) 'm s^{-1})']);box on;
pf=plot([-window:window]*TimeRes+Tshift,nanmean(SMPropAroundReturn1),'color',[0 0.45 0.74],'Linewidth',3.5);
ax.YLim=Ylim;ax.XLim=[-18 18];ax.YTick=[-0.05 0 0.05];
plot([0 0],ax.YLim,'k:','linewidth',1.2);plot(ax.XLim,[0 0],'k:','linewidth',1.2);
legend([pf,pb],{'v_f','v_b'});

edgeHist=25;
NumofPeaks=TimeRes*(Locs-window-1);
NumofPeaks(NumofPeaks>edgeHist)=[];
NumofPeaks(NumofPeaks<-edgeHist)=[];
figure(10);hold on
% hh=histogram(NumofPeaks,2*edgeHist);
[counts,edges]=histcounts(NumofPeaks,2*edgeHist);
counts=counts/size(PropAroundReturn1,1);
counts=smooth(counts,smSlip);
plot((edges(1:end-1)+edges(2:end))/2,counts,'color',[0.4 0.4 0.4],'Linewidth',3.5)
ax=gca;ax.FontSize=17;
ax.LabelFontSizeMultiplier=1.25;ax.FontName='Times New Roman';
ax.XLim=TimeRes*[-window window];
ylabel(['number of front slippage (min^{-1})']); box on;
xlabel('t - t_{rev} (min)')
ax.XLim=[-18 18];
axf=gcf;axf.Position=[518   380   517   421];
ax.XAxis.TickLength=[0 0];
figure(11);hold on
plot([-window:window]*TimeRes,smooth(nanmean(-(sign(PropAroundReturn2)-numberforBackProtrude)/2),4),'color',[0.4 0.4 0.4],'Linewidth',3.5)
ax=gca;ax.FontSize=17;
ax.LabelFontSizeMultiplier=1.3;ax.FontName='Times New Roman';
xlabel('t - t_{rev} (min)');box on;ylabel('probability of rear protruding');
ax.XLim=[-18 18];
axf=gcf;axf.Position=[518   380   517   421];
figure;hold on
TimeRes=0.5;
for i=1:size(PropAroundReturn3,1)
    plot([-window:window]*TimeRes,PropAroundReturn3(i,:),'LineWidth',0.5,'color',[0.8 0.8 0.95])
end
for i=1:10:size(PropAroundReturn4,1)
    plot([-window:window]*TimeRes,PropAroundReturn4(i,:),'LineWidth',0.5,'color',[0.95 0.8 0.8])
end
plf=plot([-window:window]*TimeRes,nanmean(PropAroundReturn3),'color',[0 0.45 0.74],'Linewidth',3.5);
plb=plot([-window:window]*TimeRes,nanmean(PropAroundReturn4),'color',[0.85 0.325 0.1],'Linewidth',3.5);
ax=gca;ax.FontSize=17;
ax.LabelFontSizeMultiplier=1.4;ax.FontName='Times New Roman';
xlabel('t - t_{rev} (min)');ylabel(['L (' char(181) 'm)']);box on;
legend([plf,plb],{'L_f','L_b'});

if Simulation==1
figure;hold on
TimeRes=0.5;
for i=1:size(PropAroundReturn5,1)
    plot([-window:window]*TimeRes,PropAroundReturn5(i,:),'LineWidth',0.5,'color',[0.8 0.8 0.95])
end
for i=1:10:size(PropAroundReturn6,1)
    plot([-window:window]*TimeRes,PropAroundReturn6(i,:),'LineWidth',0.5,'color',[0.95 0.8 0.8])
end
plf=plot([-window:window]*TimeRes+Tshift,nanmean(PropAroundReturn5),'color',[0 0.45 0.74],'Linewidth',3.5);
plb=plot([-window:window]*TimeRes+Tshift,nanmean(PropAroundReturn6),'color',[0.85 0.325 0.1],'Linewidth',3.5);
ax=gca;ax.FontSize=17;
ax.LabelFontSizeMultiplier=1.4;ax.FontName='Times New Roman';
xlabel('t - t_{rev} (min)');ylabel(['v_r (' char(181) 'm s^{-1})']);box on;
legend([plf,plb],{'v_{r,f}','v_{r,b}'});
ax.XLim=[-18 18]; ax.YLim=[0 0.08]; 
end