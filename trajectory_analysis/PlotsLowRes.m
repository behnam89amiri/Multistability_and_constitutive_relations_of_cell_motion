%% This subroutine creates plots based on the low temporal resolution data (exp/sim)
clear all
% close all
% clc

nameOfData(1).name='4_ctrl_10min';
nameOfData(2).name='1_ctrl_10min_Sim';
nameOfData(3).name='5_lat_10min';
nameOfData(4).name='2_lat_10min_Sim';
nameOfData(5).name='6_blebb_10min';
nameOfData(6).name='3_blebb_10min_Sim';
nameOfData(7).name='7_untreated_10min';

for k=1:numel(nameOfData)
    if strcmp(strrep(nameOfData(k).name,'Sim',''),nameOfData(k).name)==0
        nameOfData(k).Simulation=1;
    else
        nameOfData(k).Simulation=0;
    end
end
%% V-FN (Fig.2B) (S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% datasets
Folder='LowResDataSets';
currentFolder=pwd;
figure;hold on
cd(Folder);
Jdata=[7 2];
for j=1:1
    load([nameOfData(Jdata(j)).name 'PersistenceProps']);
    ID=DATA(:,1); T=DATA(:,2); X=DATA(:,3); FN=DATA(:,4);
    if j==1
        Bbins=[prctile(FN,[0:5:50]) prctile(FN,[50:4:90]) prctile(FN,[90:2:95]) prctile(FN,[97:1:100]) ];
    else
        Bbins=[0:5:50 55 65 80 105];FN=FN/1.1;
    end
    Bx=(Bbins(1:end-1)+Bbins(2:end))/2;
    %%% V_FN
    clear VVB VVBerr sizes
    for i=1:numel(Bbins)-1
        tmp=find(FN>=Bbins(i) & FN<Bbins(i+1));
        sizes(i)=numel(tmp);
        VVB(i)=nanmean(abs(AllV(tmp)));
        VVBerr(i)=nanstd(abs(AllV(tmp)))/sqrt(numel(tmp));
    end
    if  nameOfData(Jdata(j)).Simulation==1
        Style='-';
    else
        Style='.';
    end
    errorbar(Bx,VVB,VVBerr,Style,'linewidth',2);
    ax=gca;ax.FontSize=18;ax.LabelFontSizeMultiplier=1.1;ax.FontName='Times New Roman';
    ylabel(['v (' char(181) 'm s^{-1})'],'interpreter','tex');xlabel('B (ng cm^{-2})','interpreter','tex');box on
    ax.YLim=[0 0.01]; ax.XLim=[0 100];
end
cd(currentFolder)
figure(1);legend('Experimetns','Simulations')
%% Persistence-V (Fig.5B)- S6(bleb)
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% datasets
nameOfData(1).name='4_ctrl_10min';
nameOfData(2).name='1_ctrl_10min_Sim';
nameOfData(3).name='5_lat_10min';
nameOfData(4).name='2_lat_10min_Sim';
nameOfData(5).name='6_blebb_10min';
nameOfData(6).name='3_blebb_10min_Sim';
nameOfData(7).name='7_untreated_10min';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% datasets
Titels={'ctrl','lat','bleb'};
Folder='LowResDataSets';
currentFolder=pwd;
cd(Folder);

minEpisodeTime=11;
Colors=get(gca,'colororder');
for iii=1:2 % different conditions   
    for j=1:2 % simulation and experiment
        figure;hold on;
        load([nameOfData(2*(iii-1)+j).name 'PersistenceProps']);
        ID=DATA(:,1); T=DATA(:,2); X=DATA(:,3); FN=DATA(:,4);
        ttmp=find(EpisodePersistence<=minEpisodeTime);
        dataPers=EpisodePersistence;dataV=EpisodeV;
        dataPers(ttmp)=[];dataV(ttmp)=[];Pers=[];PersErr=[];VErr=[];
        datax=[0.0002:.001:0.009];
if iii==3
   datax=[0.0001:.001:0.011];
end
        clear Pers PersErr VErr
        for i=1:numel(datax)-1
            tmp2=find(dataV>=datax(i) & dataV<datax(i+1));
            Pers(i)=nanmean(dataPers(tmp2));
            PersErr(i)=nanstd(dataPers(tmp2))/sqrt(length(dataPers(tmp2)));
            VErr(i)=nanstd(dataV(tmp2))/sqrt(length(dataV(tmp2)));
        end
        Xdata=((datax(1:end-1)+datax(2:end))/2);
        e2=errorbar(Xdata,Pers,PersErr,'linewidth',2);
        ylabel('persistence time (min)');xlabel(['v (' char(181) 'ms^{-1})']);
        ax=gca;ax.FontSize=20;ax.LabelFontSizeMultiplier=1.1;ax.FontName='Times New Roman';
        ax.YLim=[15 180]; box on
    end
    legend('Experimetns','Simulations')
    title(Titels{iii})
end
cd(currentFolder)
