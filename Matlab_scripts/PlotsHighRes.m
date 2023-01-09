%% This subroutine make plots based on the high temporal resolution data (exp/sim)
clear all
currFold=pwd;
%% Fig 3 B (histogram of states)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA
DataName(1).name='1_ctrl_30s'; TimeRes=30;
DataName(2).name='1_ctrl_30s_Sim';TimeRes=20;
DataName(3).name='2_lat_30s'; TimeRes=30;
DataName(4).name='2_lat_30s_Sim';TimeRes=20;
DataName(5).name='3_blebb_30s'; TimeRes=30;
DataName(6).name='3_blebb_30s_Sim';TimeRes=20;
cd('StandardDataSets')
for ik=1:3
    figure;hold on;combinedCondition=[];
    for jk=1:2
        load([DataName(2*(ik-1)+jk).name 'StateEpisodes'])
        Hours(1)=sum(PhaseLength(abs(State)==1))*30/3600;
        Hours(2)=sum(PhaseLength(abs(State)==2))*30/3600;
        Hours(3)=sum(PhaseLength(abs(State)==3))*30/3600;
        Hours(4)=sum(PhaseLength(abs(State)==4))*30/3600;
        Hours=Hours/sum(Hours);
        combinedCondition=[combinedCondition ; Hours];
    end
    bar(combinedCondition')
    ax=gca;ax.XTick=[1:4];ax.XTickLabel={'SS','SO','MS','MO'};
    ax.FontSize=15;ylabel(['Fraction of measured time']) ;axis tight;
    ax.LabelFontSizeMultiplier=1.3;ax.FontName='Times New Roman';
    box on;legend('Exp','Sim')
end
cd(currFold)
%% Fig 4 (transitions)
clearvars -except currFold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA
DataName(1).name='1_ctrl_30s'; TimeRes=30;
DataName(2).name='1_ctrl_30s_Sim';TimeRes=20;
DataName(3).name='2_lat_30s'; TimeRes=30;
DataName(4).name='2_lat_30s_Sim';TimeRes=20;
DataName(5).name='3_blebb_30s'; TimeRes=30;
DataName(6).name='3_blebb_30s_Sim';TimeRes=20;
cd('StandardDataSets')
for jk=1:6
    load([DataName(jk).name])
    CellIndx=Data(:,1);
    AllCPs=Data(:,9);
    Time=Data(:,2);
    
    ChangeofPhase=[0; find(AllCPs); numel(Time)];%all detected Phases even S to O
    ChangeofCell=[find(diff(CellIndx))];%all detected Phases
    ChanngePoints=union(ChangeofPhase,ChangeofCell);
    clearvars -except  DataName jk ChanngePoints ChangeofCell currFold
    load([DataName(jk).name 'StateEpisodes'])
    State2=[State(1:end-1);State(2:end)];
    
    SS2MS=find(State2(1,:)==1 & abs(State2(2,:))==3);SS2MS=setdiff(ChanngePoints(SS2MS+1),ChangeofCell);%% 6 3
    MS2SS=find(abs(State2(1,:))==3 & State2(2,:)==1);MS2SS=setdiff(ChanngePoints(MS2SS+1),ChangeofCell);%% 2
    SO2MO=find(State2(1,:)==2 & abs(State2(2,:))==4);SO2MO=setdiff(ChanngePoints(SO2MO+1),ChangeofCell);%% 28 1 20 26
    MO2SO=find(abs(State2(1,:))==4 & State2(2,:)==2);MO2SO=setdiff(ChanngePoints(MO2SO+1),ChangeofCell);%%6 20
    
    SS2MO=find(State2(1,:)==1 & abs(State2(2,:))==4);SS2MO=setdiff(ChanngePoints(SS2MO+1),ChangeofCell);%%9
    MO2SS=find(abs(State2(1,:))==4 & State2(2,:)==1);MO2SS=setdiff(ChanngePoints(MO2SS+1),ChangeofCell);%%9 2 8
    SO2MS=find(State2(1,:)==2 & abs(State2(2,:))==3);SO2MS=setdiff(ChanngePoints(SO2MS+1),ChangeofCell);%%4 3 10
    MS2SO=find(abs(State2(1,:))==3 & State2(2,:)==2);MS2SO=setdiff(ChanngePoints(MS2SO+1),ChangeofCell);%%12
    
    SS2SO=find(State2(1,:)==1 & State2(2,:)==2);SS2SO=setdiff(ChanngePoints(SS2SO+1),ChangeofCell);%13 10
    SO2SS=find(State2(1,:)==2 & State2(2,:)==1);SO2SS=setdiff(ChanngePoints(SO2SS+1),ChangeofCell);%% 17 13
    
    MO2MSsameDir=find((State2(1,:)==4 & State2(2,:)==3) | (State2(1,:)==-4 & State2(2,:)==-3));MO2MSsameDir=setdiff(ChanngePoints(MO2MSsameDir+1),ChangeofCell);%2
    MO2MSoppDir=find((State2(1,:)==4 & State2(2,:)==-3) | (State2(1,:)==-4 & State2(2,:)==3));MO2MSoppDir=setdiff(ChanngePoints(MO2MSoppDir+1),ChangeofCell);% 6 1 3
    
    MS2MOsameDir=find((State2(1,:)==-3 & State2(2,:)==-4) | (State2(1,:)==3 & State2(2,:)==4));MS2MOsameDir=setdiff(ChanngePoints(MS2MOsameDir+1),ChangeofCell);%1
    MS2MOoppDir=find((State2(1,:)==3 & State2(2,:)==-4) | (State2(1,:)==-3 & State2(2,:)==4));MS2MOoppDir=setdiff(ChanngePoints(MS2MOoppDir+1),ChangeofCell);%3 5
    
    MO2MOoppDir=find((State2(1,:)==-4 & State2(2,:)==4) | (State2(1,:)==4 & State2(2,:)==-4));MO2MOoppDir=setdiff(ChanngePoints(MO2MOoppDir+1),ChangeofCell);%50 43 7 36
    MS2MSoppDir=find((State2(1,:)==-3 & State2(2,:)==3) | (State2(1,:)==3 & State2(2,:)==-3));MS2MSoppDir=setdiff(ChanngePoints(MS2MSoppDir+1),ChangeofCell);%1 2
    
    NumofTransitions=[numel(SS2MS)  numel(SS2MO) numel(SS2SO) nan...
        numel(SO2MO) numel(SO2MS)  numel(SO2SS) nan...
        numel(MS2SS) numel(MS2MSoppDir) numel(MS2SO) numel(MS2MOsameDir) numel(MS2MOoppDir) nan...
        numel(MO2SO)   numel(MO2MOoppDir) numel(MO2SS) numel(MO2MSsameDir) numel(MO2MSoppDir)];
    
    inds=[0 find(isnan(NumofTransitions)==1) numel(NumofTransitions)+1];
    NumofTransitionsNorm=NumofTransitions;
    dists=[1.3 2 1.1 2];
    figure;hold on
    for ii=1:numel(inds)-1
        Sums=sum(NumofTransitions(inds(ii)+1:inds(ii+1)-1));
        NumofTransitionsNorm((inds(ii)+1:inds(ii+1)-1))=NumofTransitionsNorm((inds(ii)+1:inds(ii+1)-1))/Sums;
        plot([inds(ii+1),inds(ii+1)],[0 1],'r:','linewidth',2)
        text((inds(ii)+inds(ii+1))/2 -0.5,0.95,[num2str(Sums)],'FontSize',15)
    end
    
    bar(NumofTransitionsNorm)
    ax=gca;
    ax.FontSize=15;ylabel(['Probability of transition']) ;axis tight;
    ax.LabelFontSizeMultiplier=1.3;ax.FontName='Times New Roman';
    ax.XTick=(inds(2:end)+inds(1:end-1))/2;ax.XTickLabel={'SS','SO','MS','MO'};
    ax.YLim=[0 1];
    box on
end
cd(currFold)

%% Fig 5E (back exciation duration)
clearvars -except currFold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA
DataName(1).name='1_ctrl_30s'; TimeRes=30;
DataName(2).name='2_lat_30s'; TimeRes=30;
DataName(3).name='1_ctrl_30s_Sim';TimeRes=20;
DataName(4).name='2_lat_30s_Sim';TimeRes=20;
cd('StandardDataSets')
for ik=1:2
    if ik==2
        YLIM=[5 32];
    else
    YLIM=[5 20];
    end
combinedCondition=[]; 
    for jk=1:2
            figure;hold on
        load([DataName(2*(ik-1)+jk).name 'StateEpisodes'])
        Vdata=abs(velindx);
        Vbins=[ 0.002:0.0026:0.0155];
        for j=1:numel(Vbins)-1
            tmpp1=find(Vdata<Vbins(j+1) & Vdata>=Vbins(j) );
            ExDur(j)= nanmean(AvgBackExcitationDuration(tmpp1)*1/2) ;
            ExDurer(j)= nanstd(AvgBackExcitationDuration(tmpp1)*1/2)/sqrt(numel(tmpp1)) ;
        end      
        errorbar((Vbins(1:end-1)+Vbins(2:end))/2,ExDur,ExDurer,'linewidth',1.5)
        ax=gca;ax.FontSize=18;xlabel(['v (' char(181) 'm s^{-1})']),ylabel(['Back excitation duration (min)']) ;axis tight;
        ax.LabelFontSizeMultiplier=1.2;ax.FontName='Times New Roman';box on
        ax.XLim=[Vbins(1) Vbins(end)];ax.YLim=YLIM;
        
    end
end
cd(currFold)

%% Fig 5E (front resistance)
clearvars -except currFold
% EpisodePersistenceFRet  DeltaLFRet PHASEVFRet (window before) EpisodeVFRet(last episode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA
DataName(1).name='1_ctrl_30s'; TimeRes=30;
DataName(2).name='2_lat_30s'; TimeRes=30;
DataName(3).name='1_ctrl_30s_Sim';TimeRes=20;
DataName(4).name='2_lat_30s_Sim';TimeRes=20;
cd('StandardDataSets')
for ik=1:2
   combinedCondition=[];
       if ik==1
        YLIM=[5 22];
    else
    YLIM=[3 9];
    end
    for jk=1:2
         figure;hold on;
        load([DataName(2*(ik-1)+jk).name 'FrontResLength'])
        Vdata1=EpisodeVFRet;PS=DeltaLFRet; yLabel=['Front resistance length (' char(181) 'm)'];
        if ik==2
            Vbins=[0.003:0.0012:0.0113]; % SIM%
            XLim=[Vbins(1) Vbins(end)];
        elseif ik==1
            if jk==1
                Vbins=[ 0.004:0.003:  0.02 ]; % EXP%%%%
            else
                Vbins=[ 0.004:0.003:  0.017]; % EXP%%%%
            end
            XLim=[Vbins(1) 0.019];
        end
        LfinC=[];LfinCer=[];
        for j=1:numel(Vbins)-1
            tmpp1=find(Vdata1<Vbins(j+1) & Vdata1>=Vbins(j) );
            LfinC(j)= nanmean(PS(tmpp1)) ;
            LfinCer(j)= nanstd(PS(tmpp1))/sqrt(numel(tmpp1)) ;
        end
        errorbar((Vbins(1:end-1)+Vbins(2:end))/2,LfinC,LfinCer,'linewidth',1.5)
        ax=gca;ax.FontSize=17.5;xlabel(['v (' char(181) 'm s^{-1})']);ylabel(yLabel) ;axis tight;
        ax.LabelFontSizeMultiplier=1.1;ax.FontName='Times New Roman';box on
        ax.XLim=XLim;ax.YLim=YLIM;
    end
end
cd(currFold)

