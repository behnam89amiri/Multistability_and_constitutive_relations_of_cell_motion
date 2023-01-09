%% Automatic line detection
%{
@author: Christoph Schreiber
2nd contact person: Johannes Heyn

Detects lines (e.g. fluorescently labelled protein lanes) and saves coordinates
and intensity in mat files and coordinates only in a csv table with columns
x1, x2, y1, y2 where each row signifies a line.

Overview of work flow
1. line_detection [you are here]
2. lines_track
3. combine_tracks_folders
4. combine_alltracks
5. clean_alltracks
6. analysis

For more info, see README.md

Protocol
- go to folder with pattern images. Images must be tif-files.
- define grid with step11
- detected lines will be displayed. Delete all lines-files where line
detection failed.

Parameters to change
brightnessvalue: defines brightness of pattern image, only for display
default filename is ...XY01
%}
clear all

thetastep = 0.2;
anglemaxdiff = 30;
%brightnessvalue %larger means darker
% brightnessvalue=10000; %TIRF, Cy5_bin2 without binning (config error)
% brightnessvalue = 600; %TIRF, Cy5_JH (100ms, 16bit, no binning)
%brightnessvalue = 24000; %TIRF, Cy5_JH (100ms, 16bit, no binning), 31.03.21 for 210215TIRF

%brightnessvalue = 180; %Veigel Ti, GFP Intensilight, 2.7.21

brightnessvalue = 240; %TIRF, Cy5_JH (100ms, 16bit, no binning), 10/12/21 for 211126TIRF
showAll = false; %show all intermediate steps

%% input GUI
useGUI = true;
if useGUI
    prompt = {'Start loop at position:',...1
        'Position identification tag:',...2
        'Binning factor:',...3
        'Occular lense magnification factor:',...4
        'Are the lines stepped?'}; 
    dlgtitle = 'Input parameters';
    dims = [1 35];
    definput = {'1',...1
        'XY',...2
        '2',...3
        '1',...4
        'false'};
    userInput = inputdlg(prompt, dlgtitle, dims, definput);

    %variable assignment
    startpoint = str2double(userInput{1});
    cXY = userInput{2};
    binning = str2double(userInput{3});
    mag = str2double(userInput{4});
    difFnConc = ~strcmp(userInput{5},'false');
else
    startpoint = 1; %change if you don't want to start from the first file
    binning = 2;
    mag = 1;
end

%% Get list of images

dirFolder = dir;                               %# Get the data for the current directory
dirIndex = [dirFolder.isdir];                  %# Find the index for directories
FileList = {dirFolder(~dirIndex).name}';       %# Get a list of the files
FileListIdx = contains(FileList, '.tif');
FileList = FileList(FileListIdx);

folderName = 'lines';
if ~isfolder(folderName)
    [status, msg, msgID] = mkdir(folderName);
    if status == 0
        error(msgID, msg);
    end
end


%% define Grid

%% for Unikon 200821
% 1-18; 55-72; 91-116
% step11=12.8;
% step12=55.2;
% step21=153;
% step22=147.5;
% 19-54; 73-90; 117-133
% step11=12.8;
% step12=55.2;
% step21=116;
% step22=111;

%% for Tirf 200824
% 1-7; 34-41;58-75
% step11=12.8;
% step12=55.2;
% step21=153;
% step22=147.5;
% 8-33;42-57; 76-90
% step11=12.8;
% step12=55.2;
% step21=116;
% step22=111;

%% for Unikon 200824
% 14:22; 37:44 ; 87:96
% step11=12.8;
% step12=55.2;
% step21=153;
% step22=147.5;
% 1:13; 23:36; 45:86
% step11=12.8;
% step12=55.2;
% step21=116;
% step22=111;

%% for Unikon 200826
% 20:54; 86:123; 
% step11=12.8;
% step12=55.2;
% step21=153;
% step22=147.5;
% 1:19; 55:85 ; 124:143
% step11=12.8;
% step12=55.2;
% step21=116;
% step22=111;

%% for Unikon 200828
%  
% step11=12.8;
% step12=55.2;
% step21=153;
% step22=147.5;

%% for Tirf 200831
%  16:66;75:104; 131:152
step11=12.8;
%step12=55.2;
step21=153;
step22=147.5;
% 1:15; 67:74; 105:130
% step11=12.8;
% step12=55.2;
% step21=116;
% step22=111;

%insert drop down menu where you can choose between microscope
micWhich = 'VeigelTi'; %which microscope was used?

switch micWhich
    case 'UNikon'
        micFactor = 0.99;
    case 'TIRF'
        micFactor = 1;
    case 'Nikon'
        micFactor = 1.01; %scaling factor; Nikon
    case 'VeigelTi'
        micFactor = 0.4065;
    otherwise
        error('No microscope specified.')
end
%% TIRF & UNikon 07/06/21
% step11 - lane width
% step12 - distance between lanes
% patternsize - distance between repeating lanes (=step11+step12)
step11=16.35;
%% for UNikon 211008
patternsize = 68.5*2/binning*mag * micFactor; %68 too wide, 67.
step11 = 13;
step11 = step11*2/binning*mag * micFactor; %13

%% for UNikon 211008
patternsize = 68.*2/binning*mag * micFactor; %68 too wide, 67.9
step11 = step11*2/binning*mag * micFactor; %13
patternsize = 68.5; %ca 68 +/-1
patternsize = patternsize*2/binning*mag * micFactor; 
step12 = patternsize - step11;


progressbar('Files')

for filenr=startpoint:size(FileList,1)
fname = char(FileList(filenr));
fnameL = lower(fname);
cXY = lower(cXY);
[startIdx, endIdx] = regexp(fnameL, cXY);
numberpattern = regexp(fnameL(endIdx+1:end),'\d+','match', 'once');

I = imread(fname);

if showAll
    figure('Name', strcat(fname, ': I'))
    imshow(I)
end
% BW=edge(I,'log');
% figure
% imshow(BW)
%I = imgaussfilt(I,1);
%%
BW = edge(imgaussfilt(I,1.5),'Sobel');
if showAll
    figure('Name', strcat(fname, ': BW'))
    imshow(BW)
end

BW2 = im2bw(I,graythresh(I)*1.4);
if showAll
    figure('Name', strcat(fname, ': BW2'))
    imshow(BW2)
end

% test whether binarisation worked (notoriously problematic for weak signal)
[tmpx, tmpy] = size(BW2);
if sum(BW2, 'all') == tmpx*tmpy
    BW2 = BW;
    if showAll
            disp('binarisation of pattern unsuccessful')
        figure('Name', strcat(fname, ': BW2 corrected'))
        imshow(BW2)
    end
end

[H,theta,rho] = hough(BW,'Theta', -90:thetastep:89.99);
[H2,theta2,rho2] = hough(BW2,'Theta', -90:thetastep:89.99);

if showAll
    figure('Name', strcat(fname, ': H'))
    imshow(imadjust(mat2gray(H)),[],...
           'XData',theta,...
           'YData',rho,...
           'InitialMagnification','fit');
    xlabel('\theta (degrees)')
    ylabel('\rho')
    axis on
    axis normal
    hold on
    colormap(hot)
end

%P = houghpeaks(H,300,'threshold',ceil(0.05*max(H(:))));
%x = theta(P(:,2));
%y = rho(P(:,1));




%% find angles via maxima of correlationfunction of stepfunction with H

for i=1:size(H,2)
     %maxcorrH1(i)=max(xcorr(H(:,i)/norm(H(:,i)),stepfilterfunction2(step11,step12,size(H,1))/norm(stepfilterfunction2(step11,step12,size(H,1)))));
%     maxcorrH2(i)=max(xcorr(H(:,i)/norm(H(:,i)),stepfilterfunction2(step21,step22,size(H,1))/norm(stepfilterfunction2(step21,step22,size(H,1)))));
% changed to use H2 also to detect narrow lines
    maxcorrH1(i)=max(xcorr(H2(:,i),stepfilterfunction2(step11,step12,size(H,1))/norm(stepfilterfunction2(step11,step12,size(H2,1)))));
%     maxcorrH2(i)=max(xcorr(H(:,i),stepfilterfunction2(step21,step22,size(H,1))/norm(stepfilterfunction2(step21,step22,size(H,1)))));
%     maxcorrH3(i)=max(xcorr(H(:,i),stepfilterfunction2(step22,step21,size(H,1))/norm(stepfilterfunction2(step22,step21,size(H,1)))));
    
    %maxcorrH21(i)=max(xcorr(H2(:,i),stepfilterfunction2(step11,step12,size(H2,1))/norm(stepfilterfunction2(step11,step12,size(H2,1)))));
    maxcorrH22(i)=max(xcorr(H2(:,i),stepfilterfunction2(step21,step22,size(H2,1))/norm(stepfilterfunction2(step21,step22,size(H2,1)))));
    maxcorrH23(i)=max(xcorr(H2(:,i),stepfilterfunction2(step22,step21,size(H2,1))/norm(stepfilterfunction2(step22,step21,size(H2,1)))));
    
    
    
end

if showAll
    figure('Name', strcat(fname, ': maxcorrH1'))
    hold on
    plot(maxcorrH1)
end

[M1,phi1]=max(maxcorrH1);

% if phi1+(-90-anglemaxdiff)/thetastep >= 0
%     [M2,phi2]=max(maxcorrH2(phi1+round((-90-anglemaxdiff)/thetastep):phi1+round((-90+anglemaxdiff)/thetastep)));
% elseif phi1+(-90+anglemaxdiff)/thetastep>=0
%     [M2,phi2]=max(maxcorrH2([0:phi1+round((-90+anglemaxdiff)/thetastep), size(H,2)+phi1+round((-90-anglemaxdiff)/thetastep)]));
% else
%     [M2,phi2]=max(maxcorrH2(size(H,2)+phi1+round((-90-anglemaxdiff)/thetastep):size(H,2)+phi1+round((-90+anglemaxdiff)/thetastep)));
% end
%% determination of phi2...
% if phi1+(-90-anglemaxdiff)/thetastep >= 0
%     maxcorrH2(1:phi1+round((-90-anglemaxdiff)/thetastep))=0;
%     maxcorrH2(phi1+round((-90+anglemaxdiff)/thetastep):end)=0;
%     maxcorrH3(1:phi1+round((-90-anglemaxdiff)/thetastep))=0;
%     maxcorrH3(phi1+round((-90+anglemaxdiff)/thetastep):end)=0;
% elseif phi1+(-90+anglemaxdiff)/thetastep>=0
%     maxcorrH2(phi1+round((-90+anglemaxdiff)/thetastep): size(H,2)+phi1+round((-90-anglemaxdiff)/thetastep))=0;
%     maxcorrH3(phi1+round((-90+anglemaxdiff)/thetastep): size(H,2)+phi1+round((-90-anglemaxdiff)/thetastep))=0;
% else
%     maxcorrH2(1:size(H,2)+phi1+round((-90-anglemaxdiff)/thetastep))=0;
%     maxcorrH2(size(H,2)+phi1+round((-90+anglemaxdiff)/thetastep):end)=0;
%     maxcorrH3(1:size(H,2)+phi1+round((-90-anglemaxdiff)/thetastep))=0;
%     maxcorrH3(size(H,2)+phi1+round((-90+anglemaxdiff)/thetastep):end)=0;
% end
% 
% 
% 
% [M2,phi2]=max(maxcorrH2);
% [M3,phi3]=max(maxcorrH3);
% if M3>M2
%     phi2=phi3;
%     step22h=step22;
%     step22=step21;
%     step21=step22h;
% end
%% alternative determination of phi2...
if phi1+(-90-anglemaxdiff)/thetastep >= 0
    maxcorrH22(1:phi1+round((-90-anglemaxdiff)/thetastep))=0;
    maxcorrH22(phi1+round((-90+anglemaxdiff)/thetastep):end)=0;
    maxcorrH23(1:phi1+round((-90-anglemaxdiff)/thetastep))=0;
    maxcorrH23(phi1+round((-90+anglemaxdiff)/thetastep):end)=0;
elseif phi1+(-90+anglemaxdiff)/thetastep>=0
    maxcorrH22(phi1+round((-90+anglemaxdiff)/thetastep): size(H,2)+phi1+round((-90-anglemaxdiff)/thetastep))=0;
    maxcorrH23(phi1+round((-90+anglemaxdiff)/thetastep): size(H,2)+phi1+round((-90-anglemaxdiff)/thetastep))=0;
else
    maxcorrH22(1:size(H,2)+phi1+round((-90-anglemaxdiff)/thetastep))=0;
    maxcorrH22(size(H,2)+phi1+round((-90+anglemaxdiff)/thetastep):end)=0;
    maxcorrH23(1:size(H,2)+phi1+round((-90-anglemaxdiff)/thetastep))=0;
    maxcorrH23(size(H,2)+phi1+round((-90+anglemaxdiff)/thetastep):end)=0;
end



[M2,phi2]=max(maxcorrH22);
[M3,phi3]=max(maxcorrH23);
if M3>M2
    phi2=phi3;
    step22h=step22;
    step22=step21;
    step21=step22h;
end









%%
clear M1 M2 M3
corr11=xcorr(H(:,phi1)/norm(H(:,phi1)),stepfilterfunction2(step11,step12,size(H,1))/norm(stepfilterfunction2(step11,step12,size(H,1))));
corr22=xcorr(H2(:,phi2)/norm(H2(:,phi2)),stepfilterfunction2(step21,step22,size(H2,1))/norm(stepfilterfunction2(step21,step22,size(H2,1))));

rho1= mod(find(corr11(size(H,1):end)==max(corr11(size(H,1):end)),1)-size(H,1)/2,step11+step12);
rho2= mod(find(corr22(size(H2,1):end)==max(corr22(size(H2,1):end)),1)-size(H2,1)/2,step21+step22);
%% correct angle
phi1b=pi*(phi1-1)*thetastep/180;
phi2b=pi*(phi2-1)*thetastep/180;





%% plot lines
% firstlines=myLinesCalculator2(step11,step12,phi1b+pi/2,rho1,size(I,2),size(I,1));
firstlines=linescalculator2(step11,step12,phi1b,rho1,size(I,2),size(I,1));
firstlines2=linescalculator2(step21,step22,phi2b,rho2,size(I,2),size(I,1));
figure('Name', fname);
imshow(double(I)/brightnessvalue,'InitialMagnification','fit')
hold on
for i=1:size(firstlines,1)
    plot([firstlines(i,1),firstlines(i,2)],[firstlines(i,3),firstlines(i,4)],'b')
end
for i=1:size(firstlines2,1)
    plot([firstlines2(i,1),firstlines2(i,2)],[firstlines2(i,3),firstlines2(i,4)],'b')
end


%% calculate intensity of lines ...

%% first group pairs of lines together by looking at the distance of the first lines
firstlines(isnan(firstlines(:,1)),:)=[];
firstlines2(isnan(firstlines2(:,1)),:)=[];
startline1=NaN;
if round(-((firstlines(2,3)-firstlines(2,1)*(firstlines(2,4)-firstlines(2,3))/(firstlines(2,2)-firstlines(2,1)))-...
    (firstlines(1,3)-firstlines(1,1)*(firstlines(1,4)-firstlines(1,3))/(firstlines(1,2)-firstlines(1,1))))*cos(phi1b))-step11<=1
    startline1=0;
elseif round( -((firstlines(3,3)-firstlines(3,1)*(firstlines(3,4)-firstlines(3,3))/(firstlines(3,2)-firstlines(3,1)))-...
        (firstlines(2,3)-firstlines(2,1)*(firstlines(2,4)-firstlines(2,3))/(firstlines(2,2)-firstlines(2,1))))*cos(phi1b))-step11<=1
    startline1=1;
else %quick fix. One of the conditions above should always be met
    fprintf('Could not determine startline1 for file %s', fname)
    continue
end


for i=1: floor((size(firstlines,1)-startline1)/2)
    lines(i).linexy=firstlines(i*2-1+startline1:i*2+startline1,1:4);
    lines(i).meanline=mean(lines(i).linexy);
    lines(i).width=[step11 step12];
end

%% for each pair calculate wich
for i=1:size(lines,2)
    ind=1;
    for j=1:size(firstlines2,1)
    crosslines1(1,1:2)=intersect_lines(lines(i).linexy(1,:), firstlines2(j,:));
    crosslines2(1,1:2)=intersect_lines(lines(i).linexy(2,:), firstlines2(j,:));

    if crosslines1(1,1)>=0 && crosslines1(1,1)<=size(I,2)&& crosslines1(1,2)>=0 && crosslines1(1,2)<=size(I,1) &&...
       crosslines2(1,1)>=0 && crosslines2(1,1)<=size(I,2)&& crosslines2(1,2)>=0 && crosslines2(1,2)<=size(I,1);
        lines(i).intersect1(ind,1:2)=crosslines1;
        lines(i).intersect2(ind,1:2)=crosslines2;
        lines(i).meanintersect(ind,:)=(crosslines2+crosslines1)/2;
        ind=ind+1;
    end
    
    end
    clear crosslines1 crosslines2 ind
end



for i=1:size(lines,2)
   for j=1:size(lines(i).intersect1,1)-1
%        figure
%        imshow(rectanglespec(lines(i).intersect1(j,:),lines(i).intersect1(j+1,:), lines(i).intersect2(j,:),lines(i).intersect2(j+1,:),size(I)));
      lines(i).intensitymedian(1,j)=median(I(rectanglespec2(lines(i).intersect1(j,:),lines(i).intersect1(j+1,:), lines(i).intersect2(j,:),lines(i).intersect2(j+1,:),size(I))));
      lines(i).intensitymean(1,j)=mean(I(rectanglespec2(lines(i).intersect1(j,:),lines(i).intersect1(j+1,:), lines(i).intersect2(j,:),lines(i).intersect2(j+1,:),size(I))));
      lines(i).intensitystd(1,j)=std(double(I(rectanglespec2(lines(i).intersect1(j,:),lines(i).intersect1(j+1,:), lines(i).intersect2(j,:),lines(i).intersect2(j+1,:),size(I)))));
      lines(i).length(1,j)= round(norm(lines(i).intersect1(j,:)-lines(i).intersect1(j+1,:)));
   end
   if i<size(lines,2)
      lines(i).backgrndmedian(1,1)=median(I(rectanglespec2( lines(i).linexy(2,[1,3]), lines(i).linexy(2,[2,4]), lines(i+1).linexy(1,[1,3]), lines(i+1).linexy(1,[2,4]),size(I) )));
      lines(i).backgrndmean(1,1)=mean(I(rectanglespec2( lines(i).linexy(2,[1,3]), lines(i).linexy(2,[2,4]), lines(i+1).linexy(1,[1,3]), lines(i+1).linexy(1,[2,4]),size(I) )));
      lines(i).backgrndstd(1,1)=std(double(I(rectanglespec2( lines(i).linexy(2,[1,3]), lines(i).linexy(2,[2,4]), lines(i+1).linexy(1,[1,3]), lines(i+1).linexy(1,[2,4]),size(I) ))));
   end
    
end

%% delete stripes without intersect
emptylines=NaN(size(lines)).';
for i=1:size(lines,2)
    
    emptylines(i)=isempty(lines(i).intensitymedian);
end
lines(emptylines==1)=[];
%%
%rectanglespec(lines(i).intersect1(j,:),lines(i).intersect1(j+1,:), lines(i).intersect2(j,:),lines(i).intersect2(j+1,:),size(I));

%% plot median intensity distribution

intmedian=[lines(:).intensitymedian];
intmedian1=intmedian(abs([lines(:).length]-step21)<=2);
intmedian2=intmedian(abs([lines(:).length]-step22)<=2);

intmean=[lines(:).intensitymean];
intmean1=intmean(abs([lines(:).length]-step21)<=2);
intmean2=intmean(abs([lines(:).length]-step22)<=2);

lines(1).intmedian1=intmedian1;
lines(1).intmedian2=intmedian2;


% figure
% hold on
% notBoxPlot(double(intmean1)-nanmedian([lines(:).backgrndmean]),1)
% notBoxPlot(double(intmean2)-nanmedian([lines(:).backgrndmean]),2)
%%
lines(1).filename=fname;
filename = strcat('lines', numberpattern,'.mat');
filePath = fullfile(folderName,filename);
save(filePath,'lines')

meanline = [];
for i = 1:length(lines)
    meanline = [meanline; lines(i).meanline];
end
%meanline = floor(meanline /binningFactor); %why? If I'm correct then 
%it's only important if no binning in pattern but in stack
fnameCSV = strcat('lines', numberpattern, '.csv');
CSVPath = fullfile(folderName, fnameCSV);
fid = fopen(CSVPath, 'w');
writematrix(meanline, CSVPath);
fclose(fid);

clear lines


loadfrac = (filenr-startpoint+1) / (size(FileList,1)-startpoint+1);
progressbar(loadfrac)
end

%%
scrptName = mfilename;
disp(strcat(scrptName, ' finished'))






