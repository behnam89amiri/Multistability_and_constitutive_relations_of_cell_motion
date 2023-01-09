%% Tracking of cell nuclei

%{
@author: Christoph Schreiber
2nd contact person: Johannes Heyn

Cell nuclei are tracked across a tif stack. The trajectory is represented in the
pase contrast stack.
Trajectories and information about the lines are saved in the folder 'tracks'.

Overview of work flow
1. line_detection
2. lines_track [you are here]
3. combine_tracks_folders
4. combine_alltracks
5. clean_alltracks
6. analysis

For more info, see README.md

Protocol
- go to folder with phase contrast and nucleus images
- paste lines folder.
- specify input parameters such as file name input, pixelsize, int2dens

Movies of tracked nuclei will be displayed. Check that BWtreshhold is
chosen fine so that the nuclei are segmented correctly. Track files
with all information and movies are saved. The file log.txt saves
the input parameters.
%}


clearvars
clearvars -global

%% Important constants, parameters and variables
pixelcalibration = 0.66;%size of one pixel in um [um/px]; 0.66 for no-binning at 10x
binning = 2;
mag = 1; %default: 1; if for whatever reason you used the occular lense, type in 1.5
bit16 = false;
bitDepth = 256;

%factor map intensity of pattern to surface concentration
int2dens = 1/21; %UNikon, 10x, Cy5 filter; Alexa 647 batch 02/21

%different, i.e. intersecting Fibronectin concentrations on pattern?
difFnConc = false;

%% variable definition for segmentation of nuclei
minarea=4; % minimal area of objects to be tracked in square pixels
BWthreshold= 0.01;

%% variables to filter out tracks
tmin = 30; %minimum duration of one track in timesteps
line_tolerance = 1.25; %to account for imperfections in the line detection
proximity_threshold_um = 10; % maximal distance for tracking between consecutive timesteps [um]
mincelldistance_um = 50; % minimum distance between neighbouring cells [um]

startpoint = 1; %change if you don't want to start from the first file

%% input GUI
useGUI = true;
furtherParameters = false;
if useGUI
    prompt = {'Start loop at position:',...1
        'Binning factor:',...2
        'Occular lense magnification factor:',...3
        'Factor to map intensity of pattern to surface density:',...4
        'Are the lines stepped?',...5
        'Is the image format 16-bit?',...6
        'Change more parameters (advanced)?'}; %7
    dlgtitle = 'Input parameters';
    dims = [1 35];
    definput = {'1',...
        '2',...
        '1',...
        '1/81',...
        'false',...
        'false',...
        'false'};
    userInput = inputdlg(prompt, dlgtitle, dims, definput);
    
    %variable assignment
    startpoint = str2double(userInput{1});
    binning = str2double(userInput{2});
    mag = str2double(userInput{3});
    int2dens = str2double(userInput{4});
    if isnan(int2dens)
        error('int2dens could not be converted to a number')
    end
    difFnConc = ~strcmp(userInput{5},'false');
    bit16 = ~strcmp(userInput{6},'false');
    furtherParameters = ~strcmp(userInput{7},'false');
end

if furtherParameters
    dlgtitleMore = 'Advanced parameters';
    promptMore = {'Line tolerance to account for imperfections in the line detection:',...1
        'Maximum displacement of cell between consecutive time steps (um):',...2
        'Minimum area of cells (px):',...3
        'Threshold for nucleus stack binarisation:',...4
        'Minimum distance between two cells (um):',...5
        'Minimum duration of a track (time steps):',...6 %convert to [min]
        'Pixel calibration: (um/px)'};%7
    definputMore = {'1.25',...1
        '10',...2
        '4',...3
        '0.01',...4
        '66',...5
        '30',...6 currently time steps instead of min
        '0.66'};%7
    userInputMore = inputdlg(promptMore, dlgtitleMore, dims, definputMore);
    
    %variable assignment
    line_tolerance = str2double(userInputMore{1});
    proximity_threshold_um = str2double(userInputMore{2});
    minarea = str2double(userInputMore{3});
    BWthreshold = str2double(userInputMore{4});
    mincelldistance_um = str2double(userInputMore{5});
    tmin = str2double(userInputMore{6}); %convert to [min]
    pixelcalibration = str2double(userInputMore{7});
    
end

pixelsize = pixelcalibration * binning / mag;
proximity_threshold = proximity_threshold_um * pixelsize;
mincelldistance = mincelldistance_um / pixelsize;
if bit16    bitDepth = 65536; end

patternPath = 'pattern';



%% script start
% Get the data for the current directory
dirFolder = dir;
dirIndex = [dirFolder.isdir];                  %# Find the index for directories
FileList = {dirFolder(~dirIndex).name}';       %# Get a list of the files, i.e. image stacks
FileList = char(FileList);
FileList = string(FileList);
FolderList = {dirFolder(dirIndex).name}';
%%
linesPath = 'lines';
pattern = 'pattern';
if ~isfolder(linesPath)
    if isfolder(fullfile(pattern, linesPath))
        linesPath = fullfile(pattern, linesPath);
    else
        error('Could not find lines folder');
    end
end

if ~isfolder('tracks')
    [status, msg, msgID] = mkdir('tracks');
    if status == 0
        error(msgID, msg);
    end
end

% Get the data for the lines directory
dirFolder2 = dir(linesPath);
dirIndex2 = [dirFolder2.isdir];                  %# Find the index for directories
FileList2 = {dirFolder2(~dirIndex2).name}';       %# Get a list of the files, i.e. mat files
FileList2 = char(FileList2);
FileList2 = string(FileList2);
FileListIdx = contains(FileList2, '.mat');
FileList2 = FileList2(FileListIdx);

% Get the data for the pattern directory
dirFolderPattern = dir(patternPath);
dirIndexPattern = [dirFolderPattern.isdir];
FileListPattern = {dirFolderPattern(~dirIndexPattern).name}';
FileListPattern = char(FileListPattern);
FileListPattern = string(FileListPattern);
FileListIdxPattern = contains(FileListPattern, '.tif');
FileListPattern = FileListPattern(FileListIdxPattern);

%% Variables to specify file name input
% cBF and cNucleus have to be uniqe markers to identify the correct channels

% Get the identification tag for the position, e.g. 'XY'
idx = contains(lower(FileList), '_xy');
if sum(idx) > 1
    cXY = '_xy';
else
    cXY = input('Enter identifiaction tag for the positions\n>', 's');
end

% Get the identification tags for the BF and nucleus channels
% cBF
idx = contains(FileList, '_BF');
if sum(idx) > 1
    cBF = '_BF';
else
    cBF = input("Enter identification tag for the BF channel. Type in 'none' if there is no BF channel\n> ", 's');
end

% cNucleus
idxTxs = contains(lower(FileList), '_texas');
idxMCherry = contains(lower(FileList), '_mcherry');
idxDapi = contains(lower(FileList), '_dapi');

if sum(idxTxs) > 1
    cNucleus = '_texas';
elseif sum(idxMCherry) > 1
    cNucleus = '_mcherry';
elseif sum(idxDapi) > 1
    cNucleus = '_dapi';
else
    cNucleus = input('Enter identification tag for the nucleus channel\n> ' , 's');
end

%% loop over different possitions (all line**.mat files in lines folder)
progressbar
try
    for filenr = startpoint:size(FileList2,1)
        %% load the lines data
        load(fullfile(linesPath, FileList2(filenr)));
        
        %% find the corresponding brightfield and Dapi movies
        numberlines = regexp(FileList2(filenr),'\d+','match');
        [~,isNumber] = str2num(numberlines);
        [~,isNumber2] = str2num(numberlines(end));
        if ~isNumber || ~isNumber2
            error("'numberlines' is not a number: '%s.'", numberlines);
        end
        
        [~, endIdx] = regexp(lower(FileList),cXY, 'once');
        emptyIdx = cellfun('isempty', endIdx);
        FileList = FileList(~emptyIdx);
        endIdxd = cell2mat(endIdx);
        FileListShrt = extractBetween(FileList, endIdxd+1, cellfun('length',FileList));
        
        FileListIdxNmb = regexp(FileListShrt,'\d+','match','once') == numberlines;
        FileListNmb = FileList(FileListIdxNmb); %causes problems when FileList doesn't contain the file matching numberlines
        FileListNmb = [FileListNmb; 'helperString']; %to make sure it's always a cell array
        
        if strcmpi(cBF, 'none')
            containsBF = 0;
        end
        if containsBF
            fnamebf = myfindfilename(FileListNmb,cBF);
        else
            fnamebf = 'noBFImage';
        end
        
        fnamekern = myfindfilename(FileListNmb, cNucleus);
        patternName = fullfile(patternPath, myfindfilename(FileListPattern, cXY + numberlines));
        
        if isempty(fnamebf)
            warning("No video matching 'XY%s' and '%s' could be found.", numberlines, cBF)
        elseif isempty(fnamekern)
            error("No video matching 'XY%s' and '%s' could be found.", numberlines, cNucleus)
        elseif isempty(patternName)
            error("No pattern matching 'XY%s' could be found.", numberlines)
        end
        
        % get infos about the file to be analyzed
        infokern = imfinfo(fnamekern);
        if containsBF
            infobf = imfinfo(fnamebf);
        end
        num_images = numel(infokern); % most importantly the number of images
        
        
        %% find position of cells using the nucleus movie
        %pre-allocate for xy_stack where all object positions are saved
        tiffimage = imread(fnamekern, num_images, 'Info', infokern);
        im_thres = mybinarize(tiffimage, bitDepth, BWthreshold);
        xy_stack = zeros(size(mypos(im_thres, minarea, num_images),1 )*num_images,3);
        xy_stack_row = 1;
        
        figurename = strcat('track_nucleus_XY_', numberlines);
        figurekern = figure('Name', figurename);
        
        tiffimage = zeros(infokern(1).Height, infokern(1).Width, 'uint8'); %pre-allocation for speed
        
        for k = 1:1:num_images
            %% cell detection
            tiffimage = imread(fnamekern, k, 'Info', infokern);
            im_thres = mybinarize(tiffimage, bitDepth, BWthreshold);
            
            xypos = mypos(im_thres, minarea, k);
            
            figure(figurekern);
            imshow(imadjust(tiffimage),'InitialMagnification','fit')
            hold on
            if size(xypos,1)>0
                plot(xypos(:,1),xypos(:,2),'x')
            end
            
            plotLines = true;
            if plotLines
                for i=1:size(lines,2)
                    plot([lines(i).linexy(1,1), lines(i).linexy(1,2)],...
                        [lines(i).linexy(1,3), lines(i).linexy(1,4)],'b')
                    plot([lines(i).linexy(2,1), lines(i).linexy(2,2)],...
                        [lines(i).linexy(2,3), lines(i).linexy(2,4)],'b')
                end
            end
            hold off
            %movie_kern(k) = getframe(gcf);%getframe is very slow. If not
            %saving the video, don't use getframe
            xy_stack_row2 = xy_stack_row + size(xypos,1) - 1;
            xy_stack(xy_stack_row:xy_stack_row2, :) = xypos;
            xy_stack_row = xy_stack_row2 +1 ;
        end
        myIdx = xy_stack(:,1) + xy_stack(:,2) > 0;
        xy_stack = xy_stack(myIdx,:);
        
        close(figurekern)
        
        %skip rest if no object has been detected
        if sum(xy_stack(:,1) + xy_stack(:,2)) == 0
            fprintf('No objects in XY%s', numberlines);
            continue
        end
        
        
        %% loop over each stripe
        
        for stripe=1:size(lines,2)
            
            %% find cell positions in proximity to the stripe
            if lines(stripe).meanline(1) == lines(stripe).meanline(2)
                xy_stripe = xy_stack(abs(xy_stack(:,1)-lines(stripe).meanline(1)) <= min(lines(stripe).width)*line_tolerance/2,:);
            else
                m = (lines(stripe).meanline(4)-lines(stripe).meanline(3)) / (lines(stripe).meanline(2)-lines(stripe).meanline(1));
                t = lines(stripe).meanline(3)-(m*lines(stripe).meanline(1));
                xy_stripe = xy_stack(abs(xy_stack(:,2)-((xy_stack(:,1)*m)+t)) <=...
                    min(lines(stripe).width)*line_tolerance/(2*cos(atan(m))), :);
            end
            % fix problem, xy_stripe sometimes empty !!!
            
            
            %         %% filter out when cells are close to each other
            %
            %         res = xy_stripe;
            %             for i= 1: num_images
            %                 resi= res(res(:,3)==i,:);
            %                 if size(resi,1)~=0
            %                     for j=1:size(resi)
            %                         %        dist= .sqrt((res(res(:,3)==i,1)-res(    )).^2-res(res(:,3)==i,2).^2)
            %                         dist= sqrt(((resi(:,1)-resi(j,1)).^2) +  ((resi(:,2)-resi(j,2)).^2));
            %                         if min(dist(1:end ~=j))<mincelldistance
            %                             resi(j,5)=0;
            %                         else
            %                             resi(j,5)=1;
            %                         end
            %                         clear dist
            %                     end
            %                     res(res(:,3)==i,5)=resi(:,5);
            %                 end
            %                 clear resi
            %
            %             end
            %             % delete when cells are close
            %             resd=res;
            %             clear res
            %             resd(resd(:,5)==0,:)=[];
            %             xy_stripe = resd;
            
            %% cell tracking
            if size(xy_stripe,1)>0 && sum(xy_stripe(:,3)-xy_stripe(1,3))~=0
                res=track(xy_stripe,proximity_threshold);
                tracks(stripe).alltracks=res;
                
                % filter out when cells are close to each other
                for i= 1: num_images
                    
                    resi = res(res(:,3)==i,1:4);
                    if size(resi,1)~=0
                        for j=1:size(resi)
                            dist = sqrt(((resi(:,1)-resi(j,1)).^2) +  ((resi(:,2)-resi(j,2)).^2));
                            if min(dist(1:end ~=j)) < mincelldistance
                                resi(j,5)=0;
                            else
                                resi(j,5)=1;
                            end
                            clear dist
                        end
                        res(res(:,3)==i,5) = resi(:,5);
                    end
                    clear resi
                    
                end
                
                % delete when cells are close
                resd=res;
                clear res
                resd(resd(:,5)==0,:)=[];
                
                resd(1,6)=1;
                for i=1:size(resd,1)-1
                    if resd(i,3)==resd(i+1,3)-1 & resd(i,4)==resd(i+1,4)
                        resd(i+1,6)=resd(i,6);
                    else
                        resd(i+1,6)=resd(i,6)+1;
                    end
                end
                
                % check whether remaining tracks are longer than tmin
                resd2=resd;
                k=0;
                for i= 1: resd(end,6)
                    if sum(resd(:,6)==i)<tmin
                        resd2(resd2(:,6)==i,:)=[];
                    else
                        k=k+1;
                        resd2(resd2(:,6)==i,7)=k;
                    end
                end
                
                
                % sort tracks
                resdsort=[];
                resds2=resd2;
                k=1;
                while size(resds2,1)>0
                    [a,b]=min(resds2(:,3));
                    resdsort=[resdsort;resds2(resds2(:,7)==resds2(b,7),:) ones(size(resds2(resds2(:,7)==resds2(b,7),:),1),1)*k];
                    resds2(resds2(:,7)==resds2(b,7),:)=[];
                    k=k+1;
                end
                clear k a b
                
                tracks(stripe).resdsort=resdsort;
                
                %% transform to coordinate system with stripes oriented parallel to x-axis
                % first for all tracks (could be combined with next loop)
                % in case of vertical lines rotate by 90 degree
                if lines(stripe).meanline(1) == lines(stripe).meanline(2)
                    tracks(stripe).alltracks(:,4) = tracks(stripe).alltracks(:,2);
                    tracks(stripe).alltracks(:,5) = tracks(stripe).alltracks(:,1)-lines(stripe).meanline(1);
                    if difFnConc
                        tracks(stripe).intersect(:,1) = lines(stripe).meanintersect(:,2);
                        tracks(stripe).intersect(:,2) = lines(stripe).meanintersect(:,1);
                    end
                    % otherwise calculate formula for meanline
                else
                    m =(lines(stripe).meanline(4)-lines(stripe).meanline(3))/(lines(stripe).meanline(2)-lines(stripe).meanline(1));
                    t = lines(stripe).meanline(3)-(m*lines(stripe).meanline(1));
                    tracks(stripe).alltracks(:,4:5) = tracks(stripe).alltracks(:,1:2);
                    tracks(stripe).alltracks(:,5) = tracks(stripe).alltracks(:,5)-t;
                    % rotate each point of the cell trajectory
                    for i = 1:size(tracks(stripe).alltracks,1)
                        tracks(stripe).alltracks(i,4:5)=([cos(atan(-m)) -sin(atan(-m)); sin(atan(-m)) cos(atan(-m))]*...
                            tracks(stripe).alltracks(i,4:5).').';
                    end
                    if sum(contains(fieldnames(lines(stripe)), 'meanintersect'))
                        if size(lines(stripe).meanintersect,1)>0
                            tracks(stripe).intersect = lines(stripe).meanintersect;
                            tracks(stripe).intersect(:,2) = tracks(stripe).intersect(:,2)-t;
                            % rotate the intersection points
                            for i = 1:size(lines(stripe).meanintersect,1)
                                tracks(stripe).intersect(i,1:2)=[cos(atan(-m)) -sin(atan(-m)); sin(atan(-m)) cos(atan(-m))]*...
                                    tracks(stripe).intersect(i,1:2).';
                            end
                        end
                    end
                    
                end
                if sum(contains(fieldnames(lines(stripe)), 'intersect'))
                    tracks(stripe).intersect=tracks(stripe).intersect*pixelsize;
                end
                tracks(stripe).alltracks(:,4:5)=tracks(stripe).alltracks(:,4:5)*pixelsize;
                
                
                % add median intensity to each position
                tracks(stripe).alltracks(:,6)=NaN;
                if (sum(contains(fieldnames(lines(stripe)), 'intersect'))) &...
                        (size(tracks(stripe).intersect,1)>1) & difFnConc
                    for i=1:size(tracks(stripe).intersect,1)-1
                        tracks(stripe).alltracks(tracks(stripe).alltracks(:,4)>...
                            tracks(stripe).intersect(i,1)&...
                            tracks(stripe).alltracks(:,4)<=...
                            tracks(stripe).intersect(i+1,1),6) =...
                            double(lines(stripe).intensitymedian(i)-median([lines.backgrndmedian]))*int2dens;
                    end
                else
                    tracks(stripe).alltracks(:,6) =...
                        double(median([lines(stripe).intensitymedian])-median([lines.backgrndmedian]))*int2dens;
                end
                
                
                
                
                % now only for filtered tracks
                if size(resdsort,1)>0
                    % in case of vertical lines rotate by 90 degree
                    if lines(stripe).meanline(1)==lines(stripe).meanline(2)
                        restrans=resdsort;
                        restrans(:,1)=resdsort(:,2);
                        restrans(:,2)=resdsort(:,1)-lines(stripe).meanline(1);
                        
                        
                    else % otherwise calculate formula for meanline
                        m=(lines(stripe).meanline(4)-lines(stripe).meanline(3))/...
                            (lines(stripe).meanline(2)-lines(stripe).meanline(1));
                        t=lines(stripe).meanline(3)-(m*lines(stripe).meanline(1));
                        restrans=resdsort;
                        restrans(:,2)=restrans(:,2)-t;
                        % rotate each point of the cell trajectory
                        for i = 1:size(restrans,1)
                            restrans(i,1:2)=[cos(atan(-m)) -sin(atan(-m)); sin(atan(-m)) cos(atan(-m))]*...
                                restrans(i,1:2).';
                        end
                        
                        
                        
                    end
                    restrans(:,1:2)=restrans(:,1:2)*pixelsize;
                    
                    tracks(stripe).xytrans=restrans;
                    
                    
                    
                    tracks(stripe).xytrans(:,9)=str2num(numberlines);
                    tracks(stripe).xytrans(:,10)=stripe;
                    % add median intensity to each position of filtered tracks
                    tracks(stripe).xytrans(:,11:12)=NaN;
                    if sum(contains(fieldnames(lines(stripe)), 'intersect')) && (size(tracks(stripe).intersect,1)>1)
                        for i=1:size(tracks(stripe).intersect,1)-1
                            tracks(stripe).xytrans(tracks(stripe).xytrans(:,1)>...
                                tracks(stripe).intersect(i,1)&tracks(stripe).xytrans(:,1)<=...
                                tracks(stripe).intersect(i+1,1),11) = ...
                                double(lines(stripe).intensitymedian(i)-median([lines.backgrndmedian]))*int2dens;
                        end
                    else
                        tracks(stripe).xytrans(:,11)=double(lines(stripe).intensitymedian-median([lines.backgrndmedian]))*int2dens;
                    end
                end
            end
        end
        
        %quick fix for above problem
        %inquire problem further. Do these ROIs really not display any nuclei?
        %Does it happen more often for vertical lines?
        if exist('tracks', 'var') == 0
            disp("Problem with xy_stripe for file " + fnamekern)
            disp("Variable 'tracks' does not exist")
            continue
        end
        
        %% write bf movie with tracks
        
        figurename = strcat('tracksall', numberlines);
        figurebf = figure('Name', figurename);
        figure(figurebf)
        
        patternIm = imread(patternName);
        imshow(imadjust(patternIm), 'InitialMagnification','fit')
        if containsBF
            showBFvid = true;
        else
            showBFvid = false;
        end
        
        if showBFvid == true
            for k = 1:1:num_images
                figure(figurebf);
                bfimage = imread(fnamebf, k, 'Info', infobf);
                imshow(bfimage,'InitialMagnification','fit')
                hold on
                %plot trajectories
                for l=1:size(tracks,2) %loop over all lanes
                    if numel(tracks(l).resdsort)>0
                        resd3 = tracks(l).resdsort(tracks(l).resdsort(:,3)<=k,:);
                        if size(resd3,1)~=0
                            for i=1:resd3(end,8) %loop over all tracks on lane
                                ax = gca;
                                ax.ColorOrderIndex = mod(i,7)+1;
                                plot(resd3(resd3(:,8)==i,1),resd3(resd3(:,8)==i,2))
                            end
                        end
                    end
                    clear resd3
                end
                hold off
            end
        end
        
        %show pattern with trajectories
        imshow(imadjust(patternIm), 'InitialMagnification', 'fit')
        
        if plotLines
            hold on
            for i=1:size(lines,2)
                plot([lines(i).linexy(1,1), lines(i).linexy(1,2)],...
                    [lines(i).linexy(1,3), lines(i).linexy(1,4)],'b')
                plot([lines(i).linexy(2,1), lines(i).linexy(2,2)],...
                    [lines(i).linexy(2,3), lines(i).linexy(2,4)],'b')
            end
            hold off
        end
        
        hold on
        for l=1:size(tracks,2) %loop over all lanes
            if numel(tracks(l).resdsort)>0
                resd3 = tracks(l).resdsort(:,:);
                if size(resd3,1)~=0
                    for i=1:resd3(end,8) %loop over all tracks on lane
                        
                        %observed several examples of cells that are
                        %plotted inbetween FN lanes. Problem with
                        %filtering?
                        
                        ax = gca;
                        ax.ColorOrderIndex = mod(i,7)+1;
                        plot(resd3(resd3(:,8)==i,1),resd3(resd3(:,8)==i,2))
                    end
                end
            end
            clear resd3
        end
        hold off
        %movie_bftracks(k) = getframe(gcf); %getframe is very slow. If not
        %saving the video, don't use getframe
        
    filename = strcat('tracksall', numberlines,'.mat');
    filename = fullfile('tracks', filename);
    save(filename, 'tracks', 'lines', 'int2dens', 'fnamebf')
    
    clear res resd resds2 resdsort restrans tracks
    progressbar((filenr - startpoint +1)/(size(FileList2,1) - startpoint +1))
    end

% %% make gif
% bla=figure;
% filename = 'movie_kern.gif';
% for n = 1:size(movie_bftracks,2)
%     imshow(movie_kern(n).cdata)
%     frame = getframe(bla);
%     im = frame2im(frame);
%     [A,map] = rgb2ind(im,256);
%     if n == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
%     end
% end
%



%%
catch ME
    myLog(binning, mag, int2dens, bit16, line_tolerance,...
        proximity_threshold, minarea,  mincelldistance, tmin, ...
        pixelcalibration);
    rethrow(ME)
end
    
    
    myLog(binning, mag, int2dens, bit16, line_tolerance, proximity_threshold, ...
        minarea,  mincelldistance, tmin, pixelcalibration);
    
    scrptName = mfilename;
    disp(strcat(scrptName, ' finished'))
    
    
    %% function declarations
    
    function status = myLog(binning, mag, int2dens, bit16, line_tolerance,...
        proximity_threshold, minarea,  mincelldistance, tmin, ...%BWthreshold 4th to last in list
        pixelcalibration)
    % save parameters to log
    
    fileID = fopen('log.txt', 'w');
    fprintf(fileID, 'Date: %s\n', datetime('today'));
    
    formatSpec = "Binning factor: %d\n" +...
        "Occular lense magnification factor: %.1f\n" +...
        "Factor to map intensity of pattern to surface density: %f\n" +...
        "Is the image format 16-bit? %d\n" +...
        "Line tolerance %.2f\n" +...
        "Maximum displacement of cell between consecutive time steps (um): %.3f\n" +...
        "Minimum area of cells (px): %d\n" +...    %"Threshold for nucleus stack binarisation: %.3f\n" +...%BWthreshol
        "Minimum distance between two cells (um): %.3f\n" +...
        "Minimum duration of a track (time steps): %d\n" +...
        "Pixel calibration (um/px): %.3f";
    fprintf(fileID, formatSpec, binning, mag, int2dens, bit16, line_tolerance,...
        proximity_threshold, minarea,  mincelldistance, tmin, ...%BWthreshold 4th to last in list
        pixelcalibration);
    status = fclose('all');
    end