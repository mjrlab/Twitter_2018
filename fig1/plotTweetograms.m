%2018-06-27, EL:
%   - Fig. 1B, D: make Tweetograms color-coded by trough width

close all; clear all; clc;

%% export figures?
TOEXP = 0; %1 = yes

%% load table with counties we have data for
load('datafiles/CountyTableHaveData.mat');

%% load calculations on troughs

% parameters that make up the name of the directory and file
whichDay = 'MTWThFSatSun';
randStr = 'ckvc'; %make sure this matches the unique ID!
whichBin = {2012:2013};
binType = 'annual';
binList =  {2012:2013};
geoRegion = 'FIPS';
whichYears = [2012:2013];

% EK's calculations for 1500+ counties
INDIR = '../tweetograms/tweetograms_annual_FIPS_2012_2013_sfyo';
sjlTab = load([INDIR '/' 'troughs_2018-04-09_01.13.04_sfyo.mat']);
sjlTab = sjlTab.outTable;

% merge with county info
sjlTab = join(sjlTab,cTable); 

% rename 'MTWThFSatSun' -> 'allweek'
varNames = sjlTab.Properties.VariableNames;
varNames = cellfun(@(x) strrep(x,'MTWThFSatSun','allweek'),varNames,'UniformOutput',0);
sjlTab.Properties.VariableNames = varNames;

%% plot tweetogram for Cook, IL (Fig. 1B)
mapVar = 'allweek_trWidth';
sjlTab(isnan(sjlTab.(mapVar)),:) = [];
geoName = 'Cook IL';
cartoonColors = cool(100);
fH = plotCounty(geoName, sjlTab, mapVar, whichDay, whichBin,binType,binList,...
    geoRegion,whichYears,randStr,'mapColor',0,'colorValue',cartoonColors(50,:),...
    'figSize',[0 0 1.9 2],...
    'toArrow',1);
expif(TOEXP,fH,strrep(geoName,' ',''),'fun2use','expfig',...
    'OUTDIR','exTgram','figColor','none');

%% plot tweetograms for four counties (Fig. 1D)
fH_overlay = plotCounties({ 'Suffolk NY',...
                            'Lafayette LA',...
                            'Orange CA','Wayne MI'}, ...
                    sjlTab, mapVar, whichDay,whichBin,binType,binList,...
                    geoRegion,whichYears,randStr,...
                    'figSize',[0 0 3 2]);
expif(TOEXP,fH_overlay,'overlay','fun2use','expfig',...
    'OUTDIR','exTgram','figColor','none');

%% helper funs
function tgram = loadTgram(geoCode, whichDay,whichBin,binType,binList,...
    geoRegion,whichYears,randStr)

CALCDIR = ['../tweetograms/'...
           'tweetograms' '_' binType '_' geoRegion ...
                num2str(whichYears,'_%d') '_' randStr '/' binType];

%subdirectory name is day of week
if isnumeric(whichDay)
    SUBDIR = daynum2str(whichDay);
else
    SUBDIR = whichDay;
end

%line number is the number of the bin you want (e.g., Jan or Feb if binType is monthly)
if numel(whichBin) == 1
    if isnumeric(whichBin)
        LINENUM = find(cellfun(@(x) isequal(whichBin,x), binList));
    elseif iscell(whichBin)
        LINENUM = find(cellfun(@(x) isequal(whichBin{1},x), binList));
    elseif ischar(whichBin)
        LINENUM = find(cellfun(@(x) strcmp(whichBin,x), binList));
    end
else
    if isnumeric(whichBin)
        for w=1:numel(whichBin)
            LINENUM(w) = find(cellfun(@(x) isequal(whichBin(w),x), binList));
        end
    elseif iscell(whichBin)
        for w=1:numel(whichBin)
            LINENUM(w) = find(cellfun(@(x) isequal(whichBin{w},x), binList));
        end
    end
end

%make filename
FILENAME = [CALCDIR '/' SUBDIR '/' num2str(geoCode) '.csv'];

%csv read indexes rows/columns from 0
tgram = csvread(FILENAME,LINENUM-1,0,[LINENUM-1 0 LINENUM-1 191]);

end

function hFig = plotCounty(geoName, mapTab, mapVar, whichDay,whichBin,binType,binList,...
    geoRegion,whichYears,randStr,varargin) 

%parse arguments
myParser = inputParser;
addParameter(myParser,'mapColor',1,@isnumeric);    %pick color based on the mapVar?
addParameter(myParser,'colorValue',[0.6 0.6 0.6]); %set color manually to this value
addParameter(myParser,'figSize',[0 0 2.5 2]);
addParameter(myParser,'figColor','w');
addParameter(myParser,'toArrow',0); %make arrow showing width
parse(myParser,varargin{:});
args = myParser.Results;

thisCounty = loadTgram(geoIx(geoName,mapTab), whichDay, ...
                            whichBin, binType, binList, ...
                            geoRegion,whichYears, randStr);

%need this at least 3x plotted s.t. the endpoints match up, otherwise
%smoothing messes them up; so 4x plot it here. plot the second repeat
%(24-48 hrs)
thisCounty = [thisCounty thisCounty];                       
thisCounty = smoothify(thisCounty);

tgX = [1:numel(thisCounty)]*0.25;
tgY = thisCounty;
sampleRate = 0.01;
tgXX = tgX(1):sampleRate:tgX(end);
tgYY = interp1(tgX, tgY, tgXX);

% decide if color is determined by the value of mapVar or is pre-selected
if args.mapColor
    myColor = getColor(geoName, mapTab, mapVar);
else
    myColor = args.colorValue;
end

%make figure
hFig = figure('units','inches',...
              'position',args.figSize,...
              'color',args.figColor);
hAx = axes;
if args.toArrow
    [width, xymin, xymax, xymidleft, xymidright] = ...
        fwhm([tgXX' tgYY'],'widthAt','median','mult', 0.25);
    plot(hAx,[min(tgX) max(tgX)],[xymidleft(2) xymidleft(2)],...
        'linestyle',':','color','k');
    hold on;
end

plot(tgXX,tgYY,'-','color',myColor,'linewidth',2);
hold on;

set(hAx,'xlim',[24 48],'xtick',[0:6:48],'xticklabels',[0:6:48]-24,...
    'ylim',[-0.001 0.02],'ytick',[0:0.005:0.02],...
    'yticklabels',4*[0:0.005:0.02]);
set(hAx,'fontsize',8);
xlabel(hAx, 'local time (hr)');
ylabel(hAx,'norm. activity (AU)');

if args.toArrow
    axis(hAx);
    xymidleft(1) = xymidleft(1) + 0.1; %made the arrow a little smaller so no overlap
    xymidright(1) = xymidright(1) - 0.1;
    arrow(xymidleft,xymidright,'color','m',...
        'length',3,'width',0.5,'ends',[1 2],'tipangle',30,'baseangle',90);
    
    [~, ix] = geoIx(geoName,mapTab);
    plot(mapTab.allweek_trPos(ix)+24,mapTab.allweek_trDepth(ix),...
        'marker','o','markerfacecolor','c',...
        'markeredgecolor','none','markersize',4);
end

%set(gcf,'position',args.figSize);

end

function hFig = plotCounties(manyGeoNames, mapTab, mapVar, whichDay,whichBin,binType,binList,...
    geoRegion,whichYears,randStr,varargin) 
%plot many counties

%parse arguments
myParser = inputParser;
addParameter(myParser,'mapColor',1,@isnumeric);    %pick color based on the mapVar?
addParameter(myParser,'colorValue',[0.6 0.6 0.6]); %set color manually to this value
addParameter(myParser,'figSize',[0 0 2.5 2]);
addParameter(myParser,'figColor','w');
addParameter(myParser,'toLegend',0);
addParameter(myParser,'toArrow',0); %make arrow showing width
addParameter(myParser,'toMarkWidthHeight',1); %mark horizontal dashed line showing height where width is computed
parse(myParser,varargin{:});
args = myParser.Results;


%make figure
hFig = figure('units','inches',...
              'position',args.figSize,...
              'color',args.figColor);
hAx = axes;

%plot all the counties on it
for n=1:numel(manyGeoNames)
    geoName = manyGeoNames{n};
thisCounty = loadTgram(geoIx(geoName,mapTab), whichDay, ...
                            whichBin, binType, binList, ...
                            geoRegion,whichYears, randStr);

%need this at least 3x plotted s.t. the endpoints match up, otherwise
%smoothing messes them up; so 4x plot it here. plot the second repeat
%(24-48 hrs)
thisCounty = [thisCounty thisCounty];                       
thisCounty = smoothify(thisCounty);

tgX = [1:numel(thisCounty)]*0.25;
tgY = thisCounty;
sampleRate = 0.01;
tgXX = tgX(1):sampleRate:tgX(end);
tgYY = interp1(tgX, tgY, tgXX);

if args.toMarkWidthHeight
    [width, xymin, xymax, xymidleft, xymidright] = ...
        fwhm([tgXX' tgYY'],'widthAt','median','mult', 0.25);
    if n==1
        plot(hAx,[min(tgX) max(tgX)],[xymidleft(2) xymidleft(2)],...
            'linestyle',':','color','k');
        hold on;
    end
end

% decide if color is determined by the value of mapVar or is pre-selected
if args.mapColor
    myColor = getColor(geoName, mapTab, mapVar);
else
    myColor = args.colorValue;
end

plt(n) = plot(tgXX,tgYY,'-','color',myColor,'linewidth',1);
hold on;

if args.toArrow
    axis(hAx);
    xymidleft(1) = xymidleft(1) + 0.1; %made the arrow a little smaller so no overlap
    xymidright(1) = xymidright(1) - 0.1;
    arrow(xymidleft,xymidright,'color','m',...
        'length',3,'width',0.5,'ends',[1 2],'tipangle',30,'baseangle',90);
    
    [~, ix] = geoIx(geoName,mapTab);
    plot(mapTab.allweek_trPos(ix)+24,mapTab.allweek_trDepth(ix),...
        'marker','o','markerfacecolor','c',...
        'markeredgecolor','none','markersize',4);
end
end

set(hAx,'xlim',[24 48],'xtick',[0:6:48],'xticklabels',[0:6:48]-24,...
    'ylim',[-0.001 0.02],'ytick',[0:0.005:0.02],...
    'yticklabels',4*[0:0.005:0.02]);
if args.toLegend
    hLeg = legend(hAx,plt,manyGeoNames,'location','northoutside','box','off',...
        'Orientation','horizontal');
end
set(hAx,'fontsize',8);
xlabel(hAx, 'local time (hr)');
ylabel(hAx,'norm. activity (AU)');

end

function myColor  = getColor(geoName, mapTab, mapVar)
% cmin = min(mapTab.(mapVar));
% cmax = max(mapTab.(mapVar));

%colors are trimmed to [1,99] percentile on the map
cmin = floor(prctile(mapTab.(mapVar),1));
cmax = ceil(prctile(mapTab.(mapVar), 99));

myVar = mapTab{mapTab.geoCode == ...
                      geoIx(geoName,mapTab),mapVar};
colorFrac = (myVar - cmin)/(cmax-cmin);
colorIx = floor(512*colorFrac);

colors = jet(512);
myColor = colors(colorIx,:);

end

function dayStr = daynum2str(dayNum)
%input integer 1-7 or vector of such integers and return corresponding days
%of the week as characters. 1='M', 7='Sun'.
    daysOfWeek = {'M','T','W','Th','F','Sat','Sun'};
    if numel(dayNum) == 1
        dayStr = daysOfWeek{dayNum};
    elseif numel(dayNum) == 5 & sum(dayNum == [1:5]) == 5
        dayStr = 'wkday';
    elseif numel(dayNum) == 2 & sum(dayNum == [6:7]) == 2
        dayStr = 'wkend';
    else
        dayStr = [daysOfWeek{dayNum}];
    end
end

function [geoCode, ix] = geoIx(geoName,cTable)
if ischar(geoName)
    geoCode = cTable.geoCode(cellfun(@(x) isequal(x,geoName), cTable.geoName));
    ix = find(cellfun(@(x) isequal(x,geoName), cTable.geoName));
else
    geoCode = nan;
    ix = nan;
end
end

function fixAxes(hAx)
set(hAx,'xlim',[0 24],'xtick',0:6:24,'ytick',0:0.01:0.02,'fontsize',12);
xlabel(hAx,'time of day (hr)');
ylabel(hAx,'prob. tweeting');
% hL = legend(hAx);
% set(hL,'box','off','location','southeast');

end

function xout = smoothify(xin)
%smooth array using a Gaussian smoothing kernel
    g=gausswin(10);
    g=g/sum(g);
    smoother = (@(x) conv(x,g,'same')); %2017-10-30, EL: until today had 'smooth', but this doesn't work after update
    xout = smoother(xin); %smooth array, works
end

function expif(TOEXP,fH,fName,varargin)
myParser = inputParser;
addParameter(myParser,'OUTDIR','.',@ischar);
addParameter(myParser,'ext','pdf',@ischar);
addParameter(myParser,'fun2use','saveas',@ischar);
addParameter(myParser,'nocrop','nocrop',@ischar);
addParameter(myParser,'m','-m1',@ischar);
addParameter(myParser,'figColor','none');
parse(myParser,varargin{:});
args = myParser.Results;

%export figure with handle fH and name fName if TOEXP==1
if TOEXP == 1
    mkdir(args.OUTDIR);
    fullName = [args.OUTDIR '/' getDate('yyyy-mm-dd') '_' ...
        fName '_' getDate('HH.MM.SS')];
    for f=1:numel(fH)
        set(fH,'color',args.figColor); %change figure background to 'none' on export
    end
    switch args.fun2use
        case 'saveas'
            if numel(fH) > 1
                for i=1:numel(fH)
                    saveas(fH(i),[fullName '_' num2str(i) '.' args.ext],...
                        args.ext);
                end
            else
                saveas(fH,[fullName '.' args.ext],args.ext);
            end
        case 'expfig'
            if numel(fH) > 1
                for i=1:numel(fH)
                    export_fig([fullName '_' num2str(i) '.' args.ext],...
                        ['-' args.ext],...
                        '-painters','-cmyk','-nocrop',fH(i));
                end
            else
                export_fig([fullName '.' args.ext],['-' args.ext],...
                        '-painters','-cmyk','-nocrop',args.m,fH);
            end
        case 'print'
            if numel(fH) > 1
                for i=1:numel(fH)
                    print(fH, '-dpdf', [fullName '_' num2str(i) '.' 'pdf']);
                end
            else
                print(fH, '-dpdf', [fullName '.' '.pdf']);
            end
    end
end
end

function [width, xymin, xymax, xymidleft, xymidright] = fwhm(tweetogram,varargin)
%calculate full-width at half-max for a tweetogram
%actually using 0.25x min-to-max to avoid shoulders
%pass Tweetogram in units of HOURS!

prs = inputParser;
addParameter(prs,'widthAt','fwqm',@ischar);
addParameter(prs,'mult',0.25,@isnumeric); %multiplier of height at which to calculate width
parse(prs,varargin{:});
args = prs.Results;

%grab min point
[minY,minIx] = min(tweetogram(:,2));
sampleRate = tweetogram(2,1) - tweetogram(1,1); %spacing of points, in units of HOURS

%make sure min isn't at the extremes
if minIx < 1/sampleRate || minIx > 47/sampleRate
    [minY,minIx] = min(tweetogram(floor((1/sampleRate)):floor((47/sampleRate)),2));
    minIx = minIx + floor(1/sampleRate);
end

if minIx < 24/sampleRate
    minIx = 24/sampleRate + minIx;
    minIx = floor(minIx);
end

%grab max point that preceeds min
[maxYleft,maxYleftIx] = max(tweetogram(1:minIx,2));

switch args.widthAt
    case 'fwqm'                
        %quarter-position between min and max
        %hmY = (minY + maxYleft)*0.25; %THIS IS WRONG! 2018-02-26, EL:
        hmY = minY + (maxYleft-minY)*args.mult;
    case 'fw_maxonly'
        %relative to max only
        hmY = maxYleft*args.mult;
    case 'median'
        %quarter of the median max value (precomputed)
        hmY = 0.0182*args.mult; %the value 0.0182 is for the full dataset of 1500+ counties
end

%get value of the last point to the right of maxYleft that's > hmY
midleftIx = find(tweetogram(maxYleftIx:minIx,2) > hmY,1,'last');
midleftIx = midleftIx + maxYleftIx; %started counting from left max point

%get value of first point to the right of minY that's > hmY
midrightIx = find(tweetogram(minIx:end,2) > hmY, 1, 'first');
midrightIx = midrightIx + minIx; %started counting from minimum

%outputs
width = tweetogram(midrightIx,1) - tweetogram(midleftIx,1);
xymin = tweetogram(minIx,:);
xymax = tweetogram(maxYleftIx,:);
xymidleft = tweetogram(midleftIx,:);
xymidright = tweetogram(midrightIx,:);

%detect nan's
if isempty(midrightIx) || isempty(midleftIx)
    width = nan;
    xymin = [nan nan];
    xymax = [nan nan];
    xymidleft = nan;
    xymidright = nan;
end

end


%% Recall how input variables are defined:

% %Tweetogram for avg Mon, M-F, etc.
% whichDay = {1:5,6:7}; 
% 
% %calculated on weekly, monthly, annual basis (pick one)
% binType = 'annual'; 
% 
% %which bins (e.g., {Jan, Feb, [Jan-May]})
% whichBin = {2012:2013}; 
% 
% %FIPS, State, MSA, etc.
% geoRegion = 'FIPS';
% 
% %get data from specific years (if binType = 'annual', when make sure
% %whichBins fall within min(whichYears) and max(whichYears)
% whichYears = [2012:2013]; %list of years