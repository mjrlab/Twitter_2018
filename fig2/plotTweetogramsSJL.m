%2018-06-28, EL: Fig. 2A-B
%   - modified on 2018-08-05 to change blue colors to cyan, EL

close all; clear all; clc;

%%
TOEXP = 1; %export figures? 1=yes
lineWidth = 1; %linewidth
%% load table with data on counties
load('datafiles/CountyTableHaveData.mat');

%% Tweetogram calculation details --> these determine directory name
randStr = 'fhaw';

%Tweetogram for avg Mon, M-F, etc.
whichDay = {1:5,6:7}; 

%calculated on weekly, monthly, annual basis (pick one)
binType = 'monthly'; 

%average over which bins? (e.g., {Jan, Feb, [Jan-May]})
binList = {2,3,4,10,[2 3 4 10]}; %[Feb,Mar,Apr,Oct and their sum]
whichBin = [2 3 4 10]; %used avg over 4 months in analysis

%FIPS, State, MSA, etc.
geoRegion = 'FIPS';

%get data from specific years (if binType = 'annual', when make sure
%whichBins fall within min(whichYears) and max(whichYears)
whichYears = [2012:2013]; %list of years

%% load SJL calculations based on this data
% make directory and file names
CALCDIR = ['../tweetograms/' 'tweetograms' '_' binType '_' geoRegion ...
    num2str(whichYears,'_%d') '_' randStr '/'];
CALCFILE = 'sjl_troughs_2018-03-21_14.35.15_fhaw.mat';

% load data
sjlTab = load([CALCDIR '/' CALCFILE]);
sjlTab = sjlTab.outTable;

% merge with county info
sjlTab = join(sjlTab,cTable); 

%%
binIx = cellfun(@(x) isequal(x,whichBin),binList,'UniformOutput',1);
sjlTab.whichBin = sjlTab.xcorr_wkdaywkend(:,binIx);
sjlTab = sortrows(sjlTab,'whichBin','ascend');

%% load low-sjl county
loSJLcounty = 'Orange CA';
[geoCode, ix] = geoIx(loSJLcounty,sjlTab)
head(sjlTab(ix,{'geoCode','geoName','xcorr_wkdaywkend'}))

%must put [2 3 4 10] into a cell array to get to read it in correctly
loSJL_wkday = loadTgram(geoIx(loSJLcounty,cTable), ...
    'wkday',{whichBin},binType,binList,geoRegion,whichYears,randStr);

loSJL_wkend = loadTgram(geoIx(loSJLcounty,cTable), ...
    'wkend',{whichBin},binType,binList,geoRegion,whichYears,randStr);

%% load high-sjl county
hiSJLcounty = 'Lafayette LA';
[geoCode, ix] = geoIx(hiSJLcounty,sjlTab)
head(sjlTab(ix,{'geoCode','geoName','xcorr_wkdaywkend'}))

%must put [2 3 4 10] into a cell array to get to read it in correctly
hiSJL_wkday = loadTgram(geoIx(hiSJLcounty,cTable), ...
    'wkday',{whichBin},binType,binList,geoRegion,whichYears,randStr);

hiSJL_wkend = loadTgram(geoIx(hiSJLcounty,cTable), ...
    'wkend',{whichBin},binType,binList,geoRegion,whichYears,randStr);


%% plot low SJL county
figSize = [0 0 3 2];
curves = {loSJL_wkday,loSJL_wkend};
curves = cellfun(@(x) [x x],curves, 'UniformOutput',0); %double plot again
curves = cellfun(@(x) smoothify(x), curves,'UniformOutput',0);
tt = (0:1:numel(curves{1})-1)*0.25;

fLoSJL = figure('units','inches','position',figSize,'color','w');
plot(tt,curves{1},'k','linewidth',lineWidth,'DisplayName','weekday');
hold on;
plot(tt,curves{2},'c','linewidth',lineWidth,'DisplayName','weekend');
fixAxes(gca);
title(loSJLcounty);

%% plot high SJL county
figSize = [0 0 3 2];
curves = {hiSJL_wkday,hiSJL_wkend};
curves = cellfun(@(x) [x x],curves, 'UniformOutput',0); %double plot again
curves = cellfun(@(x) smoothify(x), curves,'UniformOutput',0);
tt = (0:1:numel(curves{1})-1)*0.25;

fHiSJL = figure('units','inches','position',figSize,'color','w');
plot(tt,curves{1},'k','linewidth',lineWidth,'DisplayName','weekday');
hold on;
plot(tt,curves{2},'c','linewidth',lineWidth,'DisplayName','weekend');
fixAxes(gca);
title(hiSJLcounty);

%% save
expif(TOEXP,fLoSJL,['tgram_' strrep(loSJLcounty,' ','')],'OUTDIR','exTgrams',...
    'fun2use','expfig','ext','pdf');
expif(TOEXP,fHiSJL,['tgram_' strrep(hiSJLcounty,' ','')],'OUTDIR','exTgrams',...
    'fun2use','expfig','ext','pdf');

%% helper funs
function tgram = loadTgram(geoCode, whichDay, whichBin, binType, binList,...
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
xlabel(hAx,'local time (hr)');
ylabel(hAx,'norm. activity (AU)');
hL = legend(hAx);
set(hL,'box','off','location','southeast','fontsize',8);
set(hAx,'xlim',[24 48],'xtick',0:6:48,'xticklabels',(0:6:48)-24,...
    'ytick',[0:0.01:0.02],'yticklabel',4*[0:0.01:0.02],'fontsize',8);

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