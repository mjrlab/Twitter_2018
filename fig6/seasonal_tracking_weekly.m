%2018-06-28, EL: Figure S4
%   - plot positions of tweetogram trough minima week by week, along with
%   times of sunrise and sunset

close all; clear all; clc;

%% input parameters
% export figure?
TOEXP = 0; %(1=yes)

% which county to plot? 
myCounty = 'Suffolk NY'; % Fig. S4 shows 'Orange CA' and 'Suffolk NY'

%plot with two y axes or no?
yyax = 0; % yyax = 0 for Figure S4.A-B, yyax = 1 for Figure S4.C-D

%% load tables with data on counties & astronomical info
%loads into cTable
load('datafiles/CountyTableHaveData.mat');

%loads into data -> turn into astro
astro = load('datafiles/astroWeekly.mat');
astro = astro.data;
astro.Properties.VariableNames{'FIPS'} = 'geoCode';

equinox = load('datafiles/AllCovariates_Apr15.mat');
equinox = equinox.AllCovariates(:,{'FIPS','GeoMeteo3'});
equinox.Properties.VariableNames{'GeoMeteo3'} = 'equinoxDawn';
equinox.Properties.VariableNames{'FIPS'} = 'geoCode';

%% parameters determining the directory name with tweetogra calculations

% tweetogram calculation details --> these determine directory name
randStr = 'kekh'; %unique identifier for every set of calculations

%Tweetogram for avg Mon, M-F, etc.
whichDay = {1:5,6:7}; 

%calculated on weekly, monthly, annual basis (pick one)
binType = 'weekly'; 

%which bins (e.g., {Jan, Feb, [Jan-May]})
binList = num2cell(1:52); % all weeks

%FIPS, State, MSA, etc.
geoRegion = 'FIPS';

%get data from specific years (if binType = 'annual', when make sure
%whichBins fall within min(whichYears) and max(whichYears)
whichYears = [2013]; %list of years

%% load tweetogram trough calculations based on this data
% make directory and file names
CALCDIR = ['../tweetograms/' 'tweetograms' '_' binType '_' geoRegion ...
    num2str(whichYears,'_%d') '_' randStr '/'];
CALCFILE_wkday = 'troughs_2018-03-26_16.29.30_kekh.mat';
CALCFILE_wkend = 'troughs_wkend2018-03-28_11.43.01_kekh.mat';

% load data
sjlTab = load([CALCDIR '/' CALCFILE_wkday]);
sjlTab_wkday = sjlTab.outTable;

sjlTab = load([CALCDIR '/' CALCFILE_wkend]);
sjlTab_wkend = sjlTab.outTable;

% merge with county info & astro
sjlTab = join(sjlTab_wkday,cTable); 
sjlTab = join(sjlTab,sjlTab_wkend,'Keys','geoCode');
sjlTab = innerjoin(sjlTab,astro,'Keys','geoCode');
sjlTab = innerjoin(sjlTab,equinox,'Keys','geoCode');
tab = sjlTab;

%% Figure S4: plot weekly dawn, dusk, wkTrough positions
fWkly = figure('units','inches','position',[0 0 3.75 3.75],'color','w');
[geoCode, ix] = geoIx(myCounty,tab);
wks=1:52;
ref = 0;%9*ones(size(wks)); %tab.dawn(ix,:);

set(fWkly,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
colors = get(groot,'DefaultAxesColorOrder');

if yyax
    yyaxis left;
end
plot(wks,tab.wkday_trPos(ix,:) - ref,'.-','color',colors(1,:),'markersize',2);
hold on;
plot(wks,tab.wkend_trPos(ix,:) - ref,'.-','color',colors(2,:),'markersize',2);
set(gca,'xlim',[0 53],'ylim',[0 9],'fontsize',8);
ylabel('tweetogram trough time (hr a.m.)');

if yyax
    yyaxis right;
end
plot(wks,tab.dawn(ix,:) - ref,'k-');
hold on;
plot(wks,tab.dusk(ix,:) - ref,'k-');
plot(wks,9*ones(size(wks)) - ref,'k:');
set(gca,'xlim',[0 53],'ylim',[0 24],'xtick',[1:4:53],'ytick',0:4:24,'fontsize',8);
if yyax
    ylabel('dawn and dusk times (hr a.m.)');
end
ylabel('local time (hr)');
xlabel('weeks of the year');
title(myCounty);

%export?
expif(TOEXP,fWkly,[strrep(myCounty,' ','') '_yy' num2str(yyax)],...
    'OUTDIR','wklyCurves','fun2use','expfig');

%% helper funs
function [geoCode, ix] = geoIx(geoName,cTable)
if ischar(geoName)
    geoCode = cTable.geoCode(cellfun(@(x) isequal(x,geoName), cTable.geoName));
    ix = find(cellfun(@(x) isequal(x,geoName), cTable.geoName));
else
    geoCode = nan;
    ix = nan;
end
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
                    print(fH(i), '-dpdf', [fullName '_' num2str(i) '.' 'pdf']);
                end
            else
                print(fH, '-dpdf', [fullName '.' '.pdf']);
            end
    end
end
end