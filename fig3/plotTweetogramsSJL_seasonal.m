%2018-06-27, EL: Fig. 3

close all; clear all; clc;

%%
TOEXP = 0; %export figure? (1=yes)
lineWidth = 1; %linewidth for Fig. 3A-C
whichCounty = 'Wayne MI'; %use Wayne, MI, for example plots throughout this figure

%% load table with data on counties
load('datafiles/CountyTableHaveData.mat');

%% Tweetogram calculation details --> these determine directory name
randStr = 'wcro'; %unique identifier for every set of calculations

%Tweetogram for avg Mon, M-F, etc.
whichDay = {1:5,6:7}; 

%calculated on weekly, monthly, annual basis (pick one)
binType = 'monthly'; 

%which bins (e.g., {Jan, Feb, [Jan-May]})
binList = {1,2,3,4,5,6,7,8,9,10,11,12}; % all months

%FIPS, State, MSA, etc.
geoRegion = 'FIPS';

%get data from specific years (if binType = 'annual', when make sure
%whichBins fall within min(whichYears) and max(whichYears)
whichYears = [2012:2013]; %list of years

%% load SJL calculations based on these parameters
% make directory and file names
CALCDIR = ['../tweetograms/' 'tweetograms' '_' binType '_' geoRegion ...
    num2str(whichYears,'_%d') '_' 'hthz' '/'];
CALCFILE = 'sjl_troughs_2018-04-08_17.32.34_hthz.mat';

% load data
sjlTab = load([CALCDIR '/' CALCFILE]);
sjlTab = sjlTab.outTable;

% merge with county info
sjlTab = innerjoin(sjlTab,cTable); 

%% load winter, Wayne County
geoIx(whichCounty,cTable)
thisMo = 2;
feb_wkday = loadTgram(geoIx(whichCounty,cTable), ...
    'wkday',thisMo,binType,binList,geoRegion,whichYears,randStr);

feb_wkend = loadTgram(geoIx(whichCounty,cTable), ...
    'wkend',thisMo,binType,binList,geoRegion,whichYears,randStr);

%% load summer, Wayne County
thisMo = 8;
aug_wkday = loadTgram(geoIx(whichCounty,cTable), ...
    'wkday',thisMo,binType,binList,geoRegion,whichYears,randStr);

aug_wkend = loadTgram(geoIx(whichCounty,cTable), ...
    'wkend',thisMo,binType,binList,geoRegion,whichYears,randStr);

%% Fig. 3A: winter
figSize = [0 0 3 2];
curves = {feb_wkday,feb_wkend};
curves = cellfun(@(x) [x x],curves, 'UniformOutput',0); %double plot again
curves = cellfun(@(x) smoothify(x), curves,'UniformOutput',0);
tt = (0:1:numel(curves{1})-1)*0.25;

fOne = figure('units','inches','position',figSize,'color','w');
plot(tt,curves{1},'k','linewidth',lineWidth,'DisplayName','weekday');
hold on;
plot(tt,curves{2},'b','linewidth',lineWidth,'DisplayName','weekend');
fixAxes(gca);
title('February');

%% Fig. 3B: summer
curves = {aug_wkday,aug_wkend};
curves = cellfun(@(x) [x x],curves, 'UniformOutput',0); %double plot again
curves = cellfun(@(x) smoothify(x), curves,'UniformOutput',0);
tt = (0:1:numel(curves{1})-1)*0.25;

fTwo = figure('units','inches','position',figSize,'color','w');
plot(tt,curves{1},'k','linewidth',lineWidth,'DisplayName','weekday');
hold on;
plot(tt,curves{2},'b','linewidth',lineWidth,'DisplayName','weekend');
fixAxes(gca);
title('August');

%% save Fig. 3A-B?
expif(TOEXP,fOne,'tgram_feb','OUTDIR','exTgrams',...
    'fun2use','expfig','ext','pdf');
expif(TOEXP,fTwo,'tgram_aug','OUTDIR','exTgrams',...
    'fun2use','expfig','ext','pdf');

%% Fig. 3C: plot seasonal changes in wkday & weekend troughs
[geoCode,ix] = geoIx(whichCounty,sjlTab)

%get trough positions and SJL curve for this county
wkdayTr = sjlTab{ix,'wkday_trPos'};
wkendTr = sjlTab{ix,'wkend_trPos'};

fTroughs = figure('units','inches','position',[0 0 3.1 2],'color','w');
monthLabels = {'Jan','','Mar','','May','','July','','Sep','','Nov',''};
plot(1:12,wkdayTr,'sk-','markersize',2,'linewidth',lineWidth,...
    'DisplayName','weekday');
hold on;
plot(1:12,wkendTr,'sb-','markersize',2,'linewidth',lineWidth,...
    'DisplayName','weekend');
fixAxesTroughs(gca);
expif(TOEXP,fTroughs,'troughsEx','OUTDIR','exTgrams',...
    'fun2use','expfig','ext','pdf');

%% Fig. 3D: seasonal changes in SJL for four example counties
countyList = {'Suffolk NY','Lafayette LA','Orange CA','Wayne MI'};
%match colors with Fig. 1
countyColors = {[0.3906    1.0000    0.6094],...
                [0    0.6328    1.0000],...
                [1.0000    0.1328         0],...
                [0    0.2812    1.0000]};
fSJLmany = plotSeasonalSJL(countyList,sjlTab,'colorList',countyColors);
set(fSJLmany,'units','inches','position',[0 0 2.92 2],'color','w');
legend off;
expif(TOEXP,fSJLmany,'sjlExMany','OUTDIR','exTgrams',...
    'fun2use','expfig','ext','pdf');

%% Fig. 3E-F: plot SJL dependence on wkdays and weekends

% first, compute seasonal fluctuations in all variables
wkdayTr = sjlTab{:,'wkday_trPos'};
wkendTr = sjlTab{:,'wkend_trPos'};
sjlCurve = sjlTab{:,'xcorr_wkdaywkend'};

%in minutes
sjlTab.wkdayFluct = 60*(wkdayTr - mean(wkdayTr,2));
sjlTab.wkendFluct = 60*(wkendTr - mean(wkendTr,2));
sjlTab.sjlFluct   = 60*(sjlCurve - mean(sjlCurve,2));

%% filter out counties with small numbers of users
%load data on number of users 
popTab = load('datafiles/AllCovariates_Apr15.mat'); %loads into AllCovariates_decorrelated
popTab = popTab.AllCovariates(:,{'FIPS','Population1','Population5'}); %.AllCovariates_decorrelated;
popTab.avgUsers = (popTab.Population5.*(10.^popTab.Population1));
popTab.Properties.VariableNames{1} = 'geoCode'; %make sure name of FIPS column is right
sjlTab= innerjoin(sjlTab,popTab(:,{'geoCode','avgUsers'}),'Keys','geoCode');

%run quality checks
[sjlTab] = qcSJL(sjlTab);

%remove noisy counties (as in Fig. 4)
AVGUSERCUTOFF = 30; %min. number of avg. daily users
NOISECUTOFF = 1;    %stdev of SJL (in hrs)
sjlTab = sjlTab(sjlTab.avgUsers > AVGUSERCUTOFF & ...
                    sjlTab.sigSJL_xc_z < NOISECUTOFF,:);

%% Fig. 3E: weekday fluctuations vs SJL fluctuations
titFile='sjlFluct943';

xVar='wkdayFluct';
xlab='fluct. in weekday trough time (min)';

yVar='sjlFluct';
ylab='fluct. in social jet lag (min)';

pltTab = sjlTab;

fSJL = figure('units','inches','position',[0 0 3.35 2.26],'color','w');
xyhistogram(pltTab,xVar,yVar,...
                 'xLab',xlab,...
                 'ylab',ylab,...
                 'hFig',fSJL(1),...
                 'fontsize',8);             
                 
expif(TOEXP,fSJL,titFile,'OUTDIR','fluct','fun2use','expfig');  

%display correlation coefficient too
format long
[rho,pval] = corr(pltTab{:,xVar}(:),pltTab{:,yVar}(:),...
    'type','Pearson',...
    'rows','complete')
[rho,pval, pvalLessThan] = rhoPval(pltTab,xVar,yVar)


%% Fig. 3F: weekend fluctuations vs SJL fluctuations
xVar='wkendFluct';
xlab='fluct. in weekend trough time (min)';

yVar='sjlFluct';
ylab='fluct. in social jet lag (min)';

pltTab = sjlTab;

fSJL = figure('units','inches','position',[0 0 3.35 2.26],'color','w');
xyhistogram(pltTab,xVar,yVar,...
                 'xLab',xlab,...
                 'ylab',ylab,...
                 'hFig',fSJL(1),...
                 'fontsize',8);             
                  
expif(TOEXP,fSJL,titFile,'OUTDIR','fluct','fun2use','expfig');  

%display correlation coefficient too
[rho,pval] = corr(pltTab{:,xVar}(:),pltTab{:,yVar}(:),...
    'type','Pearson',...
    'rows','complete')

[rho,pval, pvalLessThan] = rhoPval(pltTab,xVar,yVar)

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

function hSJLfig = plotSeasonalSJL(names,sjlTab,varargin)
myParser = inputParser;
addParameter(myParser,'colorList',[]);
parse(myParser,varargin{:});
args = myParser.Results;

%assign colors now if haven't done so yet
if isempty(args.colorList)
    myColors = parula(numel(names));
    for i=1:numel(names)
        colorList{i} = myColors(i,:);
    end
else
    colorList = args.colorList;
end

hSJLfig=figure();
for i=1:numel(names)
    [geoCode, ix] = geoIx(names{i},sjlTab);
    plt(i) = plot(1:12,sjlTab{ix,'xcorr_wkdaywkend'},'s-',...
        'markersize',2,'linewidth',1,'color',colorList{i},...
        'markerfacecolor',colorList{i},'markeredgecolor','none');
    hold on;
end
legend(plt,names);
fixAxesSJL(gca);

end

function fixAxes(hAx)
xlabel(hAx,'local time (hr)');
ylabel(hAx,'prob. of tweeting');
hL = legend(hAx);
set(hL,'box','off','location','southeast','fontsize',8);
set(hAx,'xlim',[24 48],'xtick',0:6:48,'xticklabels',(0:6:48)-24,...
    'ytick',0:0.01:0.02,'yticklabels',4*[0:0.01:0.02],'fontsize',8);

end

function fixAxesTroughs(hAx)
ylabel(hAx,'Tweetogram trough time (hr)');
monthLabels = {'Jan','','Mar','','May','','July','','Sep','','Nov',''};
hrTimes = 0:1:12;
hrLabels = num2cell(hrTimes);
hrLabels = cellfun(@(x) [num2str(x) ' a.m.'],hrLabels,'UniformOutput',0);
% hL = legend(hAx);
% set(hL,'box','off','location','south',...
%     'orientation','horizontal','fontsize',8);

%make your own labels instead of legend
text(1.5,4,'weekday','color','k','fontsize',8,...
    'VerticalAlignment','bottom','HorizontalAlignment','left');
text(1.5,6.5,'weekend','color','b','fontsize',8,...
    'VerticalAlignment','top','HorizontalAlignment','left');

set(hAx,'xlim',[1 12],'xtick',1:12,'xticklabels',monthLabels,...
    'ylim',[3.5 7],'ytick',hrTimes,'yticklabels',hrLabels,'fontsize',8);

end

function fixAxesSJL(hAx)
ylabel(hAx,'social jet lag on Twitter (min)');
monthLabels = {'Jan','','Mar','','May','','July','','Sep','','Nov',''};
set(hAx,'xlim',[1 12],'xtick',1:12,'xticklabels',monthLabels,...
    'ylim',[0.25 1.75],'ytick',0:0.5:2,'yticklabels',60*(0:0.5:2),'fontsize',8);

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

function xyscatter(myTab,xVar,yVar,varargin)
%plot xyscatter, allow to pass a bunch of parameters

p = inputParser;
addParameter(p,'titStr','',@ischar);
addParameter(p,'fontSize',8,@isnumeric);
addParameter(p,'toLegend',1);
addParameter(p,'axisFormat','square');
addParameter(p,'xLab',strrep(xVar,'_',' ')); %,@ischar);
addParameter(p,'yLab',strrep(yVar,'_',' ')); %,@ischar);
addParameter(p,'hAx',[],@ishandle);
addParameter(p,'hFig',[],@ishandle);
parse(p,varargin{:});
args = p.Results;

if isempty(args.hFig)
    args.hFig = figure();
end
if isempty(args.hAx)
    args.hAx = axes;
end

figure(args.hFig);
axes(args.hAx);

%fit
myfit = fit(myTab{:,xVar}(:),myTab{:,yVar}(:),'poly1');

%plot
plot(myTab{:,xVar}(:),myTab{:,yVar}(:),'.','markersize',0.5);
hold on;
plot(myTab{:,xVar}(:),polyval([myfit.p1 myfit.p2],myTab{:,xVar}(:)),...
    'k-','linewidth',1);

xlabel(args.xLab); 
ylabel(args.yLab);
cf = corrcoef(myTab{:,xVar}(:),myTab{:,yVar}(:),'rows','complete');
if args.toLegend
%     hL = legend(['\rho=' num2str(cf(1,2),'%2.2f')]);
%     set(hL,'box','off','location','northwest','FontSize',args.fontSize);
    legTxt = ['\rho=' num2str(cf(1,2),'%2.2f')];
    xLim = get(args.hAx,'xlim');
    yLim = get(args.hAx,'ylim');
    xTxt = xLim(1) + 0.05*(xLim(2) - xLim(1));
    yTxt = yLim(1) + 0.95*(yLim(2) - yLim(1));
    text(xTxt,yTxt,legTxt,'FontSize',args.fontSize,...
        'VerticalAlignment','top');
end
set(args.hAx,'fontSize',args.fontSize);
axis(args.hAx,args.axisFormat);

end

function xyhistogram(myTab,xVar,yVar,varargin)
%plot xyscatter, allow to pass a bunch of parameters

p = inputParser;
addParameter(p,'titStr','',@ischar);
addParameter(p,'fontSize',8,@isnumeric);
addParameter(p,'toLegend',1);
addParameter(p,'axisFormat','square');
addParameter(p,'xLab',strrep(xVar,'_',' ')); %,@ischar);
addParameter(p,'yLab',strrep(yVar,'_',' ')); %,@ischar);
addParameter(p,'hAx',[],@ishandle);
addParameter(p,'hFig',[],@ishandle);
parse(p,varargin{:});
args = p.Results;

if isempty(args.hFig)
    args.hFig = figure();
end
if isempty(args.hAx)
    args.hAx = axes;
end

figure(args.hFig);
axes(args.hAx);

%fit
myfit = fit(myTab{:,xVar}(:),myTab{:,yVar}(:),'poly1');

%find min/max, trim to 1-99% 
minVal = min(prctile(myTab{:,xVar}(:),1),prctile(myTab{:,yVar}(:),1));
maxVal = min(prctile(myTab{:,xVar}(:),99),prctile(myTab{:,yVar}(:),99));
allX = myTab{:,xVar}(:);
allY = myTab{:,yVar}(:);
badX = allX <= minVal | allX >= maxVal;
badY = allY <= minVal | allY >= maxVal;
badPts = badX | badY;

%fit, excluding extremes
myfit = fit(myTab{:,xVar}(~badPts),myTab{:,yVar}(~badPts),'poly1');
cf = corrcoef(myTab{:,xVar}(~badPts),myTab{:,yVar}(~badPts),'rows','complete');

%plot
%xEdges = [minVal:10:maxVal];
xEdges = [-90:5:90];
%xTick = round([minVal:30:maxVal]);
xTick = [-90:30:90];
yEdges = xEdges;
yTick = xTick;

h = histogram2(args.hAx,myTab{:,xVar}(:),myTab{:,yVar}(:),...
    xEdges,yEdges,...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'EdgeColor',[0.9 0.9 0.9],'LineWidth',0.05);

mymap = flipud(hot(512)); %was (hot)
mymap = [1 1 1; mymap(20:end,:)];
colormap(mymap);
hold on;

colorbar(args.hAx);

%add line -- plot(...) syntax produces a jagged line! must use line()
% plot(myTab{:,xVar}(:),polyval([myfit.p1 myfit.p2],myTab{:,xVar}(:)),...
%     'k-','linewidth',1);
% hold on;
minptX = min(xEdges);
maxptX = max(xEdges);
minptY = polyval([myfit.p1 myfit.p2],minptX);
maxptY = polyval([myfit.p1 myfit.p2],maxptX);
line([minptX maxptX],[minptY maxptY],'color','k','linewidth',1);

cf = corrcoef(myTab{:,xVar}(:),myTab{:,yVar}(:),'rows','complete');
if args.toLegend
    legTxt = ['\rho=' num2str(cf(1,2),'%2.2f')];
    xLim = get(args.hAx,'xlim');
    yLim = get(args.hAx,'ylim');
    xTxt = xLim(1) + 0.05*(xLim(2) - xLim(1));
    yTxt = yLim(1) + 0.95*(yLim(2) - yLim(1));
    text(xTxt,yTxt,legTxt,'FontSize',args.fontSize,...
        'VerticalAlignment','top');
end

axis(args.hAx,args.axisFormat);

set(args.hAx,'fontSize',args.fontSize,'xtick',xTick,'yTick',yTick);
xlabel(args.xLab); 
ylabel(args.yLab);

grid off;

set(gca,'layer','top'); %this moves axis ticks, etc., above patch
set(gcf,'Renderer','opengl'); %for proper export


end

function [outTab,...
          sigSJL_xc_mu,sigSJL_xc_sig,...
          avgSJL_xc_mu,avgSJL_xc_sig] = qcSJL(SJLtable)
%run basic quality controls on dataset
%use outputs to remove outliers
      
SJLtable(:,'haveNegXc') = rowfun(@(x) sum(x < 0) > 0,SJLtable(:,'xcorr_wkdaywkend'));
SJLtable(:,'numNegXc') = rowfun(@(x) sum(x < 0),SJLtable(:,'xcorr_wkdaywkend'));
SJLtable(:,'mostNegXc') = rowfun(@(x) x(find(x == min(x),1,'first')),SJLtable(:,'xcorr_wkdaywkend'));

%how variable is this county compared to others
SJLtable.sigSJL_xc = std(SJLtable.xcorr_wkdaywkend,0,2);
sigSJL_xc_mu = mean(SJLtable.sigSJL_xc);
sigSJL_xc_sig = std(SJLtable.sigSJL_xc);

%how is the avg sjl in this county compared to others
SJLtable.avgSJL_xc = mean(SJLtable.xcorr_wkdaywkend,2);
avgSJL_xc_mu = mean(SJLtable.avgSJL_xc);
avgSJL_xc_sig = std(SJLtable.avgSJL_xc);

%convert to z scores
SJLtable.sigSJL_xc_z = (SJLtable.sigSJL_xc - sigSJL_xc_mu)./sigSJL_xc_sig;
SJLtable.avgSJL_xc_z = (SJLtable.avgSJL_xc - avgSJL_xc_mu)./avgSJL_xc_sig;

%output
outTab = SJLtable;

end

function [rho,pval, pvalLessThan] = rhoPval(myTab,xVar,yVar)
%Keep decreasing t value for given num samples until Matlab returns a
%non-zero p-value (pvalLessThan). 
%Then we know that the true p value is at least as
%significant as pvalLessThan.

[rho,pval] = corr(myTab{:,xVar}(:),myTab{:,yVar}(:),...
    'type','Pearson',...
    'rows','complete');

if pval > 0
    pvalLessThan = [];
else
    numSamples = numel(myTab{:,xVar}(:)); %no. rows in table
    tstat = rho*sqrt((numSamples-2)/(1-rho^2));
    pval2 = pval;
    t = abs(tstat);
    while pval2 == 0 && t > 0
        t = t*0.9;
        pval2 = 1-tcdf(t,numSamples);
    end
    pvalLessThan = pval2;
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