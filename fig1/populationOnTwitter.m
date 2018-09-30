%2018-06-27, EL: Fig. S1

clear all; close all; clc;
TOEXP = 0; %1 = export figures into files

%% covariates
socioData = load('datafiles/AllCovariates_Apr15.mat');
socioData = socioData.AllCovariates;

%% SJL calculations & trough widths (annual)
INDIR = '../tweetograms/tweetograms_annual_FIPS_2012_2013_sfyo';
whichDay = 'annual/MTWThFSatSun';
sjlTab = load([INDIR '/' 'troughs_2018-04-09_01.13.04_sfyo.mat']);
sjlTab = sjlTab.outTable;

% rename 'MTWThFSatSun' -> 'allweek'
varNames = sjlTab.Properties.VariableNames;
varNames = cellfun(@(x) strrep(x,'MTWThFSatSun','allweek'),varNames,'UniformOutput',0);
sjlTab.Properties.VariableNames = varNames;
sjlTab.Properties.VariableNames{'geoCode'} = 'FIPS';

%% county health rankings: % getting < 7 hrs of sleep
sleepTab_CHR = load('datafiles/sleepTable_CHR_fracSleepUnder7.mat');
sleepTab_CHR = sleepTab_CHR.sleepTab_CHR;
sleepTab_CHR.Properties.VariableNames(1:2) = {'geoCode','CHR_fracSleepUnder7'};

%% correlate with shiftwork 
mili = readtable('datafiles/NewShiftWorkSJL.xlsx');
mili = mili(:,{'FIPS','FractionShift6pm_8am','FractionShift11pm_7am'});
miliSJL = join(mili,sjlTab(:,{'FIPS','allweek_trWidth'}),'Keys','FIPS');
[rho, pval] = corr(miliSJL{:,'FractionShift6pm_8am'},miliSJL{:,'allweek_trWidth'})
[rho,pval, pvalLessThan] = rhoPval(miliSJL,'FractionShift6pm_8am','allweek_trWidth')
%writetable(miliSJL,'mili_shiftWork_SJL.xlsx');

%% merge datasets
sjlTab = innerjoin(sjlTab,socioData,'Keys','FIPS');
sjlTab = innerjoin(sjlTab,sleepTab_CHR,'LeftKeys','FIPS','RightKeys','geoCode');

sjlTab.numTwitterUsers = (sjlTab.Population5).*(10.^(sjlTab.Population1));
sjlTab.logTwitterUsers = log10(sjlTab.numTwitterUsers);

%% Fig. S1A: histogram of frac of population on Twitter
fTwitterUsersPop = figure('units','inches','position',[0 0 2.5 2],'color','w');
h = histogram(1000*sjlTab.Population5,'NumBins',20);
xlabel('Twitter users per 1000 people');
ylabel('no. counties');
set(gca,'fontsize',8);
expif(TOEXP,fTwitterUsersPop,'TwitterHistogram',...
    'OUTDIR','twitterPopPics','fun2use','expfig');  

%% correlation between num. users and population
popSize = (10.^sjlTab.Population1);
numUsers = (sjlTab.Population5.*(10.^sjlTab.Population1));
[rho, pval] = corr(popSize,numUsers,'rows','complete','type','Pearson')

rhoTab = table;
rhoTab.popSize = popSize;
rhoTab.numUsers = numUsers;
[rho,pval, pvalLessThan] = rhoPval(rhoTab,'popSize','numUsers')
     
%% Fig. S1B: does no. Twitter users correlate with Trough width?
fOne = figure('units','inches','position',[0 0 2.7 2.7]);
titFile='numUsers_trWdith';

%x,y variable names and labels
xVar='logTwitterUsers';
xlab='log10 no. Twitter users';

yVar='allweek_trWidth';
ylab='Tweetogram trough width (hr)';

pltTab = sjlTab(~isnan(sjlTab.(xVar)) & ~isnan(sjlTab.(yVar)),:);
             
%2d histogram properties
numXDivs = 30;
loX = 1;
hiX = 5;
xTick = [1:5];

numYDivs = 30;
loY = 3;
hiY = 6;
yTick = 3:6;
          
%2d histogram with the colorbar             
fThree = figure();
set(fThree,'units','inches','position',[0 0 2.7 2.7]);

xyhistogram(pltTab,xVar,yVar,...
                 'loX',loX,'hiX',hiX,...
                 'numXDivs',numXDivs,'xTick',xTick,...
                 'loY',loY,'hiY',hiY,...
                 'numYDivs',numYDivs,'yTick',yTick,...
                 'xLab',xlab,...
                 'ylab',ylab,...
                 'hFig',fOne,...
                 'fontsize',8,...
                 'toColorBar',1,...
                 'toLegend',1);
             
[rho,pval, pvalLessThan] = rhoPval(pltTab,xVar,yVar)              
expif(TOEXP,[fOne fTwo fThree],titFile,'OUTDIR','twitterPopPics','fun2use','expfig');

%% Fig. S1C: does no. Twitter users correlate with Trough position?
fOne = figure('units','inches','position',[0 0 2.7 2.7]);
titFile='numUsers_trPos';

%x,y variable names and labels
xVar='logTwitterUsers';
xlab='log10 no. Twitter users';

yVar='allweek_trPos';
ylab='Tweetogram trough time (hr)';

pltTab = sjlTab(~isnan(sjlTab.(xVar)) & ~isnan(sjlTab.(yVar)),:);
             
%2d histogram properties
numXDivs = 30;
loX = 1;
hiX = 5;
xTick = [1:5];

numYDivs = 30;
loY = 3;
hiY = 6;
yTick = 3:6;
          
%2d histogram with the colorbar             
xyhistogram(pltTab,xVar,yVar,...
                 'loX',loX,'hiX',hiX,...
                 'numXDivs',numXDivs,'xTick',xTick,...
                 'loY',loY,'hiY',hiY,...
                 'numYDivs',numYDivs,'yTick',yTick,...
                 'xLab',xlab,...
                 'ylab',ylab,...
                 'hFig',fOne,...
                 'fontsize',8,...
                 'toColorBar',1,...
                 'toLegend',1);
 
[rho,pval, pvalLessThan] = rhoPval(pltTab,xVar,yVar)             
expif(TOEXP,[fOne],titFile,'OUTDIR','twitterPopPics','fun2use','expfig');  

%% helper funs
function xyhistogram(myTab,xVar,yVar,varargin)
%plot xyscatter, allow to pass a bunch of parameters

p = inputParser;
addParameter(p,'titStr','',@ischar);
addParameter(p,'fontSize',8,@isnumeric);
addParameter(p,'toLegend',1);
addParameter(p,'toColorBar',1);
addParameter(p,'axisFormat','square');
addParameter(p,'xLab',strrep(xVar,'_',' ')); %,@ischar);
addParameter(p,'yLab',strrep(yVar,'_',' ')); %,@ischar);
addParameter(p,'hAx',[],@ishandle);
addParameter(p,'hFig',[],@ishandle);
addParameter(p,'loX',0,@isnumeric);
addParameter(p,'hiX',1,@isnumeric);
addParameter(p,'numXDivs',20,@isnumeric);
addParameter(p,'loY',0,@isnumeric);
addParameter(p,'hiY',1,@isnumeric);
addParameter(p,'numYDivs',20,@isnumeric);
addParameter(p,'xTick',0:0.2:1,@isnumeric);
addParameter(p,'yTick',0:0.2:1,@isnumeric);
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
xEdges = [args.loX:(args.hiX-args.loX)/args.numXDivs:args.hiX];
yEdges = [args.loY:(args.hiY-args.loY)/args.numYDivs:args.hiY];

h = histogram2(args.hAx,myTab{:,xVar}(:),myTab{:,yVar}(:),...
    xEdges,yEdges,...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'EdgeColor',[0.95 0.95 0.95],'LineWidth',0.01);
mymap = flipud(hot(512)); %was (hot)
mymap = [1 1 1; mymap(20:end,:)];
colormap(mymap);
hold on;

if args.toColorBar
    colorbar(args.hAx);
end

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

set(args.hAx,'fontSize',args.fontSize,'xtick',args.xTick,'yTick',args.yTick);
xlabel(args.xLab); 
ylabel(args.yLab);

grid off;

set(gca,'layer','top'); %this moves axis ticks, etc., above patch
set(gcf,'Renderer','opengl'); %for proper export


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

