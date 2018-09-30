%2018-06-28, EL: Figure 6
%   - seasonal changes in tweetogram minima, month by month
%   - linear fits to tweetogram time vs. day length

close all; clear all; clc;

% export figures? 
TOEXP = 0; %(1 = yes)

% which counties to plot as examples in panels A-B, D-E?
% must be one of the ~940 counties in this analysis (see variable 'tab' below)
whichCounty_One = 'Suffolk NY'; % used 'Suffolk NY' in text
whichCounty_Two = 'Orange CA';  % used 'Orange CA' in text

%% load tables with data on counties & astronomical info
%loads into cTable
load('datafiles/CountyTableHaveData.mat');

%astronomical data loads into data -> turn into astro
astro = load('datafiles/astroMonthly.mat');
astro = astro.data;
astro.Properties.VariableNames{'FIPS'} = 'geoCode';

%% parameters used to create folder name for SJL calculations

% Tweetogram calculation details --> these determine directory name
randStr = 'hthz'; %unique identifier for every set of calculations

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

%% load SJL calculations from this folder
% make directory and file names
CALCDIR = ['../tweetograms/' 'tweetograms' '_' binType '_' geoRegion ...
    num2str(whichYears,'_%d') '_' randStr '/'];
CALCFILE = 'sjl_troughs_2018-04-08_17.32.34_hthz.mat';

% load data
sjlTab = load([CALCDIR '/' CALCFILE]);
sjlTab = sjlTab.outTable;

% load list of FIPS used in Fig. 4 analysis, only include counties matching
% the same inclusion criteria (i.e., >30 users per day on avg, etc.)
fipsList = load('datafiles/FipsListForClust_avgUsr30_1_2018-04-18_13.13.06.mat');
fipsList = fipsList.FipsListForClust;

% merge with county info & astro
sjlTab = join(sjlTab,cTable); 
sjlTab = innerjoin(sjlTab,astro,'Keys','geoCode');
%sjlTab = innerjoin(sjlTab,equinox,'Keys','geoCode');
sjlTab = innerjoin(sjlTab,fipsList,'Keys','geoCode');
tab = sjlTab(:,{'geoCode','geoName','state','population','popdensity',...
               'latitude','longitude','timezone',...
               'dawn','dusk','daylength',...%'equinoxDawn',...
               'wkday_trPos','wkend_trPos'});

           
%% rename variables
tab.wkdayFit = tab.wkday_trPos;
tab.wkendFit = tab.wkend_trPos;
tab.Properties.VariableNames{'geoCode'} = 'FIPS'; %for compatibility with older code
tab.Properties.VariableNames{'geoName'} = 'name'; %for compatibility with older code

%% calculate times relative to dawn
%hrs after dawn
tab.wkendTroughAD = tab.wkendFit - tab.dawn;
tab.wkdayTroughAD = tab.wkdayFit - tab.dawn;

%hrs after 9 a.m.e
tab.wkendTroughAfter9 = tab.wkendFit - 9; %already in hours
tab.wkdayTroughAfter9 = tab.wkdayFit - 9;
tab.NineAM_AD = 9 - tab.dawn; %9 am in hrs after dawn

%% do linear fits or load them from file
%do you want to run new fits or load saved fit values?
%by default, load saved values
calcOrload = 'load'; 
switch calcOrload
    case 'calc'
        %details of calculations and definitions of columnnames are
        %described in fitAllQuick() function in this file
        tabTrack = fitAllQuick(tab);
        save(['allFits_dawnDusk9am_' getDate() '.mat'],'tabTrack');
    case 'load'
        load('datafiles/allFits_dawnDusk9am_2018-05-08_00.14.59.mat'); %dstLabels = {'DST','Standard','noSummer'}; dstMonths = {4:10,[1 2 11 12],[4 5 9 10]}
        tabTrack = tabTrack(ismember(tabTrack.FIPS,fipsList.geoCode),:);
end

%% parameters used to select data for calculating day length tracking slopes
% also see fitAllQuick() 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT! These parameters must match what you used in calculations in
% fitAllQuick(). Please refer to that function description for detailed
% definitions and descriptions of how variable naming works.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%categories for what time is measured relative to (dawn or 9 a.m.)
refTimeLabels = {'AD', 'After9'};      
refTrough = {'AD','After9'};                                    

%categories for which months to choose (when Daylight savings is in effect
%or not)
dstLabels = {'DST','Standard','noSummer'}; 
dstMonths = {4:10,[1 2 11 12],[4 5 9 10]};

%% Figure 6C,F: plot distributions of day length tracking slopes

% in Figure 6C, use refIx = 2, dstIx = 1
% in Figure 6F, use refIx = 2, dstIx = 3

%in local time
%which value of refTrough do you want? 
refIx = 2; %(refIx=2 refers to the second entry in refTrough, 'After9')
dstIx = 3; %(dstIx=1 refers to the first entry in dstLabels, 'DST')
trackSuffix = ['_' refTimeLabels{refIx} '_' dstLabels{dstIx}];

wkdayVar = ['wkdayTrack' trackSuffix];
wkendVar = ['wkendTrack' trackSuffix];

%print min and max slopes to make sure axis limits are set correctly
cellfun(@(x) disp([x ', min m=' num2str(min(tabTrack{:,x}.m))...
                     ', max m=' num2str(max(tabTrack{:,x}.m))]),...
                     {wkdayVar, wkendVar});

%make figure
figM = myHist(tabTrack,{wkdayVar, wkendVar},...
    {'weekday','weekend'},'m','xLab','slope m','BinLimits',[-0.8 0.8],...%[-0.6 0.6],...
    'xtick',[-0.8:0.2:0.8],...
    'FaceColor',{'auto','auto','k','k'},...
    'figSize',[0 0 2.75 2.5]);%,varargin)
set(gca,'ylim',[0 0.23]); %make all y axes comparable
hold on;

%mark dawn and dusk
line([median(tabTrack{:,'Dawn_Local_Track'}.m) ...
      median(tabTrack{:,'Dawn_Local_Track'}.m)],...
      [0 max(get(gca,'ylim'))],...
    'linewidth',2,'linestyle','--','color','k');
line([median(tabTrack{:,'Dusk_Local_Track'}.m) ...
      median(tabTrack{:,'Dusk_Local_Track'}.m)],...
      [0 max(get(gca,'ylim'))],...
    'linewidth',2,'linestyle','--','color','k');
expif(TOEXP,figM,['mSlope' trackSuffix],...
        'OUTDIR','trackingFigs_9am',...
        'fun2use','print');

% compare distributions using a t test
[ttH,ttP,ttCI,ttStats] = ttest(tabTrack{:,wkdayVar}{:,'m'},...
                               tabTrack{:,wkendVar}{:,'m'});
disp(['-log10 p value for t test: ' num2str(-log10(ttP),'%5.5f')]);                           

%% print slope values
[geoCode, ixOne] = geoIx(whichCounty_One,tabTrack)
whichCounty_One
slopeTab = tabTrack(ixOne,{'name',...
                            'wkdayTrack_AD_DST','wkendTrack_AD_DST',...
                            'wkdayTrack_After9_DST','wkendTrack_After9_DST',...
                            'wkdayTrack_After9_noSummer','wkendTrack_After9_noSummer'});
disp('wkday AD DST'); slopeTab.wkdayTrack_AD_DST.m
disp('wkend AD DST'); slopeTab.wkendTrack_AD_DST.m
disp('wkday after 9 DST'); slopeTab.wkdayTrack_After9_DST.m
disp('wkend after 9 DST'); slopeTab.wkendTrack_After9_DST.m
disp('wkday after 9 no Summer'); slopeTab.wkdayTrack_After9_noSummer.m
disp('wkend after 9 no Summer'); slopeTab.wkendTrack_After9_noSummer.m


[geoCode, ixTwo] = geoIx(whichCounty_Two,tabTrack)
whichCounty_Two
slopeTab = tabTrack(ixTwo,{'name',...
                            'wkdayTrack_AD_DST','wkendTrack_AD_DST',...
                            'wkdayTrack_After9_DST','wkendTrack_After9_DST',...
                            'wkdayTrack_After9_noSummer','wkendTrack_After9_noSummer'});
disp('wkday AD DST'); slopeTab.wkdayTrack_AD_DST.m
disp('wkend AD DST'); slopeTab.wkendTrack_AD_DST.m
disp('wkday after 9 DST'); slopeTab.wkdayTrack_After9_DST.m
disp('wkend after 9 DST'); slopeTab.wkendTrack_After9_DST.m
disp('wkday after 9 no Summer'); slopeTab.wkdayTrack_After9_noSummer.m
disp('wkend after 9 no Summer'); slopeTab.wkendTrack_After9_noSummer.m
   
%% Figure 6A-B,D-E: Plot dawn/dusk-tracking curves + dawn and 9 a.m.
% RELATIVE TO 9 a.m. IN THIS CELL

refTrough = {'AD','After9'};                                    
dstLabels = {'DST','Standard','noSummer'}; 
dstMonths = {4:10,[1 2 11 12],[4 5 9 10]};

% for Figure 6A-B, use refIx = 2, dstIx = 1
% for Figure 6D-E, use refIx = 2, dstIx = 3
refIx = 2; %2 = use local time calcs
dstIx = 1; %1 = use DST calcs
astroVar = ['daylength']; %plot daylength on x axis

% generate correct column names and look up slopes according to
% AD/After9 and DST/Std settings above
whichDays = {'wkday','wkend'};
trackSuffix = ['_' refTimeLabels{refIx} '_' dstLabels{dstIx}];
trackVar = cellfun(@(x) [x 'Trough' refTrough{refIx}],whichDays,'uniformoutput',0);
trackFitVar = cellfun(@(x) [x 'Track' trackSuffix],whichDays,'UniformOutput',0);
fitCoef = cellfun(@(x) [tabTrack.(x).m tabTrack.(x).b tabTrack.(x).rsq],trackFitVar,'UniformOutput',0);

%% plot based on settings above
xVar = astroVar;
yVar = trackVar;
ix2fit = dstMonths{dstIx};
ix2plot = 1:12;

refPoint = 'midnight'; %plot relative to mindight (local time) or dawn

fTrack(1) = plotTrack9am(tabTrack,whichCounty_One,...
                    xVar,yVar,ix2fit,ix2plot,fitCoef,...
                    'titStr',whichCounty_One,...
                    'xLab','day length (hr)',...
                    'yLab','',...
                    'refPoint',refPoint,...
                    'adj2data',0);

fTrack(2) = plotTrack9am(tabTrack,whichCounty_Two,...
                    xVar,yVar,ix2fit,ix2plot,fitCoef,...
                    'titStr',whichCounty_Two,...
                    'xLab','day length (hr)',...
                    'yLab','',...
                    'refPoint',refPoint,...
                    'adj2data',0);

expif(TOEXP,fTrack,['trackEx' trackSuffix '_ref2' refPoint],...
        'OUTDIR','trackingFigs_9am',...
        'fun2use','expfig');
    
%% map of US colored according to the day length tracking slopes
mapVarName = trackFitVar{2};
mapVarVal = tabTrack.(mapVarName).m;
mapTab = tabTrack(:,{'FIPS','name','state'});
mapTab{:,mapVarName} = mapVarVal;
mapFig = plotMap_Polygon_label(mapTab, 'FIPS', mapVarName,...
    strrep(mapVarName,'_','  '),[]);
TOEXPMAP = 1;
expif(TOEXPMAP,mapFig,['map_' 'trackEx' trackSuffix '_ref2' refPoint],...
        'OUTDIR','map',...
        'fun2use','expfig');
    
%% helper funs
function tabOut = fitAllQuick(tabIn)
%function to fit straight lines to troughTime vs daylength curves.
%tabIn must have columns {'wkendTroughAD', 'wkdayTroughAD',...
%                         'wkendTroughAfter9', 'wkdayTroughAfter9', ....
%                         'NineAM_AD', 'daylength', 'dawn','dusk'}  
%E.g., column wkendTroughAD has the monthly trough positions on weekends,
%reported relative to sunrise.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we select according to what reference we're measuring time and
% which months of the year are used in fitting:
%
%  refTrough:
%   - 'AD'  means time is measure after dawn (changes throughout
%  the year)
%   -'After9' means time is measured after 9 a.m.
%   - refTimeLabels must contain as many entries as refTrough 
%
%  dstMonths: 
%   - any vector of numbers 1-12, corresponding to months of the
%     year used in the fit 
%   - [4 5 8] = [Apr May Aug]
%   - dstLabels mus contain as many entries as dstMonhts
%
% To compare fits to different months and time reference points, 
%   you can pass a cell array of multiple options, 
%   and linear fits be computed for all pairs of possibilities. Weekdays
%   and weekends will be fit separately. If you want to only look at fits
%   in local time and 
%
% E.g., if refTrough = {'AD','After9'}, dstMonths = {4:10,7:9}:
%   - the code will do fits in both 'dawn' and '9 a.m.' timeframes
%   - and separately for Apr-Oct and also for Jul-Sep
%   - so each county will be fit in four ways 
% Fit results are stored in a tabOut, one row for each county. Each fit is
% placed in a separate column named according to refTroughLabels and
% dstLabels. In this example, if refTroughLabels = {'AD','After9'} and
% dstLabels = {'Apr2Oct','Jul2Sep'}, then the column names for weekend fits 
% will be 'wkendTrack_Apr2Oct_AD', 'wkendTrack_Jul2Sep_AD',
% 'wkendTrack_Apr2Oct_After9','wkendTrack_Jul2Sep_After9'. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%categories for what time is measured relative to (dawn or 9 a.m.)
refTimeLabels = {'AD', 'After9'};       
refTrough = {'AD','After9'};            

%categories for which months to choose (when Daylight savings is in effect
%or not)
dstLabels = {'DST','Standard','noSummer'};  %
dstMonths = {4:10,[1 2 11 12],[4 5 9 10]};  %,[1 2 11 12],1:12};

%run linear fits
for dst=1:numel(dstLabels)
    disp(['Running fit on dst ' num2str(dst) ' of ' ...
        num2str(numel(dstLabels))]);
    
    for tr=1:numel(refTimeLabels)
        disp(['Running fit on troughType ' num2str(tr) ' of ' ...
            num2str(numel(refTimeLabels))]);
        
        trackWkday = [];
        trackWkend = [];
        
        wkendTrackName = ['wkendTrack_' refTimeLabels{tr} '_' dstLabels{dst}];
        wkdayTrackName = ['wkdayTrack_' refTimeLabels{tr} '_' dstLabels{dst}];
        
        %make smaller tables for fitting
        monthIx = dstMonths{dst};
        refWkday = ['wkdayTrough' refTrough{tr}];
        refWkend = ['wkendTrough' refTrough{tr}];
        
        tab2fit = tabIn(:,{'daylength',refWkday,refWkend});
        tab2fit_small = varfun(@(x) x(:,dstMonths{dst}),tab2fit);
        tab2fit_small.Properties.VariableNames = tab2fit.Properties.VariableNames;

        %fit
        trackWkday = rowfun(@(x,y) ffun(x',y'), ...
            tab2fit_small(:,{'daylength',refWkday}));
        
        trackWkend = rowfun(@(x,y) ffun(x',y'), ...
            tab2fit_small(:,{'daylength',refWkend}));
        
        tabIn.(wkdayTrackName) = trackWkday.Var1;
        tabIn.(wkendTrackName) = trackWkend.Var1;
                
    end
end

%also calculate slope of 9 a.m. relative to dawn, slopes of dawn and dusk
%in local time
tab2fit = tabIn(:,{'daylength','NineAM_AD','dawn','dusk'});
tab2fit_small = varfun(@(x) x(:,dstMonths{dst}),tab2fit);
tab2fit_small.Properties.VariableNames = tab2fit.Properties.VariableNames;

tabIn.('NineAM_AD_Track') = rowfun(@(x,y) ffun(x',y'),...
                                tab2fit_small(:,{'daylength','NineAM_AD'}));

tabIn.('Dawn_Local_Track') = rowfun(@(x,y) ffun(x',y'),...
                                tab2fit_small(:,{'daylength','dawn'}));
                            
tabIn.('Dusk_Local_Track') = rowfun(@(x,y) ffun(x',y'),...
                                tab2fit_small(:,{'daylength','dusk'}));                                  

%a bit of clean up because the code above generates awkward nested columns 
tabIn.NineAM_AD_Track = tabIn.NineAM_AD_Track.Var1;
tabIn.Dawn_Local_Track = tabIn.Dawn_Local_Track.Var1;
tabIn.Dusk_Local_Track = tabIn.Dusk_Local_Track.Var1;

%output table
tabOut = tabIn;

    %local helper function to do linear fits for tabular data and output
    %fit results in a table
    function tabFit = ffun(x,y)
        [fitobj,gof] = fit(x,y,'poly1','normalize','off');
        tabFit = table(fitobj.p1, fitobj.p2, gof.rsquare, ...
            'VariableNames',{'m','b','rsq'});
    end

end

function [hFigUSA_Clust, shapeStruct] = ...
    plotMap_Polygon_label(ds, geoUnit, colorBy,mapTitStr,colorLims)
% 2018-01-26, EL fixed. 
%   - plot us map and color all fips as polygons

%geoUnit = FIPS, MSA --> says which regions to plot over
switch geoUnit
    % load the shape file = shapeStruct
    case 'FIPS'
        load('datafiles/mapUSA_lores_counties_shape.mat'); %pre-loaded for faster plotting
    case 'MSA'
        geoUnit = 'CBSA';
        load('datafiles/mapUSA_shape_CBSAonly.mat');
    case 'CBSA'
        load('datafiles/mapUSA_shape_CBSAonly.mat');
end
%load matlab
whoseMap = 'Matlab';
load('datafiles/mapUSA_shape_states_Matlab.mat');

%fill colors
dummy = num2cell(nan(numel(shapeStruct),1));
[shapeStruct(:).color] = deal(dummy{:});
for i=1:numel(shapeStruct)
    ix = find(ds.(geoUnit) == shapeStruct(i).(geoUnit));
    if ~isempty(ix)
        %Had a bug here before!
        %Was selecting shapeStruct(ix)! instead of shapeStruct(i)!!!
        shapeStruct(i).color = ds{ix,colorBy};
        
        %test by assigning random colors
        %shapeStruct(i).color = rand*100;
    end
end

% colors = 1000 colors from jet colormap
if isempty(colorLims)
    colorLims = [prctile(ds.(colorBy),1) ...
                 prctile(ds.(colorBy), 99)]; 
       
    %2018-06-03, EL:
    %should manually trim the values to the [1,99] percentiles. it seems like
    %Matlab assigns them to the median value otherwise...
    ds{ds.(colorBy) > colorLims(2),colorBy} = colorLims(2);
    ds{ds.(colorBy) < colorLims(1),colorBy} = colorLims(1);
end
mySymbolSpec = makesymbolspec('Polygon',...
    {'color',colorLims,... %replaced w. colorLims Mar23; before had  [min(ds.(colorBy)) max(ds.(colorBy))]trying this instead of cLims
    'FaceColor',jet(512)}); %jet by default

%get indices of Hawaii, Alaska and continental states in your dataset
S_ix_Hawaii = strcmp(ds.state,'HI');
S_ix_Alaska = strcmp(ds.state,'AK');
S_ix_ConUS = ~(S_ix_Hawaii | S_ix_Alaska);

hFigUSA_Clust = figure();
set(hFigUSA_Clust,'position',[2 2 11 8.5]);
[~, hAxUSA] = plotUSA(hFigUSA_Clust);

% mark all counties we have with empty squares
%continental USA
axes(hAxUSA(1));
hold on;
myds = ds(S_ix_ConUS,:);
drawFIPS(myds,shapeStruct,geoUnit,'none','k');
[myFIPSix,myShapeStruct] = getFIPSix(myds.(geoUnit), shapeStruct, geoUnit);
geoshow(myShapeStruct,'SymbolSpec',mySymbolSpec);
drawStates(whoseMap,'ConUS',stateShapeStruct,'none','k');

%Alaska
axes(hAxUSA(2));
myds = ds(S_ix_Alaska,:);
if numel(myds{:,1} > 0)
drawFIPS(myds,shapeStruct,geoUnit,'none','k');
[myFIPSix,myShapeStruct] = getFIPSix(myds.(geoUnit), shapeStruct, geoUnit);
geoshow(myShapeStruct,'SymbolSpec',mySymbolSpec);
end
drawStates(whoseMap,'AK',stateShapeStruct,'none','k');

%Hawaii
axes(hAxUSA(3));
myds = ds(S_ix_Hawaii,:);
if numel(myds{:,1} > 0)
    drawFIPS(myds,shapeStruct,geoUnit,'none','k');
    [myFIPSix,myShapeStruct] = getFIPSix(myds.(geoUnit), shapeStruct, geoUnit);
    % disp(numel(myShapeStruct));
    % disp({myShapeStruct.NAME});
    % disp({myShapeStruct.color});
    geoshow(myShapeStruct,'SymbolSpec',mySymbolSpec);
end
drawStates(whoseMap,'HI',stateShapeStruct,'none','k');

%label plot with cluster number
axes(hAxUSA(1));
caxis manual;
caxis(colorLims);
colormap(jet(512));
hBar = colorbar;
set(hBar,'fontsize',8);
textm(51,-100,mapTitStr,'color','k');
end

function drawFIPS(myds, shapeStruct, geoUnit, faceColor, edgeColor)
%draw FIPS geoboundaries in myds.FIPS from shapeFile loaded into 
%shapeStruct; color all FIPS according to myColor
if numel(myds{:,geoUnit} > 1)
    [myFIPSix, myShapeStruct] = getFIPSix(myds.(geoUnit), shapeStruct, geoUnit);

    geoshow(myShapeStruct,'DisplayType','Polygon',...
        'FaceColor',faceColor,'edgecolor',edgeColor,...
        'linewidth',0.1); %sually 0.01 --> too small in TIFs, make this 0.05
    hold on;
end
end

function drawStates(whoseMap, whichStates, shapeStruct, faceColor, edgeColor)
switch whoseMap
    case 'Census'
        switch whichStates
            case 'ConUS' %less that 70 (72 = Puerto Rico) and not 15 or 2
                ix = ~ismember([shapeStruct(:).FIPS],[2 15 72]);
            case 'HI' %state 15
                ix = ismember([shapeStruct(:).FIPS],[15]);
            case 'AK' %state 2
                ix = ismember([shapeStruct(:).FIPS],[2]);
        end
    case 'Matlab'
        switch whichStates
            case 'ConUS' %less that 70 (72 = Puerto Rico) and not 15 or 2
                ix = ~ismember({shapeStruct(:).Name},{'Hawaii','Alaska'});
            case 'HI' %state 15
                ix = ismember({shapeStruct(:).Name},'Hawaii');
            case 'AK' %state 2
                ix = ismember({shapeStruct(:).Name},'Alaska');
        end
end
    myShapeStruct = shapeStruct(ix);

    geoshow(myShapeStruct,'DisplayType','Polygon',...
        'FaceColor',faceColor,'edgecolor',edgeColor,...
        'linewidth',0.25);
    hold on;
end

function [myFIPSix,myShapeStruct] = getFIPSix(myFIPS, shapeStruct,geoUnit)
%return where indices of myFIPS (vector) are in structure shapeStruct,
%loaded from a shapeFile. compare myFIPS to FIPS field of shapeStruct.
%geoUnit specifies whether we care about 

if strcmp(geoUnit,'MSA')
    geoUnit = 'CBSA';
end

logicalIx = ismember([shapeStruct.(geoUnit)],myFIPS);
myFIPSix = shapeStruct(logicalIx).(geoUnit);
myShapeStruct = shapeStruct(logicalIx);

end

function [geoCode, ix] = geoIx(geoName,cTable)
if ischar(geoName)
    geoCode = cTable.FIPS(cellfun(@(x) isequal(x,geoName), cTable.name));
    ix = find(cellfun(@(x) isequal(x,geoName), cTable.name));
else
    geoCode = nan;
    ix = nan;
end
end

function fHist = myHist(myTab,varList,varNameList,subVar,varargin)
myParser = inputParser;
addParameter(myParser,'OUTDIR','.',@ischar);
addParameter(myParser,'figColor','w');
addParameter(myParser,'figSize',[0 0 2.5 2.5]);
addParameter(myParser,'BinWidth',0.04);
addParameter(myParser,'BinLimits',[0 1]);
addParameter(myParser,'yLab','frequency');
addParameter(myParser,'xLab','');
addParameter(myParser,'xtick',0:0.2:1);
addParameter(myParser,'titStr','Histogram');
addParameter(myParser,'toLegend',0);
addParameter(myParser,'fontSize',8);
addParameter(myParser,'FaceColor',{'auto'});
parse(myParser,varargin{:});
args = myParser.Results;

fHist = figure('units','inches','color',args.figColor,...
    'position',args.figSize);
if numel(args.FaceColor) < numel(varList)
    for v=1:numel(varList)
        myFaceColor{v} = 'auto'; 
    end
else
    myFaceColor = args.FaceColor;
end

for v=1:numel(varList)
    if isempty(subVar)
        plt(v) = histogram(myTab{:,varList{v}},...
            'BinLimits',args.BinLimits,...
            'BinWidth',args.BinWidth,...
            'Normalization','Probability','FaceColor',myFaceColor);
        hold on;
    else
        plt(v) = histogram(myTab{:,varList{v}}{:,subVar},...
                                   'BinLimits',args.BinLimits,...
                                   'BinWidth',args.BinWidth,...
                                   'Normalization','Probability',...
                                   'FaceColor',myFaceColor{v});
        hold on;
    end
end

if args.toLegend
hLeg = legend(plt,varNameList,'fontsize',args.fontSize,...
    'box','off','location','northwest');
end

if ~strcmp(args.titStr,'')
    title(args.titStr);
end

ylabel(args.yLab);
xlabel(args.xLab);

set(gca,'fontSize',args.fontSize,'xtick',args.xtick);

end

function fTrack = plotTrack(tab,geoName,xVar,yVarList,ix2fit,ix2plot,...
                                                      fitCoefList,varargin)
myParser = inputParser;
addParameter(myParser,'OUTDIR','.',@ischar);
addParameter(myParser,'figColor','w');
addParameter(myParser,'figSize',[0 0 2.5 2.5]);
addParameter(myParser,'yLab','');
addParameter(myParser,'xLab','');
addParameter(myParser,'xtick',1:1:24);
addParameter(myParser,'titStr','');
addParameter(myParser,'toLegend',1);
addParameter(myParser,'refPoint','dawn');
addParameter(myParser,'markDawn',1);
addParameter(myParser,'adj2data',1);
addParameter(myParser,'fontSize',8);
parse(myParser,varargin{:});
args = myParser.Results;

[geoCode, ix] = geoIx(geoName,tab);

fTrack=figure('units','inches','color',...
    args.figColor,'position',args.figSize);

colors = get(groot,'DefaultAxesColorOrder');

%make plots
for v=1:numel(yVarList)
    x2fit = tab{ix,xVar}(ix2fit);
    y2fit = tab{ix,yVarList{v}}(ix2fit); %assume all ys here are relative to dawn
    
    x2plot = tab{ix,xVar}(ix2plot);
    y2plot = tab{ix,yVarList{v}}(ix2plot); %ys relative to dawn
    
    ix2plotOnly = ~ismember(ix2plot,ix2fit);
    x2plotOnly = tab{ix,xVar}(ix2plotOnly);
    y2plotOnly = tab{ix,yVarList{v}}(ix2plotOnly); %ys relative to dawn
    
    [~,sortXix] = sort(x2plot); %sort day lengths
    [~,sortXixPonly] = sort(x2plotOnly)
    [~,sortXfitIx] = sort(x2fit);
    
    fitCoef = fitCoefList{v}; %pass only fits relative to dawn
    
    switch args.refPoint
        %plot relative to dawn (hrs after dawn)
        case 'dawn'
            plot(x2fit,y2fit,'o',...
                'markersize',3,'markerfacecolor',colors(v,:),...
                'markeredgecolor','none',...
                'linewidth',0.5);
            hold on;
            if ~isempty(ix2plot) & ~isempty(ix2plot)
                  plot(x2plotOnly(sortXixPonly),y2plotOnly(sortXixPonly),...
                      'o',...
                      'color','none',...
                      'markersize',2,'markerfacecolor','none',...
                      'markeredgecolor',colors(v,:),...
                      'linewidth',0.5);
             end
            
            fitYs = polyval([fitCoef(ix,1) fitCoef(ix,2)],x2fit')';
            plot(x2fit(sortXfitIx),fitYs(sortXfitIx),...
                '-','color',colors(v,:),'linewidth',1);
              
        %plot in local time (hrs after midnight)
        case 'midnight'
            plot(x2fit,tab{ix,'dawn'}(ix2fit)+y2fit,'o',...
                'markersize',3,'markerfacecolor',colors(v,:),...
                'markeredgecolor','none',...
                'linewidth',0.5);
            hold on;
            if ~isempty(ix2plot) & ~isempty(ix2plot)
                dawnOffset = tab{ix,'dawn'}(ix2plotOnly);
                plot(x2plotOnly(sortXixPonly),...
                    dawnOffset(sortXixPonly)+y2plotOnly(sortXixPonly),...
                    'o',...
                    'markersize',2,'markerfacecolor','none',...
                    'markeredgecolor',colors(v,:),'linewidth',0.5);
            end
            
            dawnOffset = tab{ix,'dawn'}(ix2fit);
            fitYs = polyval([fitCoef(ix,1) fitCoef(ix,2)],x2fit')';
            
            plot(x2fit(sortXfitIx),...
                dawnOffset(sortXfitIx)+fitYs(sortXfitIx),...
                '-','color',colors(v,:),'linewidth',1);           
    end
end

title(args.titStr);
xlabel(args.xLab);

%adjust ylabel to include time reference point
if strcmp(args.yLab,'')
    switch args.refPoint
        case 'dawn'
            args.yLab = 'time (hours after dawn)';
        case 'midnight'
            args.yLab = 'local time (hours)';
    end
end
ylabel(args.yLab);

%set up ylimits based on reference point
switch args.refPoint
    case 'dawn'
        yLim = [-4 4] + [-0.25 0.25];
    case 'midnight'
        yLim = [2 10] + [-0.25 0.25];
end

set(gca,'fontsize',8,'fontname','helvetica',...
                'xlim',[8 16] + [-0.25 0.25],'ylim',yLim);            
%axis('equal');            

%annotate plot with dawn, dusk, lines, etc.
if args.markDawn
    switch args.refPoint
        case 'dawn'  %mark lines relative to dawn 
        dawn = tab{ix,'dawn'} - tab{ix,'dawn'};
        dusk = tab{ix,'dusk'} - tab{ix,'dawn'};
        NineAM = 9*ones(size(tab{ix,'dawn'}))  - tab{ix,'dawn'};
        case 'midnight' %mark dawn/dusk relative to midnight
        dawn = tab{ix,'dawn'};
        dusk = tab{ix,'dusk'};
        NineAM = 9*ones(size(tab{ix,'dawn'}));
    end
    
    %adj y intercept relative to data?
    if args.adj2data
        %get min y pt of all datapoints we're plotting
        for v=1:numel(yVarList)
            y2plot = tab{ix,yVarList{v}}(ix2plot); %relative to dawn
            
            %correct relative to midnight if have to
            if strcmp(args.refPoint,'midnight')
                y2plot = y2plot + tab{ix,'dawn'};
            end
            if v==1
                minY = min(y2plot);
            else
                minY = min(minY,min(y2plot));
            end
        end
       
        dawn = dawn - min(dawn) + minY;
        dusk = dusk - min(dusk) + minY;
        NineAM = NineAM - min(NineAM) + minY;
    end
    
    %yyaxis right;
    %plot(x2plot(sortXix),dusk(sortXix),'k--','linewidth',1);
    plot(x2plot(sortXix),NineAM(sortXix),'ko-',...
        'markersize',3,...
        'markerfacecolor','k',...
        'markeredgecolor','none','linewidth',1);
    plot(x2plot(sortXix),dawn(sortXix),'k--','linewidth',1);
    
end

            
end

function fTrack = plotTrack9am(tab,geoName,xVar,yVarList,ix2fit,ix2plot,...
                                                      fitCoefList,varargin)
myParser = inputParser;
addParameter(myParser,'OUTDIR','.',@ischar);
addParameter(myParser,'figColor','w');
addParameter(myParser,'figSize',[0 0 2.5 2.5]);
addParameter(myParser,'yLab','');
addParameter(myParser,'xLab','');
addParameter(myParser,'xtick',1:1:24);
addParameter(myParser,'titStr','');
addParameter(myParser,'toLegend',1);
addParameter(myParser,'refPoint','dawn');
addParameter(myParser,'markDawn',1);
addParameter(myParser,'adj2data',1);
addParameter(myParser,'fontSize',8);
parse(myParser,varargin{:});
args = myParser.Results;

[geoCode, ix] = geoIx(geoName,tab);

fTrack=figure('units','inches','color',...
    args.figColor,'position',args.figSize);

colors = get(groot,'DefaultAxesColorOrder');

%make plots
for v=1:numel(yVarList)
    x2fit = tab{ix,xVar}(ix2fit);
    y2fit = tab{ix,yVarList{v}}(ix2fit); %rel to 9 a.m.
    
    x2plot = tab{ix,xVar}(ix2plot);
    y2plot = tab{ix,yVarList{v}}(ix2plot); %rel to 9 a.m.
    
    ix2plotOnly = ~ismember(ix2plot,ix2fit);
    x2plotOnly = tab{ix,xVar}(ix2plotOnly);
    y2plotOnly = tab{ix,yVarList{v}}(ix2plotOnly); %rel to 9 a.m.
    
    [~,sortXix] = sort(x2plot); %sort day lengths
    [~,sortXixPonly] = sort(x2plotOnly);
    [~,sortXfitIx] = sort(x2fit);
    
    fitCoef = fitCoefList{v}; %pass only fits relative to 9am!
    
    switch args.refPoint
        %plot relative to dawn (hrs after dawn)
        case 'dawn'
            dawnOffset = tab{ix,'dawn'}(ix2fit); %this is relative to midnight (local time)
            
            plot(x2fit,y2fit + (9 - dawnOffset),'o',...
                'markersize',3,'markerfacecolor',colors(v,:),...
                'markeredgecolor','none',...
                'linewidth',0.5);
            hold on;
            if ~isempty(ix2plot) & ~isempty(ix2plot)
                  dawnOffset = tab{ix,'dawn'}(ix2plotOnly);
                  plot(x2plotOnly(sortXixPonly),...
                      y2plotOnly(sortXixPonly) + (9 - dawnOffset(sortXixPonly)),...
                      'o',...
                      'color','none',...
                      'markersize',2,'markerfacecolor','none',...
                      'markeredgecolor',colors(v,:),...
                      'linewidth',0.5);
             end
            
            fitYs = polyval([fitCoef(ix,1) fitCoef(ix,2)],x2fit')';
            dawnOffset = tab{ix,'dawn'}(ix2fit);
            plot(x2fit(sortXfitIx),...
                fitYs(sortXfitIx) + (9 - dawnOffset(sortXfitIx)),...
                '-','color',colors(v,:),'linewidth',1);
              
        %plot in local time (hrs after midnight)
        case 'midnight'
            NineAMoffset = 9;
            plot(x2fit,NineAMoffset+y2fit,'o',...
                'markersize',3,'markerfacecolor',colors(v,:),...
                'markeredgecolor','none',...
                'linewidth',0.1);

            hold on;
            if ~isempty(ix2plot) & ~isempty(ix2plot)
                plot(x2plotOnly(sortXixPonly),...
                    NineAMoffset+y2plotOnly(sortXixPonly),...
                    'o',...
                    'markersize',3,'markerfacecolor','none',...
                    'markeredgecolor',colors(v,:),...
                    'linewidth',0.01);
            end
            
            NineAMoffset = 9;
            fitYs = polyval([fitCoef(ix,1) fitCoef(ix,2)],x2fit')';
            
            plot(x2fit(sortXfitIx),...
                NineAMoffset+fitYs(sortXfitIx),...
                '-','color',colors(v,:),'linewidth',1);           
    end
end

title(args.titStr);
xlabel(args.xLab);

%adjust ylabel to include time reference point
if strcmp(args.yLab,'')
    switch args.refPoint
        case 'dawn'
            args.yLab = 'time (hours after dawn)';
        case 'midnight'
            args.yLab = 'local time (hours)';
    end
end
ylabel(args.yLab);

switch args.refPoint
    case 'dawn'
        yLim = [-4 4] + [-0.25 0.25];
    case 'midnight'
        yLim = [2 10] + [-0.25 0.25];
end

set(gca,'fontsize',8,'fontname','helvetica',...
                'xlim',[8 16] + [-0.25 0.25],'ylim',yLim);            
%axis('equal');            

%annotate plot with dawn, dusk, lines, etc.
if args.markDawn
    switch args.refPoint
        case 'dawn'  %mark lines relative to dawn
            dawn = tab{ix,'dawn'} - tab{ix,'dawn'};
            dusk = tab{ix,'dusk'} - tab{ix,'dawn'};
            NineAM = 9*ones(size(tab{ix,'dawn'}))  - tab{ix,'dawn'};
        case 'midnight' %mark dawn/dusk relative to midnight
            dawn = tab{ix,'dawn'};
            dusk = tab{ix,'dusk'};
            NineAM = 9*ones(size(tab{ix,'dawn'}));
    end
    
    %adj y intercept relative to data?
    if args.adj2data
        %get min y pt of all datapoints we're plotting
        for v=1:numel(yVarList)
            y2plot = tab{ix,yVarList{v}}(ix2plot); %relative to dawn
            
            %correct relative to midnight if have to
            if strcmp(args.refPoint,'midnight')
                y2plot = y2plot + tab{ix,'dawn'};
            end
            if v==1
                minY = min(y2plot);
            else
                minY = min(minY,min(y2plot));
            end
        end
       
        dawn = dawn - min(dawn) + minY;
        dusk = dusk - min(dusk) + minY;
        NineAM = NineAM - min(NineAM) + minY;
    end
    
    %yyaxis right;
    %plot(x2plot(sortXix),dusk(sortXix),'k--','linewidth',1);
    plot(x2plot(sortXix),dawn(sortXix),'ko-',...
        'markersize',3,...
        'markerfacecolor','k',...
        'markeredgecolor','none','linewidth',1);
    plot(x2plot(sortXix),NineAM(sortXix),'k--','linewidth',1);
    
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