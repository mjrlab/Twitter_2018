%2018-08-10, EL: modified to export Fig. 1 data
%2018-06-27, EL: Fig. 1C,E

clear all;
close all;
clc;

%export figures?
TOEXP = 0; %(1 = yes)
TOEXP_TrWidth_DATA = 0; %export supporting data for this figure?

%geoLevel = 'county'; 

%% for county-level data
INDIR = '../tweetograms/tweetograms_annual_FIPS_2012_2013_sfyo';
whichDay = 'annual/MTWThFSatSun';
sjlTab = load([INDIR '/' 'troughs_2018-04-09_01.13.04_sfyo.mat']);
sjlTab = sjlTab.outTable;

% rename 'MTWThFSatSun' -> 'allweek'
varNames = sjlTab.Properties.VariableNames;
varNames = cellfun(@(x) strrep(x,'MTWThFSatSun','allweek'),varNames,'UniformOutput',0);
sjlTab.Properties.VariableNames = varNames;

%% county health rankings: % getting < 7 hrs of sleep
sleepTab_CHR = load('datafiles/sleepTable_CHR_fracSleepUnder7.mat');
sleepTab_CHR = sleepTab_CHR.sleepTab_CHR;
sleepTab_CHR.Properties.VariableNames(1:2) = {'geoCode','CHR_fracSleepUnder7'};

%% census data
census = load('datafiles/AllCovariates_Apr15.mat');
census = census.AllCovariates;
census = census(:,{'FIPS','Population1','Population2','Population5'});

%county names
nameTab = load('datafiles/countyNames.mat');
nameTab = nameTab.countyNames;

%% merge sjl table with other data
sjlTab = outerjoin(sjlTab,sleepTab_CHR,'Keys','geoCode','MergeKeys',true);
sjlTab = outerjoin(sjlTab,nameTab,'LeftKeys','geoCode','RightKeys','FIPS','MergeKeys',true);
sjlTab.Properties.VariableNames{1} = 'geoCode';
sjlTab = outerjoin(census,sjlTab,'LeftKeys','FIPS','RightKeys','geoCode');

%% export trough width data
if TOEXP_TrWidth_DATA
    TRWIDTHtable_exp =  sjlTab(~isnan(sjlTab.('allweek_trWidth')),:);
    TRWIDTHtable_exp = TRWIDTHtable_exp(:,{'FIPS','countyName',...
                                   'stateName','allweek_trWidth'});
  
    TRWIDTHtable_exp.Properties.VariableNames = {'FIPS','name','state',...
        'troughWidth_hr'};
    writetable(TRWIDTHtable_exp,['Fig1_trWidth_calcs_' getDate() '.xlsx']);
end

%% Fig. 1E: plot trough width vs CDC data on sufficient sleep
fOne = figure('units','inches','position',[0 0 2.05 2.05]);
titFile='trWidth_fracUnder7_CHR_newData';

%select table variables and give them names
xVar='CHR_fracSleepUnder7';
xlab='sufficient sleep';

yVar='allweek_trWidth';
ylab='Tweetogram trough width (hr)';

%remove rows that don't have data
pltTab = sjlTab(~isnan(sjlTab.(xVar)) & ~isnan(sjlTab.(yVar)),:);

%prevalence sufficient sleep = 1 - prevalence insufficient sleep
pltTab.(xVar) = 1-pltTab.(xVar); 
            
%histogram parameters
numXDivs = 20;
loX = 0.5;
hiX = 0.8;
xTick = [0.5:0.1:0.8];

numYDivs = 20;
loY = 2;
hiY = 7;
yTick = 2:7;
          
xyhistogram(pltTab,xVar,yVar,...
                 'loX',loX,'hiX',hiX,...
                 'numXDivs',numXDivs,'xTick',xTick,...
                 'loY',loY,'hiY',hiY,...
                 'numYDivs',numYDivs,'yTick',yTick,...
                 'xLab',xlab,...
                 'ylab',ylab,...
                 'hFig',fOne,...
                 'fontsize',8,...
                 'toColorBar',0,...
                 'toLegend',1);

expif(TOEXP,[fOne],titFile,'OUTDIR','sleepPics','fun2use','expfig');

%display correlation coefficient too
[rho,pval] = corr(pltTab.(xVar),pltTab.(yVar),'type','Pearson','rows','complete')

%% Fig. 1C: map of US colored by trough width 
mapVar = 'allweek_trWidth';
sjlTab_forMap = sjlTab(:,{'geoCode','countyName','stateName',mapVar});

%make a 'name' column by concacatenating countyName and stateName
sjlTab_forMap(:,'name') = rowfun(@(x,y) [strcat(x,y)], ...
    sjlTab_forMap(:,{'countyName','stateName'}));
sjlTab_forMap(isnan(sjlTab_forMap.(mapVar)),:) = []; %remove rows with missing values for trough width

%need to remove spaces from state names
st(:,{'stateName'}) = rowfun(@(x) strrep(x,' ',''),sjlTab_forMap(:,{'stateName'}));
sjlTab_forMap.stateName = st.stateName;
sjlTab_forMap.Properties.VariableNames = {'FIPS','countyName','state','trWidth','name'};

%make the figure
[hFigUSA, shapeStruct] = ...
    plotMap_Polygon_label(sjlTab_forMap, 'FIPS', 'trWidth','',[]);
set(hFigUSA,'units','inches','position',[0 0 5 3.5],'color','none');

%save it
expif(TOEXP,hFigUSA,'map_avgWidth_1521counties_trim1to99',...
    'OUTDIR','map','fun2use','expfig','ext','pdf','m','-m1'); 

%% Fig. 1D inset: plot empty map showing which counties are used in examples
cList = {'Suffolk NY','Lafayette LA','Orange CA','Wayne MI'};
fipsList = [];
for c=1:numel(cList)
    [geoCode, ix(c)] = geoIx(cList{c},sjlTab_forMap);
end
%subselect small tab
small_tab = sjlTab_forMap(ix,:);
%set same color limits as for big map
cLim = [floor(prctile(small_tab.trWidth,1)) ceil(prctile(small_tab.trWidth, 99))];
[hFigUSA_exCounties, shapeStruct] = ...
    plotMap_Polygon_label(small_tab, 'FIPS', 'trWidth','',cLim);
set(hFigUSA_exCounties,'units','inches','position',[0 0 5 3.5],'color','none');
expif(TOEXP,hFigUSA_exCounties,'map_exCounties',...
    'OUTDIR','map','fun2use','expfig','ext','pdf','m','-m1'); 

%% for talks: plot map with all counties in our dataset marked yellow
mapVar = 'rand';
sjlTab_forMap = sjlTab_forMap(:,{'FIPS','state'});

%need to remove spaces from state names
st = rowfun(@(x) strrep(x,' ',''),sjlTab_forMap(:,{'state'}));
sjlTab_forMap.state = st.Var1;

%set all colors to 1
sjlTab_forMap.(mapVar) = ones(size(sjlTab_forMap.FIPS));
sjlTab_forMap.Properties.VariableNames = {'FIPS','state',mapVar};

[hFigUSA_all, shapeStruct] = ...
    plotMap_Polygon_label(sjlTab_forMap, 'FIPS', mapVar,'',[0 1.5]);
set(hFigUSA_all,'units','inches','position',[0 0 5 3.5],'color','none');
colorbar off;
expif(TOEXP,hFigUSA_all,'map_allCounties_oneColor',...
    'OUTDIR','map','fun2use','expfig','ext','pdf','m','-m1');  

%% helper funs
function xyhistogram(myTab,xVar,yVar,varargin)
%plot 2D histogram, allow to pass a bunch of parameters

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

function [geoCode, ix] = geoIx(countyName,cTable)
if ischar(countyName)
    geoCode = cTable.FIPS(cellfun(@(x) isequal(x,countyName), cTable.name));
    ix = find(cellfun(@(x) isequal(x,countyName), cTable.name));
else
    geoCode = nan;
    ix = nan;
end
end

function [hFigUSA, shapeStruct] = ...
    plotMap_Polygon_label(ds, geoUnit, colorBy,mapTitStr,colorLims)
%Plot map of US and color all geoUnits according to the value of colorBy
%column in table ds.
%Inputs:
%   - ds = Matlab table or dataset containing columns 'FIPS','state', 
%          and another column colorBy containing values which will determine
%          the colors (colorBy could be any string) 
%   - colorBy = name of column containing values
%               to be converted to colors (e.g., 'trWidth')
%   - geoUnit = use 'FIPS' (could also be 'MSA' or 'CBSA'); 
%               determines which geographic entities you are plotting and
%               is used to load the appropriate shapeFile (we're focused on
%               individual counties here, so use 'FIPS')
%   - mapTitStr = set to '', currently unused; however, if the line at the
%               end of this function beginning with textm() is uncommented, 
%               the string mapTitStr will be used to title the map
%   - colorLims = vector [minVal maxVal] or empty vector [] determining the
%               max and min values of the colormap used to plot. If you
%               pass an empty vector [], the color limits will be scaled to
%               the 1st and 99th percentile of the data.
%
%Outputs:
%   - hFigUSA = handle to figure with the map; see plotUSA.m for the three
%   sets of axes defining the figure
%   - shapeStruct = shape file with geograhic coordinates of each polygon
%   used for plotting, its name and color value.

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
%load matlab map of US states
whoseMap = 'Matlab';
load('datafiles/mapUSA_shape_states_Matlab.mat');

%fill colors
dummy = num2cell(nan(numel(shapeStruct),1));
[shapeStruct(:).color] = deal(dummy{:});
for i=1:numel(shapeStruct)
    %find the row of the table corresponding to this geographic region
    ix = find(ds.(geoUnit) == shapeStruct(i).(geoUnit));
    if ~isempty(ix)
        %assign the color of this entry in the geographic shape structure
        %(shapeStruct(i).color) to the value of the mapping variable (ds{ix,colorBy})
        shapeStruct(i).color = ds{ix,colorBy};
    end
end

% set colormap limits to [1st 99th] percentile of the data
if isempty(colorLims)
    colorLims = [floor(prctile(ds.(colorBy),1)) ...
                 ceil(prctile(ds.(colorBy), 99))]; %in hours
        
    %2018-06-03, EL:
    %should manually trim the values to the [1,99] percentiles. it seems like
    %Matlab assigns them to the median value otherwise...
    ds{ds.(colorBy) > colorLims(2),colorBy} = colorLims(2);
    ds{ds.(colorBy) < colorLims(1),colorBy} = colorLims(1);
end

%make symbolspec, Matlab's way of assigning colors to geographical regions
mySymbolSpec = makesymbolspec('Polygon',...
    {'color',colorLims,...  
    'FaceColor',jet(512)}); %jet by default

%get indices of Hawaii, Alaska and continental states in your dataset
S_ix_Hawaii = strcmp(ds.state,'HI');
S_ix_Alaska = strcmp(ds.state,'AK');
S_ix_ConUS = ~(S_ix_Hawaii | S_ix_Alaska);

hFigUSA = figure();
set(hFigUSA,'position',[2 2 11 8.5]);
[~, hAxUSA] = plotUSA(hFigUSA);

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

%add colorbar
axes(hAxUSA(1));
caxis manual;
caxis(colorLims);
colormap(jet(512));
hBar = colorbar;
set(hBar,'fontsize',8);
%textm(51,-100,mapTitStr,'color','k');
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