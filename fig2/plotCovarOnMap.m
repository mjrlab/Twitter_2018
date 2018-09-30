%2018-06-27, EL: plot covariates of interest on a map of US

clear all; close all; clc;

%% input parameters
TOEXP = 1; %export figures? (1=yes)

%% list of covariates you want to plot (each on its own map)
covarList = {'GeoMeteo3'}; 
        %{'GeoMeteo1','GeoMeteo2','Population5','Age1',...
        %    'Edu8','Commute13','Health11','Health16'};

%% load sjl calculations
SJLtable = load(['../tweetograms/tweetograms_monthly_FIPS_2012_2013_dzah/'...
                   'sjl_troughs_2018-04-08_13.47.16_dzah.mat']);
sjlTab = SJLtable.outTable;
sjlTab.Properties.VariableNames{'geoCode'} = 'FIPS';

%% toss observations > 3 sigma -- all correspond to counties with low counts
sjlTab = tossNoisy(sjlTab,'xcorr_wkdaywkend');

%% load covariates
AllCovariates = load('datafiles/AllCovariates_Apr15.mat'); %loads into AllCovariates_decorrelated
AllCovariates = AllCovariates.AllCovariates; %.AllCovariates_decorrelated;
AllCovariates.Properties.VariableNames{1} = 'FIPS'; %make sure name of FIPS column is right

%% load all FIPS
geoTab = load('datafiles/CountyTableAll.mat');
geoTab = geoTab.cTable;

%% plot list of covariates on maps
geoTab = innerjoin(geoTab,AllCovariates(:,{'FIPS',covarList{:}}),...
    'LeftKeys','FIPS','RightKeys','FIPS');
geoTab = innerjoin(geoTab,sjlTab(:,{'FIPS'}),'Keys','FIPS');
%%
for cl=1:numel(covarList)
    mapTab = geoTab(~isnan(geoTab.(covarList{cl})),:);        
    [mapFig(cl), ~] = ...
        plotMap_Polygon_label(mapTab, 'FIPS', covarList{cl},'',[]);
    set(mapFig(cl),'color','w','position',[0 0 3.5 2.25]);
    expif(TOEXP,mapFig(cl),['map_' covarList{cl}],...
    'OUTDIR','./mapLinModCovars','fun2use','expfig','ext','pdf');  
end

%% helper funs
% utilities used to make figures & export them
function expif(TOEXP,fH,fName,varargin)
myParser = inputParser;
addParameter(myParser,'OUTDIR','.',@ischar);
addParameter(myParser,'ext','pdf',@ischar);
addParameter(myParser,'fun2use','expfig',@ischar);
addParameter(myParser,'nocrop','nocrop',@ischar);
addParameter(myParser,'m','-m1',@ischar);
addParameter(myParser,'renderer','-painters');
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
                        args.renderer,'-cmyk','-nocrop',args.m,fH);
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

%data cleaning, centering
function cleanTab = tossNoisy(SJLtab,responseVar)

mu = mean(SJLtab{:,responseVar}(:));
sig = std(SJLtab{:,responseVar}(:));

noisyRows = rowfun(@(x) threeSig(x,mu,sig), SJLtab(:,responseVar));
cleanTab = SJLtab(~noisyRows{:,1},:);

    function isbad = threeSig(r,mu,sig)
        %are any elements in this row more than 3 sigma away from the mean?
        %deal with one row
        badIx = (r - mu)./sig > 3 | (r - mu)./sig < -3;
        isbad = sum(badIx) > 0;
    end
% 
%     function isbad = largeStdev(r,mu,sig)
%         %does this row have a high stdev?
%         %deal with one row
%         isbad = std(r) > 40/60; %mean of distribution of these stdevs is 0.5+/-0.25
%         %isbad = isbad | min(r) < 0;
%     end
end

%map functions
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

% assign values to outliers after setting color limits!
if isempty(colorLims)
    colorLims = prctile(ds.(colorBy),[1 99]);    
%     colorLims = [floor(prctile(ds.(colorBy),1)) ...
%         ceil(prctile(ds.(colorBy), 99))];
%    colorLims = [min(ds.(colorBy)) max(ds.(colorBy))];
    
    %should manually trim the values to the [1,99] percentiles. it seems like
    %Matlab assigns them to the median value otherwise...
    ds{ds.(colorBy) > colorLims(2),colorBy} = colorLims(2);
    ds{ds.(colorBy) < colorLims(1),colorBy} = colorLims(1);
end

%fill colors
dummy = num2cell(nan(numel(shapeStruct),1));
[shapeStruct(:).color] = deal(dummy{:});
for i=1:numel(shapeStruct)
    ix = find(ds.(geoUnit) == shapeStruct(i).(geoUnit));
    if ~isempty(ix)
        %Had a bug here before!
        %Was selecting shapeStruct(ix)! instead of shapeStruct(i)!!!
        shapeStruct(i).colorVar = ds{ix,colorBy};
        
        %test by assigning random colors
        %shapeStruct(i).color = rand*100;
    end
end

mySymbolSpec = makesymbolspec('Polygon',...
    {'colorVar',colorLims,... %trying this instead of cLims
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

%make default colorbar skinnier
set(hBar,'fontsize',8);
%barTicks = 0:0.25:12;
set(hBar,'Limits',colorLims);
%set(hBar,'Limits',[min(ds.(colorBy)) max(ds.(colorBy))]);
% hSize = get(hBar,'Position');
% set(hBar,'position',[hSize(1) hSize(2) hSize(3)/2 hSize(4)]);
%textm(51,-100,mapTitStr,'color','k');
end

function drawFIPS(myds, shapeStruct, geoUnit, faceColor, edgeColor)
%draw FIPS geoboundaries in myds.FIPS from shapeFile loaded into 
%shapeStruct; color all FIPS according to myColor
if numel(myds{:,geoUnit} > 1)
    [myFIPSix, myShapeStruct] = getFIPSix(myds.(geoUnit), shapeStruct, geoUnit);

    geoshow(myShapeStruct,'DisplayType','Polygon',...
        'FaceColor',faceColor,'edgecolor','none',...%edgeColor,...
        'linewidth',0.01); %sually 0.01 --> too small in TIFs, make this 0.05
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