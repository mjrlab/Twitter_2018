%2018-06-27, EL: Code to generate Figure 4 and Figure S3
%   - clustering analysis + socioeconomic enrichment in clusters
                 
%% input parameters
clear all; close all; clc;

%which analyses would you like to run? export which figures?
%(1=yes)

% Fig. 4A, C
TOPLOT_CLUSTERS = 1; 
TOEXP_CLUSTERS = 0; %export figure?
TOEVAL_CLUSTERS = 0; %evaluate clusters using the gap statisic?

%Fig. 4B 
TOMAP = 0; % 
TOEXP_MAP = 0; %export figure? (1=yes)

% Fig. 4D
TORUN_ANOVAS = 1; 
TOEXP_BARS_MAIN = 1; %export figure? (1 = yes)

% Fig. S3
TOMAKEBARS = 0; 
TOEXP_BARS_SUP = 0; %export figure?

% sup. file with cluster assignments and activity counts?
EXP_ACTIVITY_TABLE = 0;

%% load sjl calculations
SJLtable = load(['../tweetograms/tweetograms_monthly_FIPS_2012_2013_hthz/'...
                'sjl_troughs_2018-04-08_17.32.34_hthz.mat']);
SJLtable = SJLtable.outTable;
SJLtable.tr_sjl = SJLtable.wkend_trPos - SJLtable.wkday_trPos;
SJLtable.Properties.VariableNames{1} = 'FIPS'; %geoCode = FIPS here for all 
SJLtable.Properties.VariableNames{'xcorr_wkdaywkend'} = 'xc_sjl';

%% load geo table
geoTab = load('datafiles/CountyTableHaveData.mat');
geoTab = geoTab.cTable;
SJLtable = innerjoin(geoTab,SJLtable,...
    'LeftKeys','geoCode','RightKeys','FIPS');

%% load data on number of users 
popTab = load('datafiles/AllCovariates_Apr15.mat'); %loads into AllCovariates
popTab = popTab.AllCovariates(:,{'FIPS','Population1','Population5'}); 
popTab.avgUsers = (popTab.Population5.*(10.^popTab.Population1));
popTab.Properties.VariableNames{1} = 'geoCode'; %make sure name of FIPS column is right
SJLtable = innerjoin(SJLtable,popTab(:,{'geoCode','avgUsers'}),'Keys','geoCode');

%% run quality checks
[SJLtable] = qcSJL(SJLtable);

%% remove counties with low numbers of users
AVGUSERCUTOFF = 30; %avg. number of daily users
NOISECUTOFF = 1;    %
SJLtable = SJLtable(SJLtable.avgUsers > AVGUSERCUTOFF & ...
                    SJLtable.sigSJL_xc_z < NOISECUTOFF,:);

%pick responseVariable and rename it to 'sjl'
responseVar = 'xc_sjl'; %or 'tr_sjl'
labelVar = responseVar;
if ~strcmp(responseVar,'sjl')
    SJLtable.Properties.VariableNames{responseVar} = 'sjl';
    responseVar = 'sjl';
end
      
%% subtract mean level of fluctuation
SJLtable(:,responseVar) = rowfun(@(x) x - mean(x),SJLtable(:,responseVar));
SJLtable = SJLtable(:,{'geoCode',responseVar});

%% load geo table
geoTab = load('datafiles/CountyTableHaveData.mat');
geoTab = geoTab.cTable;
geoTab = innerjoin(geoTab,SJLtable(:,{'geoCode',responseVar}),...
    'LeftKeys','geoCode','RightKeys','geoCode');
geoTab.Properties.VariableNames{1} = 'FIPS';
SJLtable.Properties.VariableNames{'geoCode'} = 'FIPS';

%% load covariates
AllCovariates = load('datafiles/allCovariates_Apr15.mat'); %loads into AllCovariates_decorrelated
AllCovariates = AllCovariates.AllCovariates; %.AllCovariates_decorrelated;
AllCovariates.Properties.VariableNames{1} = 'FIPS'; %make sure name of FIPS column is right

%% cluster SJL data
clustBy = responseVar;
SJLtable_clust = [];

[SJLtable_clust,ix] = clusterDataBy(SJLtable,clustBy,'TOORDER',0,...
    'MAXCLUST',3);
colors = hsv(max(SJLtable_clust.hierClust)); %color-code clusters

%use Gap statistic to evaluate optimal number of clusters
if TOEVAL_CLUSTERS
    eva = evalclusters(SJLtable.(clustBy),'linkage','Gap',...
        'KList',[1:10],'Distance','Euclidean')
    figure();
    plot(eva);
    xlabel('num clusters');
    ylabel('gap statistic');   
end

%make table with monthly numbers of user-days per county, as well as
%cluster assigment
if EXP_ACTIVITY_TABLE
    wkdayAct = load('datafiles/wkday_activities.csv');
    wkdayAct = table(wkdayAct(:,1),wkdayAct(:,2:end),...
        'VariableNames',{'FIPS','wkdayUserDays'});
    wkendAct = load('datafiles/wkend_activities.csv');
    wkendAct = table(wkendAct(:,1),wkendAct(:,2:end),...
        'VariableNames',{'FIPS','wkendUserDays'});
    saveTab = geoTab(:,{'FIPS','geoName','state'});
    saveTab = join(saveTab, SJLtable_clust(:,{'FIPS','hierClust'}));
    saveTab = join(saveTab,wkendAct,'Keys','FIPS');
    saveTab = join(saveTab,wkdayAct,'Keys','FIPS');
    writetable(saveTab,...
        ['clustAssign/clusters_activities_' getDate() '.xlsx']);
end

%% Figure 4A,C: plot clusters as a heatmap and also cluster avg +/- stdev
if TOPLOT_CLUSTERS
fClustHeatMap = plotHierClust(SJLtable_clust,clustBy,colors,'hierClust');
fClustAvg = plotClustAvgs(SJLtable_clust,clustBy,colors);
%%
for fc=1:numel(fClustAvg)
   set(fClustAvg(fc),'units','inches','position',[0 0 3.5 1.8]);
end
expif(TOEXP_CLUSTERS,fClustAvg,['hierClust_avg_' clustBy '_'  ...
    num2str(AVGUSERCUTOFF) 'avgUserCutoff'],...
    'OUTDIR','fig4_clust','fun2use','expfig','ext','pdf');

expif(TOEXP_CLUSTERS,fClustHeatMap,['hierClust_heatmap_' clustBy '_' ...
     num2str(AVGUSERCUTOFF) 'avgUserCutoff'],...
    'OUTDIR','fig4_clust','fun2use','print','ext','pdf');
end

%% descriptive statistics of seasonal SJL fluctuations in different clusters
% these numbers are quoted in main text
SJLtable_clust.sigSJL = 60*std(SJLtable_clust.sjl,0,2);
SJLtable_clust.maxSJL = 60*max(SJLtable_clust.sjl,[],2);
SJLtable_clust.avgSJL = 60*mean(SJLtable_clust.sjl,2);
SJLtable_clust.admSJL = 60*mean(abs(SJLtable_clust.sjl),2);
SJLtable_clust.amplSJL = 60*max(SJLtable_clust.sjl,[],2) - 60*min(SJLtable_clust.sjl,[],2);

varfun(@mean,SJLtable_clust(:,{'sigSJL','hierClust'}),...
    'GroupingVariables','hierClust');
varfun(@mean,SJLtable_clust(:,{'avgSJL','hierClust'}),...
    'GroupingVariables','hierClust');
varfun(@mean,SJLtable_clust(:,{'admSJL','hierClust'}),...
    'GroupingVariables','hierClust')
varfun(@mean,SJLtable_clust(:,{'amplSJL','hierClust'}),...
    'GroupingVariables','hierClust')

%% plot on a map -- Fig. 4B
%some geographical info from countyTable
if TOMAP 
hFigUSA = plotMap_Polygon_allClust(...
    innerjoin(geoTab,SJLtable_clust,'Keys','FIPS'),'FIPS','hierClust',...
    colors,'hierClust');
set(hFigUSA,'units','inches','position',[0 0 4 3],'color','none');
expif(TOEXP_MAP,hFigUSA,['map_hierClust' '_' ...
                         num2str(AVGUSERCUTOFF) 'avgUserCutoff'],...
    'OUTDIR','fig4_map','fun2use','expfig','ext','pdf','m','-m1');  
end

%% make table of responseVar + covariates
%myTab = [covar1,covar2,...,covarN,responseVar]
[myTab_centered,myTab_raw,y_raw,X_raw,FIPS] = ...
    makeMyTable(SJLtable_clust(:,{'FIPS','hierClust'}),...
                AllCovariates,'exclCovar',{'zzz'});
            
%remove rows that have NaN predictors
nanRows = rowfun(@(x) sum(isnan(x)) > 0,myTab_raw,'SeparateInputs',false);
myTab_raw(nanRows.Var1,:) = [];           

%% make Fig. S3, bar chars showing enrichment of socioeconomic covariates in clusters
if TOMAKEBARS == 1
%each group of covariates is plotted on one figure; 
%covariates in the same group all begin with the same text string 
%(e.g., CommuteN and CommuteM all belong to the group Commute)
groups = {'Commute','Edu','Population',...
          'Sex','Age','Race','Household',...
          'Employed','Occupation','Industry','WorkerClass',...
          'Income','GeoMeteo','Health','Religion','Politics'};

%get bonferroni corrected pvalue level
numVars = size(myTab_raw,2);
numClust = max(unique(myTab_raw.hierClust));
bonfCorrPval = 0.05/(numVars*numClust);

%make figures
for g=1:numel(groups)
    figList = [];
    [pVals, pSign, figList, hBar] = ...
        anovaBars(myTab_raw,'hierClust',groups{g},...
        'lineLevel',-log10(bonfCorrPval),...
        'figColor','w');
    expif(TOEXP_BARS_SUP,figList,['anovaP_' groups{g}],...
        'OUTDIR','./fig4_sup_socioecs_in_clusters',...
        'fun2use','expfig','ext','pdf');
end

end

%% Figure 4D -- look at cluster enrichment for manually-selected covariates 
%(same covariates as in Fig. 2)
if TORUN_ANOVAS
    %load covariate names
    handPickList = readtable('datafiles/linModelDict.xlsx','Sheet','May18_linModel');
%     handPickList = readtable('datafiles/linModelDict_v2.xlsx','Sheet','Aug5_linModel');
%     handPickList = handPickList.DescriptiveName;
    
    % join into table with cluster assignment
    avaTab = [myTab_raw(:,handPickList) myTab_raw(:,{'hierClust'})];
    
    % assign colors to clusters as above
    colors = hsv(max(avaTab.hierClust)); 
    
    % comparison type
    compareType = 'kw'; %'anova' for ANOVA; 'kw' for Kruskall-Wallis
    
    %run comparisons and plot the results
    [allAnovas,resTab] = doAnova(avaTab,'hierClust',compareType);
    
    % make Fig. 4D
    fAnova = myHistStackAnovas(avaTab,'hierClust',allAnovas,'maxFig',25);
    expif(TOEXP_BARS_MAIN,fAnova,['hierClust' '_' compareType '_by_' clustBy],...
        'OUTDIR','fig4_exampleDistr','fun2use','expfig');
    
    %save p values in Fig. 4D to text file
    if TOEXP_BARS_MAIN
        writetable(resTab,['fig4_exampleDistr/pValueTable' '_' 'hierClust' '.xlsx']);
    end
    
end

%% helper functions

%cliuster enrichment
function fAnova = myHistStackAnovas(myTab,anovaVar,allAnovas,varargin)
%make histogram plots of covariates that differ the most between clusters
%allAnovas = struct that's the output of doAnovas() fun

%parse args
myParser = inputParser;
addParameter(myParser,'toPlt',0,@isnumeric); %how covariates differ between clusters?
addParameter(myParser,'toExp',0,@isnumeric); %export figs?
addParameter(myParser,'OUTDIR','./pics',@ischar); %export where?
addParameter(myParser,'multiHypTest','Bonferroni',@ischar);
addParameter(myParser,'maxFDR',0.05,@isnumeric);
addParameter(myParser,'maxFig',10,@isnumeric);
addParameter(myParser,'colors',[1 0 0],@isnumeric);
parse(myParser,varargin{:});
args = myParser.Results;

%correct p-vals for FDR
switch args.multiHypTest
    case 'FDR'
        [FDR, Q, Pi0] = mafdr(allAnovas.p);
        diff_ix = find(Q < args.maxFDR);
    case 'Bonferroni'
        diff_ix = allAnovas.p < args.maxFDR/numel(allAnovas.p);        
end

%how many significant test are there?
numSig = sum(diff_ix);
disp(['no. significant comparisons: ' num2str(numSig)]);
[~,sort_ix] = sort(allAnovas.p,'ascend');

%how many groups to plot?
numClust = numel(unique(myTab{:,anovaVar}));

%if only two clusters (e.g., clust1 vs not-clust1)
if numClust == 2    
    fAnova=figure();
    for i=1:min(numSig,max(10,args.maxFig))
        hAx = subplot(1,10,i);
        y1ix = myTab{:,anovaVar} == 1;
        y2ix = ~y1ix;
  
        y1 = myTab{y1ix,sort_ix(i)};
        y2 = myTab{y2ix,sort_ix(i)};
        plotLRhist(fAnova,hAx, y1, y2, ...
            args.colors(1,:), args.colors(1,:),1,'k','k'); %1=plot median, 0=mean
        set(hAx,'xtick',[],'xticklabel','none');
        set(hAx,'FontSize',6);
        title(hAx,...
            [ezCovarName(myTab.Properties.VariableNames{sort_ix(i)}) 10 ...
            'p=' num2str(allAnovas.p(sort_ix(i)),'%1.1g')],...
            'fontsize',4,'VerticalAlignment','bottom');
    end
    set(fAnova,'units','inches','position',[0 0 9 3]);

%if multiple clusters (e.g., clust1-4)
else
    for i=1:min(numSig,max(8,args.maxFig))
        data = [];
        for c=1:numClust
            ix = myTab{:,anovaVar} == c;
            data{c} = myTab{ix,sort_ix(i)};
            clustName{c} = c;
        end
        varLabel = myTab.Properties.VariableNames{sort_ix(i)};
        fAnova(i) = myHistStack(data,...
            'colors',hsv(numel(data)),...
            'varLabel',varLabel,...
            'titStr',ezCovarName(varLabel),...
            'dataLabels',clustName,...
                'groupLabel','cluster',...
                'toLegend',1,'pVal',allAnovas.p(sort_ix(i)));
    end
end

end

function [allAnovas,resTab] = doAnova(myTab,groupVar,compareType)
% go through each column and do anova
% compareType = 'kw' for Kruskal-Wallis or 'anova' for ANOVA

for i=1:numel(myTab(1,:))
    if ~strcmp(myTab.Properties.VariableNames{i},groupVar)
        switch compareType
            case 'anova'
                [p,tbl,stats] = ...
                    anova1(myTab{:,i},myTab.(groupVar),'off'); %no display
                %indicator (+1 or -1) of whether value in cluster is bigger
                %than in other FIPS
                rowIx = myTab.(groupVar);
                inClust = mean(myTab{rowIx,i});
                inOther = mean(myTab{~rowIx,i});
                anovaSign = ((inClust > inOther) - 0.5)*2; 
            case 'kw'
                [p,tbl,stats] = ...
                    kruskalwallis(myTab{:,i},myTab.(groupVar),'off'); %no display
                %indicator (+1 or -1) of whether cluster median is bigger
                %than things outside the cluster
                rowIx = myTab.(groupVar);
                inClust = median(myTab{rowIx,i});
                inOther = median(myTab{~rowIx,i});
                if inClust == inOther %if medians are equal, use means
                    %anovaSign = 0;
                    inClust = mean(myTab{rowIx,i});
                    inOther = mean(myTab{~rowIx,i});
                end
                if inClust == inOther
                    anovaSign = 0;
                else
                    anovaSign = ((inClust > inOther) - 0.5)*2;
                end
                
        end
        
        allAnovas.p(i) = p;
        allAnovas.tbl{i} = tbl;
        allAnovas.stats{i} = stats;
        allAnovas.ix(i) = i;
        allAnovas.anovaSign(i) = anovaSign;
    end   
end

%turn into a table to write to file later
resTab = table;
resTab(:,'covariate') = cell2table(myTab.Properties.VariableNames(allAnovas.ix)');
resTab(:,'covariateDescr') = cell2table(cellfun(@(x) ezCovarName(x),resTab{:,'covariate'},'UniformOutput',false));
resTab.pVal = allAnovas.p(allAnovas.ix)';
resTab.ix = allAnovas.p(allAnovas.ix)';
resTab.anovaSign = allAnovas.anovaSign(allAnovas.ix)';

end

function [pVals, pSign, figList, hBar] = ...
    anovaBars(tab,clustVar,socioEcGroup,varargin)
%group = 'Commute' or 'Age' --> group all covariates starting with this 
%   string in one plot

%parse inputs
myParser = inputParser;
addParameter(myParser,'figColor','w'); %figure color? ('none' doesn't show up on screen until saved to image; use 'w' for in-Matlab analysis)
addParameter(myParser,'compareType','kw'); %Kruskall-Wallace or ANOVA?
addParameter(myParser,'toplot',1); %make charts?
addParameter(myParser,'drawLine',1); %draw horizontal line indicating pValue threshold?
addParameter(myParser,'lineLevel',[]); %where to draw the line? if [], will choose Bonferroni
parse(myParser,varargin{:});
args = myParser.Results;

%ix of non-cluster index var
otherVar = cellfun(@(x) ~strcmpi(x,clustVar),...
    tab.Properties.VariableNames);

%split up data into x (covars) and y (response) tables
yTab = tab(:,{clustVar});
xTab = tab(:,otherVar);

%which socioeconomic variables to do correlations against
if strcmp(socioEcGroup,'')
    barVarIx = logical(ones(1,numel(xTab.Properties.VariableNames)));
    titStr = 'all';
else
    barVarIx = cellfun(@(x) compStr(x,socioEcGroup),...
        xTab.Properties.VariableNames);
    titStr = socioEcGroup;
end
xTab = xTab(:,barVarIx);

%how many variables are we doing anova for? how many clusters?
[~,numVars] = size(xTab);
numClusters = numel(unique(yTab.(clustVar)));

%indicator variables for whether this FIPS is in clust n or not
for n=1:numClusters
    yTab.(['cluster' num2str(n)]) = yTab.hierClust == n;
end

%what grouping variables to do ANOVAs on
colors = hsv(numClusters); %color-code clusters
for n=1:numClusters
    anovaVars{n} = ['cluster' num2str(n)];
    anovaColors{n} = colors(n,:);
end

%compute anovas for all variables in myTab (columns) using grouping
%variable groupVar
for n=1:numClusters       
    allAnovas{n} = doAnova([xTab yTab(:,anovaVars{n})],...
        anovaVars{n},args.compareType);
end

%unpack allAnovas{p} into a numClusters x numPval matrix for bar()
for n=1:numClusters
    pVals(n,:) = allAnovas{n}.p;
    pSign(n,:) = allAnovas{n}.anovaSign;
end

%if more than 20 vars, make multiple figures
numPltsPerWin = 20;
numFig = ceil(numVars/numPltsPerWin);

%pad corrList, pList with zeros to have same number of entries for all
%groups
if numVars < numPltsPerWin
    pVals(:,end+1:numPltsPerWin) = nan;
    pSign(:,end+1:numPltsPerWin) = nan;
end

figList = [];
hBar = [];

%make figures
if args.toplot
    barVals = -log10(pVals).*pSign;
    barVals(pSign == 0) = -log10(pVals(pSign == 0)); %set the ones with same median to be positive
    yLab = ''; %'log_{10} p-value';
    yLim = [min(barVals(:)) max(barVals(:))];
    yLim = yLim + (yLim(2)-yLim(1))*0.1*[-1 1];
    %set up 3 ticks above and 3 below the 0 axis
    tickMax = max(abs(yLim));
    
    %want ticks to be round numbers
    roundNums = [1 2 5:5:100];
    try
        tickUnit = max(1,roundNums(find(roundNums > tickMax/3,1,'first')-1));
    catch err
        warning(err.message)
        tickUnit = 1; %mean(yLim);
    end
    posTick = [0:tickUnit:yLim(2)];
    negTick = fliplr([-tickUnit:-tickUnit:yLim(1)]);
    yTick = [negTick posTick];
    yTickLab = abs(yTick);
    xVals = repmat(1:numPltsPerWin,numel(barVals(:,1)),1);

%figure size
figWidth = 7.5;

for i=1:numFig
    barsRemaining = numVars - (i-1)*numPltsPerWin;
    barsInFig = (i-1)*numPltsPerWin + [1:min(numPltsPerWin,barsRemaining)];
    
    %make figure
    figList(i) = figure('units','inches',...
        'position',[0 2 0.4+figWidth*(numel(barsInFig))/(numPltsPerWin) 3.3],...
        'color',args.figColor);
    
    %draw significance level line below bar, if you want it
    if args.drawLine
        if isempty(args.lineLevel)
            lineLevel = -log10(0.05/(numVars*numClusters)); %Bonferroni corrected pvalue
        else
            lineLevel = args.lineLevel;
        end
        hL(1) = plot([0 numel(barsInFig)+1],[lineLevel lineLevel],...
            'linestyle',':','color','k');
        hold on;
        hL(2) = plot([0 numel(barsInFig)+1],-[lineLevel lineLevel],...
            'linestyle',':','color','k');
        hold on;
        
        %y limits should be bigger than significance lines
        limRng = diff(yLim);
        adj_ticks=0;
        if yLim(1) > -lineLevel
            yLim(1) = -lineLevel - 0.1*limRng;
            adj_ticks = 1;
        end
        if yLim(2) < lineLevel
            yLim(2) = lineLevel + 0.1*limRng;
            adj_ticks = 1;
        end
        if adj_ticks
            posTick = [0:tickUnit:yLim(2)];
            negTick = fliplr([-tickUnit:-tickUnit:yLim(1)]);
            yTick = [negTick posTick];
            yTickLab = abs(yTick);
        end
    end
    
    %make bars
    if numel(barsInFig) > 1
        %'grouped' only plots groups if have more than one group
        hBar{i} = bar(1:numel(barsInFig),barVals(:,barsInFig)','grouped',...
                                    'FaceColor','flat'); %the 'flat' allows to use CData below
    else
        %need to make a fake second group to get this to work --> just
        %duplicate here, will trim with axis limits below
        fakeBars = [barVals barVals];
        hBar{i} = bar(1:2,...
                      fakeBars(:,1:2)','grouped','FaceColor','flat');
    end
   
    %color each bar group according to the cluster
    for bn=1:numel(hBar{i})
        %hBar{i}(bn).FaceColor = anovaColors{bn};
        hBar{i}(bn).CData = repmat(anovaColors{bn},numel(hBar{i}(bn).XData),1);
    end
    
    for bn=1:numel(hBar{i})
        for ix=1:numel(hBar{i}(bn).XData)
            if pSign(bn,ix) == 0
                hBar{i}(bn).CData(ix,:) = [1 1 1];
            end
        end
    end
    
    %label plot
    title(titStr);
    ezNames = cellfun(@(x) ezCovarName(x,'topad',1), ...
        xTab.Properties.VariableNames,'UniformOutput',0);
    
    %format axis
    set(gca,'xtick',[1:numel(barsInFig)],...
        'xlim',[0 numel(barsInFig)] + [0.5 0.5],'xticklabels',ezNames(barsInFig),...
        'ylim',yLim,'fontsize',8);
    set(gca,'ytick',yTick,'yticklabels',yTickLab);

    ylabel(yLab);
    xtickangle(45);
    figPos = get(gcf,'position');
    
    %set position in inches so that all the figures align
    set(gca,'units','inches',...
        'outerposition',figPos,...
        'position',[0.35 1.2 figWidth*(numel(barsInFig))/(numPltsPerWin) 1.8]); %was [0.4 1.2 4.5*(numVars+1)/21 1.8]
    axPos = get(get(gcf,'children'),'position');
    axFrac = 1/max(axPos(3:4)); %1/ax length  
    set(gca,'TickLength', [0.02 0.035]*2*axFrac);
end

end
  
end

%% clustering
function [dsClustered,hierClustIx] = clusterDataBy(ds,clustBy,varargin)
%return clustered dataset, ordered by clustering variable

%parse inputs
parser = inputParser;
addParameter(parser,'TOORDER',0,@isnumeric);
addParameter(parser,'MAXCLUST',4,@isnumeric);
parse(parser,varargin{:});   

MAXCLUST = parser.Results.MAXCLUST;
TOORDER = parser.Results.TOORDER;

rng(100); %seed random number generator for reproducibility
lt = 'ward';
%MAXCLUST = 4; %MAX. NO. CLUSTERS
linkType = lt; %like ward, weighted, complete, avge
ds.hierClust = clusterdata(ds.(clustBy), 'maxclust',MAXCLUST,...
    'distance','euclidean','linkage',linkType);

%label clusters s.t. lowest label = lowest mean clustBy value in the
%cluster
avgVal = [];
for m=1:MAXCLUST
    avgVal(m) = mean(mean(ds{ds.hierClust == m,clustBy}));
    avgVal(m) = max(mean(ds{ds.hierClust == m,clustBy})); % try getting the max
end

[~,ix] = sort(avgVal);
map2newlabel = (@(x,ix) find(ix == x));
ds.hierClust_sort = zeros(size(ds.hierClust));
for m=1:MAXCLUST
    ds.hierClust_sort(ds.hierClust == m) = map2newlabel(m,ix);
end

ds = sortrows(ds,'hierClust_sort','descend');

%sort within each cluster by avg sjl
ds.avg_clustBy = mean(double(ds.(clustBy)(:,:)),2); %avg summer
if TOORDER == 1
    for m=1:MAXCLUST
        s = ds(ds.hierClust_sort == m,:);
        s = sortrows(s,'avg_clustBy','descend');
        ds(ds.hierClust_sort == m,:) = s;
    end
end

ds.hierClust = ds.hierClust_sort;
ds.hierClust_sort = [];
dsClustered = ds;
hierClustIx = ds.hierClust;
end

function hFig = plotClustAvgs(dsClustered,clustBy,colors,yTitStr)
%plot a row of subplots. each subplot shows one SJL curve from 
%the most populous county in that cluster

ds = dsClustered;
MAXCLUST = max(ds.hierClust);

minVal = [];
maxVal = [];

for i=1:MAXCLUST
    hFig(i) = figure('units','inches','position',[0 0 1.8 1.8],...
        'color','w');
    thisClust = ds.hierClust == i;
    thisDs = ds(thisClust,:);
    thisDs.stdSJL = std(thisDs.(clustBy),0,2);
    
    %avg curve
    meanCurve = mean(thisDs{:,clustBy},1);   

    %plot mean+/-stdev as a filled patch
    xs=1:12;
    ys = mean(thisDs{:,clustBy},1);
    yerr = std(thisDs{:,clustBy},0,1);
    pltYs = [ys+yerr fliplr(ys-yerr)];
    pltXs = [xs fliplr(xs)];
    fill(pltXs,pltYs,colors(i,:),'FaceAlpha',0.5,'EdgeColor','none');
    hold on;
    
    %also mark the mean with a curve
    plot(xs,meanCurve,'-',...
        'marker','none',...
        'color',colors(i,:),'linewidth',1.5);
    title(['cluster ' num2str(i)]);
    xlabel('month');
    
    if i==1
        minVal = min(mean(thisDs{:,clustBy},1) - std(thisDs{:,clustBy},0,1));
        maxVal = max(mean(thisDs{:,clustBy},1) + std(thisDs{:,clustBy},0,1));
    else
        minVal = min(minVal,...
            min(mean(thisDs{:,clustBy},1) - std(thisDs{:,clustBy},0,1)));
        maxVal = max(maxVal,...
            max(mean(thisDs{:,clustBy},1) + std(thisDs{:,clustBy},0,1)));
    end
    
    monthnames = {'Jan','Feb','Mar',...
    'Apr','May','June',...
    'July','Aug','Sep',...
    'Oct', 'Nov', 'Dec'};
    ixMoDisp  = logical(zeros(size(1:12)));
    ixMoDisp(1:3:12) = true;
    dispMonths = monthnames;
    for m=1:numel(monthnames)
        if ~ixMoDisp(m)
            dispMonths{m} = '';
        end
    end
    set(gca,'xlim',[0 13],'xtick',1:1:12,'xticklabel',...
        dispMonths,'fontsize',8);
end

%label plots, unify axis limits
for i=1:MAXCLUST
    figure(hFig(i));
    yTick = -2:0.25:2;
    for yt=1:numel(yTick)
        if mod(yt,2) == 1
            yTickLabels{yt} = num2str(yTick(yt)*60);
        else
            yTickLabels{yt} = '';
        end
    end
    
    set(gca,'ylim',[-1 1],'ytick',yTick,...
        'yticklabels',yTickLabels,'FontSize',8);
    ylabel('fluct. in SJL (min)');
end

end

function fHierClust = plotHierClust(dsClustered,clustBy,colors,titStr)
% plot output of hierarchical clustering

ds = dsClustered;
MAXCLUST = max(dsClustered.hierClust); %MAX. NO. CLUSTERS

%label left colorbar
if strcmp(clustBy,'sjl')
    leftColorBarLabel = 'fluct. in SJL (min)';
    dataLim = [min(ds{:,clustBy}(:)) max(ds{:,clustBy}(:))];
    dataLim = [prctile(ds{:,clustBy}(:),1) prctile(ds{:,clustBy}(:),99)]; %toss outliers
    cAxLim = dataLim; %[-20 140]/60;
else
    leftColorBarLabel = [clustBy ' '];
    cAxLim = [min(min(ds.(clustBy))) max(max(ds.(clustBy)))];
end

monthnames = {'Jan','Feb','Mar',...
    'Apr','May','June',...
    'July','Aug','Sep',...
    'Oct', 'Nov', 'Dec'};

%make figure
fHierClust = figure('units','inches',...
                    'position',[0 0 3.5 3],'color','w');
colormap('redbluecmap');
ax1 = axes;
hIm1 = imagesc(ds.(clustBy),cAxLim);

%link second invisible axes on top of the first one
ax2 = axes;
hIm1 = imagesc(ds.(clustBy),cAxLim);
linkaxes([ax1,ax2]);
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

set([ax1,ax2],'Position',[.2 .11 .685 .815]);
h1 = colorbar(ax1,'position',[.12 .11 .0675 .815]);
set(h1,'ylim',[cAxLim],'ytick',[-1:0.25:3],...
        'ticklabels',60*[-1:0.25:3]);
xlabel(h1, leftColorBarLabel);

h2 = colorbar(ax2,'position',[.9 .11 .0675 .815]);
sz = size(ds.(clustBy));
if sz(2) == 12 %monthly
    set(ax1,'xtick',[1:12],'xticklabel',monthnames,'ytick',[]);
    xtickangle(ax1,45);
elseif sz(2) == 51 %weekly
    set(ax1,'xtick',[4:4:48],'xticklabel',4:4:48,'ytick',[]);
    xlabel(ax1,'week no.');
else
    warning(['unrecognized number of columns in dataset! ' ...
             'neither monthly nor weekly!']);
end
xlabel(h2, 'cluster no.');

%paint appropriate regions on colormap on the right
for m=1:MAXCLUST
    clustSz(m) = sum(ds.hierClust == m);
    if m > 1
        lo(m) = hi(m-1);
        hi(m) = lo(m) + clustSz(m);
    else
        lo(m) = 1;
        hi(m) = clustSz(m);
    end
end
for m=1:MAXCLUST
    for i=lo(m):hi(m)
        mymap(i,:) = colors(m,:);
    end
end
set(h2,'colormap',mymap,'ytick',[]);

set(ax1,'fontsize',8);
set(ax2,'fontsize',8);
set(h1,'fontsize',8);
set(h2,'fontsize',8);

end

%% utilities used to make figures & export them
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

%% process strings, names of covariates
function [displayName] = ezCovarName(longName,varargin)
%Replace code name of covariate such as 'Econ1' by a human-readable one like
%'employed'. Supports exponents in polynomial models (e.g., 'Econ1^2' --> 'employed^2').

%parse inputs
myParser = inputParser;
addParameter(myParser,'topad',0);
parse(myParser,varargin{:});
args = myParser.Results;

%parse exponent
ix = find(longName == '^');
if ~isempty(ix)
    pow = str2num(longName(ix+1:end));
    longName = longName(1:ix-1);
else
    pow = 1;
end

%rename covariate
c = load('datafiles/covarNameDict.mat');
c = c.covarNameDict;
ix = find(strcmpi(c{:,1},longName));
if ~isempty(ix)
    shortName = char(c{ix,2});
else
    shortName = longName;
end

%add exponent
if pow == 1
    displayName = [shortName];
else
    displayName = [shortName '^' num2str(pow)];
end

%pad with blanks s.t. all nams have same fixed width
if args.topad 
    maxlen = max(cellfun(@(x) length(x),c{:,2}));
    if length(displayName) < maxlen
       pad = blanks(maxlen-length(displayName));
       displayName = [pad displayName];
    end
end

%print to debug
% longName
% pow
% shortName
% displayName

end

function out = compStr(test,template)
%check if the beginning of 'char' test matches the 'char' template
if numel(template) > numel(test)
    out = false;
else
    n = numel(template);
    out = strcmp(test(1:n),template);
end

end

%% data cleaning, centering
function [myTab_centered,myTab_raw,y_raw,X_raw,FIPS] = ...
    makeMyTable(SJLtable,AllCovariates,varargin)
%merge SJL calc's and covariates table by FIPS
%SJLtable=table(FIPS,responseVar1,responseVar2,...),
%AllCovariates=table(FIPS,covar1,covar2,...)
%myTab_raw=table(covar1,covar2,...,responseVar1,responseVar2,...), merged
%by FIPS. FIPS order in output variable FIPS.

%parse inputs
myParser = inputParser;
addParameter(myParser,'exclCovar',{},@iscell);
parse(myParser,varargin{:});
args = myParser.Results;

%decide which columns to exlude from linear regression
colToExclude = {}; %nothing by default
if ~isempty(args.exclCovar)
    have_ix = ismember(args.exclCovar,...
        AllCovariates.Properties.VariableNames);
    if sum(have_ix) > 0
        colToExclude = args.exclCovar{have_ix};
    end
end

%get rid of covariates you don't like
AllCovariates(:,colToExclude) = []; 

%find common FIPS indices in the two tables
[~, i_sjl, i_cov] = innerjoin(SJLtable,AllCovariates,...
    'LeftKeys','FIPS','RightKeys','FIPS');

%grab subsets of both tables
myYtab = SJLtable(i_sjl,:);
myXtab = AllCovariates(i_cov,:);

%check that FIPS match up and get rid of that column
assert(sum(myYtab.FIPS == myXtab.FIPS) == numel(myYtab.FIPS));
FIPS = myYtab.FIPS;
myYtab.FIPS = [];
myXtab.FIPS = [];

%center covariates
myXtab_centered = centerData(myXtab); %okay to center FIPS because will get rid of that column

%make table for linear model fit
myTab_centered = [myXtab_centered myYtab];
myTab_raw = [myXtab myYtab];

%store uncentred data too
y_raw = myYtab{:,:};
X_raw = myXtab{:,:};

end

function [ centeredX ] = centerData(X)
%Given input matrix X, return matrix with every column normalized to 0
%mean, unit variance. X can be numeric or table.

[nr,nc] = size(X);

if istable(X)
    centeredX = X;
    for i=1:nc
        %disp(['i=' num2str(i)]);
        centeredX{:,i} = (X{:,i} - nanmean(X{:,i}))/nanstd(X{:,i});
    end
elseif isnumeric(X)
    centeredX = zeros(size(X));
    for i=1:nc
        centeredX(:,i) = (X(:,i) - nanmean(X(:,i)))/nanstd(X(:,i));
    end
end
    
end

function [outTab,...
          sigSJL_xc_mu,sigSJL_xc_sig,...
          sigSJL_tr_mu,sigSJL_tr_sig,...
          avgSJL_xc_mu,avgSJL_xc_sig,...
          avgSJL_tr_mu,avgSJL_tr_sig] = qcSJL(SJLtable)
%run basic quality controls on dataset
%use outputs to remove outliers
      
%trim troughs
SJLtable.min_wkday_trPos = min(SJLtable.wkday_trPos,[],2);
SJLtable.min_wkend_trPos = min(SJLtable.wkend_trPos,[],2);
SJLtable.max_wkday_trPos = max(SJLtable.wkday_trPos,[],2);
SJLtable.max_wkend_trPos = max(SJLtable.wkend_trPos,[],2);

%filtering out all negative SJL (using xcorr) and then clustering on xcorr
%gives 1070 counties out of 1521 --> clustering looks O-kay, but not great
% ---
%filtering out all negative SJL (using trough positions) and then
%clustering on trough differences gives 1196 counties out of 1521 -->
%cluster O-kay into 3 groups
SJLtable(:,'haveNegXc') = rowfun(@(x) sum(x < 0) > 0,SJLtable(:,'xc_sjl'));
SJLtable(:,'haveNegTr') = rowfun(@(x) sum(x < 0) > 0,SJLtable(:,'tr_sjl'));
SJLtable(:,'numNegXc') = rowfun(@(x) sum(x < 0),SJLtable(:,'xc_sjl'));
SJLtable(:,'numNegTr') = rowfun(@(x) sum(x < 0),SJLtable(:,'tr_sjl'));
SJLtable(:,'mostNegTr') = rowfun(@(x) x(find(x == min(x),1,'first')),SJLtable(:,'tr_sjl'));
SJLtable(:,'mostNegXc') = rowfun(@(x) x(find(x == min(x),1,'first')),SJLtable(:,'xc_sjl'));

%how variable is this county compared to others
SJLtable.sigSJL_xc = std(SJLtable.xc_sjl,0,2);
sigSJL_xc_mu = mean(SJLtable.sigSJL_xc);
sigSJL_xc_sig = std(SJLtable.sigSJL_xc);

SJLtable.sigSJL_tr = std(SJLtable.tr_sjl,0,2);
sigSJL_tr_mu = mean(SJLtable.sigSJL_tr);
sigSJL_tr_sig = std(SJLtable.sigSJL_tr);

%how is the avg sjl in this county compared to others
SJLtable.avgSJL_xc = mean(SJLtable.xc_sjl,2);
avgSJL_xc_mu = mean(SJLtable.avgSJL_xc);
avgSJL_xc_sig = std(SJLtable.avgSJL_xc);

SJLtable.avgSJL_tr = mean(SJLtable.tr_sjl,2);
avgSJL_tr_mu = mean(SJLtable.avgSJL_tr);
avgSJL_tr_sig = std(SJLtable.avgSJL_tr);

%convert to z scores
SJLtable.sigSJL_xc_z = (SJLtable.sigSJL_xc - sigSJL_xc_mu)./sigSJL_xc_sig;
SJLtable.sigSJL_tr_z = (SJLtable.sigSJL_tr - sigSJL_tr_mu)./sigSJL_tr_sig;

SJLtable.avgSJL_xc_z = (SJLtable.avgSJL_xc - avgSJL_xc_mu)./avgSJL_xc_sig;
SJLtable.avgSJL_tr_z = (SJLtable.avgSJL_tr - avgSJL_tr_mu)./avgSJL_tr_sig;

sum(SJLtable.haveNegTr)
sum(SJLtable.haveNegXc)

%outptu
outTab = SJLtable;

end

%% map functions
function hFigUSA_Clust = plotMap_Polygon_allClust(ds,geoUnit,colorBy,colors,titStr)
% plot us map and color all fips as polygons

% load the county shape file
%load('mapUSA_shape.mat'); %pre-loaded for faster plotting
load('datafiles/mapUSA_lores_counties_shape.mat'); %lo res file, should be faster

whoseMap = 'Matlab';
switch whoseMap
    case 'Matlab'
        load('datafiles/mapUSA_shape_states_Matlab.mat');
    case 'Census'
        load('mapUSA_shape_states.mat');
end

% what column sets colors of counties on the map?
ds.mapClust = ds.(colorBy);

%get indices of Hawaii, Alaska and continental states in your dataset
S_ix_Hawaii = strcmp(ds.state,'HI');
S_ix_Alaska = strcmp(ds.state,'AK');
S_ix_ConUS = ~(S_ix_Hawaii | S_ix_Alaska);

hFigUSA_Clust = figure();
set(hFigUSA_Clust,'position',[2 2 11 8.5]);
[~, hAxUSA] = plotUSA(hFigUSA_Clust);

%draw colored counties
for clst=1:numel(unique(ds.hierClust))
    %continental USA
    axes(hAxUSA(1));
    hold on;
    myds = ds(S_ix_ConUS & ds.hierClust == clst,:);
    drawFIPS(myds,shapeStruct,geoUnit,colors(clst,:),'none');
    
    %Alaska
    axes(hAxUSA(2));
    myds = ds(S_ix_Alaska & ds.hierClust == clst,:);
    drawFIPS(myds,shapeStruct,geoUnit,colors(clst,:),'none');
    
    %Hawaii
    axes(hAxUSA(3));
    myds = ds(S_ix_Hawaii & ds.hierClust == clst,:);
    drawFIPS(myds,shapeStruct,geoUnit,colors(clst,:),'none');
end

% mark all counties we have with empty squares
%continental USA
axes(hAxUSA(1));
myds = ds(S_ix_ConUS,:);
drawFIPS(myds,shapeStruct,geoUnit,'none','k');
drawStates(whoseMap,'ConUS',stateShapeStruct,'none','k');

%Alaska
axes(hAxUSA(2));
myds = ds(S_ix_Alaska,:);
drawFIPS(myds,shapeStruct,geoUnit,'none','k');
drawStates(whoseMap,'AK',stateShapeStruct,'none','k');

%Hawaii
axes(hAxUSA(3));
myds = ds(S_ix_Hawaii,:);
drawFIPS(myds,shapeStruct,geoUnit,'none','k');
drawStates(whoseMap,'HI',stateShapeStruct,'none','k');

%label plot with cluster number
%axes(hAxUSA(1));
%mapTitStr = ['all FIPS. ' titStr];
%textm(51,-100,titStr,'color','k');
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

