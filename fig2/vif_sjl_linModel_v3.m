%2018-08-05, EL: modify Fig. 2 for revision, include VIFs and timezone
%analysis
%2018-06-27, EL: Fig. 2C-3, Fig. S2A

%2018-07-18, MJR: modified to add calculation of VIF

%% input parameters
clear all; close all; clc;

%% which analysis would you like to do? export figures? (1=yes)
TOMAP = 0; % make Fig. 2C? (1=yes)
TOEXPMAP = 0; %export Fig. 2C?
TOEXP = 0; %export Fig. 2D?
TOMAKE_LinModel_BARS = 1; %make Fig. 2E? bar chars of top coefficients in the linear model
TOMAKEBARS = 0; %make Fig. S2A? bar chars of correlation coef's of covars vs avgSJL
TOEXPBARS = 0;  %export Fig. 2E, Fig. S2A? 
TOEXP_coefTab = 0; %export Table S1? 
TOVIF = 1; %calculate VIF
TOTIMEZONE = 0; %time zone specific linear model
TOEXP_SJLDATA = 0; %export data used to make fig. 2?

%% load sjl calculations
SJLtable = load(['../tweetograms/tweetograms_monthly_FIPS_2012_2013_dzah/'...
                   'sjl_troughs_2018-04-08_13.47.16_dzah.mat']);
SJLtable = SJLtable.outTable;
SJLtable.fit_wkendwkday = SJLtable.wkend_trPos - SJLtable.wkday_trPos;
SJLtable.Properties.VariableNames{1} = 'FIPS'; %geoCode = FIPS here for all 

%pick which 
binList = {[2 3 4 10]};%only calculated avg SJL over [Feb,Mar,Apr,Oct]
whichBin = [2 3 4 10]; %pick this average ([Feb,Mar,Apr,Oct])
binIx = cellfun(@(x) isequal(x,whichBin),binList,'UniformOutput',1);
SJLtable.('xcorr_wkdaywkend') = SJLtable.xcorr_wkdaywkend(:,binIx);

SJLtable.Properties.VariableNames{'xcorr_wkdaywkend'} = 'avgSJL'; 
responseVar = 'avgSJL';

%% correlate with shiftwork 
mili = readtable('datafiles/NewShiftWorkSJL.xlsx');
mili = mili(:,{'FIPS','FractionShift6pm_8am','FractionShift11pm_7am'});
miliSJL = join(mili,SJLtable(:,{'FIPS','avgSJL','wkday_trWidth','wkend_trWidth'}),'Keys','FIPS');
miliSJL.sleepDebt = miliSJL.wkend_trWidth - miliSJL.wkday_trWidth;
[rhoSJL_avgSJL,pSJL_avgSJL] = corr(miliSJL.FractionShift6pm_8am,miliSJL.avgSJL)
[rhoSJL_sleepDebt,pSJL_sleepDebt] = corr(miliSJL.FractionShift6pm_8am,miliSJL.sleepDebt)
[rhoSJL_wkdayWidth,pSJL_wkdayWidth] = corr(miliSJL.FractionShift6pm_8am,miliSJL.wkday_trWidth)
[rhoSJL_wkendWidth,pSJL_wkendWidth] = corr(miliSJL.FractionShift6pm_8am,miliSJL.wkend_trWidth)

%writetable(miliSJL,'mili_shiftWork_SJL.xlsx');

%% load population data
popTab = load('datafiles/AllCovariates_Apr15.mat'); %loads into AllCovariates_decorrelated
popTab = popTab.AllCovariates(:,{'FIPS','Population1','Population5'}); %.AllCovariates_decorrelated;
popTab.avgUsers = popTab.Population5.*(10.^popTab.Population1);
popTab.Properties.VariableNames{1} = 'FIPS'; %make sure name of FIPS column is right

SJLtable = innerjoin(SJLtable,popTab(:,{'FIPS','avgUsers'}),'Keys','FIPS');

%% toss observations > 3 sigma -- all correspond to counties with low counts
SJLtable = tossNoisy(SJLtable,responseVar);

%% export SJL data
if TOEXP_SJLDATA
    cTable = load('datafiles/CountyTableHaveData.mat');
    SJLtable_exp = join(SJLtable,cTable.cTable,'LeftKeys','FIPS','RightKeys','geoCode');
    SJLtable_exp = SJLtable_exp(:,{'FIPS','geoName','state','avgSJL',...
        'wkday_trPos',...
        'wkend_trPos',...
        'wkday_trWidth',...
        'wkend_trWidth'});
    SJLtable_exp.Properties.VariableNames = {'FIPS','name','state','avg_SJL_hr',...
        'wkday_troughTime_hr_am_localTime',...
        'wkend_troughTime_hr_am_localTime',...
        'wkday_troughWidth_hr',...
        'wkend_troughWidth_hr'};
    writetable(SJLtable_exp,['SJL_calcs_' getDate() '.xlsx']);
end

%% Fig. 2D: histogram of avg SJL
hFig = responseVarDistr(SJLtable,responseVar);
expif(TOEXP,hFig,['sjlDistr_mo' sprintf('-%d', whichBin)],...
    'OUTDIR','sjlDistr','fun2use','expfig','ext','pdf');

%% sleep debt figure
sleepDebt = SJLtable.wkend_trWidth - SJLtable.wkday_trWidth;
sleepDebtFig = figure('units','inches','position',[2 2 3 2]);
hSD = histogram(sleepDebt*60,[-160:20:160]);
set(gca,'xlim',[-160 160],'xtick',-160:40:160,'fontsize',8);
xlabel({'\Delta tweetogram width (min)','wkend - wkday'});
ylabel('no. counties');
expif(TOEXP,sleepDebtFig,'sleepDebt','OUTDIR','sleepDebt',...
    'fun2use','expfig','ext','pdf');
disp(['sleep debt avg+/-stdev: '...
        num2str(60*mean(sleepDebt),'%2.2f') ...
        ' +/- ' num2str(60*std(sleepDebt),'%2.2f')]);

%% compute avg +/- SJL for different time zones
if TOTIMEZONE
    CountyData = load('./datafiles/CountyTableHaveData.mat');
    CountyData.cTable.Properties.VariableNames{1} = 'FIPS';
    CombinedTable = join(SJLtable, CountyData.cTable);
    CombinedTable = CombinedTable(:,{'avgSJL','timezone'});
        
    [tzGroups, tzGroupID] = findgroups(CombinedTable.timezone);
    avgSJL_tz = 60*splitapply(@(x) mean(x),CombinedTable.avgSJL,tzGroups);
    stdSJL_tz = 60*splitapply(@(x) std(x), CombinedTable.avgSJL,tzGroups);
    numFIPS_tz = splitapply(@(x) numel(x), CombinedTable.avgSJL,tzGroups);
    tzTabOut = table(tzGroupID,numFIPS_tz,avgSJL_tz,stdSJL_tz);
    tzTabOut.Properties.VariableNames = {'timezone','num_counties', ...
        'mean_SJL_min', 'std_SJL_min'};
    %writetable(tzTabOut,['byTimezone/sjlByTimezone_' getDate() '.xlsx']);
    
    %do anovas
    [pAnova,~,anovaStats] = anova1(60*CombinedTable.avgSJL,tzGroups);
    [compResult,muComp,hMultCompare,groupNames] = multcompare(anovaStats);
    varNames = {'tz1','tz2','lo95CI','estDiffMeans','hi95CI','pVal'};
    tzAnovaTab = array2table(compResult);
    tzAnovaTab.Properties.VariableNames = varNames;
   
    tzAnovaTabNameOne = rowfun(@(x) tzGroupID(x), tzAnovaTab(:,{'tz1'}));
    tzAnovaTabNameTwo = rowfun(@(x) tzGroupID(x), tzAnovaTab(:,{'tz2'}));
    tzAnovaTab.tz1 = tzAnovaTabNameOne.Var1;
    tzAnovaTab.tz2 = tzAnovaTabNameTwo.Var1;
    %writetable(tzAnovaTab,['byTimezone/sjlByTimezoneAnova_' getDate() '.xlsx']);
end

%% load covariates
AllCovariates = load('datafiles/AllCovariates_Apr15.mat'); %loads into AllCovariates_decorrelated
AllCovariates = AllCovariates.AllCovariates; %.AllCovariates_decorrelated;
AllCovariates.Properties.VariableNames{1} = 'FIPS'; %make sure name of FIPS column is right

%% make table of responseVar + covariates
%myTab = [covar1,covar2,...,covarN,responseVar]
[myTab_centered,myTab_raw,y_raw,X_raw,FIPS] = ...
    makeMyTable(SJLtable(:,{'FIPS',responseVar}),...
                AllCovariates,'exclCovar',{'zzz'});

%% Calculate VIF
if TOVIF == 1
    VIFTable = AllCovariates;
    dims = size(VIFTable);
    NumPredictors = dims(2) - 1; %% 1st column is FIPS
    VIFList = zeros(1,NumPredictors);
    
    %remove rows containing NaN
    nanRows = rowfun(@(x) sum(isnan(x)) > 0,VIFTable,'SeparateInputs',false);
    VIFTable(nanRows.Var1,:) = [];
    
    for j = 1%1:NumPredictors
        %note that column 1 = FIPS here
        %make table moving column j + 1 to the end
        OthersTable = VIFTable(:,2:NumPredictors+1);
        PredictorJ = OthersTable(:,j);
        OthersTable(:,j) = [];
        ModelTable = horzcat(OthersTable,PredictorJ);
        
        %linear regression
        mdl = myLinModel(ModelTable,PredictorJ,'linear','toCV',0,'truthplot',0);
        JVIF = 1 / (1 - mdl.Rsquared.Ordinary); %should be ordinary Rsq
        VIFList(j) = JVIF;
    end
    
    %% now do it with the manually curated variables
    handPickList = readtable('datafiles/linModelDict_v2.xlsx',...
                                                    'Sheet','Fig2E_dawn');
    handPickList = handPickList.DescriptiveName;
    myTab_sel = [myTab_centered(:,handPickList) myTab_centered(:,{responseVar})];

    %remove rows that have NaN predictors
    nanRows = rowfun(@(x) sum(isnan(x)) > 0,myTab_sel,'SeparateInputs',false);
    myTab_sel(nanRows.Var1,:) = [];

    handVIFTable = myTab_sel;
    dims = size(handVIFTable);
    NumPredictors = dims(2); %% 1st column is FIPS
    handVIFList = zeros(1,NumPredictors-1);
    
    %in this case the final column is <SJL>, so remove it
    for j = 1:NumPredictors-1
        %make table excluding column j + 1
        OthersTable = handVIFTable(:,1:NumPredictors-1);
        PredictorJ = OthersTable(:,j);
        OthersTable(:,j) = [];
        ModelTable = horzcat(OthersTable,PredictorJ);
        
        %linear regression
        %mdl = myLinModel(ModelTable,PredictorJ,'linear','toCV',0,'truthplot',0);
        mdl = fitlm(ModelTable,'linear');
        JVIF = 1 / (1 - mdl.Rsquared.Ordinary);
        handVIFList(j) = JVIF;
        handVIFname{j} = PredictorJ.Properties.VariableNames{1};
    end
    %%
    handVIFout = table(handVIFList','VariableNames',{'VIF'},'RowNames',handVIFname);
    handVIFout.VIF(end+1) = nan;
    handVIFout.Properties.RowNames{end} = '(Intercept)';
    
    %% pairwise correlation
    RMatrix = corrcoef(table2array(handVIFTable(:,1:NumPredictors-1)));
    figure
    imagesc(RMatrix);
    colorbar
end

%% Fig. S2A: bar chars of correlation coef's of covars vs responseVar
if TOMAKEBARS == 1
groups = {'Commute','Edu','Population',...
          'Sex','Age','Race','Household',...
          'Employed','Occupation','Industry','WorkerClass',...
          'Income','GeoMeteo','Health','Religion','Politics'};
%groups = {'Population'};
      
for g=1:numel(groups)
    figList = [];
    [corrList, pList, figList, hBar] = ...
        corrBars(myTab_raw,responseVar,groups{g},'toplot',1,'drawLine',1);
    expif(TOEXPBARS,figList,['corrCoef_' groups{g} '_mo' ...
        sprintf('-%d', whichBin)],...
        'OUTDIR','./allCorrCoef',...
        'fun2use','expfig','ext','pdf');
end

%get the list of all correlations
[corrList, pList] = ...
    corrBars(myTab_raw,responseVar,'','toplot',0);
end

%% manually select covariates
%% edited to use a covariates list that includes dawn time at equinox
handPickList = readtable('datafiles/linModelDict_v2.xlsx',...
    'Sheet','Aug5_linModel'); %'Fig2E_dawn'); %'Aug5_linModel');
handPickList = handPickList.DescriptiveName;
myTab_sel = [myTab_centered(:,handPickList) myTab_centered(:,{responseVar})];

%remove rows that have NaN predictors
nanRows = rowfun(@(x) sum(isnan(x)) > 0,myTab_sel,'SeparateInputs',false);
myTab_sel(nanRows.Var1,:) = [];

%% make linear model from pre-selected coeffs
mdl = myLinModel(myTab_sel,responseVar,'linear','toCV',0,'truthplot',0);

%% Fig. 2E: make bar chars of top coefficients in the linear model
if TOMAKE_LinModel_BARS == 1

%get significant covariates from the model
[allCoef, displayCoef, Pi0] = fdrCoeff(mdl);

%make table for supplement
saveTab = allCoef;
saveTab.name = cellfun(@(x) ezCovarName(x), saveTab.Properties.RowNames,...
                                            'UniformOutput',false);
saveTab.codeName = saveTab.Properties.RowNames;
if TOVIF
    saveTab = join(saveTab,handVIFout,'Keys','Row');
end
if TOEXP_coefTab
    writetable(saveTab,['covar_LM_table_' getDate() '.xlsx']);
end

%run Bonferroni correction and make a smaller table comprised of only those
%coefficients that meet the Bonferroni criterion (at FWER < 0.1)
bonfCorr = allCoef(allCoef.pValue < 0.1/numel(allCoef.pValue),:);
bonfCorr('(Intercept)',:) = [];

%order rows for display
[~,~,~,bonfCorr.rank] = cellfun(@(x) covarNameRank(x), ...
    bonfCorr.Properties.RowNames,'UniformOutput',0);
bonfCorr = sortrows(bonfCorr,'rank');

mdl = fitlm(myTab_sel(:,{bonfCorr.Properties.RowNames{:},'avgSJL'}),'linear');

%no need to separate into different charts by name; they're already
%'ranked' in order we want from previous command
groups = {''};
for g=1:numel(groups)
    figList = [];
    [nCoef,~] = size(allCoef);
    [figList, hBar] = lmCoefBars(bonfCorr,groups{g},'toplot',1);
    for fl=1:numel(figList)
        set(figList(fl),'units','inches','position',[0 0 3.6 2]);
    end
    if ~isempty(figList)
        expif(TOEXPBARS,figList,...
            ['lmCoef_' groups{g} '_mo' sprintf('-%d', whichBin)],...
            'OUTDIR','./linModel',...
            'fun2use','expfig','ext','pdf');
    end
end

end

%% time zone specific linear model
if TOTIMEZONE == 1
    CountyData = load('./datafiles/CountyTableHaveData.mat');
    CountyData.cTable.Properties.VariableNames{1} = 'FIPS';
    CombinedTable = join(SJLtable, CountyData.cTable);
    
    %% remove rows without desired timezone
    timezoneRows = [];
    %%timezoneRows = ~strcmp(CombinedTable.timezone,'US/Central');
    CombinedTable(timezoneRows,:) = [];
    
    [tz_myTab_centered,tz_myTab_raw,y_raw,X_raw,FIPS] = ...
    makeMyTable(CombinedTable(:,{'FIPS',responseVar}),...
                AllCovariates,'exclCovar',{'zzz'});
            
    %% manually select covariates
handPickList = readtable('datafiles/linModelDict_v2.xlsx','Sheet','Fig2E_dawn');
handPickList = handPickList.DescriptiveName;
tz_myTab_sel = [tz_myTab_centered(:,handPickList) tz_myTab_centered(:,{responseVar})];

    %% make linear model from pre-selected coeffs
tz_mdl = myLinModel(tz_myTab_sel,responseVar,'linear','toCV',0,'truthplot',0);

    
end

%% Fig. 2C: plot avg SJL on a map
if TOMAP == 1
%some geographical info from countyTable
geoTab = load('datafiles/CountyTableHaveData.mat');
geoTab = geoTab.cTable;
geoTab = innerjoin(geoTab,SJLtable(:,{'FIPS',responseVar}),...
    'LeftKeys','geoCode','RightKeys','FIPS');
geoTab.Properties.VariableNames{1} = 'FIPS';
[hFigUSA, ~] = ...
    plotMap_Polygon_label(geoTab, 'FIPS', responseVar,'',[]);
set(hFigUSA,'color','w','position',[0 0 6 4.5]);

%note: using these export settings, it looks like the colorbar isn't
%rendered in CMYK, whereas the map itself is. I couldn't figure out how to
%fix this with painters renderer (while still using cmyk), but this seems
%to sort itself out when imported into Illustrator file that's CMYK PDF.
expif(TOEXPMAP,hFigUSA,['map_avgSJL_mo' sprintf('-%d', whichBin)],...
    'OUTDIR','./map','fun2use','expfig','ext','pdf');  
end

%% helper funs

%linear model & analysis of correlations
function mdl = myLinModel(myTab,responseVar,modelSpec,varargin)

%parse args
myParser = inputParser;
addParameter(myParser,'toCV',1,@isnumeric);    %cross-validate?
addParameter(myParser,'toExpCV',0,@isnumeric); %export cross-val figures?
addParameter(myParser,'truthPlot',1,@isnumeric); %make truthplot?
addParameter(myParser,'truthPlotType','points',@ischar); %points or heatmap?
addParameter(myParser,'toExp',0,@isnumeric); %export truth plot?
addParameter(myParser,'OUTDIR','./pics',@ischar); %export where?
parse(myParser,varargin{:});
args = myParser.Results;

%run linear model with all coefficients (full model)
%by default, last column is the responseVar here
mdlDummy = fitlm(myTab,modelSpec); %to detect outliers

% %re-fit after tossing outliers ( > 3 stdevs away from the initial fit)
% outlrObs = find(abs(mdlDummy.Residuals.Standardized) > 3); 
% mdl = fitlm(myTab,modelSpec,'Exclude',outlrObs);
outlrObs = [];
mdl = mdlDummy;

disp('lin model: Rsq adj');
disp(mdl.Rsquared.Adjusted);

% run cross validation on full model
if args.toCV
numCV = 5;
[R, TEST, TRAIN,fCovars] = cvLinModel(myTab, ...
                                mdl.Formula, responseVar, ...
                                outlrObs, numCV);
disp(['cross-valid R = ' num2str(mean(R),'%2.2f') ...
    ' +/- ' num2str(std(R), '%2.2f')]);

expif(args.toExpCV,fCovars,[responseVar '_CVcovarPlot_FullModel'],...
    'OUTDIR',args.OUTDIR);
end

%plot full model
if args.truthPlot
    switch args.truthPlotType
        case 'points'
            fH_fullModel = plotModel(mdl,[],outlrObs);
            expif(args.toExp,fH_fullModel,...
                [responseVar '_fullModel'],'OUTDIR',args.OUTDIR);
        case 'heatmap'
            fH_fullModel_HM = plotModelHeatMap(mdl,myTab,[]);
            set(fH_fullModel_HM,'renderer','zbuffer');
            expif(args.toExp,fH_fullModel_HM,...
                [responseVar '_fullModelHM'],'OUTDIR',args.OUTDIR);
    end
end

end

function [R, TEST, TRAIN, fCOVARS] = cvLinModel(myTab, myFormula, ...
    ResponseVar, outlrObs, NFOLD)
%run cross-validated linear model on myTab using model in myFormula. If
%myFormula = [], use all predictors except for ResponseVar in the fit.

%toss outliers right away
if ~isempty(outlrObs)
    myTab(outlrObs,:) = [];
end

%parse how many pieces to break up data into
if isempty(NFOLD)
    NFOLD = 5; %by default
end

%break up data on 10 groups (returns random indices 1-10 for each row)
ix = crossvalind('Kfold', size(myTab,1), NFOLD);

%for each group, run linear regression on training data and get R squared on
%testing data
for i=1:NFOLD
    train_ix = ix ~= i;
    test_ix = ~train_ix;
    
    %train model
    mdl = fitlm(myTab(train_ix,:),myFormula,'ResponseVar',ResponseVar);
    
    %predict test variable
    yfit = predict(mdl,myTab(test_ix,:));
    
    %get rsq adj, save it
    R(i) = getRsq(myTab{test_ix,ResponseVar},yfit,numel(mdl.Coefficients(:,1)));
    TEST(:,i) = test_ix;
    TRAIN(:,i) = train_ix;
    
    %plot covariates
    fCOVARS(i) = plotCovariates(mdl.Coefficients,myTab(test_ix,:), mdl, 'numFig',1);
    title(['i= ' num2str(i) '. R^2_{adj, test}=' num2str(R(i),'%2.2f')]);
end

    %helper fun
    function  [rsq, rsq_adj] = getRsq(ys,linfit,numcoef)
    %helper function to compute Rsq and adj Rsq from true ys vs fit ys
        resid = ys - linfit;
        SSresid = nansum(resid.^2);
        SStotal = (length(ys)-1)*nanvar(ys); %sum of variances per pt
        rsq = 1- SSresid/SStotal;
        rsq_adj = 1 - SSresid/SStotal * (length(ys)-1)/(length(ys)-numcoef);
    end
end

function [figList, hBar] = lmCoefBars(coefTab,groupVar,varargin)

%parse inputs
myParser = inputParser;
addParameter(myParser,'toplot',1);
addParameter(myParser,'drawLine',0);
addParameter(myParser,'lineLevel',0,@isnumeric);
addParameter(myParser,'valToBar','coef',@ischar);
parse(myParser,varargin{:});
args = myParser.Results;

%which variables to group together into a single plot
if strcmp(groupVar,'')
    
    %get rid of intercept term
    if sum(cellfun(@(x) strcmp(x,'(Intercept)'),coefTab.Properties.RowNames)) > 0     
        coefTab({'(Intercept)'},:) = [];
    end
    
    %keep everything else
    barVarIx = logical(ones(1,numel(coefTab.Properties.RowNames)));
    titStr = 'all';
else
    barVarIx = cellfun(@(x) compStr(x,groupVar),...
        coefTab.Properties.RowNames);
    titStr = groupVar;
end

coefTab = coefTab(barVarIx,:)

[numVars,~] = size(coefTab);

figList = [];
hBar = [];

if args.toplot & ~isempty(coefTab)
    switch args.valToBar
        case 'coef'
            barVals = coefTab.Estimate;
            barErr = coefTab.SE;
            %yLab = '\rho (Pearson)';
            yLab = '';
            yLim = [-0.2 0.2];
            yTick = [-1:0.1:1];
            yTickLab = []; %makes it easier to position plots
            %yTickLab = {'-1','','-0.5','','0','','0.5','','1'};
        case 'p'
            barVals = -log10(coefTab.pValue);
            yLab = '-log_{10} p-value';
            yLim = [0 max(barVals)];
            yTick = [];
            yTickLab = [];
    end
    
    switch groupVar
        case ''
            %if no grouping -- plot all bars. Match ticks and size to
            %truthplot figure -- do all annotations in Illustrator for
            %expedience
            
            %first figure -- has bars but no annotations (add them in
            %manually)
            figList(1) = figure('units','inches',...
                'position',[2 2 2 2],'color','w');
            
            hBar = bar(1:numVars,barVals);
            hold on;
            hErr = errorbar(1:numVars,barVals,barErr,'k','linestyle','none');
            
            %format axis
            set(gca,'xtick',1:numVars,'xlim',[0 numVars+1],...
                'xticklabels',1:numVars,...
                'ylim',yLim,'fontsize',8);
            set(gca,'ytick',yTick,'yticklabels',yTick);
            xlabel('will - rename - bars');
            ylabel('will - rename - label');
            
            %second figure has same data + correct labels, but poor
            %formatting --> use this to add labels in Illustrator
            figList(2) = figure('units','inches',...
                'position',[0 2 2 2],'color','w');
            
            hBar = bar(1:numVars,barVals);
            hold on;
            hErr = errorbar(1:numVars,barVals,barErr,'k','linestyle','none');
            
            %format axis
            ezNames = cellfun(@(x) strrep(ezCovarName(x,'topad',1),'_',' '), ...
                coefTab.Properties.RowNames,'UniformOutput',0);
            set(gca,'xtick',1:numVars,'xlim',[0 numVars+1],'xticklabels',ezNames,...
                'ylim',yLim,'fontsize',5);
            xtickangle(45);
            ylabel(yLab);
            set(gca,'ytick',yTick,'yticklabels',yTick);
                   
        otherwise
            figList = figure('units','inches',...
                'position',[0 2 0.4+4.5*(numVars+1)/21 3.3],'color','w');
            
            hBar = bar(1:numVars,barVals);
            hold on;
            hErr = errorbar(1:numVars,barVals,barErr,'k','linestyle','none');
            
            if args.drawLine
                line([0 args.lineLevel],[numPltsPerWin args.lineLevel],...
                    'color','k');
            end
            title(titStr);
            ezNames = cellfun(@(x) strrep(ezCovarName(x,'topad',1),'_',' '), ...
                coefTab.Properties.RowNames,'UniformOutput',0);
            
            %format axis
            set(gca,'xtick',1:numVars,'xlim',[0 numVars+1],'xticklabels',ezNames,...
                'ylim',yLim,'fontsize',8);
            set(gca,'ytick',yTick,'yticklabels',yTickLab);
            
            ylabel(yLab);
            xtickangle(45);
            figPos = get(gcf,'position');
            
            %set position in inches so that all the figures align
            set(gca,'xlim',[0 numVars+1],'units','inches',...
                'outerposition',figPos,...
                'position',[0.3 1.2 4.5*(numVars+1)/21 1.8]);
            axPos = get(get(gcf,'children'),'position');
            axFrac = 1/max(axPos(3:4)); %1/ax length
            set(gca,'TickLength', [0.02 0.035]*4*axFrac);   
    end
end

end

function [allCoef, displayCoef, Pi0] = fdrCoeff(mdl,varargin)
%run FDR correction on p values, pick Q < 0.05

%parse args
myParser = inputParser;
addParameter(myParser,'toPlt',0,@isnumeric); %plot correlation against covariates?
addParameter(myParser,'toExp',0,@isnumeric); %export truth plot?
addParameter(myParser,'OUTDIR','./pics',@ischar); %export where?
parse(myParser,varargin{:});
args = myParser.Results;

%[FDR] = mafdr(mdl.Coefficients.pValue,'BHFDR',true);
%grab significant coef's
[FDR, Q, Pi0] = mafdr(mdl.Coefficients.pValue);
allCoef = mdl.Coefficients;
allCoef.FDR = FDR;

signif_coef_ix = find(FDR < 0.05);
sigCoef = allCoef(signif_coef_ix,:);
%sigCoef = sortrows(sigCoef,'pValue','ascend'); %don't sort, confusing to
%read

%how many are significant?
numTopCoefs = numel(sigCoef.Properties.RowNames);

%make formula based on significant terms
myShortFormula = [mdl.ResponseName ' ~ 1'];
for i=2:numTopCoefs 
%for i=2:min(10,numel(good_coef.Properties.RowNames)) %allow at most 10 rows
    myShortFormula = [myShortFormula ' + ' sigCoef.Properties.RowNames{i}];
end

%make display table
%topcoefs = min(10,numel(good_coef.Properties.RowNames));
displayCoef = sigCoef(1:numTopCoefs,:);
displayCoef.Properties.RowNames = cellfun(@(x) ezCovarName(x), ...
    displayCoef.Properties.RowNames,'UniformOutput',0);
fTable = makeTableFig(displayCoef);
expif(args.toExp,fTable,[mdl.ResponseName '_topCovarsTable'],'OUTDIR',args.OUTDIR)

%plot correlation against top covariates
if args.toPlt
    fH_covariates = plotCovariates(sigCoef, mdl.Variables, mdl);
    expif(args.toExp,fH_covariates,[mdl.ResponseName '_covars'],'OUTDIR',args.OUTDIR);
end

end

function [corrList, pList, figList, hBar] = ...
    corrBars(tab,responseVar,groupVar,varargin)

%parse inputs
myParser = inputParser;
addParameter(myParser,'toplot',0);
addParameter(myParser,'drawLine',0);
addParameter(myParser,'lineLevel','auto');
addParameter(myParser,'valToBar','rho',@ischar);
parse(myParser,varargin{:});
args = myParser.Results;

%ix of non-response var
otherVar = cellfun(@(x) ~strcmpi(x,responseVar),...
    tab.Properties.VariableNames);

%split up data into x (covars) and y (response) tables
yTab = tab(:,{responseVar});
xTab = tab(:,otherVar);

%which variables to do correlations against
if strcmp(groupVar,'')
    barVarIx = logical(ones(1,numel(xTab.Properties.VariableNames)));
    titStr = 'all';
else
    barVarIx = cellfun(@(x) compStr(x,groupVar),...
        xTab.Properties.VariableNames);
    titStr = groupVar;
end
xTab = xTab(:,barVarIx);

[~,numVars] = size(xTab);

%compute correlations
corrList = [];
pList = [];
for i=1:numVars
    [corrList(i), pList(i)] = ... %2018-08-06, EL: modified for spearman
        corr(yTab.(responseVar),xTab{:,i},'type','Spearman','rows','complete');
end

%if more than 20 vars, make multiple figures
numPltsPerWin = 19;
numFig = ceil(numVars/numPltsPerWin);

%pad corrList, pList with zeros to have same number of entries for all
%groups
if numVars < numPltsPerWin
    corrList(end+1:numPltsPerWin) = nan;
    pList(end+1:numPltsPerWin) = nan;
end

figList = [];
hBar = [];

%make figures
if args.toplot
    switch args.valToBar
        case 'rho'
            barVals = corrList;
            %yLab = '\rho (Pearson)';
            yLab = '';
            yLim = [-0.6 0.6];
            yTick = [-1:0.25:1];
            yTickLab = []; %makes it easier to position plots
            %yTickLab = {'-1','','-0.5','','0','','0.5','','1'};
            
            %compute line-level for Bonferroni-corrected FWER < 0.1
            numSamples = numel(tab{:,1}); %no. rows in table
            numCovarsTotal = numel(tab{1,:})-1; %total no. covariates
            alphaLevel = 0.1; % what FWER to correct to
            t_crit = tinv(1-(alphaLevel/numCovarsTotal),numSamples-2); %0.95% alpha t-val
            rho_crit = t_crit/sqrt(numSamples - 2 + t_crit^2);
            
            if strcmp(args.lineLevel,'auto')
                args.lineLevel = rho_crit;
            end
            
        case 'p'
            barVals = -log10(pList);
            yLab = '-log_{10} p-value';
            yLim = [0 max(barVals)];
            yTick = [];
            yTickLab = [];
    end
    

for i=1:numFig
    figList(i) = figure('units','inches',...
        'position',[0 2 0.4+4.5*(numVars+1)/21 3.3],'color','w');

    %draw bars
    hBar(i) = bar(1:numPltsPerWin,barVals);
    
    %draw significance line?
    if args.drawLine
        hold on;
        hL(1) = line([0 numPltsPerWin],[args.lineLevel args.lineLevel],...
            'color','k','linestyle',':');
        hL(2) = line([0 numPltsPerWin],-[args.lineLevel args.lineLevel],...
            'color','k','linestyle',':');
        uistack(hL,'bottom');
    end
    
    title(titStr);
    ezNames = cellfun(@(x) ezCovarName(x,'topad',1), ...
        xTab.Properties.VariableNames,'UniformOutput',0);
    
    %format axis
    set(gca,'xtick',1:numVars,'xlim',[0 numVars+1],'xticklabels',ezNames,...
        'ylim',yLim,'fontsize',8);
    set(gca,'ytick',yTick,'yticklabels',yTickLab);

    ylabel(yLab);
    xtickangle(45);
    figPos = get(gcf,'position');
    
    %set position in inches so that all the figures align
    set(gca,'xlim',[0 numVars+1],'units','inches',...
        'outerposition',figPos,...
        'position',[0.35 1.2 4.5*(numVars+1)/21 1.8]); %was [0.4 1.2 4.5*(numVars+1)/21 1.8]
    axPos = get(get(gcf,'children'),'position');
    axFrac = 1/max(axPos(3:4)); %1/ax length  
    set(gca,'TickLength', [0.02 0.035]*4*axFrac);
    set(gca,'layer','top'); %this moves axis ticks, etc., above patch
end

end

  
end

function  [rsq, rsq_adj] = getRsq(ys,linfit,numcoef)
%helper function to compute Rsq and adj Rsq from true ys vs fit ys
resid = ys - linfit;
SSresid = nansum(resid.^2);
SStotal = (length(ys)-1)*nanvar(ys); %sum of variances per pt
rsq = 1 - SSresid/SStotal;
rsq_adj = 1 - SSresid/SStotal * (length(ys)-1)/(length(ys)-numcoef);
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

%% graphical tools for linear model & correlation analysis
function hFig = responseVarDistr(myTab,responseVar,varargin)
%parse arguments
myParser = inputParser;
addParameter(myParser,'fH',[]);
addParameter(myParser,'xLab','social jet lag on Twitter (min)');
parse(myParser,varargin{:});
args = myParser.Results;

y = myTab{:,responseVar};

hFig = figure('units','inches','position',[0 0 2 2]);
h = histogram(y,20);

minX = min(0,min(h.BinEdges));
maxX = max(2.25,max(h.BinEdges));
xtick = minX:0.25:maxX;
xticklab = num2cell(60*xtick);
xticklab = cellfun(@num2str,xticklab,'uniformoutput',0);
for i=2:2:numel(xticklab)
    xticklab{i} = '';
end

set(gca,'xlim',[minX maxX],...
    'xtick',xtick,...
    'XTickLabel',xticklab);
xlabel(args.xLab);
ylabel('no. counties');
set(gca,'fontsize',8);


end

function modelFig = plotCovariates(covars,myTab, mdl, varargin)
%for top covariates, plot response var ~ covariate
%covars = table of covariates output by FDR calc.

%parse arguments
myParser = inputParser;
addParameter(myParser,'fH',[]);
addParameter(myParser,'outlr',[]);
addParameter(myParser,'numFig','max');
parse(myParser,varargin{:});
args = myParser.Results;
fH = args.fH;
outlr = args.outlr;

%sort covariates by pValue
covars = sortrows(covars,'pValue');
intercept_row_ix = strcmpi(covars.Properties.RowNames, '(Intercept)');
intercept_row = covars(intercept_row_ix,:); %save it for linear plots
covars(intercept_row_ix,:) = []; %get rid of intercept

%calculate how many unique covariates there are, order by p value and group
%different polynomial powers together
cnt=1;
nameList = {};
ixList = {};
powList = {};
for i=1:numel(covars.Row)
    dum = covars.Properties.RowNames{i};
    dumShort = dum(1:find(dum == '^')-1);
    dumPow = str2num(dum(find(dum == '^')+1:end));
    if isempty(dumPow)
        dumPow = 1;
        dumShort = dum; %if don't have a ^ in name, need to keep old one
    end
    %new entry that has multiple powers
    if ~isempty(dumShort) & ~ismember(dumShort,nameList) 
        nameList{cnt} = dumShort;
        ixList{cnt} = [i];
        powList{cnt} = dumPow;
        cnt=cnt+1;
    %repeat entry (second power)
    elseif ~isempty(dumShort) & ismember(dumShort,nameList)
        ix = find(cellfun(@(x) strcmp(x,dumShort),nameList));
        ixList{ix} = [ixList{ix}(:) i];
        powList{ix} = [powList{ix}(:) dumPow];
    %entry has only one power
    else
        nameList{cnt} = dum;
        ixList{cnt} = [i];
        powList{cnt} = dumPow;
        cnt=cnt+1;
    end
end

%how many figures we have -- depends on number of unique covariates
if ischar(args.numFig) & strcmp(args.numFig,'max')
    numFig = ceil(cnt/6);
elseif isnumeric(args.numFig)
    numFig = args.numFig;
end
nCovarsToPlt = min(numel(nameList),numFig*6); %six covars per fig

%parse whether passed a figure handle and the right no. of handles
if nargin > 3 && ~isempty(fH) && numel(fH) == numFig 
    modelFig = fH;
    figure(modelFig);
else
    for i=1:numFig
        modelFig(i) = figure();
    end
end

%plot covariates one by one
for i=1:nCovarsToPlt
    figure(modelFig(ceil(i/6)));
    pp(i) = subplot(2,3,mod(i-1,6)+1);
    pltVar = nameList{i};
    pow = powList{i}(:);
%     [~, ~] = coloredBins(myTab{:,pltVar}.^pow,...
%         myTab{:,mdl.ResponseName},30,gcf,subplot(2,3,i));
%     hold on;
    pp(i)=plot(myTab{:,pltVar},...
        myTab{:,mdl.ResponseName},'.');
    hold on;
    
    %plot linear fit
    minX = min(myTab{:,pltVar});
    maxX = max(myTab{:,pltVar});
    xx = [minX:0.1:maxX];
    yy = intercept_row{'(Intercept)','Estimate'};
    for j=1:numel(ixList{i})
        yy = yy + (xx.^pow(j))*...
            covars{covars.Properties.RowNames{ixList{i}(j)},'Estimate'};
    end
    plot(xx,yy,'b-');
    
    powName = [ezCovarName([pltVar '^' num2str(pow(1))])];
    for j=2:numel(pow)
        powName = [powName ' + ' ezCovarName([pltVar '^' num2str(pow(j))])];
    end
    xlabel(powName);
    ylabel(mdl.ResponseName);
end

%plot outliers
outlierStr = ['no outliers'];
if nargin > 4
    for i=1:nCovarsToPlt
        figure(modelFig(ceil(i/6)));
        subplot(2,3,mod(i-1,6)+1);
        pltVar = nameList{i};
        plot(myTab{outlr,pltVar}, ...
            myTab{outlr,mdl.ResponseName},'g.');
    end
    outlierStr = [num2str(numel(outlr)) ' outliers'];
end

%add legends (do last so get only one legend entry)
%for i=1:min(6,numel(covars.Row))
for i=1:nCovarsToPlt
    figure(modelFig(ceil(i/6)));
    subplot(2,3,mod(i-1,6)+1);
    legend boxoff;
    legend(pp(i),['p=' num2str(covars{i,4},'%0.1e')]);
end

%title
suplabel([mdl.ResponseName '. ' ...
     'top predictors.' ],'t');
 
%change size
set(modelFig,'units','inches','position',[0 0 11 6]);
end

function modelFig = plotModel(mdl, fH, outlr,varargin)
%plot true Y vs. model prediction

%parse arguments
myParser = inputParser;
addParameter(myParser,'pltPretty',0,@isnumeric); %format figure for pub (=1) or no for analysis?
addParameter(myParser,'plotType','points',@ischar);
addParameter(myParser,'toColorBar',0);
addParameter(myParser,'figSize',[2 2 2.16 2.16]);
parse(myParser,varargin{:});
args = myParser.Results;

myTab = mdl.Variables;

%parse whether passed a figure handle
if nargin > 2 && ~isempty(fH)
    modelFig = fH;
    figure(modelFig);
else
    modelFig = figure();
end

plt(2)=plot(myTab{:,mdl.ResponseName},predict(mdl,myTab),'.','markersize',2);
hold on;

%plot outliers
outlierStr = ['no outliers']; 
if nargin > 3
    plot(myTab{outlr,mdl.ResponseName},...
        predict(mdl,myTab{outlr,...
        ~strcmpi(myTab.Properties.VariableNames,mdl.ResponseName)}),'g.');
    outlierStr = [num2str(numel(outlr)) ' outliers'];
end

xlabel(mdl.ResponseName);
ylabel(['lin. model with ' num2str(numel(mdl.PredictorNames)) ' predictors']);

axis(gca,'equal');

xlim = get(gca,'Xlim');
ylim = get(gca,'Ylim');
comLim = [min(xlim(1),ylim(1)) max(xlim(2),ylim(2))];
comLim(1) = min(comLim(1),25/60); %want to get at least 30 to 120 on the plot
comLim(2) = max(comLim(2),125/60);
comTick = [comLim(1)-mod(comLim(1),1):0.5:comLim(2)];
set(gca,'xlim',comLim,'ylim',comLim,'xticklabels',comTick*60);
set(gca,'xtick',comTick','ytick',comTick,'yticklabels',comTick*60);
plt(1)=plot(get(gca,'xlim'),get(gca,'ylim'),'k-');
uistack(plt);
set(modelFig, 'units','inches','position',[2 2 2 2]);

if strcmp(args.plotType,'heatmap')
    delete(plt);
    edges = comLim(1):(5/60):comLim(2);
    h = histogram2(myTab{:,mdl.ResponseName},predict(mdl,myTab),...
        edges,edges,...
        'DisplayStyle','tile','ShowEmptyBins','on',...
        'EdgeColor',[0.95 0.95 0.95],'LineWidth',0.01);
    hold on;
    line(get(gca,'xlim'),get(gca,'ylim'),'color','k','linewidth',1);
    
    mymap = flipud(hot(512)); %was (hot)
    mymap = [1 1 1; mymap(20:end,:)];
    colormap(mymap);
    hold on;
    
    if args.toColorBar
        colorbar(gca);
    end
    hold on;
end

set(gca,'layer','top'); %this moves axis ticks, etc., above patch
set(gcf,'Renderer','opengl'); %for proper export

%nicer format
if args.pltPretty
    
if strcmp(mdl.ResponseName,'avgSJL')
    xlabel('social jet lag on Twitter (min)');
end
ylabel(['linear model']);   
    
set(modelFig,'units','inches','position',args.figSize,'color','w'); %to make it match other 2x2s
mFigAx = get(modelFig,'children');
mFigAx = mFigAx(1); %only grab first axis, second could be colorbar?
set(mFigAx,'fontsize',8);

legTxt = ['R^2_{adj}=' num2str(mdl.Rsquared.Adjusted,'%2.2f')];
xLim = get(mFigAx,'xlim');
yLim = get(mFigAx,'ylim');
xTxt = xLim(1) + 0.05*(xLim(2) - xLim(1));
yTxt = yLim(1) + 0.95*(yLim(2) - yLim(1));
text(xTxt,yTxt,legTxt,'FontSize',8,...
    'VerticalAlignment','top');

end

end

function modelFig = plotModelHeatMap(mdl,myTab,fH, outlr)
%plot true Y vs. model prediction

%parse whether passed a figure handle
if nargin > 2 && ~isempty(fH)
    modelFig = fH;
    figure(modelFig);
else
    modelFig = figure();
end

[~, ~ ] = coloredBins(myTab{:,mdl.ResponseName},predict(mdl,myTab),30,modelFig,[]);
hold on;
%plot(myTab{:,mdl.ResponseName},predict(mdl,myTab),'.');
hold on;
plot(myTab{:,mdl.ResponseName},myTab{:,mdl.ResponseName},'k-');
%coloredBins(myTab{:,mdl.ResponseName},myTab{:,mdl.ResponseName},30);

% %plot outliers
outlierStr = ['no outliers']; 
if nargin > 3
    plot(myTab{outlr,mdl.ResponseName},...
        predict(mdl,myTab{outlr,...
        ~strcmpi(myTab.Properties.VariableNames,mdl.ResponseName)}),'g.');
    outlierStr = [num2str(numel(outlr)) ' outliers'];
end

xlabel(mdl.ResponseName);
ylabel(['lin. model with ' num2str(numel(mdl.PredictorNames)) ' predictors']);
title([mdl.ResponseName '. ' outlierStr '. ' 10 ...
     num2str(numel(mdl.PredictorNames)) ' predictors.' ...
    ' R^2_{adj}=' num2str(mdl.Rsquared.Adjusted,'%2.2f') ]);
axis(gca,'equal');

xlim = get(gca,'Xlim');
ylim = get(gca,'Ylim');
comLim = [min(xlim(1),ylim(1)) max(xlim(2),ylim(2))];
comTick = [comLim(1)-mod(comLim(1),20):20:comLim(2)];
set(gca,'xlim',comLim,'ylim',comLim);
set(gca,'xtick',comTick','ytick',comTick);
set(modelFig, 'units','inches','position',[0 0 4 4.5]);
end

%% utilities used to make figures & export them
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

function fTable = makeTableFig(tabIn,fTable)
%save formatted table output as image

if nargin > 2 & ~isempty(fTable)
    figure(fTable);
else
    fTable = figure('Units','Inches','Position',[0 0 8.5 11]);
end

 % Get the table in string form.
    TString = evalc(['disp(tabIn)']);
    % Use TeX Markup for bold formatting and underscores.
    TString = strrep(TString,'<strong>','\bf');
    TString = strrep(TString,'</strong>','\rm');
    TString = strrep(TString,'_','\_');
    % Get a fixed-width font.
    FixedWidth = get(0,'FixedWidthFontName');
    % Output the table using the annotation command.
    annotation(fTable,'Textbox','String',TString,'Interpreter','Tex',...
        'FontName',FixedWidth,'FontSize',8,'units',...
        'normalized','Position',[0 0 1 1]);
end

function [hScFigOut, hScPatch] = coloredBins(x, y, nbins, hScFig, hScAx)

%parse nbins
if ~isempty(nbins)
    nbins = 30;
end

%make x and y column vectors if they're not
if size(x,2) > size(x,1)
    x = x';
end
if size(y,2) > size(y,1)
    y = y';
end
assert(numel(x) == numel(y));

%combine x and y into a matrix for hist3
z = [x y];

%get number of points per bin using hist3
n = hist3(z,[nbins nbins]); % pass in vector of num bins in [X Y] dims
n1 = n'; % n is a 2-D matrix containing the number of points per bin
n1(size(n,1)+1, size(n,2)+1) = 0; %add in one zero entry

%make x, y grid having as many points as there are bins
xb = linspace(min(z(:,1)),max(z(:,1)),size(n,1)+1);
yb = linspace(min(z(:,2)),max(z(:,2)),size(n,1)+1); %shouldn't this be size(n,2)+1 if it's not square?

%parse figure and axes arguments
if isempty(hScFig)
    hScFigOut = figure();
else
    hScFigOut = figure(hScFig);
end
if ~isempty(hScAx)
    axes(hScAx);
end

%make histogram
hScPatch = pcolor(xb,yb,n1);
hScPatch.EdgeColor = 'none';
colormap(flipud(hot));

set(gca,'layer','top'); %this moves axis ticks, etc., above patch
set(hScFigOut,'Renderer','opengl'); %for proper export
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

function [groupName, numInGroup, groupRank, covarRank] = covarNameRank(s)
%Operations on covariate names: useful to order covariates by group, etc.

%given full covariate name (e.g., Age11), return what group it belongs to
%(e.g.,Age), the number of this covariate in the group (11),
%the rank of this group in the list of covariates below (Age is
%4th), and the rank of the covariate among all covariates (4*100 + 11 =
%411, where each covariate group is given 100).


group = {'GeoMeteo',...
              'Population',...
              'Sex',...
              'Age',...
              'Race','Household',...
             'Edu','Employed','Income','Commute',...
             'Occupation','Industry','WorkerClass',...
             'Health','Religion','Politics'};

%find the name of group this belongs to
ix = cellfun(@(x) compStr(s,x),group);
s
groupName = group{ix};

%find the rank of this group among all groups (according to list above)
groupRank = find(strcmp(group,groupName));

%split 'Age11' into {'','11'}
splitName = strsplit(s,groupName);
numInGroup = str2num(splitName{2}); %

%give 100 units to each covar group + whatever rank this covariate is in
%its group; e.g., so 'Age11' = 400 + 11 = 411.
covarRank = 100*groupRank + numInGroup;

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

%% map functions
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
    colorLims = prctile(ds.(colorBy),[1 99]);    %this is in minutes, so don't round to next hour as for Fig. 1
%     colorLims = [floor(prctile(ds.(colorBy),1)) ...
%         ceil(prctile(ds.(colorBy), 99))];
%    colorLims = [min(ds.(colorBy)) max(ds.(colorBy))];

    %2018-06-03, EL:
    %should manually trim the values to the [1,99] percentiles. it seems like
    %Matlab assigns them to the median value otherwise...
%     ds{ds.(colorBy) > colorLims(2),colorBy} = colorLims(2);
%     ds{ds.(colorBy) < colorLims(1),colorBy} = colorLims(1);
end
%colormap(jet(512));
mySymbolSpec = makesymbolspec('Polygon',...
    {'color',colorLims,... %trying this instead of cLims
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
barTicks = 0:0.25:3;
set(hBar,'Limits',[min(ds.(colorBy)) max(ds.(colorBy))],...
    'Ticks',barTicks','TickLabels',60*barTicks);
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