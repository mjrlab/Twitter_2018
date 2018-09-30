%2018-08-10, EL:
%   - tweetogram trough time before and during school break?
%2018-06-27, EL: 
%    - Figure 5: plot weekly tweetogram trough times, overlay holidays
%                in the school calendar as well as public holidays

clear all; 
close all; 
clc;

%% parameters
TOEXP = 0; %export figures? (1=yes)
MARK_SCHOOL_HOLIDAYS = 0; %mark school holidays (1=yes)
MARK_PUBLIC_HOLIDAYS = 0; %mark federal/public holidays (1=yes)
TOEXP_SummerEffect = 1; %write table with avg. trough time before and during summer break?

%which county to run this analysis on?
myGeo = 'Virginia Beach VA'; %in text used: 'Orange FL', 'Henry GA' and 'Virginia Beach VA'

%% where are the raw tweetograms located?

%these variables determine the folder and file name with the data
whichDay = 'wkday';
whichBin = 1:52;
binType = 'weekly';
binList = num2cell(1:52);
geoRegion = 'FIPS';
whichYears = 2013;

%unique ID for calculations
randStr = 'kekh';

%file Tweetogram minima
INDIR = ['../tweetograms' '/' 'tweetograms' '_' binType '_' geoRegion ...
    num2str(whichYears,'_%d') '_' randStr];
SJLFILE = [INDIR '/' 'troughs_2018-03-26_16.29.30_kekh.mat'];
sjlTab = load(SJLFILE);
sjlTab = sjlTab.outTable;

%cluster assignments (for colors)
CLUST_file = load('datafiles/clustAssign_2018-08-11_19.51.01.mat'); %SJLtable_clustAssign_2018-04-20_11.25.24.mat');
sjlTab_clust = CLUST_file.SJLtable_clust;
sjlTab_clust = sjlTab_clust(:,{'FIPS','hierClust'}); %only keep FIPS and hierClust index
sjlTab_clust.Properties.VariableNames{'FIPS'} = 'geoCode';

% calendar data
CALFILE = ['datafiles/CountyTable_SchoolCalendar_10Aug2018.mat'];
calTab = load(CALFILE);
calTab = calTab.outTab;

% add columns to calTab to capture the first full week during break and
% last full school week before break
lfsw = rowfun(@(x) lastFullSchoolWeek(x), calTab(:,'lastDay1213'),...
                            'OutputVariableNames',{'lfsw1213','ffbw1213'});
calTab = [calTab lfsw];
                        
% basic county info
cTable = load('datafiles/CountyTableHaveData.mat');
cTable = cTable.cTable;

% merge tables
sjlTab = join(sjlTab,cTable,'Keys','geoCode');
sjlTab = innerjoin(sjlTab,calTab); %add school calendar
sjlTab = join(sjlTab,sjlTab_clust,'Keys','geoCode');

%table to look at
simTab = sjlTab(:,{'geoCode','geoName','state','avgUsers','hierClust'});

%colors of clusters
colors = hsv(max(sjlTab.hierClust));

%% compute avg trough time in 3 wks before and 3 wks after school start
sjlTab{:,'NumWeeks'} = 3;
for f=1:numel(sjlTab{:,1})
    sjlTab.inSummer(f) = ...
        avgOver(sjlTab{f,'ffbw1213'},... %first full break week
        sjlTab{f,'NumWeeks'},...
        sjlTab{f,'wkday_trPos'});
     sjlTab.nextSummer(f) = ...
        avgOver(sjlTab{f,'ffbw1213'}+sjlTab{f,'NumWeeks'},... %next block of weeks in summer
        sjlTab{f,'NumWeeks'},...
        sjlTab{f,'wkday_trPos'});
    sjlTab.preSummer(f) = ...
        avgOver(sjlTab{f,'lfsw1213'}-sjlTab{f,'NumWeeks'},... %last full school week
        sjlTab{f,'NumWeeks'},...
        sjlTab{f,'wkday_trPos'});
end

%export table with summer effect calculations?
if TOEXP_SummerEffect
   summerEffectTab = sjlTab(:,{'geoCode',...
                            'geoName',...
                            'hierClust',...
                            'preSummer',...
                            'inSummer',...
                            'nextSummer'});
   summerEffectTab.summerEffect = summerEffectTab.inSummer - ...
                                  summerEffectTab.preSummer;
   summerEffectTab.nextEffect = summerEffectTab.nextSummer - ...
                                  summerEffectTab.inSummer;
   summerEffectTab.Properties.VariableNames = ...
                            {'FIPS',...
                            'name',...
                            'cluster_num',...
                            'avgWkdayTrTime_3wks_beforeSummerBreak_hr_am',...
                            'avgWkdayTrTime_3wks_intoSummerBreak_hr_am',...
                            'avgWkdayTrTime_3wks_nextpartofSummerBreak_hr_am',...
                            'summerShift_hr',...
                            'nextShift_hr'}; 
   writetable(summerEffectTab,['datafiles/summerEffect_081218.xlsx']);
   disp(mean(summerEffectTab.summerShift_hr)*60);
   disp(mean(summerEffectTab.nextShift_hr)*60);
   [h, p] = ttest(summerEffectTab.nextShift_hr,summerEffectTab.summerShift_hr)
end

%% Fig. 5: plot weekly Tweetograms as heatmaps and trace minima
[myGeoCode,myGeoIx] = geoIx(myGeo,sjlTab);

% where to save data
OUTDIR = 'schoolCalendar_fig5';

% Fig. 5B-D: positions of tweetogram minima week-by-week
% along with public and school holidays 
[hImTroughs,~] = plotTroughsCalendar(myGeo,sjlTab,...
    'color',colors(sjlTab.hierClust(myGeoIx),:),...
    'holidays',MARK_PUBLIC_HOLIDAYS,...
    'schoolBars',MARK_SCHOOL_HOLIDAYS,...
    'fontsize',6,'axfontsize',8,...
    'titStr','',...
    'ylab',{'weekday tweetogram','trough time (hr a.m.)'},...
    'figsize',[0 0 7 2]);

tgram = loadTgram(geoIx(myGeo,sjlTab), 'wkday', 1:52, ...
    'weekly',num2cell(1:52),'FIPS',2013,randStr);

% Fig. 5A: 24-hr weekly tweetograms stacked as heat maps 
[fImHeatMap,hAx] = heatMapTgram(tgram,myGeo,sjlTab{myGeoIx,'wkday_trPos'},...
    'ylim',[0 53],'figSize',[0 0 7 1.5],'titStr','');
view([-90 90]);

%% export 
expif(TOEXP,hImTroughs,['clust' num2str(sjlTab.hierClust(myGeoIx)) '_' ... 
                        strrep(myGeo,' ','') '_trPos_holidays'],...
                        'OUTDIR',OUTDIR);                  
expif(TOEXP,fImHeatMap,['clust' num2str(sjlTab.hierClust(myGeoIx)) '_' ...
                        strrep(myGeo,' ','') '_heatMapTgram'],...
                        'OUTDIR',OUTDIR);

%% helper funs
function [hIm,hAx] = plotTroughsCalendar(geoCode,sjlTab,varargin)
%plot 2013 notable events

myParser = inputParser;
addParameter(myParser,'holidays',1,@isnumeric);
addParameter(myParser,'plotVar','wkday_trPos',@ischar);
addParameter(myParser,'ylab','trough time (hr a.m.)');
addParameter(myParser,'titStr','',@ischar);
addParameter(myParser,'school',0,@isnumeric);
addParameter(myParser,'schoolBars',1,@isnumeric);
addParameter(myParser,'fontsize',8,@isnumeric);
addParameter(myParser,'axfontsize',10,@isnumeric);
addParameter(myParser,'figsize',[0 0 8 3],@isnumeric);
addParameter(myParser,'color',[0    0.4470    0.7410]);
parse(myParser,varargin{:});
args = myParser.Results;

%allow to pass geoName instead of geoCode
if ischar(geoCode)
    geoix = find(cellfun(@(x) strcmp(x,geoCode),sjlTab.geoName));
else
    geoix = sjlTab.geoCode == geoCode;
end

%holidays in 2013
[school_vars, school_lab, ...
          school_vars_to_mark, ...
          holidays, holidays_lab] = holidays2013(sjlTab.geoName{geoix});

%need their weekly indices
wk_school = cellfun(@(x) ['wk_' x],school_vars,'UniformOutput',0);
wk_holidays = cellfun(@(x) ['wk_' x],holidays,'UniformOutput',0);

%plot troughs
hIm = figure('units','inches','position',args.figsize);
tr = sjlTab{geoix,args.plotVar}; %which variable to plot wk-by-wk
%tr = sjlTab.wkday_trPos(geoix,:);
plot(1:numel(tr),tr,'o-','linewidth',1,'markersize',2,...
    'color',args.color,...
    'markerfacecolor',args.color,...
    'markeredgecolor','none');
%title(sjlTab.geoName(geoix));
hAx = gca;
tickLen = get(gca,'ticklength');
set(gca,'xlim',[0 53],'xtick',1:4:53,... % 'ylim',[min(tr) max(tr)] + [-1 0.75]
    'ylim',[2.6 7],'fontsize',args.axfontsize,...
    'ticklength',tickLen*0.5);
xlabel('weeks of the year');
ylabel(args.ylab);

%mark school holidays below curve using arrows
if args.school == 1
    wk_pos = sjlTab{geoix,wk_school};
    for a=1:numel(wk_pos)
        if wk_pos(a) < 3 || wk_pos(a) > 45  
            xAlign = 'center';
        else
            xAlign = 'left';
        end
        
        if ~isnan(wk_pos(a))
            y = min(tr) - 0.05;
            x = wk_pos(a);
            arrow([x y-0.15],[x y],...
                'Length',2,'tipangle',45,'baseangle',90);
            if ismember(school_vars{a},school_vars_to_mark)
                text(x,min(tr)-0.3,school_lab{a},'fontsize',args.fontsize,...
                    'HorizontalAlignment',xAlign,'clipping','on',...
                    'VerticalAlignment','top');
            end
        end
    end
end

%mark school holidays below curve using arrows
if args.schoolBars == 1
    
    %mark first, end weeks
    sjlTab.wk_zero(:) = 0;
    sjlTab.wk_53(:) = 53;
    
    wkintervals = {{'wk_zero','wk_winterBreakEnd1213'},...
                 {'wk_febBreakStart1213','wk_febBreakEnd1213'},...
                 {'wk_springBreakStart1213','wk_springBreakEnd1213'},...
                 {'wk_lastDay1213','wk_firstDay1314'},...
                 {'wk_fallBreakStart1314','wk_fallBreakEnd1314'},...
                 {'wk_winterBreakStart1314','wk_53'}};
    for ii=1:numel(wkintervals)
        wkstart = sjlTab{geoix,wkintervals{ii}{1}};
        wkend = sjlTab{geoix,wkintervals{ii}{2}};
        %find if either isnan
        nanEntry = sum(isnan([wkstart wkend]));
        if nanEntry > 0
            continue;
        else
            y = min(tr) - 0.05;
            rectangle('position',[wkstart y-0.25 (wkend-wkstart) 0.25],...
                'Curvature',1,'facecolor','y','edgecolor','k'); %[248 222 126]/255
        end
        
    end
    
end

%mark other holidays above curve
if args.holidays == 1
    wk_pos = sjlTab{geoix,wk_holidays};
    for a=1:numel(wk_pos)
        %if isequal(holidays_lab{a}, {'Labor','Day'}) | ...
        if  isequal(holidays_lab{a}, {'Memorial','Day'}) | ...
            isequal(holidays_lab{a}, {'Columbus','Day'})  | ...
            isequal(holidays_lab{a}, {'President''s','Day'})
            
            y = max(tr)-0.75;
            texty = y + 0.45;
            xAlign = 'center';
        elseif isequal(holidays_lab{a},{'Thanks-','giving'})
            y = max(tr)-0.75;
            texty = y + 0.45;
            xAlign = 'center';
        else
            y = max(tr);
            texty = y + 0.3;
            xAlign = 'center';
        end
        
        x = wk_pos(a);
        arrow([x y+0.2],[x y+0.05],...
            'Length',2,'tipangle',45,'baseangle',90);
        
        text(x,texty,holidays_lab{a},'fontsize',args.fontsize,...
            'HorizontalAlignment',xAlign,'clipping','on',...
            'VerticalAlignment','bottom');
    end
end

set(gca,'layer','top'); %this moves axis ticks, etc., above patch

end

function [school_vars, school_lab, ...
          school_vars_to_mark, ...
          holidays, holidays_lab] = ...
                                    holidays2013(geoName)

%events occurring in 2013
school_vars = {'winterBreakEnd1213',...
               'febBreakStart1213',...
               'febBreakEnd1213',...
               'springBreakStart1213',...
               'springBreakEnd1213',...
               'lastDay1213',...
               'firstDay1314',...
               'fallBreakStart1314',...
               'fallBreakEnd1314',...
               'winterBreakStart1314'};
           
school_lab = { {'winter', 'break'},...
               {'break'},...
               {'break'},...
               {'spring', 'break'},...
               {'spring', 'break'},...
               {'last', 'day'},...
               {'first', 'day'},...
               {'fall', 'break'},...
               {'fall', 'break'},...
               {'winter', 'break'}};

school_vars_to_mark = {'winterBreakEnd1213',...
                       'febBreakEnd1213',...
                       'springBreakStart1213',...
                       'lastDay1213',...
                       'firstDay1314',...
                       'fallBreakStart1314',...
                       'winterBreakStart1314'};                  
                   
holidays = {'thanks13',...
            'ind13',...
            'dstStart13',...
            'dstEnd13',...
            'memorial13',...
            'labor13',...
            'columbus13',...
            'mlk13',...
            'pres13'};
        
holidays_lab = {{'Thanks-','giving'},...
                '4th of July',...
                'DST start',...
                'DST end',...
               {'Memorial','Day'},...
               {'Labor Day'},...
               {'Columbus','Day'},...
               {'MLK Day'},...
               {'President''s', 'Day'}};

if ~strcmp(geoName,'Henry GA')
    school_vars(2) = [];
    school_lab(2) = [];
    %school_vars_to_mark{2} = 'febBreakStart1213';
end
           
end

function [hIm,hAx] = heatMapTgram(tgram,titStr,troughs,varargin)
    myParser = inputParser;
    addParameter(myParser,'ylim',[0 53],@isnumeric);
    addParameter(myParser,'titStr','',@ischar);
    addParameter(myParser,'fontSize',8);
    addParameter(myParser,'markTrough',1);
    addParameter(myParser,'figSize',[0 0 7 2]);
    addParameter(myParser,'figColor','w');
    parse(myParser,varargin{:});
    args = myParser.Results;

    %double plot the actograms
    hIm = figure('units','inches','position',args.figSize,...
                                  'color',args.figColor);
    colormap(redbluecmap);
    imagesc(tgram);
    hold on;
    if args.markTrough
        plot(4*troughs,1:numel(troughs),':','linewidth',1,...
            'marker','none',...
            'color',[0.7 0.7 0.7]);
    end

    xticklist = [0:4*4:96 96+(16:4*4:96)];
    xticklab  = [0:4:24   (4:4:24)];

    tickLen = get(gca,'ticklength');
    
    set(gca,'xlim',[0 24*4], 'ylim', myParser.Results.ylim,...
        'xtick',xticklist,'xticklabel',xticklab,...
        'yticklabel',(1:4:2*52),...
        'ytick',1:4:2*52,'fontsize',args.fontSize,'ticklength',tickLen*0.5);
    xlabel('local time (hr a.m.)');
    ylabel('weeks of the year');
    
    hAx = gca;
    %title(titStr);
end

function [lfsw, ffbw] = lastFullSchoolWeek(lastDay)
%return lfsw (last full school week number) and ffbw (first full break
%week) number given the last day of school

if weekday(lastDay) == 1 || ... % Sun 
   weekday(lastDay) == 7        % Sat

   disp('warning: last day of school year fell on a weekend!');
   
end

%if school year ends on Friday, then this week is the last full school week
if weekday(lastDay) == 6
    lfsw = week(lastDay);
    ffbw = week(lastDay)+1;

%otherwise, this week is both break and non-break;
%rerturn week before this one and week after this one
else
    lfsw = week(lastDay)-1;
    ffbw = week(lastDay)+1;
end

end

function aw = avgOver(firstWk,NumWeeks,weekData)
%weekdata = [1x52] vec with tweetogram trough times; return the average of
%entries firstWk:firstWk+NumWeeks-1

aw = mean(weekData(firstWk:firstWk+NumWeeks-1));

end

function awfips = avgFIPSOver(geoCode,sjlTab,firstWk,NumWeeks)

%find this county
if ischar(geoCode)
    geoix = find(cellfun(@(x) strcmp(x,geoCode),sjlTab.geoName));
else
    geoix = sjlTab.geoCode == geoCode;
end

%get troughs
tr = sjlTab{geoix,'wkday_trPos'};

%get avg time
awfips = avgOver(firstWk,NumWeeks,tr);

end

function tgram = loadTgram(geoCode, whichDay, whichBin, binType, binList,...
                            geoRegion,whichYears,randStr)

CALCDIR = ['../tweetograms/'...
           'tweetograms' '_' binType '_' geoRegion ...
                num2str(whichYears,'_%d') '_' randStr '/' binType];

disp(['geoCode ' num2str(geoCode)]);            
            
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
tgram = csvread(FILENAME,LINENUM(1)-1,0,[LINENUM(1)-1 0 LINENUM(end)-1 191]);

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

function xout = smoothify(xin,varargin)
%smooth array using a Gaussian smoothing kernel

%optional parameter to pass smoothing window size
p=inputParser;
addParameter(p,'winsize',10,@isnumeric);
parse(p,varargin{:});
winsize = p.Results.winsize;

g=gausswin(winsize);
g=g/sum(g);
smoother = (@(x) conv(x,g,'same')); %2017-10-30, EL: until today had 'smooth', but this doesn't work after update
xout = smoother(xin); %smooth array, works
end

function expif(TOEXP,fH,fName,varargin)
myParser = inputParser;
addParameter(myParser,'OUTDIR','.',@ischar);
addParameter(myParser,'ext','pdf',@ischar);
addParameter(myParser,'fun2use','expfig',@ischar);
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
