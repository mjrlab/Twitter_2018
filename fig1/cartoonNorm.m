%2018-06-27, EL:
%   - make Fig. 1A
%   - make a plot showing how individual users' tweeting histories add up
%   to a normalized Tweetogram

close all; clear all; clc;

TOEXP = 0;
OUTDIR = 'cartoonPics';

%% load the tweetogram for avg day in Chicago in 2012-2013
randStr = 'ckvc';

load('datafiles/CountyTableHaveData.mat');
chicago = loadTgram(geoIx('Cook IL',cTable), 'MTWThFSatSun', 2012:2013, ...
    'annual',{2012:2013},'FIPS',2012:2013, randStr);

%smooth
chicago = smoothify(chicago);
chicago = chicago(1:96);
chicago = chicago ./ sum(chicago);

%bin hourly
bChicago = binByHr(chicago, 8);

%get prob of tweeting in any time bin
[p,binTail] = binProb(bChicago);


%% make plots of individual users tweeting patterns
N_Users = 4;
myColors = cool(N_Users);
for n=1:N_Users
    z = 1+rand*10; % #tweets/day = uniform(1,11)
    numTweets = floor(z);
    [md{n}, mnd{n}] = onePerson(numTweets,binTail);
    myNormDist = md{n};
    maxval = max(cellfun(@(x) max(x(:,2)),md,'uniformoutput',true));
    %subplot(N_Users+1,1,n);
    fUsers(n)=figure('color','none');
    %have to make your own bar function b/c Matlab's sucks == widths keep
    %changing! just make patches -- got a workaround by filling empty
    %squaares with nan's
    twf = 1:24;
    notInTwf = twf(~ismember(1:24,myNormDist(:,1)));
    bar( [myNormDist(:,1); notInTwf'],...
         [myNormDist(:,2); nan(size(notInTwf'))],...
         0.5,'FaceColor',myColors(n,:),'ShowBaseline','off');
    set(gca,'xlim',[-1 13],...
        'xtick',0:2:12,'xticklabel',[],...
        'XAxisLocation','top',...
        'ylim',[-0.25 maxval + 0.25],'ytick',[],...
        'box','on','ticklength',[0.02 0.035]);
    if n > 0
        set(gca,'xticklabel',[]);
    end
    set(fUsers(n),'units','inches','position',[0 0 1.5 0.375]);
end

expif(TOEXP,fUsers,['inds'],...
            'OUTDIR',OUTDIR,'fun2use','expfig','ext','pdf');

%% plot summed tweets
[bigDist, bigDistAgg] = addManyPeople(mnd);
%subplot(N_Users+1,1,N_Users+1);
fSummed=figure('color','none');
b = bar(bigDistAgg',0.5,'stacked','FaceColor',...
                'flat','ShowBaseline','off');
for k=1:numel(b)
    b(k).CData = myColors(k,:);
end
set(gca,'xlim',[-1 13],...
    'xtick',0:2:12,'xticklabel',[],...
    'XAxisLocation','bottom',...
    'ylim',[0 0.1+max(sum(bigDistAgg,1))],...
    'ytick',[0:0.5:0.1+max(sum(bigDistAgg,1))],...
    'yticklabel',[],...
    'box','on','ticklength',[0.02 0.035]);
set(fSummed,'units','inches','position',[0 0 1.5 0.5]);

expif(TOEXP,fSummed,['aggregate'],...
            'OUTDIR',OUTDIR,'fun2use','expfig','ext','pdf');

%%
function tgram = loadTgram(geoCode, whichDay,whichBin,binType,binList,...
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
if isnumeric(whichBin)
    LINENUM = find(cellfun(@(x) isequal(whichBin,x), binList));
elseif ischar(whichBin)
    LINENUM = find(cellfun(@(x) strcmp(whichBin,x), binList));
end

%make filename
FILENAME = [CALCDIR '/' SUBDIR '/' num2str(geoCode) '.csv'];

%csv read indexes rows/columns from 0
tgram = csvread(FILENAME,LINENUM-1,0,[LINENUM-1 0 LINENUM-1 191]);

end

function sm = smoothify(x)
g=gausswin(10);
g=g/sum(g);
sm = conv(x,g,'same');
end

function b = binByHr(x,binWidth)
%assume data is in 15-min bins, bin into 1-hr bins

b=[];
for i=1:(96/binWidth)
    start = binWidth*(i-1) + 1;
    fin = binWidth*i;
    b(i) = sum(x(start:fin));
end

end

function [p,binEnd] = binProb(bins)
%given bins, calculate probabilities of being in every bin and the range of
%values this bin occupies: e.g., if the first pin is worth 5%, it's in the
%range of [0 5], if second bin is worth 10%, it's in [5 15]. binEnd = [5
%15 ...], recording the ends of the bins in order.

p = [];
for n=1:numel(bins)
    p(n) = bins(n);
end

binEnd(1) = p(1);
for n=2:numel(bins)
    binEnd(n) = binEnd(n-1)+p(n);
end

end

function [myDist, myNormDist] = onePerson(numTweets,binTail)
%simulate one user's Tweet's over 24 hours

myRand = rand(1,numTweets); %random numbers
myBinIx = [];
for r=1:numel(myRand)
    myBinIx(r) = find(binTail >= myRand(r), 1,'first');
end

uniBins = unique(myBinIx);
myDist = [];
for u=1:numel(uniBins)
    myDist(u,:) = [uniBins(u) sum(myBinIx == uniBins(u))];
end

myNormDist = myDist;
disp(myDist);
try
    myNormDist(:,2) = myDist(:,2)./sum(myDist(:,2));
catch err
    disp(myRand);
    disp(numTweets);
    disp(err);
end
end

function [bigDist, bigDistAgg] = addManyPeople(mnd)
%mnd = {}, one cell = one person's normalized tweets

%get largest timebin value
maxBin = max(cellfun(@(x) max(x(:,1)),mnd));

bigDistAgg = zeros(numel(maxBin),maxBin);

for i=1:numel(mnd)
   bigDistAgg(i,[mnd{i}(:,1)]') = [mnd{i}(:,2)]'; 
end

bigDist = mnd{1};
for i=2:numel(mnd)
    bigDist = [bigDist; mnd{i}];
end

% %reduce by column
% uniHr = unique(bigDist(:,1));
% for u=1:numel(uniHr)
%    ix = bigDist(:,1) == uniHr(u);
%    bigDistAgg(u,1:2) = [uniHr(u) sum(bigDist(ix,2))];
% end

end

function geoCode = geoIx(geoName,cTable)
if ischar(geoName)
    geoCode = cTable.geoCode(cellfun(@(x) isequal(x,geoName), cTable.geoName));
else
    geoCode = nan;
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
