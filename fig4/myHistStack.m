%2018-03-19, EL:
%   - stack histograms on top of each other (or left-to-right), in different colors

function hFig = myHistStack(data,varargin)
%data = cell array of multiple datasets, each to be plotted as a histogram

%parse arguments
myParser = inputParser;
%addRequired(myParser,'data'); %,@(x)iscell(x));
addParameter(myParser,'colors',parula(numel(data)));
addParameter(myParser,'numBins',20);
addParameter(myParser,'markMedian',1);
addParameter(myParser,'dataLabels','');
addParameter(myParser,'varLabel','');
addParameter(myParser,'groupLabel','');
addParameter(myParser,'pVal',nan);
addParameter(myParser,'toLegend',0);
addParameter(myParser,'titStr','covariate');
addParameter(myParser,'figSize',[0 0 1.8 1.8]);
addParameter(myParser,'fontSize',8);

parse(myParser,varargin{:});
args = myParser.Results;

nd = numel(data);
myColors = args.colors;

if strcmp(args.dataLabels,'')
    for i=1:nd
        dataLabels{i} = ['data ' num2str(i)];
    end
else
    dataLabels = args.dataLabels;
end

%first, get all bin counts --> need this to set axis limits correctly
for i=1:nd
   x = data{i};
   [N{i}, edges{i}] = histcounts(x, 'Normalization', 'probability',...
       'numBins', args.numBins);
end

%get min/max bin limits and max bin value
Nmax = max(cellfun(@max,N));
edgeMin = min(cellfun(@min,edges));
edgeMax = max(cellfun(@max,edges));
e2e = edgeMax - edgeMin;

%make plots
hFig = figure('units','inches','position',args.figSize,'color','w');
hAx = axes;
for i=1:nd
   [xx,yy] = getHistXY(N{i},edges{i});
  
   yy = yy*0.35/max(yy);
   hP(i) = fill([xx fliplr(xx)],[yy -fliplr(yy)]+i,myColors(i,:));
   hold on;
   
   %mark medians with dotted lines
   if args.markMedian
        medVal = median(data{i});
        medIx = find(xx > medVal,1,'first');
        xMed = xx(medIx-1);
        yMed = yy(medIx-1);
        line([xMed xMed],[-yMed yMed]+i,...
            'linestyle','-.','color','k','linewidth',2);
   end
end

set(hAx,...
    'xlim',[edgeMin edgeMax] + 0.1*[-e2e e2e],...
    'ylim',[0.5 nd+0.5],...
    'ytick',1:nd,...
    'yticklabels', dataLabels,...
    'fontsize',args.fontSize);

xlabel(args.varLabel);
ylabel(args.groupLabel);
title(args.titStr);

%add legend with the p value for the comparison
if args.toLegend
    legTxt = ['p=' num2str(args.pVal,'%1.0e')];
    xLim = get(hAx,'xlim');
    yLim = get(hAx,'ylim');
    xTxt = xLim(1) + 0.25*(xLim(2) - xLim(1));
    yTxt = yLim(1) + 0.65*(yLim(2) - yLim(1));
    text(xTxt,yTxt,legTxt,'FontSize',args.fontSize,...
        'VerticalAlignment','top');
end

view([90 -90]);

end

function [xx,yy] = getHistXY(N,edges)
%histcounts() returns:
%   - edges = left-most bin edge, mid-points of each bin and
%   rightmost edge; 
%   - N = counts in each bin
%
%change these to return the left and right point of every bin, including 0s
%at the edges

x = edges;
y = N;

%get all the big edges (these already have leftmost and rightmost values)
for i=1:numel(x)
    xx([2*(i-1)+1 2*i]) = x(i);
end

%get corresponding y values
for i=1:numel(y)
    yy([2*(i-1)+1 2*i]) = y(i);
end
yy = [0 yy 0]; %add zeroes at left and right


end