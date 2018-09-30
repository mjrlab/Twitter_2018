function plotLRhist(hFig, hAx,y1,y2,col1,col2, meanMedian, meanCol1, meanCol2)
%Plot two vertical distributions size-by-side
%   Detailed explanation goes here

%%
TOTEST = 0;
if TOTEST == 1
    y1 = randn(1,1000);
    y2 = y1+1;
    
    hFig= figure();
    hAx1 = axes;
    hAx2 = axes;
    
    col1 = 'b';
    col2 = 'b';
    
    h1 = histogram(hAx1,y2,'DisplayStyle','bar',...
        'facecolor',col2,'facealpha',1,'edgecolor',col2,...
        'Orientation','horizontal','Normalization','probability');
    hold on;
    h2 = histogram(hAx2,y1,'DisplayStyle','stairs','edgecolor',col1,...
        'Orientation','horizontal','Normalization','probability');
    fixAxes;
end

figure(hFig);

hAx1 = hAx;
axes(hAx1);

h1 = histogram(hAx1,y1,'NumBins',20,'DisplayStyle','bar',...
    'facecolor',col1,'facealpha',1,'edgecolor',col1,...
    'Orientation','horizontal','Normalization','probability');
hold on;

if meanMedian == 1
    plot(0.9*[0 max(h1.Values)], [nanmean(y1) nanmean(y1)], ':',...
        'color',meanCol1,'linewidth',1.5);
elseif meanMedian == 0
    plot(0.9*[0 max(h1.Values)], [nanmedian(y1) nanmedian(y1)], ':',...
        'color',meanCol1,'linewidth',1.5);
end

hAx2 = axes;
axes(hAx2);
h2 = histogram(hAx2,y2,'NumBins',20,'DisplayStyle','stairs','edgecolor',col2,...
    'Orientation','horizontal','Normalization','probability');
hold on;

if meanMedian == 1
    plot(0.9*[0 max(h2.Values)], [nanmean(y2) nanmean(y2)], ':',...
        'color',meanCol2,'linewidth',1.5);
elseif meanMedian == 0
    plot(0.9*[0 max(h2.Values)], [nanmedian(y2) nanmedian(y2)], ':',...
        'color',meanCol2,'linewidth',1.5);
end


fixAxes;
%fixAxes;

    function fixAxes
        %helper function to make the axes have corresponding x,y limits
        %remove x ticks and labels off the second (phantom) axes and
        %reverse the direction of one of them
        
        hAx2.Position = hAx1.Position;
        hAx2.Visible = 'off';
        hAx2.XTick = [];
        hAx2.YTick = [];
        
        hAx2.XDir = 'reverse';
        
        hAx2.XLim = findXLims;
        hAx2.YLim = findYLims;
        
        hAx1.XLim = findXLims;
        hAx1.YLim = findYLims;
    end

    function lims = findYLims
        %helper function to find Y limits of the axes;
        %set to be 10% lower than the min value and 10% greater than the
        %max value of the greatest bin in each histogram
        
        lo = min(h1.BinLimits(1), h2.BinLimits(1));
        hi = max(h1.BinLimits(2), h2.BinLimits(2));
        dist = hi-lo;
        lo = lo - dist*0.1;
        hi = hi + dist*0.1;
        
        lims = [lo hi];
    end

    function lim = findXLims
        %helper function to find X limits;
        %set to 110% of max value of either histogram
        
        lim1 = max(h1.Values);
        lim2 = max(h2.Values);
        lim = max(lim1,lim2)*1.1;
        lim = [-lim lim];
    end

end

