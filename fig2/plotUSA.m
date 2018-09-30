function [hFig, hAx] = plotUSA(figIn)
%Plot figure of cont. US with Alaska and Hawaii below and states outlined.
%   Return hFig = figure handle; hAx = vector of three axes: 1 =
%   Continental US, 2 = Alaska, 3 = Hawaii.
%   Can pass in figure handle figIn to plot onto.

%make new figure to use passed in figure
if nargin > 0
    hFig = figure(figIn);
else
    hFig = figure();
end

%% matlab example for plotting US map in monochrome

ax = usamap('all');
set(ax, 'Visible', 'off')
states = shaperead('usastatelo', 'UseGeoCoords', true);
names = {states.Name};
indexHawaii = strcmp('Hawaii',names);
indexAlaska = strcmp('Alaska',names);
indexConus = 1:numel(states);
indexConus(indexHawaii|indexAlaska) = []; 
stateColor = [0.9 0.9 0.9];

geoshow(ax(1), states(indexConus),  'FaceColor', stateColor)
geoshow(ax(2), states(indexAlaska), 'FaceColor', stateColor)
geoshow(ax(3), states(indexHawaii), 'FaceColor', stateColor)
hold on;

%adjust axes to remove grid and longtitute/latitude labels
for k = 1:3
    setm(ax(k), 'Frame', 'off', 'Grid', 'off',...
      'ParallelLabel', 'off', 'MeridianLabel', 'off')
end

% set(ax(2),'Position',[0.24  0.24 0.3  0.3]);
% set(ax(3),'Position',[0.33  0.28 0.25  0.25]);
% set(hFig,'units','inches','position',[0 0 12 6]);

%2018-01-25, EL
%   - modified to make maps better-looking
set(ax(1),'XLim',[-2200000 2600000],'YLim',[2500000 5500000]);
set(ax(1),'Position',[0  0 1  1]);
set(ax(2),'Position',[0.01  0.04 0.3  0.3]);
set(ax(3),'Position',[0.18  0.04 0.3  0.3]);
set(hFig,'units','inches','position',[0 0 8 6]);

%tightmap;
hAx = ax;
    
end

