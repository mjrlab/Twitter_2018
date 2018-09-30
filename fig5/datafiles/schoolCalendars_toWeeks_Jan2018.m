%2018-02-27, EL:
%   - turn Excel table with school calendars into Matlab table, convert
%   dates to week numbers

clear all; clc;

TOSAVE = 0;

%% load Excel table with calendar info
tab = readtable('CountyTable_UnifiedSchoolDistrict_10Aug2018.xlsx');
tab = tab(~isnat(tab.firstDay1213),:);
% pickCounties = cellfun(@(x) ...
%     ismember(x,{'Orange FL','Virginia Beach VA', 'Henry GA'}),...
%     tab.geoName);
% tab = tab(pickCounties,:);

%% check that no calendar dates are weekends
calendarColumns = tab.Properties.VariableNames(15:end);

for i=1:numel(tab(:,1))
    wkends = weekday(tab{i,calendarColumns}) == 1 | ... %sun
             weekday(tab{i,calendarColumns}) == 7;      %sat
    if sum(wkends > 0)
        disp([tab.state{i} ' ' num2str(tab.geoCode(i)) ' ' ...
            tab.geoName{i} ' incorrect dates:']);
        disp([char(9) calendarColumns{wkends}]);
    end
end

%% add holidays and daylight savings dates 
hol.thanks13 = '28-Nov-2013';
hol.thanks14 = '27-Nov-2013';
hol.labor13 = '2-Sep-2013';
hol.labor14 = '1-Sep-2013';
hol.columbus13 = '14-Oct-2013';
hol.columbus14 = '13-Oct-2014';
hol.mlk13 = '21-Jan-13';
hol.mlk14 = '20-Jan-14';
hol.pres13 = '18-Feb-13';
hol.pres14 = '17-Feb-14';
hol.memorial13 = '27-May-2013';
hol.memorial14 = '26-May-2014';
hol.ind13 = '4-July-2013';
hol.ind14 = '4-July-2014';
hol.dstStart13 = '11-Mar-2013'; %first Mon in DST
hol.dstEnd13 = '4-Nov-2013';    %first Mon
hol.dstStart14 = '10-Mar-2014'; %first Mon in DST
hol.dstEnd14 = '3-Nov-14';      %first Mon

holidays = fieldnames(hol);

for h=1:numel(holidays)
    tab(:,holidays{h}) = {datetime(hol.(holidays{h}))};
end

%% convert calendar dates to week numbers
calendarColumns = tab.Properties.VariableNames(15:end);
for c=1:numel(calendarColumns)
    newColName = ['wk_' calendarColumns{c}];
    tab.(newColName) = week(tab.(calendarColumns{c}));
end

%% save
if TOSAVE
    outTab = tab;
    save('CountyTable_SchoolCalendar_10Aug2018.mat','outTab');
end