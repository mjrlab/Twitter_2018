clFile = readtable('clusters_activities_2018-08-06_22.50.56.xlsx');
SJLtable_clust = clFile(:,{'FIPS','hierClust'});
save(['clustAssign_' getDate() '.mat'],'SJLtable_clust');