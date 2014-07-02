function [ActivitySheet]=V2G_ActivitySheetGen(TripMatrix, chargingcodes)

ActivitySheet=cell(1,11);
ActivitySheet(1,:)={'Vehicle ID', 'State', 'Start time (hour)', ...
    'End time (hour)', 'Distance (mi)', 'Drive cycle', 'P_max (W)', ...
    'Location', 'NHTS HH Wt' , 'HHState' , 'CBSA' };
VIDs=unique(TripMatrix(:,1));
m=1;
while m<=length(VIDs)
    VehTrips=TripMatrix(VIDs(m)==TripMatrix(:,1),:);
    [Activitycell]=V2GSIM_dailytrips(VehTrips, chargingcodes);
    ActivitySheet(end+1:end+size(Activitycell,1),:)=Activitycell;
    m=m+1;
end
end