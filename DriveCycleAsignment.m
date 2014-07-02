function [DriveCycleCode]=DriveCycleAsignment(avgspd,dist)
%This function uses the distance and average speed of a trip to assign a
%drivecycle type for that trip

drivecycles = {'UDDS' 'HWFET' 'US06';...
    19.6 48.3 48.3}; %These are the three drivecycles we are currently considering and their respective average speeds according to the EPA

%This very simply selects the drivecycle by it's average speed.  If it is
%less that UDDS, then choose UDDS, if between UDDS and HWFET, then HWFET,
%if greater than HWFET, choose US06
if avgspd <= drivecycles{2,1}
    DriveCycleCode = drivecycles{1,1};
elseif avgspd <= drivecycles{2,2} & avgspd > drivecycles{2,1}
    DriveCycleCode = drivecycles{1,2};
else
    DriveCycleCode = drivecycles{1,3};
end

end