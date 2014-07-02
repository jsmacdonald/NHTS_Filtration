function [Activitycell]=V2GSIM_dailytrips(VehTrips, chargingcodes)
% This creates a cell array for a single vehicle's daily usage pattern.
% The ActivityCell will be written to make up the Activity worksheet to be
% input into V2GSIM.

%chargingcodes is a cell that has the allowable codes for charging and
%their maximum charging level at those locations.  The codes are in the
%first row of the cell. The second contains a matrix that indicates the
%probability of having a particular charging power available the first
%column is the probability (all sum to one), the second is the charging
%power of each probability in Watts).  This information will be passed ot
%the LocationDesignation function to determine the charging level.

% Activitycell is a cell that is made up of values readable in the V2GSIM
% format that will be written into an excel spreadsheet. These are
% representative of a single vehicle's daily travel.  They format of the
% rows in the cell is:
% [ Vehicle ID | State | Start time (hour)| End time (hour)| Distance (mi) |
% Drive cycle | P_max (W) | Location | NHTS HH Wt | HHState | CBSA ]
% Vehicle ID is an ID number
% State is either Charging/Parked/Driving
% Start time is the time at the start of the period in the current state hour of the day from 0 to 23.99
% End time is the end time at the current state (Start(state2)=End(state1))
% Distance is the distance traveled while driving, or "-1" if parked or
% charging
% Drive cycle is the drive cycle used for the trip, or "-1" if parked or
% charging
% P_max is the maximum power available through the charging station when
% state is charging, and "-1" for other states.
% Location is a summary phrase describing the location that a vehicle is
% parked in
% NHTS HH Wt is the weighting factor from the NHTS that the household is
% meant to represent.  These weighting factors were created to get a
% representative national sample and should be used as such.
% HHState: The state of the US the household is located in.
% CBSA: Core Based Statistical Area

% VehTrips comes in with the information in this form:
% [ NEWVEHID | STRTTIME | ENDTIME | TRPMILES | WHYTO | WHYFROM | WTHHFIN |
% WTTRDFIN | HHSTATE | HH_CBSA ]

%% Getting Started
%Check to make sure there is only one vehicle in the data
VehicleIDs=unique(VehTrips(:,1));
if length(VehicleIDs)>1
    disp('Too many vehicles in the VehTrips matrix as input into V2GSIM_dailytrips function')
    return
end

%% Initialize output
%Initialize the output and determine the first row of information using
%WHYFROM variable unless it is driving.


Activitycell = cell(1,11);
M=1; %Activitycell row counter

%% Parked befor the trips start
if VehTrips(1,2)>0 %Make sure the vehicle isnt starting the day on a trip
    Activitycell{M,1} = VehTrips(1,1); %VehicleD
    [LocActivity, charge_power, where] = LocationDesignation(VehTrips(1,6),chargingcodes);
    Activitycell{M,2} = LocActivity; %Current State
    Activitycell{M,3} = 0; %Start Time at State
    Activitycell{M,4} = VehTrips(1,2)/60; %End Time at State.  In this
    % case, it is the start time of first trip.  The units for time in the
    % VehTrips file is minutes, so this must be converted to hours.
    Activitycell{M,5} = -1; %Trip Distance
    Activitycell{M,6} = -1; %Drive cycle
    Activitycell{M,7} = charge_power; %max power available at the charging station
    Activitycell{M,8} = where; %Location
    Activitycell{M,9} = VehTrips(1,7); %The NHTS household weight for the trip
    Activitycell{M,10} = VehTrips(1,9); %The household state location
    Activitycell{M,11} = VehTrips(1,10); %The household CBSA location
    M=M+1;
end
% [ Vehicle ID | State | Start time (hour)| End time (hour)| Distance (mi) |
% Drive cycle | P_max (W) | Location | NHTS HH Wt]


%% Fill in the rest of the table
n=1; %NHTS row counter
while n<=size(VehTrips,1)
    %First while the vehicle is actually on a trip.
    Activitycell{M,1} = VehTrips(n,1); %VehicleD
    [LocActivity, charge_power, where] = LocationDesignation([],chargingcodes);
    Activitycell{M,2} = LocActivity; %Current State
    Activitycell{M,3} = VehTrips(n,2)/60; %Start Time at State
    Activitycell{M,4} = min(1439,VehTrips(n,3))/60; %End Time at State.
    Activitycell{M,5} = VehTrips(n,4); %Trip Distance
    Tripspeed = VehTrips(n,4)/(VehTrips(n,3)-VehTrips(n,2))*60; %Trip Speed in milesperhour
    [DriveCycleCode]=DriveCycleAsignment(Tripspeed,VehTrips(n,4));
    Activitycell{M,6} = DriveCycleCode; %Drive cycle
    Activitycell{M,7} = charge_power; %max power available at the charging station
    Activitycell{M,8} = where;
    Activitycell{M,9} = VehTrips(n,7); %The NHTS household weight for the trip
    Activitycell{M,10} = VehTrips(n,9); %The household state location
    Activitycell{M,11} = VehTrips(n,10); %The household CBSA location
    M=M+1;
    
    %second, while the vehicle is parked between trips
    if n<size(VehTrips,1)
        if VehTrips(n+1,2)>VehTrips(n,3) %make sure the next trip is further in time than the current trip
            Activitycell{M,1} = VehTrips(n,1); %VehicleD
            [LocActivity, charge_power, where] = LocationDesignation(VehTrips(n,5),chargingcodes);
            Activitycell{M,2} = LocActivity; %Current State
            Activitycell{M,3} = VehTrips(n,3)/60; %Start Time at State
            Activitycell{M,4} = VehTrips(n+1,2)/60; %End Time at State.  In this
            % case, it is the last minute of the day.  The units for time in the
            % VehTrips file is minutes, so this must be converted to hours.
            Activitycell{M,5} = -1; %Trip Distance    Activitycell{M,6} = -1; %Drive cycle
            Activitycell{M,6} = -1; %DriveCycleCode
            Activitycell{M,7} = charge_power; %max power available at the charging station
            Activitycell{M,8} = where;
            Activitycell{M,9} = VehTrips(n,7); %The NHTS household weight for the trip
            Activitycell{M,10} = VehTrips(n,9); %The household state location
            Activitycell{M,11} = VehTrips(n,10); %The household CBSA location
            M=M+1;
        end
    elseif n==size(VehTrips,1) & VehTrips(n,3)<1439
        Activitycell{M,1} = VehTrips(n,1); %VehicleD
        [LocActivity, charge_power, where] = LocationDesignation(VehTrips(n,5),chargingcodes);
        Activitycell{M,2} = LocActivity; %Current State
        Activitycell{M,3} = VehTrips(n,3)/60; %Start Time at State
        Activitycell{M,4} = 1439/60; %End Time at State.  In this
        % case, it is the last minute of the day.  The units for time in the
        % VehTrips file is minutes, so this must be converted to hours.
        Activitycell{M,5} = -1; %Trip Distance    Activitycell{M,6} = -1; %Drive cycle
        Activitycell{M,6} = -1; %DriveCycleCode
        Activitycell{M,7} = charge_power; %max power available at the charging station
        Activitycell{M,8} = where;
        Activitycell{M,9} = VehTrips(n,7); %The NHTS household weight for the trip
        Activitycell{M,10} = VehTrips(n,9); %The household state location
        Activitycell{M,11} = VehTrips(n,10); %The household CBSA location
        M=M+1;
    end
    n=n+1;
end
