function [LocActivity, charge_power, where]=LocationDesignation(locale,chargingcodes)
%This function uses the NHTS locations for a trip and identifies the
%whether the trip ends at a place that is just parking or a place that can
%be charged and at what level to charge

%locale is the NHTS data point that indicates where a vehicle ended it's last trip.

%chargingcodes is a cell that has the allowable codes for charging and
%their maximum charging level at those locations.  The codes are in the
%first row of the cell. The each entry of the second contains a matrix that
%indicates the probability of having a particular charging power available: 
%the first column is the probability (all sum to one), the second is the 
%charging power of each probability in Watts).

%% Initialize local variables
NHTS_Locales = {'Home' 1;...
    'Work' [10 11 12];...
    'School/Church' [20 21 22 23 24];...
    'Medical/Dental' 30;...
    'Shopping/Errands' [40 41 42];...
    'Gas Station' [43];...
    'Social/Recreational' [50 51 52 53 54];...
    'Public Building (Library, museum, park, etc.)' 55;...
    'Personal Obligations' [60:65];...
    'Transport Someone else' [70:73];...
    'Restaurant' [80 81 82];...
    'Coffee Shop' 83;...
    'Other' [97 13 14]};

AcceptableChargeList = [];
NHTS_locationcodes = [];

%creating a double with acceptable codes for NHTS locations
for l=1:size(NHTS_Locales,1)
    NHTS_locationcodes = [NHTS_locationcodes NHTS_Locales{l,2}];
end

%% Check the probability
for l=1:size(chargingcodes,2)
    tester = sum(chargingcodes{2,l}(:,1));
    if tester~=1
        display(['Error: the sum of probabilities for charging power does not equal 1 for location type ' num2str(chargingcodes{1,l})])
        return
    end
    AcceptableChargeList = [AcceptableChargeList chargingcodes{1,l}]; %Just wanted to convert the cell entries into a double so that I can compare to see if things are charging or not;
end


%% Assign the outputs
if isempty(locale)
    LocActivity = 'Driving';
    charge_power = -1;
    where = -1;
elseif ismember(locale,AcceptableChargeList)
    LocActivity = 'Charging';
    % Determine charging power
    Chargepowerprob=chargingcodes{2,AcceptableChargeList==locale};
    Chargepowerprob(:,1)=cumsum(Chargepowerprob(:,1));
    prob=rand;
    m=1;
    while prob>Chargepowerprob(m,1)
        m=m+1;
    end
    charge_power = Chargepowerprob(m,2);
    % Determine the name of the location
    m=1;
    while ~ismember(locale,NHTS_Locales{m,2})
        m=m+1;
    end
    where = NHTS_Locales{m,1};
else
    LocActivity = 'Parked';
    charge_power = -1;
    %Determine the name of the location
    if ismember(locale,NHTS_locationcodes)
        m=1;
        while ~ismember(locale,NHTS_Locales{m,2})
            
            m=m+1;
        end
        where = NHTS_Locales{m,1};
    else
        where = 'Not in NHTS Codebook';
    end
end
end