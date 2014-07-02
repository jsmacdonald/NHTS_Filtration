function [NHTS_data, NHTS_headers, numFiltered ] = NHTS_PreFilter(datayr,singledriver)
%NHTS_Filter - This function filters the raw NHTS data for either 2001 or 2009 and
%outputs set of trip matrices for use with PECM and another set to derive
%statistics of the demographics in the set.
%   Detailed explanation goes here

% "datayr" should be the data set that you would like to use.  Either 2001
% or 2009.

% "singledriver" is a binary input flag that indicates whether or not
% vehicles that were driven by multiple drivers should be removed. 1
% indicates that only vehicles with a single driver are included in the
% analysis.

% "TID" is the output trip matrix in the format required by PECM.  It is a
% struct with 5 subgroups: all, car, van, SUV, pickup.  In each of the
% subgroups there is a cell that contains the trips by day of the week
% (Sunday to Saturday, columns 1:7).  Each of these cells are a mx7 double,
% in the format: [ NEWVEHID | STRTTIME | ENDTIME | TRPMILES | WHYTO |
% WTHHFIN | WTTRDFIN ].  The STRTTIME and ENDTIME will have been changed
% from military time to sequential minutes of the day.  NEWVEHID is a
% derived variable based on the HOUSEID and VEHID NHTS variables. The
% others are NHTS variables representing the trip distance, destination,
% and weights for the household and trip.

% "DesiredOut" is the matrix of other NHTS data for statistical evaluation,
%  This matrix will include anything that was included in the first row  of
%  datadesired input variable.

% "numFiltered" indicates the number of trips and vehicles in each level of
% the filtering process.

% Sample Inputs:
% datayr=2009;
% datadesired={'R_AGE' 'TRPMILES'; -inf 0 ;inf 40};

%% Get the correct data files:
if datayr == 2001
    load NHTS2001raw.mat
    DP_col=DAYPUB_Title;
    VP_col=VEHPUB_Title;
    DP1_data=DAYPUB_Data;
    DP2_data={};
    VP_data=VEHPUB_Data;
    KEEPdataDP = {'CENSUS_D' 'CENSUS_R' 'DRVR_FLG' 'DRVRCNT' 'EDUC' 'ENDTIME' ...
        'HBHUR' 'HTPPOPDN' 'HBHRESDN' 'HHR_HISP' ...
        'HH_ONTD' 'HHR_RACE' 'HHC_MSA' 'HHFAMINC' 'HHRESP' 'HHSIZE' 'HHSTATE' 'HHSTFIPS' ...
        'HHVEHCNT' 'HOMEOWN' 'HOMETYPE' 'HOUSEID' 'HTPPOPDN' 'HTHRESDN' ...
        'LIF_CYC' 'MSACAT' 'MSASIZE' 'NONHHCNT' 'NUMADLT' 'NUMONTRP' 'PERSONID' ...
        'PRMACT' 'R_AGE' 'R_SEX' 'STRTTIME' 'TDAYDATE' 'TDTRPNUM' 'TDWKND' ...
        'TRAVDAY' 'TRIPPURP' 'TRPHHACC' 'TRPHHVEH' 'TRPMILES' 'TRPTRANS' ...
        'URBAN' 'URBRUR' 'VEHID' 'VEHTYPE' 'WHYFROM' 'WHYTO' ...
        'WHYTRP1S' 'WHYTRP90' 'WORKER' 'WRKCOUNT' 'WTTRDFIN'};
    KEEPdataVP = { 'HOUSEID' 'VEHID' 'GSCOST' 'GSTOTCST' 'GSYRGAL' 'MAKECODE' ...
        'MODLCODE' 'OD_READ1' 'OD_READ2' 'VEHAGE' 'VEHOWNMO' 'VEHYEAR' 'WHOMAIN' 'WTHHFIN'};
    clear('DAYV2PUB_Title','VEHV2PUB_Title','DAYV2PUB1_Data','DAYV2PUB2_Data','VEHV2PUB_Data')
    clear('DAYPUB_Title','VEHPUB_Title','DAYPUB_Data','VEHPUB_Data')
elseif datayr == 2009
    load NHTS2009raw.mat
    DP_col=DAYV2PUB_Title;
    VP_col=VEHV2PUB_Title;
    DP1_data=DAYV2PUB1_Data;
    DP2_data=DAYV2PUB2_Data;
    VP_data=VEHV2PUB_Data;
    KEEPdataDP = {'CENSUS_D' 'CENSUS_R' 'DRVR_FLG' 'DRVRCNT' 'EDUC' 'ENDTIME' ...
        'FLAG100' 'GASPRICE' 'HBHUR' 'HTPPOPDN' 'HBRESDN' 'HH_CBSA' 'HH_HISP' ...
        'HH_ONTD' 'HH_RACE' 'HHC_MSA' 'HHFAMINC' 'HHRESP' 'HHSIZE' 'HHSTATE' 'HHSTFIPS' ...
        'HHVEHCNT' 'HOMEOWN' 'HOMETYPE' 'HOUSEID' 'HTPPOPDN' 'HTRESDN' 'INTSTATE' ...
        'LIF_CYC' 'MSACAT' 'MSASIZE' 'NONHHCNT' 'NUMADLT' 'NUMONTRP' 'PERSONID' ...
        'PRMACT' 'R_AGE' 'R_SEX' 'STRTTIME' 'TDAYDATE' 'TDTRPNUM' 'TDWKND' ...
        'TRAVDAY' 'TRIPPURP' 'TRPACCMP' 'TRPHHACC' 'TRPHHVEH' 'TRPMILES' 'TRPTRANS' ...
        'URBAN' 'URBANSIZE' 'URBRUR' 'USEINTST' 'VEHID' 'VEHTYPE' 'WHYFROM' 'WHYTO' ...
        'WHYTRP1S' 'WHYTRP90' 'WORKER' 'WRKCOUNT' 'WTTRDFIN'};
    KEEPdataVP = { 'HOUSEID' 'VEHID' 'GSCOST' 'GSTOTCST' 'GSYRGAL' 'HYBRID' 'MAKECODE' ...
        'MODLCODE' 'OD_READ' 'VEHAGE' 'VEHOWNMO' 'VEHYEAR' 'WHOMAIN' 'WTHHFIN'};
    clear('DAYV2PUB_Title','VEHV2PUB_Title','DAYV2PUB1_Data','DAYV2PUB2_Data','VEHV2PUB_Data')
else
    disp('Error: data is only available for 2001 or 2009, please choose one of those')
    return
end
%% The data columns we have to keep:
headerDP=KEEPdataDP;
headerVP=KEEPdataVP;


%% Creating the NHTS file that is solely the data columns we need
NHTS_DP = zeros((size(DP1_data,1)+size(DP2_data,1)),1);
for m=1:size(headerDP,2)
    if ismember(headerDP{1,m},DP_col);
        col = strcmp(headerDP{1,m},DP_col);
        if size(DP1_data,2)==size(DP2_data,2)
            NHTS_DP(:,m)=[DP1_data(:,col) ; DP2_data(:,col)];
        else
            NHTS_DP(:,m)=[DP1_data(:,col)];
        end
    else
        disp([headerDP{1,m} ' is not a header in the ' num2str(datayr) ' NHTS DAYPUB file.'])
    end
end

NHTS_VP = zeros((size(VP_data,1)),1);
for m=1:size(headerVP,2)
    if ismember(headerVP{1,m},VP_col);
        col = strcmp(headerVP{1,m},VP_col);
        NHTS_VP(:,m)=VP_data(:,col);
    else
        disp([headerVP{1,m} ' is not a header in the ' num2str(datayr) ' NHTS VEHPUB file.'])
    end
end

clear('m','col','DP1_data','DP2_data','DP_col','VP_data','VP_col')
%% Now for filtering phase 1
%  This is the filtering that has nothing to do with the user inputs, but
%  must be done for use in PECM.
%  The only variables we are carrying at this point is:
%  NHTS_DP, NHTS_VP, headerDP, headerVP, mustdataDP, mustdataVP, datasize,
%  datayr, and datadesired.

%initialize the variable that tracks the number of vehicles that have been
%filtered.
numFiltered = {'Filtered by:' 'TRIPS Left' 'Vehicles Left';...
    'Unfiltered' size(NHTS_DP,1) size(NHTS_VP,1)};
f=2+1;

% First purge some of the trips data that can be purged independently of
% creating the unique vehicle ID.

% Let's purge the TRPTRANS variable.  Remove anything that is not traveled
% by a car, van, suv, or pickup (TRPTRANS 1-4)
header='TRPTRANS';
col = strcmp(header,headerDP);
NHTS_DP = NHTS_DP(NHTS_DP(:,col)>=1 & NHTS_DP(:,col)<= 4,:);
numFiltered(f,:)={header size(NHTS_DP,1) size(NHTS_VP,1)};
f=f+1;
clear('col','header')

% Remove every trip that is no reported by the person driving (because we don’t want to know
% when you went as a passenger because we have the driver in here as well DRVR_FLAG =1
header='DRVR_FLG';
col = strcmp(header,headerDP);
NHTS_DP = NHTS_DP(NHTS_DP(:,col)==1,:);
numFiltered(f,:)={header size(NHTS_DP,1) size(NHTS_VP,1)};
f=f+1;
clear('col','header')

% Now Let's purge VEHID for both the NHTS_VP and NHTS_DP.  Any VEHID <=0
% must be removed.
% Remove vehicle IDs == -1 (or keep all things ~= -1)
header = 'VEHID';
colDP = strcmp(header,headerDP);
colVP = strcmp(header,headerVP);
NHTS_DP = NHTS_DP(NHTS_DP(:,colDP) ~= -1,:);
NHTS_VP = NHTS_VP(NHTS_VP(:,colVP) ~= -1,:);
numFiltered(f,:)={header size(NHTS_DP,1) size(NHTS_VP,1)};
f=f+1;
clear('colDP','colVP','header')

% We now need to append the vehicle information in the VEHPUB file onto the Trip info
% so we need to create the unique vehicle IDs.  To do this we must first
% make sure there is agreement between vehicles that are in the vehpub file
% and the daypub file.
% Create a unique “in house” household ID
header = 'HOUSEID';
colVP = strcmp(header,headerVP);
colDP = strcmp(header,headerDP);
[HHID,indx,indx2] = unique(NHTS_VP(:,colVP));
TotalVehicles=length(HHID);
HHID(:,2)=[(10^floor(log10(TotalVehicles))+1):(10^floor(log10(TotalVehicles))+TotalVehicles)]';%Creates a unique, sequential HHID with the minimum number of digits required

% Remove any vehicles reporting a HHID that is not reported in VEHPUB
[TF, indx3] = ismember(NHTS_DP(:,colDP),HHID(:,1));
NHTS_DP = NHTS_DP(TF,:);
numFiltered(f,:)={header size(NHTS_DP,1) size(NHTS_VP,1)};
f=f+1;
% Create a unique vehicle ID (NEWVEHID) from a concatenation of our "in house" HHID and the VEHID number in NHTS
headerVP(1,size(headerVP,2)+1)={'NEWVEHID'};
headerDP(1,size(headerDP,2)+1)={'NEWVEHID'};
NHTS_VP(:,size(headerVP,2)) = HHID(indx2,2)*100+NHTS_VP(:,strcmp('VEHID',headerVP));
NHTS_DP(:,size(headerDP,2)) = HHID(indx3,2)*100+NHTS_DP(:,strcmp('VEHID',headerDP));
clear('header','colVP','colDP','HHID','indx','indx2','indx3','TF')

%Append the Vehicle File data (that is not included in the DAYPUB file) to
%the trip data
CurrentDP=size(headerDP,2);
headerDP=[headerDP headerVP(1,3:(size(headerVP,2)-1))];
header='NEWVEHID';
colVP = strcmp(header,headerVP);
colDP = strcmp(header,headerDP(1,1:CurrentDP));
[TF,indx]=ismember(NHTS_DP(:,colDP),NHTS_VP(:,colVP));
%remove the vehicles in the daypub file that are not in the VEHPUB file
NHTS_DP=NHTS_DP(TF,:);
%concatenate the two together
NHTS_DP = [NHTS_DP NHTS_VP(indx(TF),3:(size(headerVP,2)-1))];
numFiltered(f,:)={header size(NHTS_DP,1) size(unique(NHTS_DP(:,strcmp('NEWVEHID',headerDP))),1)};
f=f+1;
clear('header','colVP','colDP','indx','TF','CurrentDP','NHTS_VP','headerVP')

%For the future, we are going to need a column identifier for the NHTS_DP
%for the NEWVEHID
colVID = strcmp('NEWVEHID',headerDP);

% Identify start times ('STRTTIME') < 0
header = 'STRTTIME';
colS = strcmp(header,headerDP);
badData = NHTS_DP(NHTS_DP(:,colS) < 0,colVID);
% Filter data to keep only vehicles that are not in badStartData
NHTS_DP = NHTS_DP(~ismember(NHTS_DP(:,colVID),badData),:);
% Add filtering Information
numFiltered(f,:)={header size(NHTS_DP,1) size(unique(NHTS_DP(:,colVID)),1)};
f=f+1;
clear('header','badData')

% % Identify end times < 0 ('ENDTIME')
header = 'ENDTIME';
colE = strcmp(header,headerDP);
badData = NHTS_DP(NHTS_DP(:,colE) < 0,colVID);
% Filter data to keep only vehicles that are not in badEndData
NHTS_DP = NHTS_DP(~ismember(NHTS_DP(:,colVID),badData),:);
% Add filtering Information
numFiltered(f,:)={header size(NHTS_DP,1) size(unique(NHTS_DP(:,colVID)),1)};
f=f+1;
% Change the start and end times to minutes instead of military time.
NHTS_DP(:,colS)=floor(NHTS_DP(:,colS)/100)*60+(NHTS_DP(:,colS)-100*floor(NHTS_DP(:,colS)/100));
NHTS_DP(:,colE)=floor(NHTS_DP(:,colE)/100)*60+(NHTS_DP(:,colE)-100*floor(NHTS_DP(:,colE)/100));
clear('header','badData')

% Define maximum speed (120 mph)
MaxSpd = 120;
%Define max drive time (16 hr )
MaxDrvTm = 16;
% Identify trips ('TRPMILES') of 0 miles or less, and remove trips that are longer than MaxSpd*MaxDrvTm
header = 'TRPMILES';
colTD = strcmp(header,headerDP);
badData = NHTS_DP(NHTS_DP(:,colTD) <= 0 | NHTS_DP(:,colTD)>= MaxSpd*MaxDrvTm,colVID);
% Filter data to keep only vehicles that are not in badTrvDist
NHTS_DP = NHTS_DP(~ismember(NHTS_DP(:,colVID),badData),:);
% Add filtering Information
numFiltered(f,:)={header size(NHTS_DP,1) size(unique(NHTS_DP(:,colVID)),1)};
f=f+1;
clear('header','badData','MaxDrvTm')

% % Identify average vehicle speed for trip 60 [min/hour] * distance [mile]/ (end time - start time)[min]
header = 'TripSpeed';
avgTripSpeed = 60*NHTS_DP(:,colTD)./(NHTS_DP(:,colE)-NHTS_DP(:,colS));
%
% Identify speeds of 0 or less, and remove speeds that are greater than MaxSpd
badSpd = NHTS_DP(avgTripSpeed(:,1)<=0 | avgTripSpeed(:,1)> MaxSpd,colVID);

% Filter data to keep only vehicles that are not in badTrvDist
NHTS_DP = NHTS_DP(~ismember(NHTS_DP(:,colVID),badSpd),:);
% Add filtering Information
numFiltered(f,:)={header size(NHTS_DP,1) size(unique(NHTS_DP(:,colVID)),1)};
f=f+1;
clear('header','badSpd','avgTripSpeed')
%% Filtering Phase 2 - Singledriver filter
if singledriver ==1
    [NHTS_DP] = Phase2_filter(NHTS_DP,colVID,headerDP);
    % Add filtering Information
    header='SingleDriver';
    numFiltered(f,:)={header size(NHTS_DP,1) size(unique(NHTS_DP(:,colVID)),1)};
    f=f+1;
end
clear('header','singledriver')
%% Construct the outputs and DesiredOut matrices
NHTS_data=NHTS_DP;
NHTS_headers=headerDP;
end


function [NHTS_DP] = Phase2_filter(NHTS_DP,colVID,headerDP)
colP=strcmp('PERSONID',headerDP);
% we need to use PERSONID and NEWVEHID to determine if there is more
% than one person driving a vehicle.  This could take a while cause the
% only way I can think to do it is by looping
VehicleIDs = unique(NHTS_DP(:,colVID));
badData=[];
for l = 1:size(VehicleIDs,1)
    Vehicle=NHTS_DP(NHTS_DP(:,colVID)==VehicleIDs(l),:);
    if size(unique(Vehicle(:,colP)),1)>1
        badData = [badData ; VehicleIDs(l)];
    end
end
% Filter data to keep only vehicles that are not in badData
NHTS_DP = NHTS_DP(~ismember(NHTS_DP(:,colVID),badData),:);
end