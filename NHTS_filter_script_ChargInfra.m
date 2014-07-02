clear
close all
clc

%Basic Filter

%% Inputs 
NHTS = 2009;
singledriver=0; %Flag to indicate whether or not we are using data that has been filtered such that vehicles are only driven by a single driver
distmin = 1; %[miles]
distmax = inf; %[miles]
vehtype = 1; %1-5 => {'all','car','van','SUV','pickup'}
NHTS_State = 5; %State number (alphabetical includes DC, we are interested in california = 5)
StateTxt = 'CA';

Level1Charge = 120*12; 
Level2Charge = 240*30; %This is the most common Level 2 charge available in most stations, they are current limited to 30A at 240V.
FastDCCharge = 50000; %50kW Fast DC Charging

chargingcodes = {1 10 11 12 ;...%The NHTS codes for home (1) and work (10-12)
    [1 Level1Charge; 0 Level2Charge]...
    [0 Level1Charge; 1 Level2Charge]...
    [0 Level1Charge; 1 Level2Charge]...
    [0 Level1Charge; 1 Level2Charge]}; 
%chargingcodes is a cell that has the allowable codes for charging and
%their maximum charging level at those locations.  The codes are in the
%first row of the cell. The second contains a matrix that indicates the
%probability of having a particular charging power available the first
%column is the probability (all sum to one), the second is the charging
%power of each probability in Watts).


% BattSize = 24; % Battery size in kW of the vehicles used
% RoadConsume = 0.24; % rate of energy consumption when vehicle is in a trip in kW per mile
% MaxBattUt = 0.8; % Maximum battery utilization
% seas = 4; %Season of interest (if the NHTS_filter is run with datadesired(:,1) as the input this is ignored) 
% timestep = 5; %combining the NHTS into timesteps of this many minutes
% daytype = 0; %Zero for weekdays, 1 for weekends
% FltSize =  100;  %This is the number of vehicles we want in our ficitious fleet.
 
% Inputs as a function of reading the csv file (FleetGenInputs.csv)
% [num,txt]=xlsread('FleetGenInputs.csv');
% NHTS = num(1);
% singledriver=num(2); %Flag to indicate whether or not we are using data that has been filtered such that vehicles are only driven by a single driver
% distmin = num(3); %[miles]
% BattSize = num(4); % Battery size in kW of the vehicles used
% RoadConsume = num(5); % rate of energy consumption when vehicle is in a trip in kW per mile
% MaxBattUt = num(6); % Maximum battery utilization
% seas = num(7); %Season of interest (if the NHTS_filter is run with datadesired(:,1) as the input this is ignored)
% NHTS_State = num(8); %State number (alphabetical includes DC, we are interested in california = 5)
% StateTxt = txt{9,2};
% timestep = num(10); %combining the NHTS into timesteps of this many minutes
% daytype = num(11); %Zero for weekdays, 1 for weekends
% FltSize =  num(12);  %This is the number of vehicles we want in our ficitious fleet.

% These are the NHTS month callouts for each season, they will be input into 
% season = {'Winter' 'Spring' 'Summer' 'Fall';...
%     [200812 200901 200902] [200803 200903 200804 200904 200805] 200806:200808 200809:200811};

%% Data Desired for filtration
%specify the datadesired variables, remember the max is included but the minimum isn't

% datadesired = {'HHSTATE' 'TDAYDATE'; NHTS_State season{2,seas} ; NHTS_State []}; %California is the fifth state in an alphabetical list, and that is what we are looking for
% datadesired = {'HHC_MSA' ; [7362 0680 1620 2840 4940 5170 6690 7120 8120
% 7460] ; []}; %All PG&E MSA's, note that the only MSA that is actually in
% the NHTS of this list is the SF Bay area (7362)

% datadesired = {'HHC_MSA' ; [7362] ; [7362]}; %San Francisco Bay Area
% notefilter = 'MSA-SFBAY';
% 
% datadesired = {'HHC_MSA' ; [6922] ; [6922]}; %Sacramento Area
% notefilter = 'MSA-SAC';
% 
% datadesired = {'HHC_MSA' ; [4472] ; [4472]}; %Los Angeles Area
% notefilter = 'MSA-LA';

% datadesired = {'HHC_MSA' ; [7320] ; [7320]}; %San Diego Area
% notefilter = 'MSA-SD';
% 
% datadesired = {'HHSTATE' ; NHTS_State ; NHTS_State}; %All CA
% notefilter = StateTxt;    

%%All CA that is not in one of the big MSA's
% datadesired = {'HHSTATE' 'HHC_MSA' ; NHTS_State [0 -1] ; NHTS_State []};
% notefilter = [StateTxt '-smallMSA'];      
% 
%All CA in NHTS that is not in LA or San Diego
% datadesired = {'HHSTATE' 'HHC_MSA' ; NHTS_State [0 -1 6922 7362] ; NHTS_State []};
% notefilter = [StateTxt '-noLAorSD']; 
% 
%Let's get the whole NHTS (cheating by calling out a variable with no
%bounds)
datadesired = {'HHC_MSA' ; -inf ; inf};
notefilter = 'whole_enchilada'; 


%% Perform Filtration
% Filter out the important data
[ TID , DesiredOut , numFiltered ] = NHTS_Filter(NHTS, datadesired, singledriver);

% Filter out by total daily travel
% distmax = BattSize*MaxBattUt/RoadConsume; % Determine the maximum distance a vehicle could travel in a day, given it's battery life
% distmax = inf;
[TID,removed,numFiltered] = TID_distfilter(distmax,distmin,TID,numFiltered);

% Things to remember about the TID
% It is broken up by weekend (Sunday to Saturday, columns 1:7)
% It is a struct: TID.vehicletype => vehicletype = {'all','car','van','SUV','pickup'}; %indicates which sub-group of vehicles in the TID to use (NHTS classifications)
%TID colums: [ NEWVEHID | STRTTIME | ENDTIME | TRPMILES | WHYTO | WHYFROM | WTHHFIN | WTTRDFIN | HHSTFIPS | HH_CBSA ]

%% Pick the correct vehicle type
vehicletypes = {'all','car','van','SUV','pickup'};
eval(['TripMatrices = TID.' vehicletypes{vehtype} ';'])

%% Create V2GSim initiation files
% One for each day of the week (daywk)
Sheetname1 = 'Activity';
Sheetname2 = 'Initialization';
for daywk=1:7 %Sunday to Saturday
    TripMatrix=TripMatrices{daywk};
    TripMatrix=sortrows(TripMatrix,2);
    TripMatrix=sortrows(TripMatrix,1);
    filenm = ['V2GSimInit\initfile_NTHStype-' vehicletypes{vehtype}  '_' ...
        notefilter '_dayweek-' num2str(daywk) '.xlsx'];
    [ActivitySheet]=V2G_ActivitySheetGen(TripMatrix, chargingcodes);
    xlswrite(filenm,ActivitySheet,Sheetname1);
end

% %Combine all of the weekdays
% TripMatrix = [];
% 
% for daywk = 2:6
%     TripMatrix=[TripMatrix ; TripMatrices{daywk}];
% end
% 
% TripMatrix=sortrows(TripMatrix,2);
% TripMatrix=sortrows(TripMatrix,1);


% filenm = ['V2GSimInit\initfile_NTHStype-' vehicletypes{vehtype}  '_' ...
%     notefilter '_Weekdays.xlsx'];
% [ActivitySheet]=V2G_ActivitySheetGen(TripMatrix, chargingcodes);
% xlswrite(filenm,ActivitySheet,Sheetname1);
% 
% %Let's do weekends
% TripMatrix = [TripMatrices{1} ; TripMatrices{1}];
% TripMatrix=sortrows(TripMatrix,2);
% TripMatrix=sortrows(TripMatrix,1);

% filenm = ['V2GSimInit\initfile_NTHStype-' vehicletypes{vehtype}  '_' ...
%     notefilter '_Weekends.xlsx'];
% [ActivitySheet]=V2G_ActivitySheetGen(TripMatrix, chargingcodes);
% xlswrite(filenm,ActivitySheet,Sheetname1);

% save(['FilteredData-' num2str(NHTS) '-' note '-' season{1,seas}],'TID','DesiredOut','numFiltered')


miles = TID.all{1}(:,4);
for m = 2:length(TID.all)
    miles = [miles ; TID.all{m}(:,4)];
end
AvTravel = sum(miles)/numFiltered{end,3};
AvDist = sum(miles)/numFiltered{end,2};

%% Prefilter script when Necessary
% datayr=NHTS;
% [NHTS_data, NHTS_headers, numFiltered ] = NHTS_PreFilter(datayr,singledriver);
% save(['NHTS' num2str(datayr) 'prefiltered_1drv-' num2str(singledriver)],...
%     'NHTS_data','NHTS_headers','numFiltered')

% %% Code for Generating availability 
% 
% 
% [EVTrav, EVDist, EVHome, EVWork]= EVavailability(TID, vehicletype{1},timestep, daytype);
% 
% %% Code for Generating our fictitious fleet
% 
% [FltTrav, FltDist, FltHome,FltWork]=FleetGenerator(EVTrav, EVDist, EVHome, EVWork,FltSize);
% 
% %% Using the distance traveled to determine an Energy Consumption Matrix
% 
% FltConsume=FltDist*RoadConsume; %energy consumed at every time step for our fleet
% 
% % We may want to change this so that it reflects different consumption
% % rates based on how fast we are going (we could do just two, one for
% % highway driving and one for regular road, we could even have a weighted
% % average of the two based on the distance: if you have traveled more than
% % 4 miles in five minutes, then use highway, if less than 1.5 use the
% % regular, and anything between is a weighted average).  We would need to
% % implement a similar scheme in our filtration, but the code should be
% % nearly copy/paste.  
% 
% %% Sanity Check!  
% % Let's Graph the inputs and outputs to see how well our fleet represents
% % the average.
% 
% toplot=[sum(EVWork,1)/size(EVWork,1);sum(EVHome,1)/size(EVHome,1);sum(EVTrav,1)/size(EVTrav,1)];
% toplot2=[sum(FltWork,1)/size(FltWork,1);sum(FltHome,1)/size(FltHome,1);sum(FltTrav,1)/size(FltTrav,1)];
% subplot(2,1,1)
% plot(1:288,toplot*100)
% title('NHTS data')
% ylabel('Percent of Vehicles')
% axis([1 288 0 100])
% 
% subplot(2,1,2)
% plot(1:288,toplot2*100)
% title('Generated Fleet')
% ylabel('Percent of Vehicles')
% axis([1 288 0 100])
% legend('Work','Home','Traveling')
% 
% %% Copying the Data into relevant csv files
% 
% filenm = 'initfile_NTHS.xlsx';
% Sheetname1 = 'Activity';
% Sheetname2 = 'Initialization';
% csvwrite(fn_home,FltHome')
% csvwrite(fn_work,FltWork')
% csvwrite(fn_consume,FltConsume')