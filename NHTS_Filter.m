function [ TID , DesiredOut , numFiltered ] = NHTS_Filter(datayr, datadesired, singledriver)
%NHTS_Filter - This function filters the raw NHTS data for either 2001 or 2009 and
%outputs set of trip matrices for use with PECM and another set to derive
%statistics of the demographics in the set.
%   Detailed explanation goes here

% "datayr" should be the data set that you would like to use.  Either 2001
% or 2009.

% "datadesired" should be a cell in the form: row 1 contains the exact
% string matches to the column headers in the NHTS dataset that you would
% like to examine. The "DesiredOut" result will be in the same order will
% be in the same order. The next two rows of the cell contain filtration
% data for each of the chosen variables.Row 2 contains either a minimum
% value or a set of codes. The set will replace the min/max filtration. Row
% 3 contains a maximum value.  If there is no minimum, then the code should
% be set to -inf, if no max ==> inf.

% "TID" is the output trip matrix in the format required by PECM.  It is a
% struct with 5 subgroups: all, car, van, SUV, pickup.  In each of the
% subgroups there is a cell that contains the trips by day of the week
% (Sunday to Saturday, columns 1:7).  Each of these cells are a mx8 double,
% in the format: [ NEWVEHID | STRTTIME | ENDTIME | TRPMILES | WHYTO |
% WHYFROM | WTHHFIN | WTTRDFIN ].  The STRTTIME and ENDTIME will have been changed
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


%% load the prefiltered data files:
load(['NHTS' num2str(datayr) 'prefiltered_1drv-' num2str(singledriver) '.mat'])
% What is saved in these data files is the 'numFiltered' variable (a cell 
% containing first column is the column header filteres, second is the 
% remaining trips after filter, and the third is the remaining vehicles), 
% the 'NHTS_data' variable (containing the raw data) and the 'NHTS_headers'
% variable, which is a cell

%%Prep some of the header columns and the f for numFiltered
f = size(numFiltered,1)+1;
datasize = size(datadesired,2); % This is the number of things that the data will be filtered by
colVID = strcmp('NEWVEHID',NHTS_headers);
colS = strcmp('STRTTIME',NHTS_headers);
colE = strcmp('ENDTIME',NHTS_headers);
colTD = strcmp('TRPMILES',NHTS_headers);
colWh2 = strcmp('WHYTO',NHTS_headers);
colWhF = strcmp('WHYFROM',NHTS_headers);
colWgH = strcmp('WTHHFIN',NHTS_headers);
colWgT = strcmp('WTTRDFIN',NHTS_headers);
colDoW = strcmp('TRAVDAY',NHTS_headers);
colTT = strcmp('TRPTRANS',NHTS_headers);
colST = strcmp('HHSTFIPS',NHTS_headers);
colCBSA = strcmp('HH_CBSA',NHTS_headers);

%% Filtering Phase 3 - User specified filtering
for d=1:datasize
    if length(datadesired{2,d})==1
        % min max filter
        minLim = datadesired{2,d};
        maxLim = datadesired{3,d};
        header = datadesired{1,d};
        col = strcmp(header,NHTS_headers);
        if minLim==maxLim %the instance when there is a single value for acceptable data
            badData = NHTS_data(~(NHTS_data(:,col)==minLim),colVID);
        else
            badData = NHTS_data(NHTS_data(:,col) <= minLim | NHTS_data(:,col)> maxLim,colVID);
        end
        % Filter data to keep only vehicles that are not in badData
        NHTS_data = NHTS_data(~ismember(NHTS_data(:,colVID),badData),:);
        % Add filtering Information
        numFiltered(f,:)={header size(NHTS_data,1) size(unique(NHTS_data(:,colVID)),1)};
        f=f+1;
        clear('header','badData','col','maxLim','minLim')
    else
        % Limit based on specific values provided
        values = datadesired{2,d};
        header = datadesired{1,d};
        col = strcmp(header,NHTS_headers);
        badData = NHTS_data(~(ismember(NHTS_data(:,col),values)),colVID);
        % Filter data to keep only vehicles that are not in badData
        NHTS_data = NHTS_data(~ismember(NHTS_data(:,colVID),badData),:);
        % Add filtering Information
        numFiltered(f,:)={header size(NHTS_data,1) size(unique(NHTS_data(:,colVID)),1)};
        f=f+1;
        clear('header','badData','col','values')
    end
end
clear('d')
%% Construct the TID and DesiredOut matrices

% Create matrix holding properly formatted info (per PECM)
% tid = [NEWVEHID STRTTIME ENDTIME TRIPMILES WHYTO WTHHFIN 'WTTRDFIN' HHSTFIPS HH_CBSA]
tid = [NHTS_data(:,colVID) NHTS_data(:,colS) NHTS_data(:,colE) NHTS_data(:,colTD) ...
    NHTS_data(:,colWh2) NHTS_data(:,colWhF) NHTS_data(:,colWgH)  NHTS_data(:,colWgT) NHTS_data(:,colST) NHTS_data(:,colCBSA)];
% Create a structure housing all the data
% DaysTID = {SuTID,MTID,TTID,WTID,ThTID,FTID,SaTID};
TID.all={};
TID.car={};
TID.van={};
TID.SUV={};
TID.pickup={};
for day = 1:7
    TID.all{1,day}=tid((NHTS_data(:,colDoW)==day),:);
    TID.car{1,day} =tid((NHTS_data(:,colTT)==1 & NHTS_data(:,colDoW)==day),:);
    TID.van{1,day}=tid((NHTS_data(:,colTT)==2 & NHTS_data(:,colDoW)==day),:);
    TID.SUV{1,day}=tid((NHTS_data(:,colTT)==3 & NHTS_data(:,colDoW)==day),:);
    TID.pickup{1,day}=tid((NHTS_data(:,colTT)==4 & NHTS_data(:,colDoW)==day),:);
end

DesiredOut.all={};
DesiredOut.car={};
DesiredOut.van={};
DesiredOut.SUV={};
DesiredOut.pickup={};
for m=1:datasize
    header=datadesired{1,m};
    col = strcmp(header,NHTS_headers);
    for day = 1:7
        DesiredOut.all{1,day}(:,m)=NHTS_data((NHTS_data(:,colDoW)==day),col);
        DesiredOut.car{1,day}(:,m) =NHTS_data((NHTS_data(:,colTT)==1 & NHTS_data(:,colDoW)==day),col);
        DesiredOut.van{1,day}(:,m)=NHTS_data((NHTS_data(:,colTT)==2 & NHTS_data(:,colDoW)==day),col);
        DesiredOut.SUV{1,day}(:,m)=NHTS_data((NHTS_data(:,colTT)==3 & NHTS_data(:,colDoW)==day),col);
        DesiredOut.pickup{1,day}(:,m)=NHTS_data((NHTS_data(:,colTT)==4 & NHTS_data(:,colDoW)==day),col);
    end
end
end

function [NHTS_data] = Phase2_filter(NHTS_data,colVID,NHTS_headers)
colP=strcmp('PERSONID',NHTS_headers);
% we need to use PERSONID and NEWVEHID to determine if there is more
% than one person driving a vehicle.  This could take a while cause the
% only way I can think to do it is by looping
VehicleIDs = unique(NHTS_data(:,colVID));
badData=[];
for l = 1:size(VehicleIDs,1)
    Vehicle=NHTS_data(NHTS_data(:,colVID)==VehicleIDs(l),:);
    if size(unique(Vehicle(:,colP)),1)>1
        badData = [badData ; VehicleIDs(l)];
    end
end
% Filter data to keep only vehicles that are not in badData
NHTS_data = NHTS_data(~ismember(NHTS_data(:,colVID),badData),:);
end