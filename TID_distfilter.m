function [TID_L,removed,numFiltered] = TID_distfilter(distmax,distmin,TID,numFiltered)
%
%limit total daily travel for vehicles in each day of Trip Matrix,

subgroup = {'all' 'car' 'van' 'SUV' 'pickup'};
for l=1:length(subgroup);
    eval(['DaysTID=TID.' subgroup{l} ';'])
    daysofweek = length(DaysTID);
    removed=zeros(daysofweek,2);
    %this is just right, logical indexing, not too slow and it works
    for day = 1:1:daysofweek
        TID2 = DaysTID{day};
        VID=zeros(length(unique(TID2(:,1))),2);
        VID(:,1)=unique(TID2(:,1));
        for m = 1:length(VID(:,1)) %all
            %         VID(m,2)=sum(TID(ismember(TID(:,1),VID(m,1)),4));
            VID(m,2)=sum(TID2((TID2(:,1)==VID(m,1)),4));
        end
        BadVID = VID(VID(:,2)>distmax,:);
        RBVID = VID(VID(:,2)<=distmin,:);
        removed(day,:)=[length(BadVID(:,1)) , length(RBVID(:,1))];
        BadVID = [BadVID(:,1) ; RBVID(:,1)];
        newTID = TID2(~(ismember(TID2(:,1),BadVID(:,1))),:);
        newDaysTID{day,1}= newTID;
    end
    
    eval(['TID_L.' subgroup{l} '=newDaysTID'';'])
end

FTID=TID_L.all;
numFiltered{end+1,1}='DailyTravel';
Vehs = 0;
Trips = 0;
for m = 1:length(FTID)
    Trips = Trips+size(FTID{1,m},1);
    Vehs = Vehs + size(unique(FTID{1,m}(:,1)),1);
end
numFiltered{end,2}=Trips;
numFiltered{end,3}=Vehs;
