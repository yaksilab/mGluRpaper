function [freezing_percentage,freezing_time]=me_freezing_fromFP(X,Y,D,parameters) 
%% parameters needs to be = [dist to freeze in mm , sec to freeze in sec]
% D = diff(dt)
dist_to_freeze=parameters(2);
second_to_freeze=parameters(1);
for j=1 %%loop over sessions    
    for i=1  %%loop over fish
        tempd= [];tempd= D*1000; %time needs to be in ms, you have it in sec, so multiple by 1000
        tempx= [];
        tempy= [];
        tempx= X;
        tempy= Y;    
        deltax=[];
        deltay=[];
        dist_tot=[];
        deltax=diff(tempx);
        deltay=diff(tempy);
        deltax=fillmissing(deltax,'next');
        deltay=fillmissing(deltay,'next');
         % the code from fabrizio gave some slightly different results for
        % the distance, use this one instead
            for k=2:size(deltax,1)-1
               dis(1,k-1) = pdist(cat(2,deltax([k,k+1]),deltay([k,k+1])));
            end
        dist_tot=dis;  % mine is already in mm/3; %conversion pix to mm
        tempd=tempd(1:size(dist_tot,2)); %distance and tempd need to be the same
     dist_cumulative=floor(cumsum(dist_tot)./dist_to_freeze);
     dist_cum=diff(dist_cumulative);    
     moved=[];moved=[find(dist_cum>0)];
        deltat=[];
        deltat(1)=sum(tempd(1:moved(1))) ;
    
      for t= 1: max(size(moved))
         if(t< max(size(moved)))
          deltat(t)=sum(tempd((moved(t)+1):moved(t+1))) ;
         else
          deltat(t)=sum(tempd((moved(t)+1):moved(end))) ;
         end
      end
    freezing_time(j,i)=sum(deltat(deltat>=(second_to_freeze*1000)));%% convert second_to_freeze in ms
    freezing_percentage(j,i)=freezing_time(j,i)./sum(deltat);
    end
end
end


%i am still testing it now, but it should be easy to use ! you throw in X, Y, D(is the time between frames), 
%and info_vector must be a vector with the size of the arena in pixel for each fish(use squared roi so you can just give one dimension) 
%and arena_in_mm the same vector but with the dimension in mm ! 
%how confusing is it ? XD