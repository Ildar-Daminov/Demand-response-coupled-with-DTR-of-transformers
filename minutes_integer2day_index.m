function [day_index]=minutes_integer2day_index(minutes_integer) 
%% Purpose of this function
% This function converts the interval duration in minutes to corresponding 
% day index

% SYNTAX 
% Input: 
% minutes_integer - vector of minutes 

% Output: 
% day_index - vector of day indexes 
%% Data processing 
% Converting to minute intervals
for i=1:length(minutes_integer)
    
    if isreal(minutes_integer(i)) && rem(minutes_integer(i),1)==0 % previous day
        % Find the start and end of the day
        start=fix(minutes_integer(i))*1440-1439;
        finish=fix(minutes_integer(i))*1440;
        day_index(:,i)=[start:1:finish];
    elseif isreal(minutes_integer(i)) && ~(rem(minutes_integer(i),1)==0)% next day
        % Find the start and end of the day
        start=fix(minutes_integer(i))*1440+1;
        finish=start+1439;
        day_index(:,i)=[start:1:finish];
    else 
        error ('Check for cycle in minutes_integer2day_index. It seems not a real number in minutes_integer ')
    end % end of if condition 
    
end % end of for cycle 

end % end of the function 