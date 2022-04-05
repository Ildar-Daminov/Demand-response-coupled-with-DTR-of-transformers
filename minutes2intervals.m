function [intervals]=minutes2intervals(minutes)
%% Goal of this function 
%  This function converts minutes values (where transformer constraints 
%  are violated) to intevals [strart end]

% SYNTAX

% Input:
% minutes - a vector of minutes 

% Output:
% intervals - an array of intervals 

%% Data processing 

if length(minutes)>1440 % 1440 minutes 
    % Find consecutive intervals
    consecutive_intervals = diff(minutes);
    
    % Find the end of intervals 
    finish_intervals=find(consecutive_intervals>1);
    
    if ~(finish_intervals==0)
        
        % Find the start of intervals  
        for i=1:length(finish_intervals)
            if i==1 % index 
                start_intervals=1;
            else % i>1
                start_intervals(i,1)=finish_intervals(i-1)+1;
            end
        end % end of for cycle
        
        % Calculate the interval [start end]
        start_intervals(end+1,1)=finish_intervals(i)+1;
        finish_intervals(end+1,1)=length(minutes);
        intervals=[start_intervals finish_intervals];
        
    else %(finish_intervals==0)
        intervals=[1 525600]; % whole year in minutes [start end]
    end % end of "if ~(finish_intervals==0)"
    
else % length(minutes)<=1440
    intervals=[minutes(1) minutes(end)];
    
end % end of "if length(minutes)>1440 "

end % end of function