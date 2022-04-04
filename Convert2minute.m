function [AMB]=Convert2minute(AMB)
%% Goal of this function
% This function converts the data in hours to the date in minutes 

% Syntax 
% Input: 
% AMB - a profile of ambient temperature (as an example) in hours

% Output :
% AMB - a profile of ambient temperature (as an example)in minutes

%%
t=60; % 60 min per  hour

% number of hours
nhours=fix(length(AMB)*t/t);

% intermediate variable: ambient temperature in minutes
AMB_minutes=zeros(length(AMB)*t,1);

for i=1:nhours % for each hour
    % Set the mean ambient temperature
    AMB_minutes(((i-1)*t+1):i*t)=AMB(i);
end

% Change a variable name
AMB=AMB_minutes;


end % end of function
