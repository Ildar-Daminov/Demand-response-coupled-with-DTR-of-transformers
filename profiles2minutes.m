function [minutes,intervals] = profiles2minutes(PUL,theta_h,theta_0,AEQ)
%% Goal of the  function
% This function his function analyzes the power and temperature profiles 
% and shows minutes and intervals where constraints are violated

% Syntax

% Input:
% PUL - time-series profile of transformer loadings in pu 
% theta_h - time-series profile of transformer's hot spot temperature in °C 
% theta_0 - time-series profile of transformer top-oil temperature in °C 
% AEQ - ageing f transformer insulation 

% Output:
% minutes - vector of minutes where transformer constraints are violated
% intervals - array of intervals where transformer constraints are violated

%% Processing the input data

% Analyzing the temperature and loading profiles which limits are exceeded.
index_hst=find(theta_h>120);
index_tot=find(theta_0>105);
index_pul=find(PUL>1.5);
% index_aeq=AEQ>1;

% Attention!!! It is better to import the limits (add more inputs) and not
% use the numbers inside of function! global variable as alternative?

% Little trick. Dividing by 60 we get hours and dividing by 24 we get a day
% number and remainder helps us to choose the horizon in minutes

minutes_integer_hst=index_hst/60/24;
minutes_integer_pul=index_pul/60/24;
minutes_integer_tot=index_tot/60/24;

% Find minutes vector for each constraints where they were violated 

% For hot spot temperature
if ~(length(minutes_integer_hst)==0)
    [day_index_hst]=minutes_integer2day_index(minutes_integer_hst);
    minutes_theta_h=unique(day_index_hst);
end % end of if

% For transformer loadings
if ~(length(minutes_integer_pul)==0)
    [day_index_pul]=minutes_integer2day_index(minutes_integer_pul);
    minutes_PUL=unique(day_index_pul);
end % end of if

% For top-oil temperature 
if ~(length(minutes_integer_tot)==0)
    [day_index_tot]=minutes_integer2day_index(minutes_integer_tot);
    minutes_theta_0=unique(day_index_tot);
end % end of if

% Combine minutes where different constraints (temperature, current) are
% violated
if ~(length(minutes_integer_hst)==0)&&~(length(minutes_integer_pul)==0)&&~(length(minutes_integer_tot)==0)
    minutes=[minutes_theta_h; minutes_theta_0; minutes_PUL];
elseif ~(length(minutes_integer_hst)==0)&&~(length(minutes_integer_pul)==0)&&(length(minutes_integer_tot)==0)
    minutes=[minutes_theta_h; minutes_PUL];
elseif ~(length(minutes_integer_hst)==0)&&(length(minutes_integer_pul)==0)&&~(length(minutes_integer_tot)==0)
    minutes=[minutes_theta_h; minutes_theta_0];
elseif (length(minutes_integer_hst)==0)&&~(length(minutes_integer_pul)==0)&&~(length(minutes_integer_tot)==0)
    minutes=[minutes_theta_0; minutes_PUL];
elseif ~(length(minutes_integer_hst)==0)&&(length(minutes_integer_pul)==0)&&(length(minutes_integer_tot)==0)
    minutes=[minutes_theta_h];
elseif (length(minutes_integer_hst)==0)&&~(length(minutes_integer_pul)==0)&&(length(minutes_integer_tot)==0)
    minutes=[minutes_PUL];
elseif (length(minutes_integer_hst)==0)&&(length(minutes_integer_pul)==0)&&~(length(minutes_integer_tot)==0)
    minutes=[minutes_theta_0];
end % end of if

minutes=unique(minutes); % it will be needed when combining minutes from different limits

minutes=sort(minutes); % sort in acsending order

% Convert minutes to intervals 
[intervals]=minutes2intervals(minutes);

end % end of function