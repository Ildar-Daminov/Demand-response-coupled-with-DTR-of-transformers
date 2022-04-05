function [Pt,tt]=Convert2hours(P1,t)
%% Goal of this function
% This function converts the data in minutes to the data in hours 

% Syntax 
% Input: 
% AMB - a profile of ambient temperature (as an example) in minutes
% t - always 60 (minutes) in this article

% Output :
% AMB - a profile of ambient temperature (as an example)in hours

%% Converting the data to hour resolution

N=length(P1); % find the length of vector

% if ~(N==1440)
%     error('длина суточного графика должна быть равна 1440 минут');
% end;

if ~(t==5||t==15||t==30||t==60||t==3||t==2||t==6||t==12||t==20||t==4||t==10||t==8)
    error('Problem with t. Check it');
end

n=fix(N/t); % number of hours 

%tt=((1:n)-1)*t+1;
tt=(1:n)*t; % resolution in hours 
tt=tt';

% Create a empty vector  
Pt=zeros(n,1);

% Convert to hour resolution using the mean values
for i=1:n
    Pt(i)=mean(P1(((i-1)*t+1):i*t));
end % end of for cycle 

end % end of function
