function [intervals]=minutes2intervals(minutes)
if length(minutes)>1440
    consecutive_intervals = diff(minutes);
    finish_intervals=find(consecutive_intervals>1);
    if ~(finish_intervals==0)
        for i=1:length(finish_intervals)
            if i==1
                start_intervals=1;
            else %i>1
                start_intervals(i,1)=finish_intervals(i-1)+1;
            end
        end
        start_intervals(end+1,1)=finish_intervals(i)+1;
        finish_intervals(end+1,1)=length(minutes);
        intervals=[start_intervals finish_intervals];
    else %(finish_intervals==0)
        intervals=[1 525600];
    end
else % length(minutes)<=1440
    intervals=[minutes(1) minutes(end)];
end
end