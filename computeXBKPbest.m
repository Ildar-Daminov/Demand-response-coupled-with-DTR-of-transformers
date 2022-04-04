% Function that computes the best breakpoints for simple PWL
function [Xbkp]=computeXBKPbest(X,Y, Nbkp, plotoption)

if nargin ==3
    plotoption=0;
end

% Enumerate all the combinaisons
Npt=length(X);                        % Number of samples
IdVect=linspace(1, Npt, Npt);         % total number of samples (=Npt)
CombXbkp=combnk(IdVect, Nbkp) ;       % Matlab combinaison function
CombXbkp=CombXbkp( find(CombXbkp(:,1)==1 & CombXbkp(:,end)==Npt) , :) ;  % keep only bkp 1 and 10 (start/end)

% Loop over the number of combinaisons
BestErr=1e10;   % initialize best fitting error
for c=1:length(CombXbkp(:,1))
   Err=0 ;      % initialize fitting error
   for bkp=1:Nbkp-1
       %  extract bkpoints info and slope
       id0=CombXbkp(c, bkp);
       id1=CombXbkp(c, bkp+1);
       X0=X(id0);
       X1=X(id1);
       Y0=Y(id0);
       Y1=Y(id1);
       m=(Y1-Y0)/(X1-X0);
       % compute the err in the interval
       for pt=id0:id1-1
          Err=Err + ( Y0 + m*(X(pt)-X0) - Y(pt))^2; 
       end      
   end
   Err=Err + ( Y0 + m*(X(end)-X0) - Y(end))^2; % last point
   
   % Store the best combinaison
   Err=sqrt(1/Npt*Err);
   if Err < BestErr
      BestErr=Err;
      BestXbkp=c;
   end
end

% Return Best Solution
Xbkp=CombXbkp(BestXbkp, :);

% Plot 
if plotoption==1
   disp([ 'Fitting Error = ', num2str(BestErr) ] );
   Ymod=zeros(Npt,1);
   for bkp=1:Nbkp-1
       %  extract bkpoints info and slope
       id0=Xbkp(bkp);
       id1=Xbkp(bkp+1);
       X0=X(id0);
       X1=X(id1);
       Y0=Y(id0);
       Y1=Y(id1);
       m=(Y1-Y0)/(X1-X0);
       % compute the err in the interval
       for pt=id0:id1-1
          Ymod(pt)= Y0 + m*(X(pt)-X0);
       end      
   end
   Ymod(end)= Y0 + m*(X(end)-X0); % last point
   
   figure()
   plot(X, Y, 'linewidth', 3 , 'marker' , 'o' )
   xlabel('X')
   ylabel('Y')
   title('Fitting Results')
   grid on
   hold on
   plot(X, Ymod, 'linewidth', 2, 'linestyle', '--' )
end

end