function [PUL_optim,flex_KW,flex_KWh,time]=fmincon_time_test(LOADprofile,AMB,x0)

close all
%% Decoupled formulation with transformer constraints
% Convert data into hour format
[LOADprofile,~]=create_load(LOADprofile,60);% Convert transformer loadings
[AMB,~]=create_load(AMB,60);% Convert ambient temperature with function "mean"
HorizonD = length(AMB); % horizon in hours
HorizonT=length(AMB)*60;% horizon in minutes
Nominal_rating=500; %KVA



%% Similation Parameters (to set by the user)
%-------------------------------------------
SOCflex_min=0;   %battery minimum SOC in %
SOCflex_max=1;  %battery maximum SOC in %
SOCflex_init=1;  %battery initial SOC in %
SOCflex_end=0;  %battery initial SOC in %

% Transformer limits
% % current_limit=1.5; % current limit,pu
theta_h_limit=120; % winding hot spot tempearature,degC
theta_0_limit=105; % top oil temperature,degC

% Build simulation profiles
%-------------------------------------------
Pload=LOADprofile;


%% Optimization problem and the optimization variables
% Create a minimization problem
Flex_minimization=optimproblem('ObjectiveSense','minimize');

% Variables from LP problem fromulation
Pflex_max=optimvar('Pflex_max',1);
Eflex_max=optimvar('Eflex_max',1);
Pflex=optimvar('Pflex',HorizonD); % battery discharge in kW
Ptr=optimvar('Ptr',HorizonD,'LowerBound',0); % transformer
%% Objective function
Flex_min = Pflex_max + Eflex_max;   % Electricity bill penalty for peak + price for purchased energy

% showexpr(Flex_min)

% set objective
Flex_minimization.Objective = Flex_min;

%% Constraints
% type distribution_transformer
[HST_max,TOT_max,AEQ,~,~]=fcn2optimexpr(@distribution_transformer,Ptr/Nominal_rating,AMB);
Flex_minimization.Constraints.Ageing=AEQ<=1;
Flex_minimization.Constraints.Hst=HST_max<=120;
Flex_minimization.Constraints.Tot=TOT_max<=105;
Flex_minimization.Constraints.current=Ptr<=1.5*Nominal_rating;


Flex_minimization.Constraints.balance=Pflex + Ptr == Pload ;% Balance
Flex_minimization.Constraints.Pflex1=Pflex>=0;
Flex_minimization.Constraints.Pflex2=Pflex<=Pflex_max;

% SOC_calculation
for t=1:HorizonD
    eval(['Flex_minimization.Constraints.t',num2str(t),' = SOCflex_init*Eflex_max - sum(Pflex(1:',num2str(t),')) <=SOCflex_max*Eflex_max;']);
    eval(['Flex_minimization.Constraints.T',num2str(t),' = SOCflex_init*Eflex_max - sum(Pflex(1:',num2str(t),')) >=SOCflex_min*Eflex_max;']);
end
Flex_minimization.Constraints.SOCbat_end= SOCflex_init*Eflex_max - sum(Pflex)>= SOCflex_end*Eflex_max ;


%     showproblem(Flex_minimization)
%% OPtimal solution
% options for the optimization algorithm, here we set the max time it can run for
% call the optimization solver to find the best solution
% [sol,TotalCost,exitflag,output] = solve(Flex_minimization,x0,'Options',options);


options=optimset('disp','iter','LargeScale','off','TolFun',.001,'MaxIter',10000000,'MaxFunEvals',100000000);
% options.Algorithm='sqp-legacy';
tic
[sol,fval,~,output] = solve(Flex_minimization,x0,'Options',options);
time=toc
PUL_optim=PUL_to_1min(sol.Ptr/Nominal_rating,60);
flex_KW=sol.Pflex_max;
flex_KWh=sol.Eflex_max;

%% Plot/Display Results
%------------------------
disp(['Flex Power (kW) = ', num2str(sol.Pflex_max)])
disp(['Flex Energy (kWh) = ', num2str(sol.Eflex_max)])
% figure()
% subplot(2,3,1), stairs(0:(HorizonD-1),[ Pload sol.Ptr ], 'linewidth', 2)
% hold on
% subplot(2,3,1), stairs(0:(HorizonD-1), ones(HorizonD,1)*1.5*Nominal_rating, 'color', 'k', 'linewidth', 2 , 'linestyle', '--')
% 
% xlabel('Time (h)')
% ylabel('Power (kW)')
% legend('No Flex', 'With Flex', 'Limit')
% title('Power profiles')
% grid on
% 
% SOCflex=zeros(HorizonD+1,1);
% SOCflex(1)=SOCflex_init;
% for t=1:HorizonD
%     SOCflex(t+1)=SOCflex_init - sum(sol.Pflex(1:t))/sol.Eflex_max;
% end
% 
% subplot(2,3,2), stairs(0:(HorizonD),SOCflex*100, 'linewidth', 2)
% 
% xlabel('Time (h)')
% ylabel('SOCflex (%)')
% title('Flex SOC ')
% grid on
% subplot(2,3,3), stairs(0:(HorizonD-1),sol.Pflex , 'linewidth', 2 )
% xlabel('Time (h)')
% ylabel('Power (kW)')
% title('Flex Power')
% grid on
% 
% [theta_h,theta_0,AEQ_optim]=distribution_transformer_short(sol.Ptr/Nominal_rating,AMB);
% [theta_h_NoFlex,theta_0_NoFlex,AEQ]=distribution_transformer_short(Pload/Nominal_rating,AMB);
% 
% subplot(2,3,4), plot((0:(HorizonT-1))/60,[ theta_0_NoFlex  theta_0],  'linewidth', 2)
% hold on
% subplot(2,3,4), plot(0:(HorizonT-1)/60, ones(HorizonD,1)*theta_0_limit, 'color', 'k', 'linewidth', 2 , 'linestyle', '--')
% xlabel('Time (h)')
% ylabel('T (deg)')
% legend('No Flex', 'With Flex', 'Limit')
% title('Oil Temperature')
% grid on
% subplot(2,3,5), plot((0:(HorizonT-1))/60,[ theta_h_NoFlex  theta_h],  'linewidth', 2)
% hold on
% subplot(2,3,5), plot(0:(HorizonT-1)/60, ones(HorizonD,1)*theta_h_limit, 'color', 'k', 'linewidth', 2 , 'linestyle', '--')
% subplot(2,3,5), plot(0:(HorizonT-1)/60, ones(HorizonD,1)*98, 'color', 'k', 'linewidth', 1 , 'linestyle', '--')
% 
% xlabel('Time (h)')
% ylabel('T (deg)')
% legend('No Flex', 'With Flex', 'Limit', 'Design')
% title('Hot spot Temperature')
% grid on
% subplot(2,3,6), stairs(0:(HorizonD-1),[AMB ], 'linewidth', 2)
% xlabel('Time (h)')
% ylabel('T (deg)')
% title('Ambient Temperature')
% grid on

end
