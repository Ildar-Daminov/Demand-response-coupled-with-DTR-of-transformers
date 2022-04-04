function [PUL_optim,Pflex_max,Eflex_max,theta_h_optim,theta_0_optim,time]=Linearilized_optim(LOADprofile,AMB,SOCflex_init,SOCflex_end)


%% Load Simulation profiles
%--------------------------
[LOADprofile,~]=create_load(LOADprofile,60);% Convert transformer loadings
% [AMB,~]=create_load(AMB,60);% Convert ambient temperature with function "mean"
dT=1; % time step of 1 minute
HorizonD = length(LOADprofile); % horizon in hours
HorizonT=length(LOADprofile)*60;% horizon in minutes
theta_a=AMB;    % Change of variable for ambiant temperature

%% Similation Parameters (to set by the user)
%-------------------------------------------
SOCflex_min=0;   %battery minimum SOC in %
SOCflex_max=1;  %battery maximum SOC in %
% SOCflex_init=1;  %battery initial SOC in %
% SOCflex_end=0;  %battery initial SOC in %

% Transformer limits
current_limit=1.5; % current limit,pu
theta_h_limit = 120; % winding hot spot tempearature,degC
theta_0_limit= 105; % top oil temperature,degC

% Transformer thermal characterisitcs
Nominal_rating=500; %KVA
delta_theta_or = 55;     % Top-oil (in tank) temperature rise in steady state at rated losses (no-load losses + load losses),K
delta_theta_hr = 23;     % Hot-spot-to-top-oil (in tank) gradient at rated current, K
tao_0 = 180;             % Average oil time constant, min
tao_w = 4;               % Winding time constant, min
R = 5;                   % Ratio of load losses at rated current to no-load losses
x = 0.8;                 % Exponential power of total losses versus top-oil (in tank) temperature rise (oil exponent)
y = 1.6;                 % Exponential power of current versus winding temperature rise (winding exponent)
k11 = 1;                 % Thermal model constant
k21 = 1;                 % Thermal model constant
k22 = 2;                 % Thermal model constant

% Build simulation profiles
Pload=Nominal_rating*LOADprofile;

%% Extract PWL parameters for the transformer temp. model -- DO NOT CHANGE
%--------------------------------------------------------------------------
disp('Compute PWL coefficients...')
% Generate/Fit non linear/convex functions in the temperature model
Npt=20;     % Number of samples for model fitting
Npwl=6;     % Number of pwl segments
Npwl_AEQ=12;
plotoption=0;   % plot fitting results - 0 : no  - 1 : yes
X=linspace(0, 1.5, Npt);  % X=Pgd/NominalRating
Y_f1=((1+X.^2*R)/(1+R)).^x*delta_theta_or;
Y_f2=X.^y*delta_theta_hr;
X_AEQ=0:5:120;
Y_AEQ=2.^((X_AEQ-98)/6);
% [Xbkp_f1]=computeXBKPbest(X, Y_f1, Npwl+1, plotoption);
% [Xbkp_f2]=computeXBKPbest(X, Y_f2, Npwl+1, plotoption);
Xbkp_f1=[ 1     4     7    10    13    16    20];
Xbkp_f2=[ 1     3     6     9    12    16    20];
%[Xbkp_AEQ]=computeXBKPbest(X_AEQ, Y_AEQ, Npwl_AEQ+1, plotoption);
Xbkp_AEQ=[1    11    14    16    17    18    19    20    21    22    23    24    25];
% Define parameters for PWL functions/constraints
Xk_f1=X(Xbkp_f1);
Xk_f2=X(Xbkp_f2);
Y_f1=Y_f1(Xbkp_f1);
Y_f2=Y_f2(Xbkp_f2);
%%%%%%%%%%%%%%%%%%
Xk_AEQ=X_AEQ(Xbkp_AEQ);
Y_AEQ=Y_AEQ(Xbkp_AEQ);
% Xk_AEQ=0:5:120;
% Y_AEQ=2.^((Xk_AEQ-98)/6);
% Npwl_AEQ=length(Xk_AEQ)-1;
%%%%%%%%%%%%%%%%%%
Y0_f1=Y_f1(1);
Y0_f2=Y_f2(1);
Y0_AEQ=Y_AEQ(1);
Ak_f1=[];
Ak_f2=[];
Ak_AEQ=[];
Xk_f1_max=[];
Xk_f2_max=[];
Xk_AEQ_max=[];
for k=2:Npwl+1
    Ak_f1(k-1) = (Y_f1(k)-Y_f1(k-1))/(Xk_f1(k)-Xk_f1(k-1));
    Ak_f2(k-1) = (Y_f2(k)-Y_f2(k-1))/(Xk_f2(k)-Xk_f2(k-1));
    Xk_f1_max(k-1) = Xk_f1(k)-Xk_f1(k-1);
    Xk_f2_max(k-1) = Xk_f2(k)-Xk_f2(k-1);
end
for k=2:Npwl_AEQ+1
    Ak_AEQ(k-1) = (Y_AEQ(k)-Y_AEQ(k-1))/(Xk_AEQ(k)-Xk_AEQ(k-1));
    Xk_AEQ_max(k-1) = Xk_AEQ(k)-Xk_AEQ(k-1);
end


%% Write / Solve Problem
%---------------------------------------

% Variables
Pflex=sdpvar(HorizonD, 1,'full');      % battery charge in kW
Ptr=sdpvar(HorizonD, 1,'full');      % energy imported from grid in kW  (only import  --- export/feed in forbiden)
Pflex_max=sdpvar(1,1);
Eflex_max=sdpvar(1,1);

% Transformer temporal variables
Xk_f1=sdpvar(HorizonT, Npwl,'full');        % for PWL of functions f1
Xk_f2=sdpvar(HorizonT, Npwl,'full');        % for PWL of functions f2
delta_theta_h1=sdpvar(HorizonT, 1,'full'); % component 1 of winding hot spot temperature rise,K
delta_theta_h2=sdpvar(HorizonT, 1,'full'); % Component 2 of winding hot spot temperature rise,K
theta_h=sdpvar(HorizonT, 1,'full');        % Winding hot spot temperature,degC
% theta_h_max=sdpvar(1, 1,'full');        % Winding hot spot temperature,degC
theta_0=sdpvar(HorizonT, 1,'full');        % Top oil temperature, degC
Xk_AEQ=sdpvar(HorizonT, Npwl_AEQ,'full');        % for PWL of functions f2

% Objective
disp('Write Objective Function...')
alfa_PWL=1e-3;
Objective =  (Pflex_max + Eflex_max);
Objective = Objective + alfa_PWL*( Y0_f1*HorizonT + sum(sum(repmat(Ak_f1,HorizonT,1).*Xk_f1)));
Objective = Objective + alfa_PWL*( Y0_f2*HorizonT + sum(sum(repmat(Ak_f2,HorizonT,1).*Xk_f2)));

% Constraints
disp('Write Dispatch Constraints...')
Constraints=[];
Constraints=[Constraints ; Pflex <= Pflex_max ];
Constraints=[Constraints ; Pflex >= - Pflex_max ];
Constraints=[Constraints ; Ptr >= 0 ];
Constraints = [Constraints , Ptr <= current_limit*Nominal_rating  ];

%Constraints=[Constraints ; Pflex == 0 ];
Constraints=[Constraints ; Pflex + Ptr  == Pload ]; % power balance gen=load
for t=1:HorizonD
    Constraints=[ Constraints ; SOCflex_init*Eflex_max - sum(Pflex(1:t)) <=SOCflex_max*Eflex_max];
    Constraints=[ Constraints ; SOCflex_init*Eflex_max - sum(Pflex(1:t)) >=SOCflex_min*Eflex_max];
end
Constraints=[ Constraints ;  SOCflex_init*Eflex_max - sum(Pflex)>= SOCflex_end*Eflex_max ];

% PWL
disp('Write PWL Constraints...')
for k=1:Npwl
    Constraints = [Constraints, Xk_f1(:,k) <= Xk_f1_max(k) ];
    Constraints = [Constraints, Xk_f1(:,k) >= 0  ];
    Constraints = [Constraints, Xk_f2(:,k) <= Xk_f2_max(k) ];
    Constraints = [Constraints, Xk_f2(:,k) >= 0  ];
end
for k=1:Npwl_AEQ
    Constraints = [Constraints, Xk_AEQ(:,k) <= Xk_AEQ_max(k) ];
    Constraints = [Constraints, Xk_AEQ(:,k) >= 0  ];
end
Mat_Td_Tt=zeros(HorizonT, HorizonD);
for t=1:HorizonD
    Mat_Td_Tt((HorizonT/HorizonD*(t-1)+1):HorizonT/HorizonD*t, t)=1;
end
Constraints = [Constraints , sum(Xk_f1,2)*Nominal_rating == Mat_Td_Tt*Ptr  ];
Constraints = [Constraints,  sum(Xk_f2,2)*Nominal_rating == Mat_Td_Tt*Ptr  ];
Constraints = [Constraints,  sum(Xk_AEQ,2) == theta_h  ];
% Transformer constraints
disp('Write Transformer Constraints...')
%% Write MAT
MAT_V_theta0=sparse(HorizonT-1,HorizonT);
MAT_V_deltah1=sparse(HorizonT-1,HorizonT);
MAT_V_deltah2=sparse(HorizonT-1,HorizonT);
for row=1:HorizonT-1
    MAT_V_theta0(row, row) =-1+dT/(k11*tao_0);
    MAT_V_theta0(row, row+1) =1;
    MAT_V_deltah1(row, row) =-1+dT/(k22*tao_w);
    MAT_V_deltah1(row, row+1) =1;
    MAT_V_deltah2(row, row) =-1+dT/(k22*tao_0);
    MAT_V_deltah2(row, row+1) =1;
end
%%
% Initialisation
Constraints=[ Constraints ;  theta_0(1) == Y0_f1 + sum(Ak_f1.*Xk_f1(1,:)) + theta_a(1)];
Constraints=[ Constraints ;  delta_theta_h1(1) == k21*(Y0_f2 + sum(Ak_f2.*Xk_f2(1,:)))];
Constraints=[ Constraints ;  delta_theta_h2(1) == (k21-1)*(Y0_f2 + sum(Ak_f2.*Xk_f2(1,:)))];

    
% MAT formulation
Constraints=[ Constraints ;  MAT_V_theta0*theta_0 == (dT/(k11*tao_0))*( Y0_f1 + sum(repmat(Ak_f1,HorizonT-1,1).*Xk_f1(2:end,:),2) + theta_a(2:end)) ];
Constraints=[ Constraints ;  MAT_V_deltah1*delta_theta_h1 == dT/(k22*tao_w)*(k21*(Y0_f2 + sum(repmat(Ak_f2,HorizonT-1,1).*Xk_f2(2:end,:),2)))];
Constraints=[ Constraints ;  MAT_V_deltah2*delta_theta_h2 == dT/(k22*tao_0)*((k21-1)*(Y0_f2 + sum(repmat(Ak_f2,HorizonT-1,1).*Xk_f2(2:end,:),2)))];

%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  You can active/disactivate the
%  constraint to consider the AEQ or not
Constraints=[ Constraints ; sum(sum(Y0_AEQ + sum(repmat(Ak_AEQ, HorizonT,1).*Xk_AEQ)))/HorizonT <=1 ];

Constraints=[ Constraints ; theta_h == theta_0 + delta_theta_h1 - delta_theta_h2];
Constraints=[ Constraints ; theta_0 >= 0 ]; % Top oil temperature constraint
Constraints=[ Constraints ; theta_h >= 0 ]; % Hot spot temperature constraint
% Constraints=[ Constraints ; theta_h_max >= theta_h ]; % Hot spot temperature constraint
Constraints=[ Constraints ; theta_0 <= theta_0_limit ]; % Top oil temperature constraint
Constraints=[ Constraints ; theta_h <= theta_h_limit ]; % Hot spot temperature constraint

%% Solve & Extract Results
%-------------------------
% ops = sdpsettings('solver', 'cplex','savesolveroutput',1 ,'verbose',1);
ops = sdpsettings('solver', 'linprog','savesolveroutput',1 ,'verbose',1);
tic
sol=optimize(Constraints,Objective, ops);
time=toc

[model] = export(Constraints ,Objective );
[Nvar, Ncol]=size(model.A);
disp(sol.info)
disp(['Solving Time (s): ', num2str(sol.solvertime)])
Pflex=value(Pflex);
Ptr=value(Ptr);
PUL_optim=Convert2minute(Ptr/Nominal_rating);
Pflex_max=value(Pflex_max);
Eflex_max=value(Eflex_max);
Xk_AEQ=value(Xk_AEQ);
Xk_f1=value(Xk_f1);
Xk_f2=value(Xk_f2);
delta_theta_h1=value(delta_theta_h1);
delta_theta_h2=value(delta_theta_h2);
theta_h=value(theta_h);
theta_h_optim=theta_h;
theta_0_optim=theta_0;
theta_0=value(theta_0);
SOCflex=zeros(HorizonD+1,1);
SOCflex(1)=SOCflex_init;
for t=1:HorizonD
    SOCflex(t+1)=SOCflex_init - sum(Pflex(1:t))/Eflex_max;
end

%% PWL validation with transformer temperature model
%--------------------------------------------------------------------
[~,~,AEQ_ref,theta_h_ref,theta_0_ref]=distribution_transformer(Mat_Td_Tt*Ptr/Nominal_rating,theta_a);
[~,~,AEQ,theta_h_NoFlex,theta_0_NoFlex,]=distribution_transformer(Mat_Td_Tt*Pload/Nominal_rating,theta_a);
Y_f1_pwl=zeros(HorizonT,1);
Y_f2_pwl=zeros(HorizonT,1);
Y_f1_ref=zeros(HorizonT,1);
Y_f2_ref=zeros(HorizonT,1);
Ptr_resample=Mat_Td_Tt*Ptr;
for t=1:HorizonT
    Y_f1_pwl(t)=Y0_f1 + sum(Ak_f1.*Xk_f1(t,:)) ;
    Y_f2_pwl(t)=Y0_f2 + sum(Ak_f2.*Xk_f2(t,:)) ;
    Y_f1_ref(t)=((1+ (Ptr_resample(t)/Nominal_rating).^2*R)/(1+R))^x*delta_theta_or;
    Y_f2_ref(t)=(Ptr_resample(t)/Nominal_rating)^y*delta_theta_hr;
end
% figure()
% subplot(2,2,1), plot((0:(HorizonT-1))/60,[Y_f1_pwl Y_f1_ref ])
% subplot(2,2,2), plot((0:(HorizonT-1))/60,[Y_f2_pwl Y_f2_ref ])
disp(['PWL Error (nRMSE) for teta_o (%) = ', num2str( sqrt(1/HorizonT*sum((theta_0-theta_0_ref).^2))/mean(theta_0_ref) *100  ) ] )
disp(['PWL Error (nRMSE) for teta_h (%) = ', num2str( sqrt(1/HorizonT*sum((theta_h-theta_h_ref).^2))/mean(theta_h_ref) *100  ) ] )
disp('==============================================')
disp(['AEQ teta_h_ref (%) = ', num2str( sum(2.^((theta_h_ref-98)/6))/HorizonT ) ] );
disp(['AEQ teta_h_OPTI (%) = ', num2str( sum(2.^((theta_h-98)/6))/HorizonT ) ] );
% disp(['AEQ PWL (%) = ', num2str(sum(sum(Y0_AEQ + sum(repmat(Ak_AEQ, HorizonT,1).*Xk_AEQ)))/HorizonT ) ] );


%% Plot/Display Results
%------------------------
disp(['Flex Power (kW) = ', num2str(Pflex_max)])
disp(['Flex Energy (kWh) = ', num2str(Eflex_max)])
figure()
subplot(2,3,1), stairs(0:(HorizonD-1),[ Pload Ptr ], 'linewidth', 2)
xlabel('Time (h)')
ylabel('Power (kW)')
legend('No Flex', 'With Flex')
title('Power profiles')
grid on
subplot(2,3,2), stairs(0:(HorizonD),SOCflex*100, 'linewidth', 2)
xlabel('Time (h)')
ylabel('SOCflex (%)')
title('Flex SOC ')
grid on
subplot(2,3,3), stairs(0:(HorizonD-1),Pflex , 'linewidth', 2 )
xlabel('Time (h)')
ylabel('Power (kW)')
title('Flex Power')
grid on
subplot(2,3,4), plot((0:(HorizonT-1))/60,[ theta_0_NoFlex  theta_0],  'linewidth', 2)
hold on
subplot(2,3,4), plot(0:(HorizonT-1)/60, ones(HorizonD,1)*theta_0_limit, 'color', 'k', 'linewidth', 2 , 'linestyle', '--')
xlabel('Time (h)')
ylabel('T (deg)')
legend('No Flex', 'With Flex', 'Limit')
title('Oil Temperature')
grid on
subplot(2,3,5), plot((0:(HorizonT-1))/60,[ theta_h_NoFlex  theta_h],  'linewidth', 2)
hold on
subplot(2,3,5), plot(0:(HorizonT-1)/60, ones(HorizonD,1)*theta_h_limit, 'color', 'k', 'linewidth', 2 , 'linestyle', '--')
xlabel('Time (h)')
ylabel('T (deg)')
legend('No Flex', 'With Flex', 'Limit')
title('Hot spot Temperature')
grid on
subplot(2,3,6), plot((0:(HorizonT-1))/60,[theta_a ], 'linewidth', 2)
xlabel('Time (h)')
ylabel('T (deg)')
title('Ambient Temperature')
grid on
end