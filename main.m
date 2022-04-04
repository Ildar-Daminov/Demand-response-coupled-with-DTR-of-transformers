clc
clear all
close all

%% Goal of the script
% This scripts reproduces the Figures from the article [1]:

% If you use this code, please cite this article:
% [1] Ildar Daminov, Rémy Rigo-Mariani, Raphaël Caire, Anton Prokhorov,
% Marie-Cécile Alvarez-Hérault, "Demand Response Coupled with Dynamic
% Thermal Rating for Increased Transformer Reserve and Lifetime" in Energies
% 2021, https://doi.org/10.3390/en14051378

% Other articles on this topic are available:
% https://www.researchgate.net/profile/Ildar-Daminov-2

% Note that the figures generated in this script and those given in the
% article may differ as latter had been additionally redrawn
% for a publication.

% Each section (Plotting the Figure X) is independent from each other. So
% you may launch the entire script (using the button "Run") to get all
% figures at one moment or you may launch a special section (using the
% button "Run Section" at the top)to get a specific figure

% Execution time of entire script ≈ 15 min

tic
%% Plotting the Figure 1
% Figure name: Scope of the paper with regard to the literature survey.

% Figure 1 in the article [1] was ploted without using MATLAB

%% Plotting the Figure 2
% Figure name: Case study—(a) outdoor secondary substation; (b) hourly
% load in kilovolt ampere (kVA) and monthly maximum ambient temperature
% in Grenoble, France

% Figure 2 a in the article [1] was ploted without using MATLAB

% Figure2 (b) is plotted as follows

% Load aggregated power rpofiles of houses simulated via the application
% House Load % https://fr.mathworks.com/matlabcentral/fileexchange/63375-house-load-electricity
% for the profile of ambient temperature in grenoble France (see AMB)
load('Aggregated_load_profile_100_houses.mat') % in W!
Load_agg=Load_agg/1000; % kVA, cosphi=1

load('Ambient_temperature_Grenoble.mat') % in °C
% Note that ambient temperature was processed.

% Prepare a time vector
t1 = datetime(2019,1,1,0,0,0,'Format','HH:SS');
t2 = datetime(2019,12,31,23,59,0,'Format','HH:SS');
time = t1:minutes(1):t2; time=time';

% Create figure
figure('InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized');

% Create axes
axes1 = axes('Position',...
    [0.102707749766573 0.11 0.782212885154062 0.858973747016706]);
hold(axes1,'on');
colororder([0 0 0]);

% Activate the left side of the axes
yyaxis(axes1,'left');

% Create plot
plot(time,Load_agg,'DisplayName','Hour load, kVA','LineWidth',1,...
    'Color',[0 0.4470 0.7410]);

% Create ylabel
ylabel('Transformer load, kVA');

% Preserve the Y-limits of the axes
% ylim(axes1,[0 450]);

% Set the remaining axes properties
set(axes1,'YColor',[0 0.447058823529412 0.741176470588235]);

% Activate the right side of the axes
yyaxis(axes1,'right');

% Create plot
plot(time,AMB,'DisplayName','Month maximum of ambient temperature,°C',...
    'LineWidth',3,'Color',[0.8500 0.3250 0.0980]);

% Create ylabel
ylabel('Ambient temperature, °C');

% Preserve the Y-limits of the axes
% ylim(axes1,[15 40]);

% Set the remaining axes properties
set(axes1,'YColor',[0.85 0.325 0.098]);

% Preserve the Z-limits of the axes
% zlim(axes1,[-1 0]);

box(axes1,'on');

hold(axes1,'off');

% Set the remaining axes properties
set(axes1,'FontSize',20,'GridColor',[0 0.447058823529412 0.741176470588235]);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.369894400614767 0.505554305650581 0.336848965572814 0.0943818897383089],...
    'FontSize',18,...
    'EdgeColor',[1 1 1]);

%% Plotting the Figure 3
% Figure name: Preliminary thermal results for yearly simulation: load
% growing from 1% to 100% of nominal rating (reserve in p.u.)

clc;clear all % clear a command window and a workspace

% Load aggregated power rpofiles of houses simulated via the application
% House Load % https://fr.mathworks.com/matlabcentral/fileexchange/63375-house-load-electricity
% for the profile of ambient temperature in grenoble France (see AMB)
load('Aggregated_load_profile_100_houses.mat') % in W!
Load_agg=Load_agg/1000; % kVA, cosphi=1

load('Ambient_temperature_Grenoble.mat') % in °C
% Note that ambient temperature was processed.

% Prepapre the data for thermal simulations of distribution transformer
Nominal_rating=500; % kVA
PUL_init=Load_agg/Nominal_rating; % Initial Per Unit Load (PUL) 

Reserve_margins=0.01:0.01:1; % load growth pu = reserve values

% Perform thermal simulations 
for given_reserve=1:length(Reserve_margins)

    % Select the reserve for study
    studied_reserve=Reserve_margins(given_reserve);

    % Reconstruct the load profile in pu for studied_reserve
    PUL=PUL_init+studied_reserve;

    % Find HST, TOT and AEQ 
    [HST_max,TOT_max,AEQ,~,~]=distribution_transformer(PUL,AMB);

    % Save results 
    Results.HST_max(given_reserve,1)=HST_max;
    Results.TOT_max(given_reserve,1)=TOT_max;
    Results.AEQ(given_reserve,1)=AEQ;
    Results.Load_max(given_reserve,1)=max(PUL);

end % end of "for given_reserve=1:length(Reserve_margins)"

% Create figure
figure1 = figure('InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized');

% Create axes
axes1 = axes('Position',[0.1 0.1 0.636363636363636 0.8]);
hold(axes1,'on');

subplot(3,1,1) % create a top figure out of three

% Plot temperatures 
plot1 = plot(Reserve_margins,[Results.TOT_max Results.HST_max],'LineWidth',3);
set(plot1(1),'DisplayName','Maximal TOT, °C','LineStyle','--',...
    'Color',[0.847058832645416 0.321568638086319 0.0941176488995552]);
set(plot1(2),'DisplayName','Maximal HST, °C','Color',...
    [0.847058832645416 0.321568638086319 0.0941176488995552]);

% Create datatip
datatip(plot1(1),'DataIndex',74,'Location','southeast');

% Create datatip
datatip(plot1(2),'DataIndex',49,'Location','northwest');

% Create datatip
datatip(plot1(2),'DataIndex',29,'Location','northwest');

% Create ylabel
ylabel('Maximal transformer temperature, °C');

% Create xlabel
xlabel('Reserve, pu');

% Show the legend
legend ('show')

subplot(3,1,2) % create a middle figure out of three
% Create line
line1=line(Reserve_margins,Results.Load_max,'DisplayName','Maximal load, pu','LineWidth',3,...
    'Color',[0.494117647409439 0.184313729405403 0.556862771511078]);

% Create datatip
datatip(line1,'DataIndex',64);

% Create ylabel
ylabel('Maximal transformer loading,pu');
ylim([0.8,2])

% Create xlabel
xlabel('Reserve, pu');

% Show the legend
legend ('show')

subplot(3,1,3) % create a bottom figure out of three

% Create plot
plot2 = plot(Reserve_margins,Results.AEQ,'DisplayName','Equivalent annual ageing, pu',...
    'LineWidth',3,...
    'Color',[0.301960796117783 0.745098054409027 0.933333337306976]);

% Create datatip
datatip(plot2,'DataIndex',71,'Location','northwest');

% Create ylabel
ylabel('Equivalent annual ageing, pu');

% Create xlabel
xlabel('Reserve, pu');

% Show the legend
legend ('show')

%% Plotting the Figure 4
% Figure name: Number of days where thermal limits are reached.

clc;clear all % clear a command window and a workspace

% Load aggregated power rpofiles of houses simulated via the application
% House Load % https://fr.mathworks.com/matlabcentral/fileexchange/63375-house-load-electricity
% for the profile of ambient temperature in grenoble France (see AMB)
load('Aggregated_load_profile_100_houses.mat') % in W!
Load_agg=Load_agg/1000; % kVA, cosphi=1

load('Ambient_temperature_Grenoble.mat') % in °C
% Note that ambient temperature was processed. Only monthly  maximums 
% are kept over the duration of the whole month. This is done to make
% thermal simulations more conservative. 

% Prepapre the data for thermal simulations of distribution transformer
Nominal_rating=500; % kVA
PUL_init=Load_agg/Nominal_rating; % Initial Per Unit Load (PUL) 
Reserve_margins=0.01:0.01:1.63; % load growth pu = reserve values

% % Perform thermal simulations 
% for given_reserve=1:length(Reserve_margins)
% 
%     % Select the reserve for study
%     studied_reserve=Reserve_margins(given_reserve);
% 
%     % Reconstruct the load profile in pu for studied_reserve
%     PUL=PUL_init+studied_reserve;
% 
%     % Find HST, TOT and AEQ 
%     [HST_max,TOT_max,AEQ,HST,TOT]=distribution_transformer(PUL,AMB);
% 
%     % Save results 
%     Results.HST_max(given_reserve,1)=HST_max;
%     Results.TOT_max(given_reserve,1)=TOT_max;
%     Results.AEQ(given_reserve,1)=AEQ;
%     Results.Load_max(given_reserve,1)=max(PUL);
%     Results.PUL{given_reserve,1}=PUL;
%     Results.HST{given_reserve,1}=HST;
%     Results.TOT{given_reserve,1}=TOT;
% 
%     if HST_max>120 || TOT_max>105 || max(PUL)>1.5
%         [minutes,intervals] = profiles2minutes(PUL,HST,TOT,AEQ);
%         Results.minutes{given_reserve,1}=minutes;
%         Results.intervals{given_reserve,1}=intervals;
% 
%     else % no thermal or load violations
% 
%         Results.minutes{given_reserve,1}=0;
%         Results.intervals{given_reserve,1}=0;
%     end 
% end % end of "for given_reserve=1:length(Reserve_margins)"

load('all_intervals.mat')
Results.intervals=all_intervals;

% Find the longest interval
for given_reserve=1:length(Reserve_margins)
    interm_array=Results.intervals{given_reserve,1};
    if interm_array==0
        % The longest interval is zero
        the_longest_interval(given_reserve,1)=0;

        % Number of days is zero
        days_number(given_reserve,1)=0;

    else % if interm_array is not zero

        % Calculate the interval lengths 
        for j=1:length(interm_array(:,1))
            interval_lengths(j)=interm_array(j,2)-interm_array(j,1)+1;
        end

        % Find the longest interval
        the_longest_interval(given_reserve,1)=max(interval_lengths);

        % Days number with thermal violations
        days_number(given_reserve,1)=sum(interval_lengths/1440);

        % Delete the variable interval_length
        interval_lengths=[];

    end % if interm_array==0

end % given_reserve=1:length(Reserve_margins)

% Convert to number of days
the_longest_interval=the_longest_interval/1440;

% Create figure
figure1 = figure('InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized');

plot1=plot(Reserve_margins,[days_number,the_longest_interval],'linewidth',3);
xlabel('Reserve, pu')
ylabel('The number of days')

set(plot1(1),'DisplayName','Days number','LineStyle','-');
set(plot1(2),'DisplayName','The longest interval','LineStyle','--');

legend('show')
%% Plotting the Figure 5
% Figure name: . The procedure for finding the needed Demand Response (DR) 
% volume to interconnect the studied reserve

% Figure 5 in the article [1] was ploted without using MATLAB

% UPD (04/04/22): the MATLAB code of this algorithm will be uploaded during 
% next two weeks as it requires us doublechecking before sharing it on
% GitHub

%% Plotting the Figure 6
% Figure name: Piece-wise linearization (PWL) process: (a) functions 
% fitting and breakpoints; (b) mathematical formulation for function f

clc;clear all % clear a command window and a workspace

% Some thermal characteristics of ONAN distribution transformer per IEC60076-7
delta_theta_or = 55;     % Top-oil (in tank) temperature rise in steady state at rated losses (no-load losses + load losses),K
delta_theta_hr = 23;     % Hot-spot-to-top-oil (in tank) gradient at rated current, K
R = 5;                   % Ratio of load losses at rated current to no-load losses
x = 0.8;                 % Exponential power of total losses versus top-oil (in tank) temperature rise (oil exponent)
y = 1.6;                 % Exponential power of current versus winding temperature rise (winding exponent)


% Generate/Fit non linear/convex functions in the temperature model
Npt=20;     % Number of samples for model fitting
Npwl=6;     % Number of PWL segments
Npwl_AEQ=12; % Number of PWL segments for ageing curve 
plotoption=1;   % plot fitting results - 0 : no  - 1 : yes
X=linspace(0, 1.5, Npt);  % X=Pgd/NominalRating

% y-axis functions: f1 and f2 (temperature-related curves)
Y_f1=((1+X.^2*R)/(1+R)).^x*delta_theta_or;
Y_f2=X.^y*delta_theta_hr;

% x-axis value for AEQ curve 
X_AEQ=0:5:120;

% y-axis function for AEQ
Y_AEQ=2.^((X_AEQ-98)/6);

% Compute the best breakpoints for simple PWL
[Xbkp_f1]=computeXBKPbest(X, Y_f1, Npwl+1, plotoption);
[Xbkp_f2]=computeXBKPbest(X, Y_f2, Npwl+1, plotoption);
[Xbkp_AEQ]=computeXBKPbest(X_AEQ, Y_AEQ, Npwl_AEQ+1, 0);

% Plotting happens inside of computeXBKPbest.m

%-----------------next code for info how to calculate Ak_f1,Ak_f2 ---------
% Define parameters for PWL functions/constraints
Xk_f1=X(Xbkp_f1);
Xk_f2=X(Xbkp_f2);
Y_f1=Y_f1(Xbkp_f1);
Y_f2=Y_f2(Xbkp_f2);
Xk_AEQ=X_AEQ(Xbkp_AEQ);
Y_AEQ=Y_AEQ(Xbkp_AEQ);

%%%%%%%%%%%%%%%%%%
Y0_f1=Y_f1(1);
Y0_f2=Y_f2(1);
Y0_AEQ=Y_AEQ(1);

% Create empty variables 
Ak_f1=[];
Ak_f2=[];
Ak_AEQ=[];
Xk_f1_max=[];
Xk_f2_max=[];
Xk_AEQ_max=[];

% Calculate Ak values 
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


%% Plotting the Figure 7
% Figure name: 7. PWL performances for different numbers of blocks—(a) 
% top oil temperature; (b) hot spot temperature

clc;clear all % clear a command window and a workspace

load('fig_nRMSE.mat') % load data of precalculated figure 
% We will upload the code for calculating fig_nRMSE.mat later 

% create a figure
figure('InvertHardcopy','off','Color',[1 1 1],'DefaultAxesFontSize',14);

%---------------------------Figure 7 (a)----------------------------------
subplot(1,2,1)
plot1=plot(x_axis,[ref_TOT',TOT{1,1}',TOT{2,1}',TOT{3,1}']);
set(plot1(1),'DisplayName','ref','LineStyle','-','Color','g');
set(plot1(4),'DisplayName','C=2 -9.5%','LineStyle','-','Color',[17 17 17]/255);
set(plot1(3),'DisplayName','C=4 -2.1%','LineStyle','-','Color','k');
set(plot1(2),'DisplayName','C=8 -0.5%','LineStyle','--','Color','k');

ylabel('Top-oil temperature,°C')
xlabel('Time, hours')

legend('Location','northwest')
%---------------------------Figure 7 (b)----------------------------------
subplot(1,2,2)
plot2=plot(x_axis,[ref_HST',HST{1,1}',HST{2,1}',HST{3,1}']);
set(plot2(1),'DisplayName','ref','LineStyle','-','Color','g');
set(plot2(4),'DisplayName','C=2 -9.5%','LineStyle','-','Color',[17 17 17]/255);
set(plot2(3),'DisplayName','C=4 -2.1%','LineStyle','-','Color','k');
set(plot2(2),'DisplayName','C=8 -0.5%','LineStyle','--','Color','k');

ylabel('Hot-spot temperature,°C')
xlabel('Time, hours')

legend('Location','northwest')


%% Plotting the Figure 8
% Figure name: PWL process: (a) functions fitting and breakpoints for fAEQ; 
% (b) aging function—fAEQ error versus θh error

clc;clear all % clear a command window and a workspace

%---------------------------Figure 8 (a)----------------------------------

% Number of PWL segments for ageing curve 
Npwl_AEQ=12;

% plot fitting results - 0 : no  - 1 : yes
plotoption=1; 

% x-axis value for AEQ curve 
X_AEQ=0:5:120;

% y-axis function for AEQ
Y_AEQ=2.^((X_AEQ-98)/6);

% Plot figure 
[~]=computeXBKPbest(X_AEQ, Y_AEQ, Npwl_AEQ+1, plotoption);


%---------------------------Figure 8 (b)----------------------------------
% Select the three hot spot temperature (HST) as a reference 
HST_ref=[60,80,120]; % °C

% Calculate the aging equivalent for reference HST
AEQ_ref=2.^((HST_ref-98)/6); 

% Create the initial deviations:  +-5% from 100 %
Deviations=[95:1:105]';

% Calculate absolute values of HST for +-5% deviations
for i=1:length(HST_ref)
    HST_deviations(i,:)=HST_ref(i)*Deviations/100;
end

% Calculate the corresponding AEQ deviations
AEQ_deviations=2.^((HST_deviations-98)/6);

% Convert AEQ_deviations from pu to % (relative to 100%)
for i=1:length(AEQ_deviations)
    AEQ_deviations(:,i)=AEQ_deviations(:,i)./AEQ_ref'*100;
end

% Create figure
figure('InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized',...
    'DefaultAxesFontSize',14);

Relative_error=Deviations-100;
Relative_error_AEQ=AEQ_deviations-100;
plot1=plot(Relative_error,Relative_error_AEQ);
xlabel('Relative error HST, %')
ylabel('Relative error AEQ, %')

set(plot1(1),'DisplayName','HST=60 °C - AEQ=0.01 pu','LineStyle','-',...
    'LineWidth',3,'Color',[17 17 17]/255);
set(plot1(2),'DisplayName','HST=80 °C - AEQ=0.12 pu','LineStyle','--',...
    'LineWidth',3,'Color','k');
set(plot1(3),'DisplayName','HST=120 °C- AEQ=12.69 pu','LineStyle','-',...
    'LineWidth',3,'Color','g');

% Show the legend
legend('show')

%% Plotting the Figure 9
% Figure name: Simulation with single and multiple time sets: 
% (a) DR flexibility power; (b) transformer temperatures.

% The code is in prepation. It will be uploaded soon

%% Plotting the Figure 10
% Figure name: Validation run with “energy shifting”: 
% (a) transformer loadings; (b) temperatures

clc;clear all % clear a command window and a workspace


load LOADprofile.mat % normalized load profile in pu

% Convert
[LOADprofile]=Convert2minute(LOADprofile);

% Load ambient temperature 
load('initial_data.mat','AMB');
AMB=AMB+15; % make an ambient temperature hotter for 15 degc


% Prepare the intial conditions
ageing_constraint_on=0; % consider the transformer ageing 1- yes 0 - no

% Select the "energy shifting" mode 
SOCflex_init=0.5;  % battery initial SOC in 50 %
SOCflex_end=0.5;  % battery initial SOC in 50 %
 

% Do a validation run
[PUL_optim,Pflex_max,Eflex_max,theta_h_optim,theta_0_optim]=...
    validation_run(LOADprofile,AMB,SOCflex_init,SOCflex_end,...
    ageing_constraint_on);

% Note that additional (to article [1]) figures appear 

disp('-------------------------------------Attention to figure 10!---------------------------------')
disp('   The flexible power (kW) and energy (kWh)  is a bit different than in the article Energies [1]')
disp('   This is because in this version we use the internal MATLAB solver linprog and not a cplex (external solver) as in [1].') 
disp('This allows users, who does not have CPLEX, to run the MATLAB code. So we decided to keep the version which can be used for larger number of users')
disp('          Anyway, this small difference does not change the article conclusions')
disp('-------------------------------------Attention to figure 10!----------------------------------')
%% Plotting the Figure 11
% Figure name: Validation run with “energy shedding”: 
% (a) transformer loadings; (b) temperatures

clc;clear all % clear a command window and a workspace


load LOADprofile.mat % normalized load profile in pu

% Convert load from 24-hours format to 1440-minute format
[LOADprofile]=Convert2minute(LOADprofile);

% Load ambient temperature 
load('initial_data.mat','AMB');
AMB=AMB+15; % make an ambient temperature hotter for 15 degc


% Prepare the intial conditions
ageing_constraint_on=0; % consider the transformer ageing 1- yes 0 - no

% Select the "energy shedding" mode 
SOCflex_init=1;  % battery initial SOC in 100 %
SOCflex_end=0;  % battery initial SOC in 0 %

% Do a validation run
[PUL_optim,Pflex_max,Eflex_max,theta_h_optim,theta_0_optim]=...
    validation_run(LOADprofile,AMB,SOCflex_init,SOCflex_end,...
    ageing_constraint_on);

% Note that additional (to article [1]) figures appear
disp('-------------------------------------Attention to figure 11!---------------------------------')
disp('   The flexible power (kW) and energy (kWh)  is a bit different than in the article Energies [1]')
disp('   This is because in this version we use the internal MATLAB solver linprog and not a cplex (external solver) as in [1].') 
disp('This allows users, who does not have CPLEX, to run the MATLAB code. So we decided to keep the version which can be used for larger number of users')
disp('          Anyway, this small difference does not change the article conclusions')
disp('-------------------------------------Attention to figure 11!----------------------------------')

%% Plotting the Figure 12
% Figure name: Validation run with “energy shedding” + ageing constraint: 
% (a) transformer loadings; (b) temperatures

clc;clear all % clear a command window and a workspace


load LOADprofile.mat % normalized load profile in pu

% Convert
[LOADprofile]=Convert2minute(LOADprofile);

% Load ambient temperature 
load('initial_data.mat','AMB');
AMB=AMB+15; % make an ambient temperature hotter for 15 degc


% Prepare the intial conditions
ageing_constraint_on=1; % consider the transformer ageing 1- yes 0 - no

% Select the "energy shedding" mode 
SOCflex_init=1;  % battery initial SOC in 100 %
SOCflex_end=0;  % battery initial SOC in 0 %

% Do a validation run
[PUL_optim,Pflex_max,Eflex_max,theta_h_optim,theta_0_optim]=...
    validation_run(LOADprofile,AMB,SOCflex_init,SOCflex_end,ageing_constraint_on);

% Note that additional (to article [1]) figures appear 

disp('-------------------------------------Attention to figure 12!---------------------------------')
disp('   The flexible power (kW) and energy (kWh)  is a bit different than in the article Energies [1]')
disp('This is because in this version we use the internal MATLAB solver "linprog" and not a CPLEX (external solver) as in [1].') 
disp('This allows users, who does not have CPLEX, to run the MATLAB code. So we decided to keep the version which can be used for larger number of users')
disp('          Anyway, this small difference does not change the article conclusions')
disp('-------------------------------------Attention to figure 12----------------------------------')

%% Plotting the Figure 13
% Figure name: One-week profiles in January, comparison of the base case 
% and Dynamic Thermal Rating (DTR)/DR in “energy shedding mode”: 
% (a) transformer loading; (b) hot-spot temperature

clc;clear all % clear a command window and a workspace

% Load results for 75% reserve (+existing loading 86%)
load('result_AEQ_100.mat')
Base_load=PUL_init;

% Find annual profile of transformer with reserve =  initial PUL +reserve
Base_load_reserve=Base_load+result.headroom(end);

% Find optimized annual profile of transformer i.e. with DR
Optimal_loading=PUL_optim;

% Caclcaulte thermal regime of distrbution transformer 
[~,~,AEQ_init,theta_h_init,theta_0_init]=distribution_transformer(Base_load,AMB);
[~,~,AEQ,theta_h,theta_0]=distribution_transformer(Base_load_reserve,AMB);
[~,~,AEQ_optim,theta_h_optim,theta_0_optim]=distribution_transformer(Optimal_loading,AMB);

% Prepare a time vector
t1 = datetime(2019,1,1,0,0,0,'Format','HH:SS');
t2 = datetime(2019,12,31,23,59,0,'Format','HH:SS');
t = [t1:minutes(1):t2]';

% Create figure
figure('WindowState','maximized','DefaultAxesFontSize',14);

% Plot transformer loadings 
subplot(2,1,1)
plot1=plot(t,[Base_load,Base_load_reserve,Optimal_loading]);
set(plot1(1),'DisplayName','Base load','LineStyle','-',...
    'LineWidth',2,'Color',[17 17 17]/255);
set(plot1(2),'DisplayName','Base load+reserve 75%','LineStyle',':',...
    'LineWidth',2,'Color','k');
set(plot1(3),'DisplayName','Optimal loading','LineStyle','-',...
    'LineWidth',2,'Color','g');
xlabel('Time')
ylabel('Transformer loading, pu')
ylim([0,1.6]);
xlim([t(1441) t(192*60)])
legend ('show')

% Plot transformer temperatures 
subplot(2,1,2)
plot2=plot(t,[theta_h_init,theta_h,theta_h_optim]);
set(plot2(1),'DisplayName','Base load','LineStyle','-',...
    'LineWidth',2,'Color',[17 17 17]/255);
set(plot2(2),'DisplayName','Base load+reserve 75%','LineStyle',':',...
    'LineWidth',2,'Color','k');
set(plot2(3),'DisplayName','Optimal loading','LineStyle','-',...
    'LineWidth',2,'Color','g');

ylabel('Hot spot temperature, °C')
xlabel('Time')
ylim([0,140]);
xlim([t(1441) t(192*60)]) % Select one week from January 2

%% Plotting the Figure 14
% Figure name:  Obtained results for different reserve margins and DR in 
% “energy-shedding” mode: (a) DR rated power; (b) DR power share compared 
% to the added load; (c) DR rated energy; (d) DR energy share compared to 
% the total energy of load

clc;clear all % clear a command window and a workspace

% Transformer rating
Nominal_rating=500; %kVA

% Step of power connected 
Power_connected_delta=25; % kW cosphi=1

% energy for connected constant load 25 kW (5%)
Energy_connected_delta=Power_connected_delta*8760; % kWh per year 

%-----------------------AEQ≤1pu Ptr≤1.5pu θh≤120℃ θo≤105℃---------------
% Load precalculated results
load('result_AEQ_100.mat')

% extract the vector of studied reserve in pu 
reserve=result.headroom(2:end)*100;

% Reserve in kW
reserve_kW=Nominal_rating*reserve/100;

% Extract flexibility metrics: power in kW and energy in kWh
for i=1:length(reserve)
    P_flex(i,1)=result.flex_KW(i+1);
    E_flex(i,1)=result.flex_KWh(i+1);    
end
% Load aggregated load profile in W 
load('Aggregated_load_profile_100_houses.mat')  

% Convert load to kW and to hour resolution
Load_agg_kW=Convert2hours(Load_agg,60)/1000;

% Estimate the existing energy as a sum (but better as integral)
Energy_existing=sum(Load_agg_kW); %Energy in kWh

% Find a peak load 
Load_agg_peak=max(Load_agg_kW);

% Create a figure 
figure()

% Plot the flexibility in kW as a function of reserve margin
subplot(2,2,1), plot(reserve,P_flex, 'linewidth', 2)
xlabel('Reserve (%)')
ylabel('Power (kW)')
title('Power shedding, kW')
legend(' AEQ(kW)')
grid on

% Plot power shedding share in %
subplot(2,2,2), plot(reserve,(P_flex./...
    (Load_agg_peak+reserve/100*Nominal_rating))*100, 'linewidth', 2)
grid on
ylabel('Power shedding share (%)')
xlabel('Reserve (%)')
title('Power shedding/Power peak (%)')
legend('AEQ ≤1pu+Ptr≤ 1.5pu +θh≤120℃+θo≤105℃')

% plot values of energy shedding
subplot(2,2,3), plot(reserve,E_flex, 'linewidth', 2)
ylabel('Energy (kWh)') 
title('Energy shedding(kWh)')
legend('AEQ ≤1pu+Ptr≤ 1.5pu +θh≤120℃+θo≤105℃')
xlabel('Reserve (%)')
set(gca, 'YScale', 'log')
grid on

% Plot energy shedding in % 
subplot(2,2,4) 
Energy_total=Energy_existing+Energy_connected_delta.*reserve/5; %  5 = 5% of Power_connected_delta 
plot(reserve,(E_flex./(Energy_total))*100, 'linewidth', 2)
title('Energy shedding/Energy total(%)')
legend('AEQ≤1pu+Ptr≤ 1.5pu +θh≤120℃+θo≤105℃')
grid on
ylabel('The Energy shedding share(%)')
xlabel('Reserve (%)')
hold on
%------------------------Ptr≤ 1.5pu θh≤120℃ θo≤105℃--------------------
clear
% Load precalculated results for given formulation
load('result_PUL_100.mat')

% Redefine the transformer rating
Nominal_rating=500;

% Step of power connected 
Power_connected_delta=25; % kW

% Calculate the annual energy connected 
Energy_connected_delta=Power_connected_delta*8760; % KWh per year for 25 KW connected (5%)

% Extract the reserve margins 
reserve=result.headroom(2:end)*100;

% Convert into kW
reserve_kW=Nominal_rating*reserve/100;

% Extract flexibility metrics: power in kW and energy in kWh
for i=1:length(reserve)
    P_flex(i,1)=result.flex_KW(i+1);
    E_flex(i,1)=result.flex_KWh(i+1);    
end

% Load aggregated load profile in W 
load('Aggregated_load_profile_100_houses.mat')

% Convert to kW and transform into hour resolution
Load_agg_kW=Convert2hours(Load_agg,60)/1000;

% Estimate the energy of load
Energy_existing=sum(Load_agg_kW); %Energy in kWh

% Find a peak load 
Load_agg_peak=max(Load_agg_kW);

% Plot the flexibility in kW as a function of reserve margin
subplot(2,2,1),hold on
plot(reserve,P_flex, 'linewidth', 2)
xlabel('Reserve (%)')
ylabel('Power (kW)')
legend('AEQ ≤1pu+Ptr≤ 1.5pu +θh≤120℃+θo≤105℃','Ptr≤ 1.5pu +θh≤120℃+θo≤105℃')
grid on

% Plot power shedding share in %
subplot(2,2,2), hold on
plot(reserve,(P_flex./(Load_agg_peak+reserve/100*Nominal_rating))*100, 'linewidth', 2)
grid on
ylabel('Power shedding share (%)')
xlabel('Reserve (%)')

legend('Power shedding/Power peak (%) AEQ','Power shedding/Power peak (%) PUL')

% Plot energy shedding 
subplot(2,2,3), hold on
plot(reserve,E_flex, 'linewidth', 2)
ylabel('Energy (kWh)') 
legend('AEQ ≤1pu+Ptr≤ 1.5pu +θh≤120℃+θo≤105℃','Ptr≤ 1.5pu +θh≤120℃+θo≤105℃')
xlabel('Reserve (%)')

set(gca, 'YScale', 'log')
grid on

% Plot energy shedding in % 
subplot(2,2,4), hold on
Energy_total=Energy_existing+Energy_connected_delta.*reserve/5;
plot(reserve,(E_flex./(Energy_total))*100, 'linewidth', 2)
legend('AEQ ≤1pu+Ptr≤ 1.5pu +θh≤120℃+θo≤105℃','Ptr≤ 1.5pu +θh≤120℃+θo≤105℃')
grid on
ylabel('The Energy shedding share(%)')
xlabel('Reserve (%)')

%-------------------------% θh≤120℃ θo≤105℃------------------------------
clear

% Load precalculated results for given formulation
load('result_temp_100.mat')

% Redefine the transformer rating
Nominal_rating=500;

% Step of power connected 
Power_connected_delta=25; % kW

% Calculate the annual energy connected 
Energy_connected_delta=Power_connected_delta*8760; % KWh per year for 25 KW connected (5%)

% Extract the reserve margins 
reserve=result.headroom(2:end)*100;

% Convert into kW
reserve_kW=Nominal_rating*reserve/100;

% Extract flexibility metrics: power in kW and energy in kWh
for i=1:length(reserve)
    P_flex(i,1)=result.flex_KW(i+1);
    E_flex(i,1)=result.flex_KWh(i+1);    
end

% Load aggregated load profile in W 
load('Aggregated_load_profile_100_houses.mat')

% Convert to kW and transform into hour resolution
Load_agg_kW=Convert2hours(Load_agg,60)/1000;

% Estimate the energy of load
Energy_existing=sum(Load_agg_kW); %Energy in kWh

% Find a peak load 
Load_agg_peak=max(Load_agg_kW);

% Plot the flexibility in kW as a function of reserve margin
subplot(2,2,1),hold on
plot(reserve,P_flex, 'linewidth', 2)
legend('AEQ ≤1pu+Ptr≤ 1.5pu +θh≤120℃+θo≤105℃','Ptr≤ 1.5pu +θh≤120℃+θo≤105℃','θh≤120℃+θo≤105℃')

grid on

% Plot power shedding share in %
subplot(2,2,2), hold on
plot(reserve,(P_flex./(Load_agg_peak+reserve/100*Nominal_rating))*100, 'linewidth', 2)
grid on
legend('AEQ ≤1pu+Ptr≤ 1.5pu +θh≤120℃+θo≤105℃','Ptr≤ 1.5pu +θh≤120℃+θo≤105℃','θh≤120℃+θo≤105℃')


% Plot energy shedding 
subplot(2,2,3), hold on
plot(reserve,E_flex, 'linewidth', 2)
legend('AEQ ≤1pu+Ptr≤ 1.5pu +θh≤120℃+θo≤105℃','Ptr≤ 1.5pu +θh≤120℃+θo≤105℃','θh≤120℃+θo≤105℃')


set(gca, 'YScale', 'log')
grid on

% Plot energy shedding in % 
subplot(2,2,4), hold on
Energy_total=Energy_existing+Energy_connected_delta.*reserve/5;
plot(reserve,(E_flex./(Energy_total))*100, 'linewidth', 2)
grid on
legend('AEQ ≤1pu+Ptr≤ 1.5pu +θh≤120℃+θo≤105℃','Ptr≤ 1.5pu +θh≤120℃+θo≤105℃','θh≤120℃+θo≤105℃')


%-----------------------AEQ≤1pu Ptr≤1.5pu θh≤98℃ θo≤105℃----------------
clear
% Load precalculated results for given formulation
load('result98_AEQ_100.mat')
Nominal_rating=500;

% Step of power connected 
Power_connected_delta=25;

% Calculate the annual energy connected 
Energy_connected_delta=Power_connected_delta*8760; % KWh per year for 25 KW connected (5%)

% Extract the reserve margins 
reserve=result.headroom(2:end)*100;

% Convert into kW
reserve_kW=Nominal_rating*reserve/100;

% Extract flexibility metrics: power in kW and energy in kWh
for i=1:length(reserve)
    P_flex(i,1)=result.flex_KW(i+1);
    E_flex(i,1)=result.flex_KWh(i+1);    
end

% Load aggregated load profile in W 
load('Aggregated_load_profile_100_houses.mat')

% Convert to kW and transform into hour resolution
Load_agg_kW=Convert2hours(Load_agg,60)/1000;

% Estimate the energy of load
Energy_existing=sum(Load_agg_kW); %Energy in kWh

% Find a peak load 
Load_agg_peak=max(Load_agg_kW);

% Plot the flexibility in kW as a function of reserve margin
subplot(2,2,1),hold on
plot(reserve,P_flex, 'linewidth', 2)
xlabel('Reserve (%)')
ylabel('Power (kW)')
legend('AEQ ≤1pu+Ptr≤ 1.5pu +θh≤120℃+θo≤105℃',...
    'Ptr≤ 1.5pu +θh≤120℃+θo≤105℃','θh≤120℃+θo≤105℃',...
    'Ptr≤1.5pu +θh≤98℃+θo≤105℃')

grid on

% Plot power shedding share in %
subplot(2,2,2), hold on
plot(reserve,(P_flex./(Load_agg_peak+reserve/100*Nominal_rating))*100,...
    'linewidth', 2)
grid on
ylabel('Power shedding share (%)')
xlabel('Reserve (%)')
legend('AEQ ≤1pu+Ptr≤ 1.5pu +θh≤120℃+θo≤105℃',...
    'Ptr≤ 1.5pu +θt≤120℃+θo≤105℃','θh≤120℃+θo≤105℃',...
    'Ptr≤1.5pu +θh≤98℃+θo≤105℃')


% Plot energy shedding 
subplot(2,2,3), hold on
plot(reserve,E_flex, 'linewidth', 2)
ylabel('Energy (kWh)') 
xlabel('Reserve (%)')
legend('AEQ ≤1pu+Ptr≤ 1.5pu +θh≤120℃+θo≤105℃',...
    'Ptr≤ 1.5pu +θh≤120℃+θo≤105℃','θh≤120℃+θo≤105℃',...
    'Ptr≤1.5pu +θh≤98℃+θo≤105℃')

set(gca, 'YScale', 'log')
grid on

% Plot energy shedding in % 
subplot(2,2,4), hold on
Energy_total=Energy_existing+Energy_connected_delta.*reserve/5;
plot(reserve,(E_flex./(Energy_total))*100, 'linewidth', 2)
grid on
ylabel('The Energy shedding share(%)')
xlabel('Reserve (%)')
set(gca, 'YScale', 'log')
legend('AEQ ≤1pu+Ptr≤ 1.5pu +θh≤120℃+θo≤105℃',...
    'Ptr≤ 1.5pu +θh≤120℃+θo≤105℃','θh≤120℃+θo≤105℃',...
    'Ptr≤1.5pu +θh≤98℃+θo≤105℃')
%% Plotting the Figure 15
% Figure name: Power shedding over the year: 75% reserve in 
% “energy-shedding” mode

clc;clear all % clear a command window and a workspace

% Load results for 75% reserve (+existing loading 86%)
load('result_AEQ_100.mat')

% Find the power profile with reserve 75%
Profile_reserve=PUL_init+result.headroom(end);

% Find the profile reconstructed with DR
Profile_reconstructed=PUL_optim;

% Power sheeding 
Power_shedding=(Profile_reserve-Profile_reconstructed)*500;

% Create figure
figure('DefaultAxesFontSize',14); 

% index of times when DR is applied 
Minutes_with_DR=find(Power_shedding>1);

% Find the DR activation (kW)
DR_values=Power_shedding(Minutes_with_DR);

% Plot histogram
h=histogram(DR_values);

xlabel('Power schedding, kW')
ylim([0 1800]) % limiting the y-axis

ytix = get(gca, 'YTick'); ytix(end+1)=1800; 
set(gca, 'YTick',ytix, 'YTickLabel',round(ytix/60)); % convert to hours
ylabel('DR activation, hours')

%% Plotting the Figure 16
% Figure name: Yearly temperature duration curves: 75% reserve in 
% “energy-shedding” mode leading to less than 1% energy curtailment

clc;clear all % clear a command window and a workspace

%-----------------------AEQ≤1pu Ptr≤1.5pu θh≤120℃ θo≤105℃---------------

% Load precalculated results 
load('result_AEQ_100.mat')

% Reconstruct the power profile with reserve 75% (without power shedding)
Profile_reserve=PUL_init+result.headroom(end);

% Reconstruct the power profile with reserve 75% (with power shedding)
Profile_reconstructed=PUL_optim;

%  Find the duration curve for Profile_reserve
Duration=sort(Profile_reserve, 'descend');

%  Find the duration curve for Profile_reconstructed
Duration_opt=sort(Profile_reconstructed, 'descend');

% Thermal study for the case without DR
[~,~,AEQ,theta_h,theta_0]=distribution_transformer(Profile_reserve,AMB);

% Thermal study for the case with DR
[~,~,AEQopt,theta_hopt,theta_0opt]=distribution_transformer(Profile_reconstructed,AMB);

% Find the duration curve of hot spot temperature for the case without DR
Duration=sort(theta_h, 'descend');

% Find the duration curve of hot spot temperature for the case with DR
Duration_opt=sort(theta_hopt, 'descend');

% Find the duration values (x-axis)
Duration_values=[1:525600]/525600*100;
Duration_values=Duration_values';

% Create figure
figure('DefaultAxesFontSize',14); 

% create smaller axes in top right, and plot on it
plot(Duration_values,Duration,'linewidth',3)
hold on
plot(Duration_values,Duration_opt,'linewidth',3)

ylabel('Hot spot temperature,°C')
xlabel ('Duration,%')
clear
%-----------------------Ptr≤1.5pu θh≤120℃ θo≤105℃-----------------------
% Load precalculated results 
load('result_PUL_100.mat')

% Reconstruct the power profile with reserve 75% (without power shedding)
Profile_reserve=PUL_init+result.headroom(end);

% Reconstruct the power profile with reserve 75% (with power shedding)
Profile_reconstructed=PUL_optim;

%  Find the duration curve for Profile_reconstructed
Duration_opt=sort(Profile_reconstructed, 'descend');

% Thermal study for the case without DR
[~,~,AEQ,theta_h,theta_0]=distribution_transformer(Profile_reserve,AMB);

% Thermal study for the case with DR
[~,~,AEQopt,theta_hopt,theta_0opt]=distribution_transformer(Profile_reconstructed,AMB);

% Find the duration curve of hot spot temperature for the case without DR
Duration=sort(theta_h, 'descend');

% Find the duration curve of hot spot temperature for the case with DR
Duration_opt=sort(theta_hopt, 'descend');

% Find the duration values (x-axis)
Duration_values=[1:525600]/525600*100;
Duration_values=Duration_values';

% create smaller axes in top right, and plot on it
hold on
plot(Duration_values,Duration_opt,'linewidth',3)

clear
%-----------------------θh≤120℃ θo≤105℃---------------------------------
% Load precalculated results 
load('result_temp_100.mat')

% Reconstruct the power profile with reserve 75% (without power shedding)
Profile_reserve=PUL_init+result.headroom(end);

% Reconstruct the power profile with reserve 75% (with power shedding)
Profile_reconstructed=PUL_optim;

% Find the duration curve of hot spot temperature for the case with DR
Duration_opt=sort(Profile_reconstructed, 'descend');

% Find the duration curve of initial load profile 
Duration_init=sort(PUL_init, 'descend');

% Thermal study for the case without DR
[~,~,AEQ,theta_h,theta_0]=distribution_transformer(Profile_reserve,AMB);

% Thermal study for the case with DR
[~,~,AEQopt,theta_hopt,theta_0opt]=distribution_transformer(Profile_reconstructed,AMB);

% Find the duration curve of hot spot temperature for the case without DR
Duration=sort(theta_h, 'descend');

% Find the duration curve of hot spot temperature for the case with DR
Duration_opt=sort(theta_hopt, 'descend');

% Find the duration values (x-axis)
Duration_values=[1:525600]/525600*100;
Duration_values=Duration_values';

% Thermal study for initial load profile 
[~,~,AEQ_init,theta_h_init,theta_0_init]=distribution_transformer(PUL_init,AMB);

% Find the duration curve of hot spot temperature for the initial load profile
Duration_init=sort(theta_h_init, 'descend');

% create smaller axes in top right, and plot on it
hold on
plot(Duration_values,Duration_opt,'linewidth',3)

plot(Duration_values,Duration_init,'linewidth',3)
Limit98=linspace(98,98, length(Duration_values))';
Limit120=linspace(120,120, length(Duration_values))';
plot(Duration_values,Limit98,'g')
plot(Duration_values,Limit120,'r')
legend('Not optimized: Existing θh+reserve','Optimized with AEQ ≤1pu... Ptr≤ 1.5pu  θh≤120℃ θo≤105℃ ','Optimized with Ptr≤ 1.5pu θh≤120℃ θo≤105℃ ','Optimized with θh≤120℃+θo≤105℃ ', 'Existing θh','Design temperature','Temperature limit' )

%% Plotting the Figure 17
% Figure name: Results: DR modeled with or without payback effect/energy 
% conservation: (a) DR power in kW; (b) DR energy in kWh

clc;clear all % clear a command window and a workspace

% Load precalculated results for 'energy shedding' mode
load('result_AEQ_100.mat')

% Extract studied reserve margins 
reserve=result.headroom(2:end)*100;

% Convert reserve margins into kW
reserve_kW=Nominal_rating*reserve/100;

% Extract flexibility metrics: power in kW and energy in kWh
for i=1:length(reserve)
    P_flex(i,1)=result.flex_KW(i+1);
    E_flex(i,1)=result.flex_KWh(i+1);    
end

% Load aggregated load profile in W 
load('Aggregated_load_profile_100_houses.mat')  

% Convert Load_agg to kW and into hour resolution 
Load_agg_kW=Convert2hours(Load_agg,60)/1000;

% Calculate the energy of load as a sum (but better as integral)
Energy_existing=sum(Load_agg_kW); %Energy in kWh

% Find the peak load
Load_agg_peak=max(Load_agg_kW);

% Create the figure 
figure()

% Plot power shedding for 'energy shedding' mode 
subplot(1,2,1), plot(reserve,P_flex, 'linewidth', 2)
xlabel('Reserve (%)')
ylabel('Power (kW)')
title('Power shedding, kW')
hold on

grid on

% Plot energy shedding for 'energy shedding' mode 
subplot(1,2,2), plot(reserve,E_flex, 'linewidth', 2)
ylabel('Energy (kWh)') 
title('Energy shedding(kWh)')
xlabel('Reserve (%)')
set(gca, 'YScale', 'log')

grid on

hold on

% Load precalculated results for 'energy shifting' mode
load('result_AEQ_50_60.mat')

% Extract the studied reserve margins 
reserve=result.headroom(2:end)*100;

% Create new variables 
P_flex=[];
E_flex=[];

% Extract flexibility metrics: power in kW and energy in kWh
for i=1:length(reserve)
    P_flex(i,1)=result.flex_KW(i+1);
    E_flex(i,1)=result.flex_KWh(i+1);    
end

% Plot power shedding for 'energy shifting' mode 
subplot(1,2,1)
plot(reserve,P_flex, 'linewidth', 2)

% Plot energy shedding for 'energy shifting' mode 
subplot(1,2,2)
plot(reserve,E_flex, 'linewidth', 2)
legend('Energy shedding', 'Energy shifting')

%% Plotting the Figure 18
% Figure name: Computation times for the nonlinear formulation (the blue 
% line) and for the suggested PWL formulation (the green line). 
% Both formulations were tested for the “energy-shedding” case (SOCDR0 =
% 100% and SOCDRt=t = 0%)

clc;clear all % clear a command window and a workspace

% Horizon of optimization problem in minutes
length_in_min=[1 2 3 4 5 6 7 10 14  20 25 45]*1440;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATTENTION! The code commented below is not approved yet. We will check it 
% before april 10 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve a nonlinear problem
% for i=1:length(length_in_min)
%     try
%         LOADprofile=linspace(600,600,length_in_min(i))';
%         AMB=linspace(20,20,length_in_min(i))';
%         x0.Eflex_max=1000;
%         x0.Pflex_max=553;
%         x0.Pflex=ones(length(LOADprofile)/60,1)*553;
%         x0.Ptr=ones(length(LOADprofile)/60,1)*500;
%         [PUL_optim,flex_KW,flex_KWh,time]=fmincon_time_test(LOADprofile,AMB,x0);
%         Time_solve_fmincon(i,1)=time;
%     catch
%         error('Error in for cycle')
%     end
% end
% 
% % Solve a linearized problem
% for i=1:length(length_in_min)
%     try
%         LOADprofile=linspace(600,600,length_in_min(i))';
%         AMB=linspace(20,20,length_in_min(i))';
%         [PUL_optim,Pflex_max,Eflex_max,theta_h_optim,theta_0_optim,time]=Linearilized_optim(LOADprofile,AMB);
%         Time_solve_PWL(i,1)=time;
%     catch
%         error('Error in for cycle')
%     end
% end

load('time_test_result_5.mat', 'Time_solve')
Time_solve_fmincon=Time_solve(1:12);

load('time_test_result_linear.mat')
Time_solve_PWL=Time_solve;

% Create figure
figure('InvertHardcopy','off','Color',[1 1 1]);

% Create axes
axes1 = axes('Position',...
    [0.125900217884299 0.117271760784138 0.753003542472616 0.849672489082969]);
hold(axes1,'on');

% Create plot
plot(length_in_min/1440,Time_solve_fmincon,'DisplayName','Nonlinear formulation: fmincon','MarkerSize',12,...
    'Marker','o',...
    'LineWidth',3);

% Create plot
plot(length_in_min/1440,Time_solve_PWL,'DisplayName','Suggested linear formulation: PWL',...
    'MarkerSize',12,...
    'Marker','^',...
    'LineWidth',3,...
    'Color',[0.466666668653488 0.674509823322296 0.18823529779911]);

% Create ylabel
ylabel('Solving time, min');

% Create xlabel
xlabel('Horizon of optimization problem, days');

box(axes1,'on');
hold(axes1,'off');

% Set the remaining axes properties
set(axes1,'FontSize',20);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.147014922290036 0.730501802396799 0.435499412955561 0.105193664322437],...
    'FontSize',20,...
    'EdgeColor',[1 1 1]);

% Elapsed time of script 
Elapsed_time=toc;