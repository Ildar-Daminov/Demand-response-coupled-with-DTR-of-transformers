# Demand Response Coupled with Dynamic Thermal Rating for Increased Transformer Reserve and Lifetime
<img align="left" alt="Coding" width="275" src="https://www.i3upgrade.eu/files/2021/07/logo-journal-energies.png">

This repository shares the MATLAB code and data for the research article :\
I. Daminov, R. Rigo-Mariani, R. Caire, A. Prokhorov, M-C Alvarez-Herault, [“Demand Response Coupled with Dynamic Thermal Rating for Increased Transformer Reserve and Lifetime”](https://doi.org/10.3390/en14051378) in Energies (IF: 2.702, Q2), 2021

## Article's abstract
(1) Background: This paper proposes a strategy coupling Demand Response Program with Dynamic Thermal Rating to ensure a transformer reserve for the load connection. This solution is an alternative to expensive grid reinforcements. (2) Methods: The proposed methodology firstly considers the N-1 mode under strict assumptions on load and ambient temperature and then identifies critical periods of the year when transformer constraints are violated. For each critical period, the integrated management/sizing problem is solved in YALMIP to find the minimal Demand Response needed to ensure a load connection. However, due to the nonlinear thermal model of transformers, the optimization problem becomes intractable at long periods. To overcome this problem, a validated piece-wise linearization is applied here. (3) Results: It is possible to increase reserve margins significantly compared to conventional approaches. These high reserve margins could be achieved for relatively small Demand Response volumes. For instance, a reserve margin of 75% (of transformer nominal rating) can be ensured if only 1% of the annual energy is curtailed. Moreover, the maximal amplitude of Demand Response (in kW) should be activated only 2–3 h during a year. (4) Conclusions: Improvements for combining Demand Response with Dynamic Thermal Rating are suggested. Results could be used to develop consumer connection agreements with variable network access. 

## How to run a code 
There are two ways how you may run this code:
  
I. Launching all calculations at once. This will reproduce all figures in the article but it would take 5 minutes:
1. Copy this repository to your computer
2. Extract YALMIP.zip to any place you wish. 
3. Add YALMIP to Set Path. 
   Go to Matlab's Home tab -> Click on "Set Path" -> Click on "Add with Subfolders" and Choose the folder where YALMIP is located -> Click on Select the folder
4. Open the script main.m
5. Launch the script "main.m" by clicking on the button "Run" (usually located at the top of MATLAB bar).\
As alternative, you may type ```main``` 
in Command Window to launch the entire script. 


II. Launching the specific section of the code to reproduce the particular figure: 
1. Copy this repository to your computer 
2. Extract YALMIP.zip to any place you wish. 
3. Add YALMIP to Set Path. 
   Go to Matlab's Home tab -> Click on "Set Path" -> Click on "Add with Subfolders" and Choose the folder where YALMIP is located -> Click on Select the folder
4. Open the script main.m 
5. Find the section (Plotting the Figure XX) corresponding to the Figure you would like to reproduce. 
6. Put the cursor at any place of this section and click on the button "Run Section" (usually located at the top of MATLAB bar)

Attention! Some code uses [fcn2optimexpr](https://fr.mathworks.com/help/optim/ug/fcn2optimexpr.html), which becomes available since the version MATLAB 2019a! For previous MATLAB version [fcn2optimexpr](https://fr.mathworks.com/help/optim/ug/fcn2optimexpr.html), as far as we know, does not work


## Files description
Main script:
* main.m - the principal script which launches all calculations
  
Additional functions: 
* computeXBKPbest.m - this function computes the best breakpoints for simple PWL piecewise linearization
* Convert2hours.m - this function converts from 1-minute resolution to hours resolution
* Convert2minute.m - this function converts from hours resolution to 1-minute resolution
* distribution_transformer.m - a thermal model of distribution transformer (up to 2.5 MVA) per the loading guide IEC 60076-7
* fmincon_time_test.m - this function formulates and then solves a nonlinear optimization problem (using fmincon)
* minutes_integer2day_index.m - this function converts the interval duration in minutes to corresponding day index 
* minutes2intervals.m  - this function transforms minutes values (where transformer constraints are violated) to intevals 
* profiles2minutes.m - this function analyzes the power and temperature profiles and shows minutes and intervals where constraints are violated
* validation_run.m - this function performs the validation runs (using MILP optimization problem) 

More details are given inside of each functions and script "main.m"

Initial data:
* Aggregated_load_profile_100_houses.mat - the load profile in W of 100 houses 
* all_intervals.mat - precalulated results 
* Ambient_temperature_Grenoble.mat - annual ambient temperature in Grenoble, France ([weather data](https://www.meteoblue.com/en/historyplus) provided by [meteoblue](https://www.meteoblue.com/)) 
* fig_nRMSE.mat - precalulated results
* initial_data.mat - daily profile of transfomer loading and ambient temperature 
* LOADprofile.mat - transformer load profile 
* result_AEQ_50_60.mat - precalulated results for energy shifting mode**  
* result_AEQ_100.mat - precalulated results for energy shedding mode** (considering all constraints : ageing , current and temperatures 120°C)
* result_PUL_100.mat - precalulated results for energy shedding mode (considering current and temperature constraints)
* result_temp_100.mat - precalulated results for energy shedding mode (considering only temperature  constraints)
* result98_AEQ_100.mat - precalulated results for energy shedding mode (considering ageing and design temperature 98 °C constraints)
* time_test_result_5.mat - precalulated results of time tests (of solving  the nonlinear problem via fmincon)
* time_test_result_linear.mat - precalulated results of time tests (of solving the linearized problem via linprog)

** energy shifting mode: SOC at the beginnnig =50 % & SOC at the end 50%; energy shedding mode: SOC at the beginnnig =100 % & SOC at the end 0%


## How to cite this article 
Ildar Daminov, Rémy Rigo-Mariani, Raphael Caire, Anton Prokhorov, Marie-Cécile Alvarez-Herault, "Demand Response coupled with Dynamic Thermal Rating for increased transformer reserve and lifetime." Energies 14.5 (2021): 1378. https://doi.org/10.3390/en14051378

## More about DTR of power transformers 
<img align="left" alt="Coding" width="250" src="https://sun9-19.userapi.com/impg/3dcwjraHJPNgrxtWv7gEjZTQkvv5T0BttTDwVg/e9rt2Xs8Y5A.jpg?size=763x1080&quality=95&sign=7c57483971f31f7009fbcdce5aafd97e&type=album">This paper is a part of PhD thesis "Dynamic Thermal Rating of Power Transformers: Modelling, Concepts, and Application case". The full text of PhD thesis is available on [Researchgate](https://www.researchgate.net/publication/363383515_Dynamic_Thermal_Rating_of_Power_Transformers_Modelling_Concepts_and_Application_case) or [HAL theses](https://tel.archives-ouvertes.fr/tel-03772184). Other GitHub repositories on DTR of power transformers:
* Article: Assessment of dynamic transformer rating, considering current and temperature limitations. [GitHub repository](https://github.com/Ildar-Daminov/Assessment_Dynamic_Thermal_Rating_of_Transformers)
* Article: Energy limit of oil-immersed transformers: A concept and its application in different climate conditions. [GitHub repository](https://github.com/Ildar-Daminov/Energy-limit-of-power-transformer)
* Conference paper: Optimal ageing limit of oil-immersed transformers in flexible power systems [GitHub repository](https://github.com/Ildar-Daminov/MATLAB-code-for-CIRED-paper)
* Conference paper: Application of dynamic transformer ratings to increase the reserve of primary substations for new load interconnection. [GitHub repository](https://github.com/Ildar-Daminov/Reserve-capacity-of-transformer-for-load-connection)
* Conference paper: Receding horizon algorithm for dynamic transformer rating and its application for real-time economic dispatch. [GitHub repository](https://github.com/Ildar-Daminov/Receding-horizon-algorithm-for-dynamic-transformer-rating)
