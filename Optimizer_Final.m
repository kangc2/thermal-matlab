clc; close all; clear
%% Optimizer Starts Here: Step One: Define parameters we want to optimize
% PCM Material Parameters
% rhoL, Lh, f 
% Array of the parameters we want to optimize: x
x_PCM = [];
% Set initial Values - focused on just PCM parameters
% rhoL
x_PCM(1) = 760;
% Lh
x_PCM(2) = 260;
% Tm
x_PCM(3) = 15;

% Geometric Parameters
% L1, b1, t
x_geo = [];
% L1
x_geo(1) = 1;
%b1
x_geo(2) = .15;
% t
x_geo(3) = 0.01;

%% Step 2: Run fmincon PCM Material
% options = optimoptions('fmincon', 'MaxFunctionEvaluations', 50000, 'MaxIterations', 3000)
[xopt_PCM, fval_PCM, exitflag_PCM, output_PCM] = fmincon(@objective_PCM, x_PCM, [], [], [], [], [], [], @constraint_PCM, [])

[c0_PCM, ceq0_PCM] = constraint_PCM(x_PCM)
[c1_PCM, ceq1_PCM] = constraint_PCM(xopt_PCM)

% Results: Eff = 6.9675%
% rhoL = 800, Lh = 210, Tm = 18

%% Step 2: Run fmincon (GEOMETRIC)

% options = optimoptions('fmincon', 'MaxFunctionEvaluations', 50000, 'MaxIterations', 3000)
[xopt_geo, fval_geo, exitflag_geo, output_geo] = fmincon(@objective_geo, x_geo, [], [], [], [], [], [], @constraint_geo, [])

[c0_geo, ceq0_geo] = constraint_geo(x_geo)
[c1_geo, ceq1_geo] = constraint_geo(xopt_geo)

% Results: Eff = 14.9440%
% L1 = 1.3, b1 = .2, t = .01
%% ENGINE STRUCTURE - CAN BE USE TO TEST OPTIMAL PARAMETERS FOUND
% Thermal Engine Struct

%Environmental Parameters
T = 29; % Working Temperature [degC]
engine.T = T;
engine.Thigh = T; % High Temperature
engine.Tlow = 5; % Low Temperature [degC]
engine.Po = 0.101; % Initial Pressure/ambient pressure = atmospheric pressure [MPa]

% -> Geometry & Material Properties of Cylinder
engine.L1 = 1.3; % Length of the cylinder [m]
%engine.a1 = 7.7e-02; % Internal diameter of the cylinder [m]
engine.b1 = 0.15; % External diameter of the cylinder [m]
engine.t = 0.01; %wall thickness [m]
engine.a1 = engine.b1 - 2*engine.t; % Internal diameter of the cylinder [m]
engine.Ey = 105e03; % Young's Modulus [MPa]
engine.v = 0.33; % Poisson's ratio

% -> Material Properties of PCM
engine.delta_rho = 70;
engine.rhoL = 773; % Density of PCM - Liquid phase [kg/m3]
engine.rhoS = engine.rhoL + engine.delta_rho; %TTTTHIIISSSSSSSSS
engine.f = 0.6557; %Volume fraction of PCM
engine.mPCM = (engine.f*pi*engine.L1*( (engine.a1 / 2)^2 ))/(1/engine.rhoS);
engine.Tm = 18.2; % Melting Temperature of PCM [degC]

% -> Heat Transfer Properties
engine.csd = 1.64; % Specific heat in solid state [kJ/kg K]
engine.cld = 2.09; % Specific heat in liquid state [kJ/kg K]
engine.Lh = 236; % Latent heat of fusion [kJ/kg]

% -> Material Properties of Hydraulic Fluid
engine.voH = 1/1000; % Specific volume of hydraulic fluid [m3/kg]
engine.v1H = engine.voH; % Specific volume of hydraulic fluid at State 1 [m3/kg]

% -> B and C values of PCM (P) and HF (H)
engine.CP = 2.66e-04;
engine.BP = 102.12;
engine.BH = ( 2672.9 + (15.97*T) - (0.166*(T^2) ) )*1e-01;
engine.CH = 0.3150*engine.voH;

% -> Accumulator Parameters
engine.V1N = 2e-03; % Initial Volume of N2 gas
engine.ar = 6.573 / 100; % Volume fraction of residual air
%% Add in Optimal Parameters
engine1 = engine;
engine1.rhoL = xopt_PCM(1);
engine1.Lh = xopt_PCM(2);
engine1.Tm = xopt_PCM(3);
engine1.rhoS = engine1.rhoL + engine1.delta_rho; %TTTTHIIISSSSSSSSS
engine1.mPCM = (engine1.f*pi*engine1.L1*( (engine1.a1 / 2)^2 ))/(1/engine1.rhoS);

% Find Efficiency with engine structure
Eff_PCM = findEfficiency(engine1)

%%
engine2 = engine;
engine2.L1 = xopt_geo(1); % Length of the cylinder [m]
engine2.b1 = xopt_geo(2); % External diameter of the cylinder [m]
engine2.t = xopt_geo(3); %wall thickness [m]
engine2.a1 = engine2.b1 - 2*engine2.t; % Internal diameter of the cylinder [m]
engine2.mPCM = (engine2.f*pi*engine2.L1*( (engine2.a1 / 2)^2 ))/(1/engine2.rhoS);

Eff_geo = findEfficiency(engine2)

%% Optimizer Functions

% Define objective function for optimization

% This function is for PCM Material Parameters
function obj = objective_PCM(x)
% add engine parameters    
    %Environmental Parameters
    T = 29; % Working Temperature [degC]
    engine.T = T;
    engine.Thigh = T; % High Temperature
    engine.Tlow = 5; % Low Temperature [degC]
    engine.Po = 0.101; % Initial Pressure/ambient pressure = atmospheric pressure [MPa]

    % -> Geometry & Material Properties of Cylinder
    engine.L1 = 1.3; % Length of the cylinder [m]
    engine.b1 = 0.15; % External diameter of the cylinder [m]
    engine.t = 0.01; %wall thickness [m]
    engine.a1 = engine.b1 - 2*engine.t; % Internal diameter of the cylinder [m]
    engine.Ey = 105e03; % Young's Modulus [MPa]
    engine.v = 0.33; % Poisson's ratio
    
    % -> Material Properties of PCM
    engine.delta_rho = 70;
    engine.f = 0.6557; %Volume fraction of PCM
    
    % -> Heat Transfer Properties
    engine.csd = 1.64; % Specific heat in solid state [kJ/kg K]
    engine.cld = 2.09; % Specific heat in liquid state [kJ/kg K]
    engine.Lh = 236; % Latent heat of fusion [kJ/kg]
    
    % -> Material Properties of Hydraulic Fluid
    engine.voH = 1/1000; % Specific volume of hydraulic fluid [m3/kg]
    engine.v1H = engine.voH; % Specific volume of hydraulic fluid at State 1 [m3/kg]
    
    % -> B and C values of PCM (P) and HF (H)
    engine.CP = 2.66e-04;
    engine.BP = 102.12;
    engine.BH = ( 2672.9 + (15.97*T) - (0.166*(T^2) ) )*1e-01;
    engine.CH = 0.3150*engine.voH;
    
    % -> Accumulator Parameters
    engine.V1N = 2e-03; % Initial Volume of N2 gas
    engine.ar = 6.573 / 100; % Volume fraction of residual air

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add in the x parameters into engine

    %%% ADD ARRAY %%%
    engine.rhoL = x(1);
    engine.Lh = x(2);
    engine.Tm = x(3);
    % Need to update these values with new 'x' values
    engine.rhoS = engine.rhoL + engine.delta_rho; %TTTTHIIISSSSSSSSS
    engine.mPCM = (engine.f*pi*engine.L1*( (engine.a1 / 2)^2 ))/(1/engine.rhoS);

    % Find objective: max efficiency
    obj = -findEfficiency(engine);
end

% Define constraint for optimization
% c(x) <= 0
%ceq(x) = 0
function [c, ceq] = constraint_PCM(x)
    c(1) = x(1) - 800; % rhoL < 800
    c(2) = 720 - x(1); % rhoL > 720
    c(3) = x(2) - 270; % Lh < 270
    c(4) = 210 - x(2); % Lh < 210
    
    ceq = x(3) - 18; % Tm
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function is for Geometric Parameters
function obj = objective_geo(x)
% add engine parameters    
    %Environmental Parameters
    T = 29; % Working Temperature [degC]
    engine.T = T;
    engine.Thigh = T; % High Temperature
    engine.Tlow = 5; % Low Temperature [degC]
    engine.Po = 0.101; % Initial Pressure/ambient pressure = atmospheric pressure [MPa]

    % -> Geometry & Material Properties of Cylinder
    engine.Ey = 105e03; % Young's Modulus [MPa]
    engine.v = 0.33; % Poisson's ratio
    
    % -> Material Properties of PCM
    engine.delta_rho = 70;
    engine.f = 0.6557; %Volume fraction of PCM
    engine.rhoL = 773;
    engine.Lh = 236;
    engine.Tm = 18.2;
    engine.rhoS = engine.rhoL + engine.delta_rho; %TTTTHIIISSSSSSSSS
    
    % -> Heat Transfer Properties
    engine.csd = 1.64; % Specific heat in solid state [kJ/kg K]
    engine.cld = 2.09; % Specific heat in liquid state [kJ/kg K]
    engine.Lh = 236; % Latent heat of fusion [kJ/kg]
    
    % -> Material Properties of Hydraulic Fluid
    engine.voH = 1/1000; % Specific volume of hydraulic fluid [m3/kg]
    engine.v1H = engine.voH; % Specific volume of hydraulic fluid at State 1 [m3/kg]
    
    % -> B and C values of PCM (P) and HF (H)
    engine.CP = 2.66e-04;
    engine.BP = 102.12;
    engine.BH = ( 2672.9 + (15.97*T) - (0.166*(T^2) ) )*1e-01;
    engine.CH = 0.3150*engine.voH;
    
    % -> Accumulator Parameters
    engine.V1N = 2e-03; % Initial Volume of N2 gas
    engine.ar = 6.573 / 100; % Volume fraction of residual air

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add in the x parameters into engine

    %%% ADD ARRAY %%%
    engine.L1 = x(1); % Length of the cylinder [m]
    engine.b1 = x(2); % External diameter of the cylinder [m]
    engine.t = x(3); %wall thickness [m]

    % Need to update due to new 'x'
    engine.a1 = engine.b1 - 2*engine.t; % Internal diameter of the cylinder [m]
    engine.mPCM = (engine.f*pi*engine.L1*( (engine.a1 / 2)^2 ))/(1/engine.rhoS);

    % Find objective: max efficiency
    obj = -findEfficiency(engine);
end

% Define constraint for optimization
% c(x) <= 0
%ceq(x) = 0
function [c, ceq] = constraint_geo(x)
    c(1) = x(1) - 1.3; % L1 < 1.3
    c(2) = .3 - x(1); % L1 > .3
    c(3) = x(2) - .2; % b1 < 0.1
    c(4) = .1 - x(2); % b1 < .2
    
    ceq = x(3) - 0.01; % t
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Efficiency Equation with engine structure
% 1. Finding Specific Volume of PCM [vP]
function Eff = findEfficiency(engine)
    % a. Find specific volume of PCM under ambient pressure Po [voP]
    voP = ( (1.0307e03 - ( 1.2596*(engine.T + 273.15) )  + ...
        (1.8186e-3* (engine.T + 273.15)^2) -(1.9555e-6* (engine.T + 273.15)^3) ) )^-1;
    % b. Find the volume of PCM (liquid state) [v1P]
    v1P = 1/engine.rhoS; % Specific volume of PCM in liquid state
    
    
    % 3. Find inner volume of cylinder, Volume air and volume fraction of PCM, volume/mass of HF 
%     V = pi*engine.L1*( (engine.a1 / 2)^2 ); % Inner volume of Cylinder
    V = pi*engine.L1*( (engine.b1 / 2)^2 - (engine.a1 / 2)^2 ); % Inner volume of Cylinder

    V1A = engine.ar*V; % Volume of residual air
    %     V1H = ( V*(1 - engine.f) ) - V1A; % volume of HF at state 1
    mH = ( (V*(1 - engine.f) ) - V1A) / engine.v1H; % mass of HF
    
    
    % 5. Finding Max Pressure [P2]
    
    F = @(P)((pi / 4)*(engine.L1*( ( (2*engine.a1) + ...
            ( ( (P - engine.Po)*engine.a1*(1 - engine.v^2) ) / engine.Ey)*( ( (engine.b1^2 + engine.a1^2) / (engine.b1^2 - engine.a1^2) ) +...
        (engine.v / (1 - engine.v) ) ))*...
            ( ( (P - engine.Po)*engine.a1*(1 - engine.v^2) ) / engine.Ey)*( ( (engine.b1^2 + engine.a1^2) / (engine.b1^2 - engine.a1^2) ) + ...
        (engine.v / (1 - engine.v) ) )) ))... % delta_V1
    - ...
    (( engine.mPCM*((1.3e-03 - (2.66e-04*log10( 1 + ( (P - engine.Po) / 102.12) ) )) - v1P) ) + ...
    ( mH*((engine.voH - (engine.CH*log10(1 + ( (P - engine.Po) / engine.BH) ) )) - engine.voH) ) + (((V1A*engine.Po) / P) - V1A)); %delta_V2
    
    % Options: sets tolerance of function close to 0 (1-e14) and displays the
    % iteration, this could help with the optimization
    options = optimoptions('fsolve','TolFun',1e-14);
    
    %solves for P, same answer as the engine.P2 with the current finding
    %Pressure function
    P2 = fsolve(F,5,options);
    
    % 4. Finding change of Inner Diameter of Tube Under Internal Pressure
    delta_a1 = ( ( (P2 - engine.Po)*engine.a1*(1 - engine.v^2) ) / engine.Ey)*( ( (engine.b1^2 + engine.a1^2) / ...
        (engine.b1^2 - engine.a1^2) ) + (engine.v / (1 - engine.v) ) );
    delta_V1 = (pi / 4)*(engine.L1*( ( (2*engine.a1) + delta_a1)*delta_a1) );
    
    
    % 6. Finding Efficiency [Eff] 
    % a. Find Pre-charged pressure in acculimlator [Pa]
    Pa = (P2 / engine.V1N)*(delta_V1 + engine.V1N - V1A*( (engine.Po / P2) - 1) ...
       - (V*engine.f / v1P)*(voP - engine.CP*log10(1 + ((P2 - engine.Po) / engine.BP) ) - v1P) + ...
       ((V*(1 - engine.f) - V1A) / engine.v1H)*engine.CH*log10(1 + ((P2 - engine.Po) / engine.BH) ));  
    
    % b. find the total energy stored [Est]
    Est = -Pa*1e6*engine.V1N*log(1 - (engine.mPCM / engine.V1N)*((1 / engine.rhoL) - (1 / engine.rhoS)) );
    
    
    % c. Find total thermal energy absorbed by PCM in 1 melting/freezing cycle
    % [Qin]
    Qin = engine.mPCM*engine.csd*(engine.Tm - engine.Tlow) + engine.mPCM*engine.Lh + engine.mPCM*engine.cld*(engine.Thigh - engine.Tm);
    % d. find the theorectical Efficiency % [Eff] 
    
    Eff = Est / (Qin*1e3) * 100;
end