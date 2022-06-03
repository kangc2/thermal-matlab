clc; close all; clear
%% Optimizer Starts Here: Step One: Define parameters we want to optimize
% rhoL, Lh, f 
% Array of the parameters we want to optimize: x
x_PCM = [];
% Set initial Values - focused on just PCM parameters
% rhoL
x_PCM(1) = 720;
% Lh
x_PCM(2) = 280;
% Tm
x_PCM(3) = 18; 

%% Step 2: Run in fmincon
% options = optimoptions('fmincon', 'MaxFunctionEvaluations', 50000, 'MaxIterations', 3000)
[xopt, fval, exitflag, output] = fmincon(@objective_PCM, x_PCM, [], [], [], [], [], [], @constraint_PCM, [])

[c, ceq] = constraint_PCM(x_PCM)
[c1, ceq1] = constraint_PCM(xopt)

% Results: Eff = 
% High as possible: 
% Inbetween: 
% Low as Possible: 

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
engine.b1 = 0.15; % External diameter of the cylinder [m]
engine.t = 0.01; %wall thickness [m]
engine.a1 = engine.b1 - 2*engine.t; % Internal diameter of the cylinder [m]
engine.Ey = 105e03; % Young's Modulus [MPa]
engine.v = 0.33; % Poisson's ratio

% -> Material Properties of PCM
engine.delta_rho = 70;
engine.rhoL = 800; % Density of PCM - Liquid phase [kg/m3]
engine.rhoS = engine.rhoL + engine.delta_rho;
engine.f = 0.6557; %Volume fraction of PCM
engine.mPCM = (engine.f*pi*engine.L1*( (engine.a1 / 2)^2 ))/(1/engine.rhoS);
engine.Tm = 18; % Melting Temperature of PCM [degC]

% -> Heat Transfer Properties
engine.csd = 1.64; % Specific heat in solid state [kJ/kg K]
engine.cld = 2.09; % Specific heat in liquid state [kJ/kg K]
engine.Lh = 210; % Latent heat of fusion [kJ/kg]

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

% Find Efficiency with engine structure
Eff = findEfficiency(engine)

%% Optimizer Functions

% Define objective function for optimization
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

    %%% ADD ARRAY x %%%
    engine.rhoL = x(1);
    engine.Lh = x(2);
    engine.Tm = x(3);
    % Need to update these values with new 'x' values
    engine.rhoS = engine.rhoL + engine.delta_rho; %TTTTHIIISSSSSSSSS
    engine.mPCM = (engine.f*pi*engine.L1*( (engine.a1 / 2)^2 ))/(1/engine.rhoS);
    engine
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


% Find Efficiency with engine structure
% 1. Finding Specific Volume of PCM [vP]
function Eff = findEfficiency(engine)
    % a. Find specific volume of PCM under ambient pressure Po [voP]
    voP = ( (1.0307e03 - ( 1.2596*(engine.T + 273.15) )  + ...
        (1.8186e-3* (engine.T + 273.15)^2) -(1.9555e-6* (engine.T + 273.15)^3) ) )^-1;
    % b. Find the volume of PCM (liquid state) [v1P]
    v1P = 1/engine.rhoS; % Specific volume of PCM in liquid state
   
    % 3. Find inner volume of cylinder, Volume air and volume fraction of PCM, volume/mass of HF 
    V = pi*engine.L1*( (engine.b1 / 2)^2 - (engine.a1 / 2)^2 ); % Inner volume of Cylinder
    V1A = engine.ar*V; % Volume of residual air
    mH = ( (V*(1 - engine.f) ) - V1A) / engine.v1H; % mass of HF
        
    % 5. Finding Max Pressure [P2]
    F = @(P)...         % delta_V1
        ((pi / 4)*(engine.L1*( ( (2*engine.a1) + ...
            ( ( (P - engine.Po)*engine.a1*(1 - engine.v^2) ) / engine.Ey)*( ( (engine.b1^2 + engine.a1^2) / (engine.b1^2 - engine.a1^2) ) +...
        (engine.v / (1 - engine.v) ) ))*...
            ( ( (P - engine.Po)*engine.a1*(1 - engine.v^2) ) / engine.Ey)*( ( (engine.b1^2 + engine.a1^2) / (engine.b1^2 - engine.a1^2) ) + ...
        (engine.v / (1 - engine.v) ) )) ))...x
        - ... %delta_V2
        (( engine.mPCM*((1.3e-03 - (2.66e-04*log10( 1 + ( (P - engine.Po) / 102.12) ) )) - v1P) ) + ...
        ( mH*((engine.voH - (engine.CH*log10(1 + ( (P - engine.Po) / engine.BH) ) )) - engine.voH) ) + (((V1A*engine.Po) / P) - V1A));
        
    % Options: sets tolerance of function close to 0 (1-e14) and displays the
    % iteration, this could help with the optimization
    options = optimoptions('fsolve','TolFun',1e-14);
    %solves for P2, P = P2 in Function F(P)
    P2 = fsolve(F,5,options);
    
    % 4. Finding change of Inner Diameter of Tube Under Internal Pressure
    delta_a1 = ( ( (P2 - engine.Po)*engine.a1*(1 - engine.v^2) ) / engine.Ey)*( ( (engine.b1^2 + engine.a1^2) / ...
        (engine.b1^2 - engine.a1^2) ) + (engine.v / (1 - engine.v) ) );
    delta_V1 = (pi / 4)*(engine.L1*( ( (2*engine.a1) + delta_a1)*delta_a1) );
    
    % 6a. Find Pre-charged pressure in acculimlator [Pa]
    Pa = (P2 / engine.V1N)*(delta_V1 + engine.V1N - V1A*( (engine.Po / P2) - 1) ...
       - (V*engine.f / v1P)*(voP - engine.CP*log10(1 + ((P2 - engine.Po) / engine.BP) ) - v1P) + ...
       ((V*(1 - engine.f) - V1A) / engine.v1H)*engine.CH*log10(1 + ((P2 - engine.Po) / engine.BH) ));  
    
    % 6b. find the total energy stored [Est]
    Est = -Pa*1e6*engine.V1N*log(1 - (engine.mPCM / engine.V1N)*((1 / engine.rhoL) - (1 / engine.rhoS)) );
    
    
    % 6c. Find total thermal energy absorbed by PCM in 1 melting/freezing cycle
    % [Qin]
    Qin = engine.mPCM*engine.csd*(engine.Tm - engine.Tlow) + engine.mPCM*engine.Lh + engine.mPCM*engine.cld*(engine.Thigh - engine.Tm);
    % 6d. find the theorectical Efficiency % [Eff] 
    Eff = Est / (Qin*1e3) * 100;
end