clc; close all; clear

%% Define 'engine' structure
%Environmental Parameters
% Variables of the environment, not specific to the engine
T = 29; % Working Temperature [degC]
engine.T = T;
engine.Thigh = T; % High Temperature
engine.Tlow = 5; % Low Temperature [degC]
engine.Po = 0.101; % Initial Pressure/ambient pressure = atmospheric pressure [MPa]

% Thermal Engine from our Thermal_Engine code
% The baseline of our model
% Note: make sure the field name is right (voH vs. VoH)
engine.name = 'Original';
engine.working_fluid = 'Water';
engine.hull_material = 'TA2M (Titanium alloy)';

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
% engine.rhoS = 864; % Density of PCM - Solid phase [kg/m3] 
engine.rhoS = engine.rhoL + engine.delta_rho; %TTTTHIIISSSSSSSSS
engine.f = 0.6557; %Volume fraction of PCM
%engine.mPCM = 3.456; % Mass of PCM [kg]
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


%% Find Efficiency with engine structure
Eff = findEfficiency(engine)

%% Optimizer Starts Here
% Parameters we want to optimizeL
% rhoL, Lh, f, 
% Array of the parameters we want to optimize: x
x = [];
% rhoL
x(1) = 730;
% Lh
x(2) = 240;
% f
x(3) = 0.65;
% Step 2: Run in findEfficiency3
Eff2 = optEfficiency(x)

%%
% Step 3: Run in fmincon
x0 = inputs;

x0(5) = .8; %L1
x0(6) = .09; % b1
x0(16) = 270; %Lh
x0(25) = 0.50; %f
x0(10) = 720; %rhoL
x0(13) = 17; %Tm
% options = optimoptions('fmincon', 'MaxFunctionEvaluations', 50000, 'MaxIterations', 3000)

[xopt, fval, exitflag, output] = fmincon(@objective, x0, [], [], [], [], [], [], @constraint, [])

[c, ceq] = constraint(x0)
effOpt = findEfficiency(xopt)
[c1, ceq1] = constraint(xopt)

% Results: Eff = 
% High as possible: 
% Inbetween: 
% Low as Possible: 



%% Optimizer Functions

% Define objective function for optimization
function obj = objective(x)
    obj = -findEfficiency(x);
end

% Define constraint for optimization
% c(x) <= 0
%ceq(x) = 0
function [c, ceq] = constraint(x)
    c(1) = x(16) - 270; % Lh < 270
    c(2) = 210 - x(16); % Lh < 210
    c(3) = x(10) - 800; % rhoL < 800
    c(4) = 720 - x(10); % rhoL > 720
    c(5) = x(12) - 0.75; % f < 0.75
    c(6) = 0.25 - x(12); % f > 0.25
    %c(7) = 5.5 - x(13); %Tm > 5.5
    %c(8) = x(13) - 20; %Tm > 20
end



% Wrapper Function: put x array into engine structure
function newEngine = addInX(x, engine)
    engine.rhoL = x(1);
    engine.Lh = x(2);
    engine.f = x(3);
    
    % Need to update these values with new 'x' values
    engine.rhoS = engine.rhoL + engine.delta_rho; %TTTTHIIISSSSSSSSS
    engine.mPCM = (engine.f*pi*engine.L1*( (engine.a1 / 2)^2 ))/(1/engine.rhoS);


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
    V = pi*engine.L1*( (engine.a1 / 2)^2 ); % Inner volume of Cylinder
    V1A = engine.ar*V; % Volume of residual air
    V1H = ( V*(1 - engine.f) ) - V1A; % volume of HF at state 1
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