%Check Parameters with STRUCTURES
clc; close all; clear
%% Constants: Environmental Parameters
% Variables of the environment

T = 29; % Working Temperature [degC]
Thigh = T; % High Temperature
Tlow = 5; % Low Temperature [degC]
Po = 0.101; % Initial Pressure/ambient pressure = atmospheric pressure [MPa]

%% Parameters: Declare Variables and intial guesses
% Variables can be changed by the optimizer

% -> Intro names
engine.name = "Original";
engine.working_fluid = "Water";
engine.hull_material = "TA2M (Titanium alloy)";

% -> Geometry Properties of Cylinder
engine.L1 = 1.31; % Length of the cylinder [m]
engine.a1 = 7.7e-02; % Internal diameter of the cylinder [m]
engine.b1 = 8.6e-02; % External diameter of the cylinder [m]

% -> Hull Material Properties
engine.Ey = 105e03; % Young's Modulus [MPa]
engine.v = 0.33; % Poisson's ratio

% -> Material Properties of PCM
engine.rhoS = 864; % Density of PCM - Solid phase [kg/m3]
engine.rhoL = 773; % Density of PCM - Liquid phase [kg/m3]
engine.mPCM = 3.456; % Mass of PCM [kg]
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
%% Find Efficiency: get efficiency using all of our inputs
% adds a new structure field Eff2, same as Eff, but using function that
% takes in all the inputs at once

% engine.Eff2 = findEfficiency2(T, Tlow, Thigh, Po, engine);

engine % prints out the engine structure

%% Optimization: 
% Step 1: Put all of the engine and environmental inputs into one array

inputs = [T Thigh Tlow Po];
fields = string(fieldnames(engine));

lenEngine = length(fields);
for i = 1:1:lenEngine
    if i > 3 % skip the 3 naming fields
        inputs = [ inputs engine.(fields(i))];
    end
end
format long
inputs
%%
% Step 2: Run in findEfficiency3
Eff = findEfficiency3(inputs)


%%
% Step 3: Run in fmincon
x0 = inputs;
x0(5) = 1.3;
xopt = fmincon(@objective, x0, [], [], [], [], [], [], @constraint, []);
xopt(5)
effOpt = findEfficiency3(xopt)
%% Functions
% Define objective function for optimization
function obj = objective(x)
    obj = -findEfficiency3(x);
end

% Define constraint for optimization
function [c, ceq] = constraint(x)
%     [V,V1A, f, V1H, mH] = findVolume(a1, L1, ar, VPCM, v1H)
%     [V,V1A, f, V1H, mH] = findVolume(x(1), x(2), x(3), x(4), x(5));
    c = x(5) - 1.3;
    ceq = [];
end

% OPTIMIZER FUNCTION FOR EFFICIENCY
function Eff = findEfficiency3(x)  
    % 1. Finding Specific Volume of PCM [vP]
    voP = ( (1.0307e03 - ( 1.2596*(x(1) + 273.15) )  + ...
        (1.8186e-3* (x(1) + 273.15)^2) -(1.9555e-6* (x(1) + 273.15)^3) ) )^-1;
    v1P = 1 / x(10); % Specific volume of PCM in liquid state
    VPCM = x(12)*v1P;


    % 3. Finding Volume Of Cylinder
    V = pi*x(5)*( (x(6) / 2)^2 ); % Inner volume of Cylinder
    
    V1A = x(24)*V; % Volume of residual air
    rPCM = VPCM / V;
    VH1 = ( V*(1 - rPCM) ) - V1A;
    mH = ( (V*(1 - rPCM) ) - V1A) / x(18);
    f = rPCM;

    F = @(P)((pi / 4)*(x(5)*( ( (2*x(6)) + ...
        ( ( (P - x(4))*x(6)*(1 - x(9)^2) ) / x(8)) ...
            *( ( (x(7)^2 + x(6)^2) / (x(7)^2 - x(6)^2) ) +...
        (x(9) / (1 - x(9)) ) ))*...
            ( ( (P - x(4))*x(6)*(1 - x(9)^2) ) / x(8)) ...
            *( ( (x(7)^2 + x(6)^2) / (x(7)^2 - x(6)^2) ) + ...
        (x(9) / (1 - x(9)) ) )) ))... % delta_V1
        - ...
        (( x(12)*((1.3e-03 - (2.66e-04*log10( 1 + ( (P - x(4)) / 102.12) ) )) - v1P) ) + ...
        ( mH*((x(17) - (x(22)*log10(1 + ( (P - x(4)) / x(21)) ) )) - x(17)) ) ...
        + (((V1A*x(4)) / P) - V1A)); %delta_V2

    % Options: sets tolerance of function close to 0 (1-e14) and displays the
    % iteration, this could help with the optimization
    options = optimoptions('fsolve','Display','iter','TolFun',1e-14);
    
    %solves for P, same answer as the engine.P2 with the current finding
    %Pressure function
    P2 = fsolve(F,5,options);
    
    % 6. Finding Efficiency [Eff] 
    delta_a1 = ( ( (P2 - x(4))*x(6)*(1 - x(9)^2) ) / x(8))*( ( (x(7)^2 + x(6)^2) / (x(7)^2 - x(6)^2) ) + (x(9) / (1 - x(9)) ) );
    delta_V1 = (pi / 4)*(x(5)*( ( (2*x(6)) + delta_a1)*delta_a1) );
    
    Pa = (P2 / x(23))*(delta_V1 + x(23) - V1A*( (x(4) / P2) - 1) ...
       - (V*f / v1P)*(voP - x(19)*log10(1 + ((P2 - x(4)) / x(20)) ) - v1P) + ((V*(1 - f) - V1A) / x(18)) ...
       *x(22)*log10(1 + ((P2 - x(4)) / x(21)) ));
    
    Qin = x(12)*x(14)*(x(13) - x(3)) + x(12)*x(16) + x(12)*x(15)*(x(2) - x(13));
    Est = -Pa*1e6*x(23)*log(1 - (x(12) / x(23))*((1 / x(11)) - (1 / x(10))) );
    Eff = Est / (Qin*1e3) * 100;

end