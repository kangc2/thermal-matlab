%Check Parameters with STRUCTURES
clc; close all; clear
%% Environmental Parameters
% Variables of the environment, not specific to the engine
T = 29; % Working Temperature [degC]
Thigh = T; % High Temperature
Tlow = 5; % Low Temperature [degC]
Po = 0.101; % Initial Pressure/ambient pressure = atmospheric pressure [MPa]
%% Define Intial Structure: get all inputs in an 'engine' structure
% Thermal Engine from our Thermal_Engine code
% The baseline of our model
% Note: make sure the field name is right (voH vs. VoH)
engine.name = 'Original';
engine.working_fluid = 'Water';
engine.hull_material = 'TA2M (Titanium alloy)';

% -> Geometry & Material Properties of Cylinder
engine.L1 = 1.31; % Length of the cylinder [m]
engine.a1 = 7.7e-02; % Internal diameter of the cylinder [m]
engine.b1 = 8.6e-02; % External diameter of the cylinder [m]
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
%% Define Outputs in 'engine' structure using Functions
% 1 Find specific volume of PCM under ambient pressure Po, volume and mass of PCM [voP]
[engine.voP, engine.v1P, engine.VPCM] = specificVolPCM(T, engine.rhoS, engine.mPCM);
% 3. Find inner volume of cylinder, Volume fraction of air and PCM 
[engine.V,engine.V1A, engine.f, engine.V1H, engine.mH] = findVolume(engine.a1, engine.L1, engine.ar, engine.VPCM, engine.v1H);

% 5. Finding Max Pressure [P2]
engine.P2 = findPressure(Po, engine.a1, engine.v, engine.Ey, engine.b1, engine.V1A, ...
    engine.mPCM, engine.mH, engine.voH, engine.L1, engine.CH, engine.BH, engine.v1P);

% 4. Finding change of Inner Diameter of Tube Under Internal Pressure
engine.delta_V1 = findChangeInnerVolume(Po, engine.a1, engine.b1, engine.L1, engine.v, engine.Ey, engine.P2);

% 6a. Find Pre-charged pressure in acculimlator
engine.Pa = findPa(Po, engine.P2, engine.V1N, engine.delta_V1, engine.V1A, engine.V, engine.f, engine.v1P, ...
    engine.voP, engine.CP, engine.BP, engine.CH, engine.BH, engine.v1H);
% 6b. find the total energy stored
engine.Est = findEst(engine.Pa, engine.V1N, engine.mPCM, engine.rhoL, engine.rhoS);
% 6c. Find total thermal energy absorbed by PCM in 1 melting/freezing cycle
engine.Qin = findQin(Tlow, Thigh, engine.Tm, engine.mPCM, engine.csd, engine.Lh, engine.cld);
% 6d. find the theorectical Efficiency % 
engine.Eff = findEfficiency(engine.Est, engine.Qin);  

%% Find Efficiency: get efficiency using all of our inputs
% adds a new structure field Eff2, same as Eff, but using function that
% takes in all the inputs at once

engine.Eff2 = findEfficiency2(T, Tlow, Thigh, Po, engine);

engine % prints out the engine structure

%% Loops thru every parameter and returns a structured data for each parameter
%create a new structure 'parameter' to add in all our of tests into structures

% Parameters to change: a new structure called input_range
% includes all the parameters we want to test via structures
% Field name, lowerbound, upperbound, delta, original value
input_range.L1 = [1, 3, engine.L1];
input_range.b1 = [0.08, 0.1, engine.b1];
input_range.Lh = [0.49, 0.69, engine.Lh];
input_range.Tm = [18, 37, engine.Tm];
input_range.rhoS = [700, 1000, engine.rhoS];
input_range.rhoL = [600, 900, engine.rhoL];
input_range.mPCM = [3, 5, engine.mPCM];
fields = string(fieldnames(input_range));

delta = 21; % number of steps between upper and lower bounds, 
% ^^NOTE: we can put this into the structure field list, so we can have different values of delta for each parameter

original = engine; % to reconvert the engine back to its original values when the second for loop finishes

% For loop that goes thru every variable we want to change
for j = 1:length(fields) % For each parameter we want to change
    input = input_range.(fields(j)); % get the information from each field in the structure
    l = input(1); % lower bound value
    u = input(2); % upper bound value
    
    j % check at which variable is passing thru

    newParam = []; % this is all the changed values for each variable goes
    newEff = []; % this is where each efficiency value is for each variable
   
    for i= 0:delta
       % Find what the value is for the changing variable
        z = ((u - l) / delta) * i;
        original.(fields(j)) = z + l; % Changed value for a variable
        
        % Find Efficiency value
        % each changed value is represented by p(#) instead of engine.parameter
        answer = findEfficiency2(T, Tlow, Thigh, Po, original);
        % add all changed values into a row
        newParam = [newParam, original.(fields(j))];
        % add all the efficiency values into one row
        newEff = [newEff, answer];
    end
    
   % adds into a new 'parameter' structure
   % New structure: parameter
   % fields: name, parameter, efficiency
   parameter(j).name = fields(j);
   parameter(j).param = newParam;
   parameter(j).efficiencies = newEff;

   original.(fields(j)) = input(3); % reverts the value back to its original value
end


%% Functions
% 1. Finding Specific Volume of PCM [vP]
function [voP, v1P, VPCM] = specificVolPCM(T,rhoS, mPCM)
    % a. Find specific volume of PCM under ambient pressure Po [voP]
    voP = ( (1.0307e03 - ( 1.2596*(T + 273.15) )  + ...
        (1.8186e-3* (T + 273.15)^2) -(1.9555e-6* (T + 273.15)^3) ) )^-1;
    % b. Find the volume of PCM (liquid state) [v1P]
    v1P = 1/rhoS; % Specific volume of PCM in liquid state
    VPCM = mPCM*v1P;
end

% 3. Find inner volume of cylinder, Volume fraction of air and PCM 
function [V,V1A, f,V1H, mH] = findVolume(a1, L1, ar, VPCM, v1H)
    V = pi*L1*( (a1 / 2)^2 ); % Inner volume of Cylinder
    V1A = ar*V; % Volume of residual air
    f = VPCM / V; % volume fraction of PCM
    V1H = ( V*(1 - f) ) - V1A; % volume of HF at state 1
    mH = ( (V*(1 - f) ) - V1A) / v1H; % mass of HF
end

% 5. Finding Max Pressure [P2]
function P2 = findPressure(Po, a1, v, Ey, b1, V1A, mPCM, mH, voH, L1, CH, BH, v1P)
    % creates a symbolic variable P to plug into equations and solve for it
    syms P; 

    % Change in inner diameter of cylinder at pressure P
    delta_a1 = ( ( (P - Po)*a1*(1 - v^2) ) / Ey)*( ( (b1^2 + a1^2) / (b1^2 - a1^2) ) + (v / (1 - v) ) );
    
    % First Equation of change in inner volume of cylinder at pressure P,
    % based on geometry of cylinder
    delta_V1 = (pi / 4)*(L1*( ( (2*a1) + delta_a1)*delta_a1) );
   
    vP = 1.3e-03 - (2.66e-04*log10( 1 + ( (P - Po) / 102.12) ) ); % specific volume of PCM
    vH = voH - (CH*log10(1 + ( (P - Po) / BH) ) ); % specific volume of HF
    VA = (V1A*Po) / P; % Volume of residual air
   
    % Second Equation of change in inner volume of cylinder at pressure P,
    % based on HF and PCM mass and volume properties
    delta_V2 = ( mPCM*(vP - v1P) ) + ( mH*(vH - voH) ) + (VA - V1A); 

    % P2 is the 1st instance where delta_V1 - delta_V2 == 0
    P2 = vpasolve(delta_V1 - delta_V2 == 0, P);
    P2 = double(P2); % Convert to a numerical value with precision

end

% 4. Finding change of Inner Diameter of Tube Under Internal Pressure
function delta_V1 = findChangeInnerVolume(Po, a1, b1, L1, v, Ey, P2)
    delta_a1 = ( ( (P2 - Po)*a1*(1 - v^2) ) / Ey)*( ( (b1^2 + a1^2) / ...
        (b1^2 - a1^2) ) + (v / (1 - v) ) );
    delta_V1 = (pi / 4)*(L1*( ( (2*a1) + delta_a1)*delta_a1) );
end

% 6. Finding Efficiency [Eff] 
% a. Find Pre-charged pressure in acculimlator [Pa]
function Pa = findPa(Po, P2, V1N, delta_V1, V1A, V, f, v1P, voP, CP, BP, CH, BH, v1H)
  Pa = (P2 / V1N)*(delta_V1 + V1N - V1A*( (Po / P2) - 1) ...
       - (V*f / v1P)*(voP - CP*log10(1 + ((P2 - Po) / BP) ) - v1P) + ...
       ((V*(1 - f) - V1A) / v1H)*CH*log10(1 + ((P2 - Po) / BH) ));  
end

% b. find the total energy stored [Est]
function Est = findEst(Pa, V1N, mPCM, rhoL, rhoS)
    Est = -Pa*1e6*V1N*log(1 - (mPCM / V1N)*((1 / rhoL) - (1 / rhoS)) );
end

% c. Find total thermal energy absorbed by PCM in 1 melting/freezing cycle
% [Qin]
function Qin = findQin(Tlow, Thigh, Tm, mPCM, csd, Lh, cld)
    Qin = mPCM*csd*(Tm - Tlow) + mPCM*Lh + mPCM*cld*(Thigh - Tm);
end

% d. find the theorectical Efficiency % [Eff] 
function Eff = findEfficiency(Est, Qin)  
    Eff = Est / (Qin*1e3) * 100;
end

% This takes in all inputs and finds the Efficiency using all functions
% above
function Eff = findEfficiency2(T, Tlow, Thigh, Po, engine)  
    [voP, v1P, VPCM] = specificVolPCM(T, engine.rhoS, engine.mPCM);
    [V,V1A, f, V1H, mH] = findVolume(engine.a1, engine.L1, engine.ar, VPCM, engine.v1H);
    P2 = findPressure(Po, engine.a1, engine.v, engine.Ey, engine.b1, V1A, engine.mPCM, mH, engine.voH, ...
        engine.L1, engine.CH, engine.BH, v1P);
    delta_V1 = findChangeInnerVolume(Po, engine.a1, engine.b1, engine.L1, engine.v, engine.Ey, P2);
    Pa = findPa(Po, P2, engine.V1N, delta_V1, V1A, V, f, v1P, voP, engine.CP, engine.BP, engine.CH, engine.BH, engine.v1H);
    Est = findEst(Pa, engine.V1N, engine.mPCM, engine.rhoL, engine.rhoS);
    Qin = findQin(Tlow, Thigh, engine.Tm, engine.mPCM, engine.csd, engine.Lh, engine.cld);
    Eff = Est / (Qin*1e3)*100;
end

%Pressure to Stress Equations
function [sigma_tan, sigma_rad, sigma_long] = PtoStress(P2, a1, b1)
    % sigma_tan is the tangential stress
    sigma_tan = P2*(((b1/2)^2 + (a1/2)^2) / ((b1/2)^2 - (a1/2)^2));
    
    %sigma_rad is the radial stress
    sigma_rad = -P2;
    
    %sigma_long is the longitudinal stress on the ends
    sigma_long = (P2*(a1/2)^2)/((b1/2)^2 - (a1/2)^2);
end
