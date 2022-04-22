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

% engine % prints out the engine structure
%%
s = struct('f1', [1, 3, 2, 4]);
a=[];
fields = string(fieldnames(engine));

lenS=length(fields);
for i=1:1:lenS
    if i > 3
        a = [ a engine.(fields(i))];
    end
end

%% Optimizer Testing
% intial guess
a(1) = 1; % Length of the cylinder [m]
% Call solver to minimize the objective function given the constraint

% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% problem.options = options;
% problem.solver = 'fmincon';
% problem.objective = objective();
% problem.x0 = [0,0];

x0 = [a, T, Tlow, Thigh, Po];
xopt = fmincon(@objective, x0, [], [], [], [], [], [], @constraint, [])
engine.L1 = xopt;
effOpt = findEfficiency2(T, Tlow, Thigh, Po, engine)

%% Functions
% 1. Finding Specific Volume of PCM [vP]
function [voP, v1P, VPCM] = specificVolPCM(T, rhoS, mPCM)
    % a. Find specific volume of PCM under ambient pressure Po [voP]
    voP = ( (1.0307e03 - ( 1.2596*(T + 273.15) )  + ...
        (1.8186e-3* (T + 273.15)^2) -(1.9555e-6* (T + 273.15)^3) ) )^-1;
    % b. Find the volume of PCM (liquid state) [v1P]
    v1P = 1/rhoS; % Specific volume of PCM in liquid state
    VPCM = mPCM*v1P;
end

% 3. Find inner volume of cylinder, Volume fraction of air and PCM 
function [V,V1A, f, V1H, mH] = findVolume(a1, L1, ar, VPCM, v1H)
    V = pi*L1*( (a1 / 2)^2 ); % Inner volume of Cylinder
    V1A = ar*V; % Volume of residual air
    f = VPCM / V; % volume fraction of PCM
    V1H = ( V*(1 - f) ) - V1A; % volume of HF at state 1
    mH = ( (V*(1 - f) ) - V1A) / v1H; % mass of HF
end

% % 5. Finding Max Pressure [P2]
% function P2 = findPressure(Po, engine)
% 
%     F = @(P)((pi / 4)*(engine.L1*( ( (2*engine.a1) + ...
%                 ( ( (P - Po)*engine.a1*(1 - engine.v^2) ) / engine.Ey)*( ( (engine.b1^2 + engine.a1^2) / (engine.b1^2 - engine.a1^2) ) +...
%             (engine.v / (1 - engine.v) ) ))*...
%                 ( ( (P - Po)*engine.a1*(1 - engine.v^2) ) / engine.Ey)*( ( (engine.b1^2 + engine.a1^2) / (engine.b1^2 - engine.a1^2) ) + ...
%             (engine.v / (1 - engine.v) ) )) ))... % delta_V1
%         - ...
%         (( engine.mPCM*((1.3e-03 - (2.66e-04*log10( 1 + ( (P - Po) / 102.12) ) )) - engine.v1P) ) + ...
%         ( engine.mH*((engine.voH - (engine.CH*log10(1 + ( (P - Po) / engine.BH) ) )) - engine.voH) ) + (((engine.V1A*Po) / P) - engine.V1A)); %delta_V2
%     
%         % P2 is the 1st instance where delta_V1 - delta_V2 == 0
%     options = optimoptions('fsolve','Display','iter','TolFun',1e-14);
%     
%     %solves for P, same answer as the engine.P2 with the current finding
%     %Pressure function
%     P2 = fsolve(F,5,options);
% 
% end
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
    P2 = findPressure(Po, engine);
    delta_V1 = findChangeInnerVolume(Po, engine.a1, engine.b1, engine.L1, engine.v, engine.Ey, P2);
    Pa = findPa(Po, P2, engine.V1N, delta_V1, V1A, V, f, v1P, voP, engine.CP, engine.BP, engine.CH, engine.BH, engine.v1H);
    Est = findEst(Pa, engine.V1N, engine.mPCM, engine.rhoL, engine.rhoS);
    Qin = findQin(Tlow, Thigh, engine.Tm, engine.mPCM, engine.csd, engine.Lh, engine.cld);
    Eff = Est / (Qin*1e3)*100;
end

function Eff = findEfficiency4(a)  
    [voP, v1P, VPCM] = specificVolPCM(a(35), a(6), a(8));
    [V,V1A, f, V1H, mH] = findVolume(a(2), a(1), a(20), VPCM, a(14));
    P2 = findPressure(a(38), a(2), a(5), a(4), a(3), V1A, a(8), mH, a(13), a(1), a(18), a(17), v1P);
    delta_V1 = findChangeInnerVolume(a(38), a(2), a(3), a(1), a(5), a(4), P2);
%     Pa = findPa(Po, P2, engine.V1N, delta_V1, V1A, V, f, v1P, voP, engine.CP, engine.BP, engine.CH, engine.BH, engine.v1H);
%     Est = findEst(Pa, engine.V1N, engine.mPCM, engine.rhoL, engine.rhoS);
%     Qin = findQin(Tlow, Thigh, engine.Tm, engine.mPCM, engine.csd, engine.Lh, engine.cld);
%     Eff = Est / (Qin*1e3)*100;
    Pa = findPa(a(38), P2, a(19), delta_V1, V1A, V, f, v1P, voP, a(15), a(16), a(18), a(17), a(14));
    Est = findEst(Pa, a(19), a(8), a(7), a(6));
    Qin = findQin(a(36), a(37), a(9), a(8), a(10), a(12), a(11));
    Eff = Est / (Qin*1e3)*100;
end


%Pressure to Stress Equations
function [sigma_tan, sigma_rad, sigma_long] = pressuretoStress(P2, a1, b1)
    % sigma_tan is the tangential stress
    sigma_tan = P2*(((b1/2)^2 + (a1/2)^2) / ((b1/2)^2 - (a1/2)^2));
    
    %sigma_rad is the radial stress
    sigma_rad = -P2;
    
    %sigma_long is the longitudinal stress on the ends
    sigma_long = (P2*(a1/2)^2)/((b1/2)^2 - (a1/2)^2);
end

function obj = objective(x0)
%     obj = findEfficiency4(T, Tlow, Thigh, Po, a);
    obj = findEfficiency4(x0);

end

function [c, ceq] = constraint(engine)
    c = engine.L1 - 3;
    ceq = [];
end
