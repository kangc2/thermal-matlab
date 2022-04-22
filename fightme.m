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
inputs
%%
% Step 2: Run in findEfficiency3
Eff = findEfficiency3(inputs)


%%
% Step 3: Run in fmincon
x0 = inputs;
x0(5) = 0.8;
x0(6) = 7.7e-2;

xopt = fmincon(@objective, x0, [], [], [], [], [], [], @constraint, []);
xopt(5)
effOpt = findEfficiency3(xopt)
%% Functions: Optimization

% Define objective function for optimization
function obj = objective(x)
    obj = -findEfficiency3(x);
end

% Define constraint for optimization
function [c, ceq] = constraint(x)
    c(1) = x(5) - 1.3; % L1
    c(2) = 0.5 - x(5); % L1
    c(3) = x(7) - 1e-01; % b1
    c(4) = 8e-2 - x(7); % mPCM
    c(5) = x(12) - 5; % mPCM
    c(6) = 18 - x(13); % Tm

    ceq(1) = x(1) - 29; % T
    ceq(2) = x(2) - 29; % Thigh
    ceq(3) = x(3) - 5; % Tlow
    ceq(4) = x(4) - 0.101; % Po
    ceq(5) = x(8) - 105e3; % Ey
    ceq(6) = x(9) - 0.33; % v
    ceq(7) = x(10) - 773; % rhoS
    ceq(8) = x(11) - 864; % rhoL
    ceq(9) = x(6) - 7.7e-2; % a1
    ceq(10) = x(14) - 7.7e-2; % csd
    ceq(11) = x(15) - 7.7e-2; % cld
    ceq(12) = x(16) - 18.2; % Lh
    ceq(13) = x(17) - 1/1000; % voH
    ceq(14) = x(18) - 1/1000; % a1
    ceq(15) = x(19) - 2.66e-4; % CP
    ceq(16) = x(20) - 102.12; % BP
    ceq(17) = x(21) - 2.996424e2; % BH
    ceq(18) = x(22) - (1/1000)*0.3150; % CH
    ceq(19) = x(23) - 2e-3; % V1N
    ceq(20) = x(24) - 6.573/100; % ar


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
function [V,V1A, f, V1H, mH] = findVolume(a1, L1, ar, VPCM, v1H)
    V = pi*L1*( (a1 / 2)^2 ); % Inner volume of Cylinder
    V1A = ar*V; % Volume of residual air
    f = VPCM / V; % volume fraction of PCM
    V1H = ( V*(1 - f) ) - V1A; % volume of HF at state 1
    mH = ( (V*(1 - f) ) - V1A) / v1H; % mass of HF
end

% 5. Finding Max Pressure [P2]
function P2 = findPressure(Po, a1, v, Ey, b1, V1A, mPCM, mH, voH, L1, CH, BH, v1P)
%     F = @(P)((pi / 4)*(engine.L1*( ( (2*engine.a1) + ...
%             ( ( (P - Po)*engine.a1*(1 - engine.v^2) ) / engine.Ey)*( ( (engine.b1^2 + engine.a1^2) / (engine.b1^2 - engine.a1^2) ) +...
%         (engine.v / (1 - engine.v) ) ))*...
%             ( ( (P - Po)*engine.a1*(1 - engine.v^2) ) / engine.Ey)*( ( (engine.b1^2 + engine.a1^2) / (engine.b1^2 - engine.a1^2) ) + ...
%         (engine.v / (1 - engine.v) ) )) ))... % delta_V1
%     - ...
%     (( engine.mPCM*((1.3e-03 - (2.66e-04*log10( 1 + ( (P - Po) / 102.12) ) )) - engine.v1P) ) + ...
%     ( engine.mH*((engine.voH - (engine.CH*log10(1 + ( (P - Po) / engine.BH) ) )) - engine.voH) ) + (((engine.V1A*Po) / P) - engine.V1A)); %delta_V2
% 
%     % Options: sets tolerance of function close to 0 (1-e14) and displays the
%     % iteration, this could help with the optimization
%     options = optimoptions('fsolve','Display','iter','TolFun',1e-14);
%     
%     %solves for P, same answer as the engine.P2 with the current finding
%     %Pressure function
%     P2 = fsolve(F,5,options);
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
function [sigma_tan, sigma_rad, sigma_long] = pressuretoStress(P2, a1, b1)
    % sigma_tan is the tangential stress
    sigma_tan = P2*(((b1/2)^2 + (a1/2)^2) / ((b1/2)^2 - (a1/2)^2));
    
    %sigma_rad is the radial stress
    sigma_rad = -P2;
    
    %sigma_long is the longitudinal stress on the ends
    sigma_long = (P2*(a1/2)^2)/((b1/2)^2 - (a1/2)^2);
end
