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
engine.b1 = 0.15; % External diameter of the cylinder [m]
t = 0.035; % wall thickness [m]
engine.a1 = engine.b1 - 2*t; % Internal diameter of the cylinder [m]

% -> Hull Material Properties
engine.Ey = 105e03; % Young's Modulus [MPa]
engine.v = 0.33; % Poisson's ratio

% -> Material Properties of PCM
engine.rhoS = 864; % Density of PCM - Solid phase [kg/m3]
engine.rhoL = 773; % Density of PCM - Liquid phase [kg/m3]

f = 0.6557; %Volume fraction of PCM
engine.mPCM = (f*pi*engine.L1*( (engine.a1 / 2)^2 ))/(1/engine.rhoS);

% engine.mPCM = 3.456; % Mass of PCM [kg]
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

% Inputs List
% x(1) = T
% x(2) = Thigh
% x(3) = Tlow
% x(4) = Po
% x(5) = L1
% x(6) = b1
% x(7) = a1
% x(8) = Ey
% x(9) = v
% x(10) = rhoS
% x(11) = rhoL
% x(12) = mPCM
% x(13) = Tm
% x(14) = csd
% x(15) = cld
% x(16) = Lh
% x(17) = voH
% x(18) = v1H
% x(19) = CP
% x(20) = BP
% x(21) = BH
% x(22) = CH
% x(23) = V1N
% x(24) = ar
%%
% Step 2: Run in findEfficiency3
Eff = findEfficiency(inputs)

%%
% Step 3: Run in fmincon
x0 = inputs;
x0(12) = 3;
x0(5) = .8;
x0(6) = .12;
x(16) = 240;
x(10) = 800;
xopt = fmincon(@objective, x0, [], [], [], [], [], [], @constraint, []);
effOpt = findEfficiency(xopt)

%% Functions: Optimization

% Define objective function for optimization
function obj = objective(x)
    obj = -findEfficiency(x);
end

% Define constraint for optimization
function [c, ceq] = constraint(x)
    c(1) = x(5) - 1.3; % L1 < 1.3
    c(2) = 0.5 - x(5); % L1 > 0.5
    c(3) = x(6) - 0.15; % b1 < 0.15
    c(4) = 0.09 - x(6); % b1 > 0.08 
    c(5) = x(12) - 5; % mPCM < 15
    c(6) = x(16) - 270; % Lh < 270
    c(7) = 210 - x(16); % Lh < 210
    c(8) = x(10) - 1200; % rhoS < 1200
    c(9) = 750 - x(10); % rhoS > 750

    ceq(1) = x(1) - 29; % T = 29
    ceq(2) = x(2) - 29; % Thigh = 29
    ceq(3) = x(3) - 5; % Tlow = 5
    ceq(4) = x(4) - 0.101; % Po = 0.101
    ceq(5) = x(8) - 105e3; % Ey = 105e3
    ceq(6) = x(9) - 0.33; % v = 0.33
    ceq(7) = x(14) - 1.64; % csd = 7.7e2 
    ceq(8) = x(15) - 2.09; % cld = 7.7e2
    ceq(9) = x(17) - 1/1000; % voH = 1/1000
    ceq(10) = x(18) - 1/1000; % v1H = 1/1000
    ceq(11) = x(19) - 2.66e-4; % CP = 2.66e-4
    ceq(12) = x(20) - 102.12; % BP = 102.12
    ceq(13) = x(21) - 2.996424e2; % BH  = 2.99e2
    ceq(14) = x(22) - (1/1000)*0.3150; % CH = (1/1000)*0.315
    ceq(15) = x(23) - 2e-3; % V1N = 2e-3
    ceq(16) = x(24) - 6.573/100; % ar = 6.573/100
    ceq(17) = x(13) - 18.2; % Tm = 18.2
    ceq(18) = x(7) - (x(6) - 2*.035); %a1 = b1 - 2*t

end

% OPTIMIZER FUNCTION FOR EFFICIENCY
function Eff = findEfficiency(x)  
    % 1. Finding Specific Volume of PCM [vP]
    voP = ( (1.0307e03 - ( 1.2596*(x(1) + 273.15) )  + ...
        (1.8186e-3* (x(1) + 273.15)^2) -(1.9555e-6* (x(1) + 273.15)^3) ) )^-1;
    v1P = 1 / x(10); % Specific volume of PCM in liquid state
    VPCM = x(12)*v1P;


    % 3. Finding Volume Of Cylinder
    V = pi*x(5)*( (  x(7) / 2)^2 ); % Inner volume of Cylinder
    
    V1A = x(24)*V; % Volume of residual air
    rPCM = VPCM / V;
    VH1 = ( V*(1 - rPCM) ) - V1A;
    mH = ( (V*(1 - rPCM) ) - V1A) / x(18);
    f = rPCM;

        F = @(P)((pi / 4)*(x(5)*( ( (2*  x(7)) + ...
        ( ( (P - x(4))*  x(7)*(1 - x(9)^2) ) / x(8)) ...
            *( ( ( x(6)^2 +  x(7)^2) / ( x(6)^2 -   x(7)^2) ) +...
        (x(9) / (1 - x(9)) ) ))*...
            ( ( (P - x(4))*  x(7)*(1 - x(9)^2) ) / x(8)) ...
            *( ( ( x(6)^2 +  x(7)^2) / ( x(6)^2 -   x(7)^2) ) + ...
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
    delta_a1 = ( ( (P2 - x(4))*  x(7)*(1 - x(9)^2) ) / x(8))*( ( ( x(6)^2 +   x(7)^2) / ( x(6)^2 -   x(7)^2) ) + (x(9) / (1 - x(9)) ) );
    delta_V1 = (pi / 4)*(x(5)*( ( (2*  x(7)) + delta_a1)*delta_a1) );
    
    Pa = (P2 / x(23))*(delta_V1 + x(23) - V1A*( (x(4) / P2) - 1) ...
       - (V*f / v1P)*(voP - x(19)*log10(1 + ((P2 - x(4)) / x(20)) ) - v1P) + ((V*(1 - f) - V1A) / x(18)) ...
       *x(22)*log10(1 + ((P2 - x(4)) / x(21)) ));
    
    Qin = x(12)*x(14)*(x(13) - x(3)) + x(12)*x(16) + x(12)*x(15)*(x(2) - x(13));
    Est = -Pa*1e6*x(23)*log(1 - (x(12) / x(23))*((1 / x(11)) - (1 / x(10))) );
    Eff = Est / (Qin*1e3) * 100;
end

function stress = PtoStress(x)
    % sigma_tan is the tangential stress
    sigma_tan = P2*(((  x(6)/2)^2 + ( x(7)/2)^2) / ((  x(6)/2)^2 - ( x(7)/2)^2));
    
    %sigma_rad is the radial stress
    sigma_rad = -P2;
    
    %sigma_long is the longitudinal stress on the ends
    sigma_long = (P2*( x(7)/2)^2)/((  x(6)/2)^2 - ( x(7)/2)^2);
    stress = max([sigma_tan, sigma_rad, sigma_long]);
end
