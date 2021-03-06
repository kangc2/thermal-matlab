%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%From Wang's An Investigation of PCM Based Ocean Thermal Energy Harvesting
% Chapter 3.4: The Deformation of the cylindrical tubes 
% under internal pressure loading using linear elasticity theory

% PNNL's Mathmatica Script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
%%
clc; clear; close all
%% Data Input
% PCM Material = Hexadecane (C16H34)
% Working fluid = Water
% Cylinder Material = TA2M (Titanium alloy)

T = 29; % Working Temperature [degC]
Thigh = T; % High Temperature
Tm = 18.2; % Melting Temperature of PCM [degC]
Tlow = 5; % Low Temperature [degC]
csd = 1.64; % Specific heat in solid state [kJ/kg K]
cld = 2.09; % Specific heat in liquid state [kJ/kg K]
Lh = 236; % Latent heat of fusion [kJ/kg]
Ey = 105e03; % Young's Modulus [MPa]
v = 0.33; % Poisson's ratio
L1 = 1.31; % Length of the cylinder [m]
a1 = 7.7e-02; % Internal diameter of the cylinder [m]
b1 = 8.6e-02; % External diameter of the cylinder [m]
Po = 0.101; % Initial Pressure/ambient pressure = atmospheric pressure [MPa]
rhoS = 864; % Density of PCM - Solid phase [kg/m3]
rhoL = 773; % Density of PCM - Liquid phase [kg/m3]
mPCM = 3.456; % Mass of PCM [kg]
ar = 6.573/100; % Volume fraction of residual air
V1N = 2e-03; % Initial Volume of N2 gas

%% Calculations
% 1. Finding Specific Volume of PCM [vP]
voP = ( (1.0307e03 - ( 1.2596*(T + 273.15) )  + ...
    (1.8186e-3* (T + 273.15)^2) -(1.9555e-6* (T + 273.15)^3) ) )^-1;
CP = 2.66e-04;
BP = 102.12;
v1P = 1/rhoS; % Specific volume of PCM in liquid state
VPCM = mPCM*v1P;

% 2. Finding Specific Volume of Hydraulic Fluid [vH]
BH = ( 2672.9 + (15.97*T) - (0.166*(T^2) ) )*1e-01;
voH = 1e-03;
v1H = voH;
CH = 0.3150*voH;

% 3. Finding Volume Of Cylinder
V = pi*L1*( (a1 / 2)^2 ); % Inner volume of Cylinder

V1A = ar*V; % Volume of residual air
rPCM = VPCM / V;
VH1 = ( V*(1 - rPCM) ) - V1A;
mH = ( (V*(1 - rPCM) ) - V1A) / v1H;
f = rPCM;
% 4. Finding change of Inner Diameter of Tube Under Internal Pressure
syms P;
delta_a1 = ( ( (P - Po)*a1*(1 - v^2) ) / Ey)*( ( (b1^2 + a1^2) / (b1^2 - a1^2) ) + (v / (1 - v) ) );
delta_V1 = (pi / 4)*(L1*( ( (2*a1) + delta_a1)*delta_a1) );
vP = 1.3e-03 - (2.66e-04*log10( 1 + ( (P-Po) / 102.12) ) );
v1H = voH;
vH = voH - (CH*log10(1 + ( (P - Po) / BH) ) );
VA = (V1A*Po) / P;
delta_V2 = ( mPCM*(vP - v1P) ) + ( mH*(vH - voH) ) + (VA - V1A);
% 5. Finding Max Pressure [P2]
P2 = vpasolve(delta_V1 - delta_V2 == 0, P); % P2 is the 1st instance where delta_V1 - delta_V2 == 0
P = double(P2); % Convert to a numerical value with precision
% 6. Finding Efficiency [Eff] 
delta_a1 = ( ( (P - Po)*a1*(1 - v^2) ) / Ey)*( ( (b1^2 + a1^2) / (b1^2 - a1^2) ) + (v / (1 - v) ) );
delta_V1 = (pi / 4)*(L1*( ( (2*a1) + delta_a1)*delta_a1) );

Pa = (P2 / V1N)*(delta_V1 + V1N - V1A*( (Po / P2) - 1) ...
   - (V*f / v1P)*(voP - CP*log10(1 + ((P2 - Po) / BP) ) - v1P) + ((V*(1 - f) - V1A) / v1H) ...
   *CH*log10(1 + ((P2 - Po) / BH) ))

Qin = mPCM*csd*(Tm - Tlow) + mPCM*Lh + mPCM*cld*(Thigh - Tm)
Est = -Pa*1e6*V1N*log(1 - (mPCM / V1N)*((1 / rhoL) - (1 / rhoS)) )
Eff = Est / (Qin*1e3) * 100;
%7. Finding Energy Output [W] (kJ)
W = Eff*Qin*mPCM; 

%% Outputs
V
P2
Pa
f = f*100
Eff
