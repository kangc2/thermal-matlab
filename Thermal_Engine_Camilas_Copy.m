%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%From Wang's An Investigation of PCM Based Ocean Thermal Energy Harvesting
% Chapter 3.4: The Deformation of the cylindrical tubes 
% under internal pressure loading using linear elasticity theory

% PNNL's Mathmatica Script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Notes: 

% Updated 2/10: Camila
% 1. Used vpasolve() to find P
% been playing around with the formatting/documenting

% Original 2/7: Rivan
% 1. Wrote up Data Input & Calculations Sections

   
%%
clc; clear; close all
%% Data Input
% PCM Material = Hexadecane (C16H34)
% Working fluid = Water
% Cylinder Material = TA2M (Titanium alloy)

T = 29; %[degC]
Thigh = T; % High Temperature
Tm = 18.2; % Melting Temperature of PCM [degC]
Tlow = 5; %L ow Temperature [degC]
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
voP = ( (1.0307e03 - ( 1.2596*(T + 273.15) ) ) + ...
    (1.8186e03*( (T + 273.15)^2) )-(1.9555e-06*( (T + 273.15)^3) ) )^-1;
CP = 2.66e-04;
BP = 102.12;
v1P = 1/rhoS; % Specific volume of PCM in liquid state
BH = ( 2672.9 + (15.97*T) - (0.166*(T^2) ) )*1e-01;
voH = 1e-03;
v1H = voH;
CH = 0.3150*voH;
V = pi*L1*( (a1 / 2)^2 ); % Inner volume of Cylinder
VPCM = mPCM*v1P;
V1A = ar*V; % Volume of residual air
rPCM = VPCM / V;
VH1 = ( V*(1 - rPCM) ) - V1A;
mH = ( (V*(1 - rPCM) ) - V1A) / v1H;

syms P;

delta_a1 = ( ( (P - Po)*a1*(1 - v^2) ) / Ey)*( ( (b1^2 + a1^2) / (b1^2 - a1^2) ) + (v / (1 - v) ) );
delta_V1 = (pi / 4)*(L1*( ( (2*a1) + delta_a1)*delta_a1) );
vP = 1.3e-03 - (2.66e-04*log10( 1 + ( (P-Po) / 102.12) ) );
v1H = voH;
vH = voH - (CH*log10(1 + ( (P - Po) / BH) ) );
VA = (V1A*Po) / P;
delta_V2 = ( mPCM*(vP - v1P) ) + ( mH*(vH - voH) ) + (VA - V1A);

vpasolve(delta_V1 - delta_V2 == 0, P)
   
