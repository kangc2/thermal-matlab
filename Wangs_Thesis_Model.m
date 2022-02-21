%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%From Wang's An Investigation of PCM Based Ocean Thermal Energy Harvesting
% Chapter 2: Harvesting Environmental Thermal Energy % 
% Using Solid/Liquid Phase Change Materials

% PNNL's Mathmatica Script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clc; close all; clear
%% Data Inputs
T = 20; % [degC]
Ey = 190e3; % Young's Modulus
v = 0.265; % Poisson's Ratio
L1 = 0.52; % Length of the cylinder of tube 1 [m]
L2 = 0.1; % Internal diameter of the cylinder of tube 2 [m]
a1 = 2.12e-2; % Internal diameter of the cylinder 1 [m]
a2 = 1.02e-2; % Internal diameter of the cylinder 2 [m]
b1 = 2.54e-2; % External diameter of the cylinder 1 [m]
b2 = 1.26e-2; % External diameter of the cylinder 2 [m]
Po = 0.101; % Initial Pressure/ambient pressure = atmospheric pressure [MPa]
mp = 0.07; % Mass of PCM [kg]
mH = 0.1091; % Mass of hydraulic fluid [kg]
v1P = 1.16e-3; % specific volume of PCM at state 1
V1A = 6.9078e-6; % volume of residual air at state 1
V1N = 8e-05; % Initial Volume of N2 gas
f = 0.6652;

%% Calculations
% Finding properties of PCM [subscript P]
voP = ( (1.0307e03 - ( 1.2596*(T + 273.15) )  + ...
    (1.8186e-3* (T + 273.15)^2) -(1.9555e-6* (T + 273.15)^3) ) )^-1;
BP = 762.8 - 4.805*(T - 79.4) + 0.0116*(T - 79.4)^2;
CP = 0.2058*voP;

% Finding properties of hydraulic oil [subscript H]
voH = 1e-03;
BH = (2672.9 + 15.97*T - 0.166*T^2)*10^-1;
% Other BH values inputted in Mathematica Script
% BH = (2668 + 10.867*T - 0.3111*T^2 + 1.778*10e-3*T^3)10^-1
% BH = (2670.8 + 19.9*T - 0.26*T^2)*10^-1
% BH = 2689.81 + 20.233*T^2 + 1.38*10^-3*T^3
CH = 0.3150*voH;

% Volume of Cylinder
V1 = pi*L1*(a1 / 2)^2;
V2 = pi*L2*(a2 / 2)^2;
V = V1 +V2;

%Symbolic variable P defined
syms P;

% Delta_a and delta_V
delta_a1 = ( ( (P - Po)*a1*(1 - v^2) ) / Ey)*( ( (b1^2 + a1^2) / (b1^2 - a1^2) ) + (v / (1 - v) ) );
delta_a2 = ( ( (P - Po)*a2*(1 - v^2) ) / Ey)*( ( (b2^2 + a2^2) / (b2^2 - a2^2) ) + (v / (1 - v) ) );
delta_V1 = (pi / 4)*(L1*(2*a1 + delta_a1)*delta_a1 + 2*L2*(2*a2 + delta_a2)*delta_a2);

vP = voP - (2.68e-04*log10( 1 + ( (P-Po) / 108.91) ) );
v1H = voH;
vH = voH - (CH*log10(1 + ( (P - Po) / BH) ) );
VA = (V1A*Po) / P;
delta_V2 = ( mp*(vP - v1P) ) + ( mH*(vH - voH) ) + (VA - V1A);

P2 = vpasolve(delta_V1 - delta_V2 == 0, P); % P2 is the 1st instance where delta_V1 - delta_V2 == 0

P = double(P2); % Convert to a numerical value with precision
delta_a1 = ( ( (P - Po)*a1*(1 - v^2) ) / Ey)*( ( (b1^2 + a1^2) / (b1^2 - a1^2) ) + (v / (1 - v) ) );
delta_V1 = (pi / 4)*(L1*( ( (2*a1) + delta_a1)*delta_a1) );
vP = voP - (2.68e-04*log10( 1 + ( (P-Po) / 108.91) ) );
v1H = voH;
vH = voH - (CH*log10(1 + ( (P - Po) / BH) ) );
VA = (V1A*Po) / P;
delta_V2 = ( mp*(vP - v1P) ) + ( mH*(vH - voH) ) + (VA - V1A);

Pa = (P2 / V1N)*(delta_V2 + V1N - V1A*( (Po / P2) - 1) ...
   - (V*f / v1P)*(voP - CP*log10(1 + ((P2 - Po) / BP) ) - v1P) ...
   + ((V*(1 - f) - V1A) / v1H)*CH*log10(1 + ((P2 - Po) / BH) ));

%% Outputs
V
Pa
P2