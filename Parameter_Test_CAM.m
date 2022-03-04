
clc; clear; close all

% PCM Material = Hexadecane (C16H34)
% Working fluid = Water
% Cylinder Material = TA2M (Titanium alloy)

% Environmental
T = 29; % Working Temperature [degC]
Thigh = T; % High Temperature
Tlow = 5; % Low Temperature [degC]
Po = 0.101; % Initial Pressure/ambient pressure = atmospheric pressure [MPa]

% Geometry & Material Properties of Cylinder
L1 = 1.31; % Length of the cylinder [m]
a1 = 7.7e-02; % Internal diameter of the cylinder [m]
b1 = 8.6e-02; % External diameter of the cylinder [m]

Ey = 105e03; % Young's Modulus [MPa]
v = 0.33; % Poisson's ratio

% Material Properties

% -> Material Properties of PCM
rhoS = 864; % Density of PCM - Solid phase [kg/m3]
rhoL = 773; % Density of PCM - Liquid phase [kg/m3]
mPCM = 3.456; % Mass of PCM [kg]
Tm = 18.2; % Melting Temperature of PCM [degC]

% -> Heat Transfer Properties
csd = 1.64; % Specific heat in solid state [kJ/kg K]
cld = 2.09; % Specific heat in liquid state [kJ/kg K]
Lh = 236; % Latent heat of fusion [kJ/kg]


% -> Material Properties of Hydraulic Fluid
voH = 1/1000; % Specific volume of hydraulic fluid [m3/kg]
v1H = voH; % Specific volume of hydraulic fluid at State 1 [m3/kg]

% Accumulator Parameters
V1N = 2e-03; % Initial Volume of N2 gas

% Misc.
ar = 6.573 / 100; % Volume fraction of residual air
%% Changing mPCM
mPCM = 3:0.1:5; % Mass of PCM [kg]
for i = 1:numel(mPCM)
    % 1. Finding Specific Volume of PCM [vP]
    voP = ( (1.0307e03 - ( 1.2596*(T + 273.15) )  + ...
        (1.8186e-3* (T + 273.15)^2) -(1.9555e-6* (T + 273.15)^3) ) )^-1;
    CP = 2.66e-04;
    BP = 102.12;
    v1P = 1 / rhoS; % Specific volume of PCM in liquid(solid??) state
    VPCM(i) = mPCM(i)*v1P;
    
    % 2. Finding Specific Volume of Hydraulic Fluid [vH]
    BH = ( 2672.9 + (15.97*T) - (0.166*(T^2) ) )*1e-01;
    CH = 0.3150*voH;
    
    % 3. Finding Volume Of Cylinder
    V = pi*L1*( (a1 / 2)^2 ); % Inner volume of Cylinder [m^3]
    V1A = ar*V; % Volume of residual air
    rPCM(i) = VPCM(i) / V;
    VH1 = ( V*(1 - rPCM(i)) ) - V1A;
    mH = ( (V*(1 - rPCM(i)) ) - V1A) / v1H;
    f(i) = rPCM(i);
    
    % 4. Finding change of Inner Diameter of Tube Under Internal Pressure
    syms P;
    delta_a1 = ( ( (P - Po)*a1*(1 - v^2) ) / Ey)*( ( (b1^2 + a1^2) / (b1^2 - a1^2) ) + (v / (1 - v) ) );
    delta_V1 = (pi / 4)*(L1*( ( (2*a1) + delta_a1)*delta_a1) );
    vP = 1.3e-03 - (2.66e-04*log10( 1 + ( (P - Po) / 102.12) ) );
    vH = voH - (CH*log10(1 + ( (P - Po) / BH) ) );
    VA = (V1A*Po) / P;
    delta_V2 = ( mPCM(i)*(vP - v1P) ) + ( mH*(vH - voH) ) + (VA - V1A);
    
    % 5. Finding Max Pressure [P2]
    P2(i) = vpasolve(delta_V1 - delta_V2 == 0, P); % P2 is the 1st instance where delta_V1 - delta_V2 == 0
%     P2(i) = double(P2(i)); % Convert to a numerical value with precision
    
    % 6. Finding Efficiency [Eff] 
    delta_a1 = ( ( (P2(i) - Po)*a1*(1 - v^2) ) / Ey)*( ( (b1^2 + a1^2) / (b1^2 - a1^2) ) + (v / (1 - v) ) );
    delta_V1 = (pi / 4)*(L1*( ( (2*a1) + delta_a1)*delta_a1) );
    Pa(i) = (P2(i) / V1N)*(delta_V1 + V1N - V1A*( (Po / P2(i)) - 1) ...
       - (V*f(i) / v1P)*(voP - CP*log10(1 + ((P2(i) - Po) / BP) ) - v1P) + ((V*(1 - f(i)) - V1A) / v1H)*CH* ...
       log10(1 + ((P2(i) - Po) / BH) ));
    
    Qin(i) = mPCM(i)*csd*(Tm - Tlow) + mPCM(i)*Lh + mPCM(i)*cld*(Thigh - Tm);
    Est(i) = -Pa(i)*1e6*V1N*log(1 - (mPCM(i) / V1N)*((1 / rhoL) - (1 / rhoS)) );
    Eff(i) = Est(i) / (Qin(i)*1e3) * 100;
end
% Turns from sym to double array
P2 = double(P2);
Eff = double(Eff);
Pa = double(Pa);
Qin = double(Qin);

%Create Text file Change_mPCM with outputs Eff, P2, Pa, VPCM, rPCM, Qin
writematrix([Eff;P2;Pa;VPCM; rPCM; Qin]','Change_mPCM.txt','Delimiter','tab')
mPCM = 3.456; % Mass of PCM [kg]

%% Changing rhoL
rhoL = 600:25:900; % Mass of PCM [kg]
for i = 1:numel(rhoL)
    % 1. Finding Specific Volume of PCM [vP]
    voP = ( (1.0307e03 - ( 1.2596*(T + 273.15) )  + ...
        (1.8186e-3* (T + 273.15)^2) -(1.9555e-6* (T + 273.15)^3) ) )^-1;
    CP = 2.66e-04;
    BP = 102.12;
    v1P = 1 / rhoS; % Specific volume of PCM in liquid(solid??) state
    VPCM = mPCM*v1P;
    
    % 2. Finding Specific Volume of Hydraulic Fluid [vH]
    BH = ( 2672.9 + (15.97*T) - (0.166*(T^2) ) )*1e-01;
    CH = 0.3150*voH;
    
    % 3. Finding Volume Of Cylinder
    V = pi*L1*( (a1 / 2)^2 ); % Inner volume of Cylinder [m^3]
    V1A = ar*V; % Volume of residual air
    rPCM = VPCM / V;
    VH1 = ( V*(1 - rPCM) ) - V1A;
    mH = ( (V*(1 - rPCM) ) - V1A) / v1H;
    f = rPCM;
    
    % 4. Finding change of Inner Diameter of Tube Under Internal Pressure
    syms P;
    delta_a1 = ( ( (P - Po)*a1*(1 - v^2) ) / Ey)*( ( (b1^2 + a1^2) / (b1^2 - a1^2) ) + (v / (1 - v) ) );
    delta_V1 = (pi / 4)*(L1*( ( (2*a1) + delta_a1)*delta_a1) );
    vP = 1.3e-03 - (2.66e-04*log10( 1 + ( (P - Po) / 102.12) ) );
    vH = voH - (CH*log10(1 + ( (P - Po) / BH) ) );
    VA = (V1A*Po) / P;
    delta_V2 = ( mPCM*(vP - v1P) ) + ( mH*(vH - voH) ) + (VA - V1A);
    
    % 5. Finding Max Pressure [P2]
    P2(i) = vpasolve(delta_V1 - delta_V2 == 0, P); % P2 is the 1st instance where delta_V1 - delta_V2 == 0
%     P2(i) = double(P2(i)); % Convert to a numerical value with precision
    
    % 6. Finding Efficiency [Eff] 
    delta_a1 = ( ( (P2(i) - Po)*a1*(1 - v^2) ) / Ey)*( ( (b1^2 + a1^2) / (b1^2 - a1^2) ) + (v / (1 - v) ) );
    delta_V1 = (pi / 4)*(L1*( ( (2*a1) + delta_a1)*delta_a1) );
    Pa(i) = (P2(i) / V1N)*(delta_V1 + V1N - V1A*( (Po / P2(i)) - 1) ...
       - (V*f / v1P)*(voP - CP*log10(1 + ((P2(i) - Po) / BP) ) - v1P) + ((V*(1 - f) - V1A) / v1H)*CH* ...
       log10(1 + ((P2(i) - Po) / BH) ));
    
    Qin = mPCM*csd*(Tm - Tlow) + mPCM*Lh + mPCM*cld*(Thigh - Tm);
    Est(i) = -Pa(i)*1e6*V1N*log(1 - (mPCM / V1N)*((1 / rhoL(i)) - (1 / rhoS)) );
    Eff(i) = Est(i) / (Qin*1e3) * 100;
end
% Turns from sym to double array
P2 = double(P2);
Eff = double(Eff);
Pa = double(Pa);
Est = double(Est);
%Create Text file Change_rhoS with outputs Eff, P2, Pa, Est
writematrix([Eff;P2;Pa;Est]','Change_rhoS.txt','Delimiter','tab')
rhoL = 773; % Density of PCM - Liquid phase [kg/m3]

%% Changing rhoS
rhoS = 600:25:1000; % Mass of PCM [kg]
for i = 1:numel(rhoS)
    % 1. Finding Specific Volume of PCM [vP]
    voP = ( (1.0307e03 - ( 1.2596*(T + 273.15) )  + ...
        (1.8186e-3* (T + 273.15)^2) -(1.9555e-6* (T + 273.15)^3) ) )^-1;
    CP = 2.66e-04;
    BP = 102.12;
    v1P(i) = 1 / rhoS(i); % Specific volume of PCM in liquid(solid??) state
    VPCM(i) = mPCM*v1P(i);
    
    % 2. Finding Specific Volume of Hydraulic Fluid [vH]
    BH = ( 2672.9 + (15.97*T) - (0.166*(T^2) ) )*1e-01;
    CH = 0.3150*voH;
    
    % 3. Finding Volume Of Cylinder
    V = pi*L1*( (a1 / 2)^2 ); % Inner volume of Cylinder [m^3]
    V1A = ar*V; % Volume of residual air
    rPCM(i) = VPCM(i) / V;
    VH1 = ( V*(1 - rPCM(i)) ) - V1A;
    mH = ( (V*(1 - rPCM(i)) ) - V1A) / v1H;
    f = rPCM(i);
    
    % 4. Finding change of Inner Diameter of Tube Under Internal Pressure
    syms P;
    delta_a1 = ( ( (P - Po)*a1*(1 - v^2) ) / Ey)*( ( (b1^2 + a1^2) / (b1^2 - a1^2) ) + (v / (1 - v) ) );
    delta_V1 = (pi / 4)*(L1*( ( (2*a1) + delta_a1)*delta_a1) );
    vP = 1.3e-03 - (2.66e-04*log10( 1 + ( (P - Po) / 102.12) ) );
    vH = voH - (CH*log10(1 + ( (P - Po) / BH) ) );
    VA = (V1A*Po) / P;
    delta_V2 = ( mPCM*(vP - v1P(i)) ) + ( mH*(vH - voH) ) + (VA - V1A);
    
    % 5. Finding Max Pressure [P2]
    P2(i) = vpasolve(delta_V1 - delta_V2 == 0, P); % P2 is the 1st instance where delta_V1 - delta_V2 == 0
%     P2(i) = double(P2(i)); % Convert to a numerical value with precision
    
    % 6. Finding Efficiency [Eff] 
    delta_a1 = ( ( (P2(i) - Po)*a1*(1 - v^2) ) / Ey)*( ( (b1^2 + a1^2) / (b1^2 - a1^2) ) + (v / (1 - v) ) );
    delta_V1 = (pi / 4)*(L1*( ( (2*a1) + delta_a1)*delta_a1) );
    Pa(i) = (P2(i) / V1N)*(delta_V1 + V1N - V1A*( (Po / P2(i)) - 1) ...
       - (V*f / v1P(i))*(voP - CP*log10(1 + ((P2(i) - Po) / BP) ) - v1P(i)) + ((V*(1 - f) - V1A) / v1H)*CH* ...
       log10(1 + ((P2(i) - Po) / BH) ));
    
    Qin = mPCM*csd*(Tm - Tlow) + mPCM*Lh + mPCM*cld*(Thigh - Tm);
    Est(i) = -Pa(i)*1e6*V1N*log(1 - (mPCM / V1N)*((1 / rhoL) - (1 / rhoS(i))) );
    Eff(i) = Est(i) / (Qin*1e3) * 100;
end

% Turns from sym to double array
P2 = double(P2);
Eff = double(Eff);
Pa = double(Pa);
v1P = double(v1P);
% Create Text file Change_mPCM with outputs Eff, P2, Pa, v1P
writematrix([Eff;P2;Pa;v1P],'Change_rhoS.txt','Delimiter','tab')
rhoS = 864; % Density of PCM - Solid phase [kg/m3]