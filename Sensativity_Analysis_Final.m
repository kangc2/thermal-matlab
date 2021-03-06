% Sensativity Analysis of MATLAB Model of our PCM engine
clc; close all; clear

%% Define Intial Struct 'engine': has all the baseline parameters of our intial engine
% Thermal Engine from our Thermal_Engine code
% The baseline of our model
% Note: make sure the field name is right (voH vs. VoH)
engine.name = 'Original';
engine.working_fluid = 'Water';
engine.hull_material = 'TA2M (Titanium alloy)'; 
% Variables of the environment, not specific to the engine
T = 29; % Working Temperature [degC]
engine.T = T;
engine.Thigh = T; % High Temperature
engine.Tlow = 5; % Low Temperature [degC]
engine.Po = 0.101; % Initial Pressure/ambient pressure = atmospheric pressure [MPa]

% -> Geometry & Material Properties of Cylinder
engine.L1 = 1.3; % Length of the cylinder [m]
engine.b1 = 0.15; % External diameter of the cylinder [m]
engine.t = 0.01; %wall thickness [m]
engine.a1 = engine.b1 - 2*engine.t; % Internal diameter of the cylinder [m]

engine.Ey = 105e03; % Young's Modulus [MPa]
engine.v = 0.33; % Poisson's ratio
engine.yield = 340; % Compressive/Tensile Yield Stress (from google) / Factor of Safety [MPa]

% -> Material Properties of PCM
engine.delta_rho = 70;
engine.rhoL = 773; % Density of PCM - Liquid phase [kg/m3]
engine.rhoS = engine.rhoL + engine.delta_rho; %TTTTHIIISSSSSSSSS

engine.f = 0.6557; %Volume fraction of PCM
engine.mPCM = (engine.f*pi*engine.L1*( (engine.a1 / 2)^2 ))/(1/engine.rhoS); % Mass of PCM [kg]
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
% mL -> m^3 = 1e-6
engine.desiredWork = 6.096*248*1e-3; %Work needed in float = P*delta_V [MPa*m^3*1e3 = kJ]

%% Define Outputs in 'engine' using Functions
% 1 Find specific volume of PCM under ambient pressure Po, volume and mass of PCM [voP]
[engine.voP, engine.v1P, engine.VPCM] = specificVolPCM(engine);
% 3. Find inner volume of cylinder, Volume fraction of air and PCM 
[engine.V,engine.V1A, engine.V1H, engine.mH] = findVolume(engine);
% 5. Finding Max Pressure [P2]
engine.P2 = findPressure(engine);
% 4. Finding change of Inner Diameter of Tube Under Internal Pressure
engine.delta_V1 = findChangeInnerVolume(engine);
% 6a. Find Pre-charged pressure in acculimlator
engine.Pa = findPa(engine);
% 6b. find the total energy stored
engine.Est = findEst(engine);
% 6c. Find total thermal energy absorbed by PCM in 1 melting/freezing cycle
engine.Qin = findQin(engine);
% 6d. find the theorectical Efficiency % 
engine.Eff = findEfficiency(engine);  
engine.Work = engine.Est / 1000; % Est from [J] to [kJ]

%% Find Efficiency: get efficiency using all of our inputs
% adds a new struct field Eff2, same as Eff, but using function that
% takes in all the inputs at once
engine.Eff2 = findEfficiency2(engine);

% Finds the stress of the pcm hull using the Pressure to Stress function
engine.stress = PtoStress(engine);
engine.FoS = engine.stress / (engine.yield); % Factor of Safety

engine % prints engine struct

%% Sensativity Analysis starts here: Loops thru every parameter and returns a structured data for each parameter
%create a new struct 'parameter' to add in all our of tests into structures

% vvv CAN EDIT SENSATIVITY ANALYSIS HERE vvv
% Add in what parameters we want to change and vary in a new struct
% 'input_range'
% variable name = 'input_range.VARIABLE_NAME' (i.e. 'L1', 'b1')
% Fields: [lowerbound value, upperbound value, original value, step size between bounds (delta)]
input_range.L1 = [0.2, 1.3, engine.L1, 21];
input_range.b1 = [0.12, 0.2, engine.b1, 21];
input_range.Lh = [210, 280, engine.Lh, 21];
input_range.Tm = [5.5, 20, engine.Tm, 21];
input_range.rhoL = [720, 800, engine.rhoL, 21];
input_range.f = [0, 1, engine.f, 21];
% ^^^ can comment out or add in any new variables to check thru here ^^^

% Get all the fields in struct to prepare for the for loop (sensativity
% analysis)
fields = string(fieldnames(input_range));

original = engine; % defining a new engine to match the original 'engine' structure
% for reconverting the values in the loop back to its original values when the second for loop finishes


% This For loop goes thru every variable we want to change in the range
for j = 1:length(fields) % loops through each parameter we want to change
    % Defines which parameter we are at
    input = input_range.(fields(j)); % get the information from each field in the structure
    l = input(1); % lower bound value
    u = input(2); % upper bound value
    delta = input(4); % step size of each value
    
    j % prints which number of the loop we are at

    newParam = []; % this is all the changed values for each variable goes
    newEff = []; % this is where each efficiency value is for each variable
   
    % This for loop goes through each new value of that parameter and
    % stores all the data into a new structure called 'parameter'
    for i= 0:delta
       % Find what the value is for the changing variable
        z = ((u - l) / delta) * i;
        
        original.(fields(j)) = z + l; % Changed value for a variable
        
        original.a1 = original.b1 - 2*original.t; %Updating thickness/interior diam.
        original.rhoS = original.rhoL + original.delta_rho;
        original.mPCM = (original.f*pi*original.L1*( (original.a1 / 2)^2 ))/(1/original.rhoS);
      
        % Find Efficiency value
        % each changed value is represented by p(#) instead of engine.parameter
        answer = findEfficiency2(original);
        % add all changed values into a row
        newParam = [newParam, original.(fields(j))];
        % add all the efficiency values into one row
        newEff = [newEff, answer];
    end
    
   % Now add into a new  struct called 'parameter'
   % New structure: parameter
   % fields: name, parameter, efficiency
   parameter(j).name = fields(j);
   parameter(j).param = newParam;
   parameter(j).efficiencies = newEff;
   
   original = engine; % to reconvert all the parameter changes back to its original values when the second for loop finishes
end

% COMMENTED OUT: Prints out all the data into an excel file
% writetable(struct2table(parameter), 'someexcelfile.xlsx')


%% This prints out plots of all the parameters looped through above
for i = 1:length(parameter)
    figure
    plot(parameter(i).param, parameter(i).efficiencies)
    xlabel(parameter(i).name)
    ylabel('Efficiency [%]')
end

%% Functions
% 1. Finding Specific Volume of PCM [vP]
function [voP, v1P, VPCM] = specificVolPCM(engine)
    % a. Find specific volume of PCM under ambient pressure Po [voP]
    voP = ( (1.0307e03 - ( 1.2596*(engine.T + 273.15) )  + ...
        (1.8186e-3* (engine.T + 273.15)^2) -(1.9555e-6* (engine.T + 273.15)^3) ) )^-1;
    % b. Find the volume of PCM (liquid state) [v1P]
    v1P = 1/engine.rhoS; % Specific volume of PCM in liquid state
    VPCM = engine.mPCM*v1P;
end

% 3. Find inner volume of cylinder, Volume fraction of air and PCM 
function [V, V1A, V1H, mH] = findVolume(engine)
    V = pi*engine.L1*( (engine.a1 / 2)^2 ); % Inner volume of Cylinder

    V1A = engine.ar*V; % Volume of residual air
    V1H = ( V*(1 - engine.f) ) - V1A; % volume of HF at state 1
    mH = ( (V*(1 - engine.f) ) - V1A) / engine.v1H; % mass of HF
end

% 5. Finding Max Pressure [P2]
function P2 = findPressure(engine)
    F = @(P)((pi / 4)*(engine.L1*( ( (2*engine.a1) + ...
            ( ( (P - engine.Po)*engine.a1*(1 - engine.v^2) ) / engine.Ey)*( ( (engine.b1^2 + engine.a1^2) / (engine.b1^2 - engine.a1^2) ) +...
        (engine.v / (1 - engine.v) ) ))*...
            ( ( (P - engine.Po)*engine.a1*(1 - engine.v^2) ) / engine.Ey)*( ( (engine.b1^2 + engine.a1^2) / (engine.b1^2 - engine.a1^2) ) + ...
        (engine.v / (1 - engine.v) ) )) ))... % delta_V1
    - ...
    (( engine.mPCM*((1.3e-03 - (2.66e-04*log10( 1 + ( (P - engine.Po) / 102.12) ) )) - engine.v1P) ) + ...
    ( engine.mH*((engine.voH - (engine.CH*log10(1 + ( (P - engine.Po) / engine.BH) ) )) - engine.voH) ) + (((engine.V1A*engine.Po) / P) - engine.V1A)); %delta_V2

    % Options: sets tolerance of function close to 0 (1-e14) and displays the
    % iteration, this could help with the optimization
    options = optimoptions('fsolve','TolFun',1e-14);
    
    %solves for P, same answer as the engine.P2 with the current finding
    %Pressure function
    P2 = fsolve(F,5,options);

end

% 4. Finding change of Inner Diameter of Tube Under Internal Pressure
function delta_V1 = findChangeInnerVolume(engine)
    delta_a1 = ( ( (engine.P2 - engine.Po)*engine.a1*(1 - engine.v^2) ) / engine.Ey)*( ( (engine.b1^2 + engine.a1^2) / ...
        (engine.b1^2 - engine.a1^2) ) + (engine.v / (1 - engine.v) ) );
    delta_V1 = (pi / 4)*(engine.L1*( ( (2*engine.a1) + delta_a1)*delta_a1) );
end

% 6. Finding Efficiency [Eff] 
% a. Find Pre-charged pressure in acculimlator [MPa]
function Pa = findPa(engine)
  Pa = (engine.P2 / engine.V1N)*(engine.delta_V1 + engine.V1N - engine.V1A*( (engine.Po / engine.P2) - 1) ...
       - (engine.V*engine.f / engine.v1P)*(engine.voP - engine.CP*log10(1 + ((engine.P2 - engine.Po) / engine.BP) ) - engine.v1P) + ...
       ((engine.V*(1 - engine.f) - engine.V1A) / engine.v1H)*engine.CH*log10(1 + ((engine.P2 - engine.Po) / engine.BH) ));  
end

% b. find the total energy stored [Est]
function Est = findEst(engine)
    Est = -engine.Pa*1e6*engine.V1N*log(1 - (engine.mPCM / engine.V1N)*((1 / engine.rhoL) - (1 /engine. rhoS)) );
end

% c. Find total thermal energy absorbed by PCM in 1 melting/freezing cycle
% [Qin]
function Qin = findQin(engine)
    Qin = engine.mPCM*engine.csd*(engine.Tm - engine.Tlow) + engine.mPCM*engine.Lh + engine.mPCM*engine.cld*(engine.Thigh - engine.Tm);
end

% d. find the theorectical Efficiency % [Eff] 
function Eff = findEfficiency(engine)  
    Eff = engine.Est / (engine.Qin*1e3) * 100;
end

% This takes in all inputs and finds the Efficiency using all functions
% above
function Eff = findEfficiency2(engine)  
    [engine.voP, engine.v1P, engine.VPCM] = specificVolPCM(engine);
    [engine.V,engine.V1A, engine.V1H, engine.mH] = findVolume(engine);
    engine.P2 = findPressure(engine);
    engine.delta_V1 = findChangeInnerVolume(engine);
    engine.Pa = findPa(engine);
    engine.Est = findEst(engine);
    engine.Qin = findQin(engine);
    Eff = findEfficiency(engine);
end

%Pressure to Stress Equations
function stress_vm = PtoStress(engine)
    % sigma_tan is the tangential (axial) stress
    stress_tan = engine.P2*(((engine.b1/2)^2 + (engine.a1/2)^2) / ((engine.b1/2)^2 - (engine.a1/2)^2));
    
    %sigma_rad is the radial stress
    stress_rad = -engine.P2;
    
    %sigma_long is the longitudinal (hoop) stress on the ends
    stress_long = (engine.P2*(engine.a1/2)^2)/((engine.b1/2)^2 - (engine.a1/2)^2);

    % von Mises Stress
    stress_vm = sqrt(((stress_tan - stress_rad)^2 + (stress_rad - stress_long)^2 ...
    + (stress_long - stress_tan)^2) / 2);
end

