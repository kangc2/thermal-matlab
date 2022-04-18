clear; close all; clc
%%
% Set intial guess values for box dimensions
lengthGuess = 1;
widthGuess = 1;
heightGuess = 1;

% Load Guess values into array
x0 = [lengthGuess, widthGuess, heightGuess];

% Call solver to minimize the objective function given the constraint
xopt = fmincon(@objective, x0, [], [], [], [], [], [], @constraint, [])
volumeOpt = calcVolume(xopt)
surfaceAreaOpt = calcSurface(xopt)
%% functions
% Define function to calculate volume of box
function volume = calcVolume(x)
    length = x(1);
    width = x(2);
    height = x(3);
    volume = length*width*height;
end

% Define function to calculate surface area of box
function surfaceArea = calcSurface(x)
    length = x(1);
    width = x(2);
    height = x(3);
    surfaceArea = 2*length*width + 2*length*height +...
        2*height*width;
end

% Define objective function for optimization
function obj = objective(x)
    obj = -calcVolume(x);
end

% Define constraint for optimization
function [c, ceq] = constraint(x)
    c = calcSurface(x) - 10;
    ceq = [];
end