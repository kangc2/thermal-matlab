clear all
close all
clc
% Names of the Y axis ticks
names={'ABS';
'Adhesive';
'Copper';
'Glass fibre';
'Melamine';
'Paper';
'Polyethylene';
'Polypropylene';
'Polystyrene';
'Polyvinylchloride';
'Steel';
'Printed Wiring Board';
'Electricity (coal)';
'Electricity (bio)';
'Truck Transport';
'Ocean Transport';
'Train Transport'};
% Import data here
CO2=rand([length(names) 1]);
% Sensitivity percentage to calculate with
% 0.2 corresponds to a plus or minus 20% 
% sensitivity analysis
sensitivity=0.2;
% Initialize low and high matricies for speed
CO2_low=zeros(1,length(names));
CO2_high=zeros(1,length(names));
CO2_low_sum=zeros(1,length(names));
CO2_high_sum=zeros(1,length(names));
% Calculate the base change due to a single factor at a time
for i=1:length(names)
    CO2_low=CO2;
    CO2_high=CO2;
    CO2_low(i)=CO2_low(i)*(1-sensitivity);
    CO2_high(i)=CO2_high(i)*(1+sensitivity);
    CO2_low_sum(i)=sum(CO2_low);
    CO2_high_sum(i)=sum(CO2_high);
end
% The base value is where the y axis is centered
CO2_base_value=sum(CO2);
% Sort the values based on the lower change
% Sort the higher values and the names arrays
%    using the same indices
[CO2_low_sum ind]=sort(CO2_low_sum,'descend');
CO2_high_sum=CO2_high_sum(ind);
names_CO2=names(ind);
% Create a figure and plot the low and high horizontally
figure
h = barh(CO2_high_sum);
hold on
barh(CO2_low_sum,'r')
bh = get(h,'BaseLine');
set(bh,'BaseValue',CO2_base_value);
title('Sensitivities')
set(gca,'yticklabel',names)
set(gca,'Ytick',[1:length(names)],'YTickLabel',[1:length(names)])
set(gca,'yticklabel',names_CO2)
xlabel('CO_2 (kg/dry)')
saveas(gcf,'CO2.png')
