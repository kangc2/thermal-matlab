clear all
close all
clc
% Names of the Y axis ticks
names={'Length (L1)';
'Thickness (b1)';
'Latent Heat (Lh)';
'Melting Temp (Tm)';
'Solid Density PCM';
'Liquid Density PCM';
'PCM Mass'};

% Import data here
CO2=rand([length(names) 1]);

% Sensitivity percentage to calculate with 0.2 corresponds to 
% a plus or minus 20% sensitivity analysis
sensitivity=0.2;

% Initialize low and high matricies for speed
low=zeros(1,length(names));
high=zeros(1,length(names));
low_sum=zeros(1,length(names));
high_sum=zeros(1,length(names));

% Calculate the base change due to a single factor at a time
for i=1:length(names)
    low=CO2;
    high=CO2;
    low(i)=low(i)*(1-sensitivity);
    high(i)=high(i)*(1+sensitivity);
    low_sum(i)=sum(low);
    high_sum(i)=sum(high);
end

% The base value is where the y axis is centered
base_value=sum(CO2);

% Sort the values based on the lower change
% Sort the higher values and the names arrays
%    using the same indices
[low_sum ind]=sort(low_sum,'descend');
high_sum=high_sum(ind);
names_plot=names(ind);

% Create a figure and plot the low and high horizontally
figure
h = barh(high_sum);
hold on
barh(low_sum,'r')
bh = get(h,'BaseLine');
set(bh,'BaseValue',base_value);
title('Sensitivities')
set(gca,'yticklabel',names)
set(gca,'Ytick',[1:length(names)],'YTickLabel',[1:length(names)])
set(gca,'yticklabel',names_plot)
xlabel('CO_2 (kg/dry)')
saveas(gcf,'CO2.png')
