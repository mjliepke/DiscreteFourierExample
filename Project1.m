%Matthew Liepke, MA 441H Project 1
%This script is aimed to read data from a file, and make a cos/full Fourier
%series from the data

%Read data & define constants
filename = 'rocketProbe-1.txt'; %this works for simple delimited data
xIndex = 1; yIndex = 2;

data = rot90(importdata(filename));

%If the x-data is in MM-DD_YYYY format do the following:
%data(xIndex,:) = linspace(0,length(data)-1,length(data));

halfPeriodCos = data(xIndex,length(data(1,:)))-data(xIndex,1); %assuming this data is the length of the half-period
halfPeriodFull = halfPeriodCos/2;
order = 1000;

%find the coefficients for the series
[cosCoeff,fullSinCoeff,fullCosCoeff] = deal(zeros(1,order));

cosCoeff(1) = A_o(data(xIndex,:),data(yIndex,:), halfPeriodCos);
fullCosCoeff(1) = cosCoeff(1);
fullSinCoeff(1) = 0;

for i=1:order-1
    cosCoeff(i+1) = A_n_even(data(xIndex,:),data(yIndex,:),halfPeriodCos,i);
    fullCosCoeff(i+1) = A_n(data(xIndex,:),data(yIndex,:),halfPeriodFull,i);
    fullSinCoeff(i+1) = B_n(data(xIndex,:),data(yIndex,:),halfPeriodFull,i);
end

%Make a representation of the calculated series to plot (numericially)
cosSeries = zeros(size(data(yIndex,:)));
sinSeries = zeros(size(data(yIndex,:)));
fullSeries = zeros(size(data(yIndex,:)));
for i=1:order
    cosSeries = cosSeries(1,:) + cosCoeff(i)*cos((i-1)*pi*data(xIndex,:)/halfPeriodCos);
    fullSeries = fullSeries(1,:) + fullCosCoeff(i)*cos((i-1)*pi*data(xIndex,:)/halfPeriodFull)...
        +fullSinCoeff(i)*sin((i-1)*pi*data(xIndex,:)/halfPeriodFull);
end

%% %Plot it
tiledlayout(3,1);nexttile;
plot(data(xIndex,:),data(yIndex,:));
title('Original Data');

nexttile;
hold on
plot(data(xIndex,:),data(yIndex,:));
plot(data(xIndex,:),cosSeries(1,:));
title(append('Fourier Cosine Series, order ',num2str(order)));

nexttile;
hold on
plot(data(xIndex,:),data(yIndex,:));
plot(data(xIndex,:),fullSeries(1,:));
title(append('Fourier Full Series, order ',num2str(order)));

 
%% %Functions for the Fourier Series - All the math is done here
function A_o = A_o(dataX, dataY, halfPeriod)
    A_o = (1/halfPeriod)*trapz(dataX,dataY);
end

function A_n_even = A_n_even(dataX, dataY, halfPeriod, n)
    A_n_even = (2/(halfPeriod))*trapz(dataX,dataY.*cos(n*pi*dataX/halfPeriod));
end

function A_n = A_n(dataX, dataY, halfPeriod, n)
    A_n = (1/(halfPeriod))*trapz(dataX,dataY.*cos(n*pi*dataX/halfPeriod));
end

function B_n = B_n(dataX, dataY, halfPeriod, n)
    B_n = (1/halfPeriod)*trapz(dataX,dataY.*sin(n*pi*dataX/halfPeriod));
end
