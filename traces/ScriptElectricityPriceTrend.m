% Reset
clc; clear; close all;    
% Add class paths
addpath('../functions');

%% Initialize variables.
filename = 'C:\Users\NhatTan\Dropbox\TCO\code\traces\electricity_prices\price_0029.76_-095.36_tx.houston.csv';
delimiter = ',';

%% Format string for each line of text:
%   column1: text (%s)
%	column2: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
datetimes = dataArray{:, 1};
rawPrices = dataArray{:, 2};

%% cut off the peak prices
xData = 1:24*365*3;
prices= rawPrices(xData);
prices = min(prices, 250);
prices = max(prices, 0);
priceLen = length(prices);


figure(1);
isFit = true;
plot(xData,prices);
if isFit 
    [fitresult, gof] = createPolyFit(1:priceLen, prices);
    hold on;
    plot(xData, fitresult(xData), '-k', 'linewidth', 1.5);    
end
%% monthly price
oneMonth = 30*24;
monthLen =  floor(priceLen/oneMonth);
tempPrices = reshape(prices(1:oneMonth* monthLen), oneMonth, monthLen);
monthlyPrices = mean(tempPrices);
figure(2);
plot(monthlyPrices);
xData = 1:length(monthlyPrices);
if isFit 
    [fitresult, gof] = createPolyFit(1:monthLen, monthlyPrices);
    hold on;
    plot(xData, fitresult(xData), '-k', 'linewidth', 1.5);    
end
%% Average daily prices
oneDay = 24;
dayLen =  floor(priceLen/oneDay);
tempPrices = reshape(prices(1:oneDay* dayLen), oneDay, dayLen);
dailyPrices = mean(tempPrices);
figure(3);
plot(dailyPrices);
xData = 1:length(dailyPrices);
if isFit 
    [fitresult, gof] = createPolyFit(1:dayLen, dailyPrices);
    hold on;
    plot(xData, fitresult(xData), '-k', 'linewidth', 1.5);    
end
%% Average weekly prices
oneWeek = 24*7;
weekLen =  floor(priceLen/oneWeek);
tempPrices = reshape(prices(1:oneWeek* weekLen), oneWeek, weekLen);
weeklyPrices = mean(tempPrices);
figure(4);
plot(weeklyPrices);
xData = 1:length(weeklyPrices);
if isFit 
    [fitresult, gof] = createPolyFit(1:weekLen, weeklyPrices);
    hold on;
    plot(xData, fitresult(xData), '-k', 'linewidth', 1.5);    
end