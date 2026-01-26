%Authors: Alex Godbout
%Date: 1/25/26
%Purpose: Complete Lab 1 Part 2 Task 1
%Inputs: Collected Data from Truss Demo Including Force and Displacment
%       Measurments
%Outputs: Plots of Data along with linear Regression Error Analysis

clear; clc; close all

%Load in Data
Case1 = ReadData("Data/Case 1 data");
Case2 = ReadData("Data/Case 2 data");
Case3 = ReadData("Data/Case 3 data");
Data = {Case1,Case2,Case3};
%Plots of Data
plotData(2001,"F0",Data, [0,65,0,35],{"N","Kg"})
plotData(2002,"F1",Data, [0,65,0,35],{"N","Kg"})
plotData(2003,"F2",Data, [0,65,0,35],{"N","Kg"})
plotData(2004,"F3D",Data, [0,65,-90,-10],{"N","Kg"})
plotData(2005,"LVDT",Data, [0,65,0,.07],{"m","Kg"})

RSquared = calculateUncertainty(Data)


function plotData(figNum, dataElement, data, axisLimits,units)
    figure(figNum);
    hold on;
    colors = {"r","g","b"};
    labels = {"Load in Center","Single Load (Unknown Location)","Two Loads (Symmetic Placement)"};
    axis(axisLimits)
    title(sprintf("%s (%s) vs Load Magnitude (%s)", dataElement,units{1},units{2}))
    xlabel(sprintf("Load Maginitude (%s)",units{2}))
    ylabel(sprintf("%s (%s)", dataElement,units{1}))
    grid on;
    x = linspace(0,65,10);
    for i = 1:numel(data)
        regression = polyfit(data{i}.LoadingCase(:),data{i}.(dataElement)(:),1);
        y = regression(2) +  x*regression(1);
        plot(x,y,'Color',colors{i}, 'HandleVisibility','off')
        scatter(data{i}.LoadingCase(:),data{i}.(dataElement)(:),36, colors{i}, 'filled',DisplayName = labels{i});
    end
    legend('Location', 'northwest')
end

function RSquared = calculateUncertainty(data)
    RSquared = struct();
    
    fieldNames = fieldnames(data{1});
    dataFields = fieldNames(2:5);

    for f = 1:numel(dataFields)
        fieldName = dataFields{f};
        RSquared.(fieldName) = struct();
    
        for i = 1:numel(data)
    
            x = data{i}.LoadingCase(:);
            y = data{i}.(fieldName)(:);
    
            Stats = fitlm(x,y);
            caseName = sprintf("Case%d", i);
            RSquared.(fieldName).(caseName) = Stats.Rsquared.Adjusted;
        end
    end
end
