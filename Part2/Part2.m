%Authors: Alex Godbout
%Date: 1/25/26
%Purpose: Complete Lab 1 Part 2 Task 1
%Inputs: Collected Data from Truss Demo Including Force and Displacment
%       Measurments
%Outputs: Plots of Data along with linear Regression Error Analysis

clear; clc; close all

%Constants
R = 3/16; % inches
r = 2/16; % inches
L =  4*3.28084*12; %inches
E = 10*10^6; %psi
c = .125*3.28084*12; %Inches
I = pi/4*(R^4-r^4); % inches^4
I_x = 4*I+4*(R^2-r^2)*pi*c^2; %I aint calculating these units
A = (R^2-r^2)*pi; %Inches^2

%Load in Data
Case1 = ReadData("Data/Case 1 data");
Case2 = ReadData("Data/Case 2 data");
Case3 = ReadData("Data/Case 3 data");
Data = {Case1,Case2,Case3};
%Plots of Data
plotData(2001,"F0",Data, [0,65,0,35],{"lb","lb"})
plotData(2002,"F1",Data, [0,65,0,35],{"lb","lb"})
plotData(2003,"F2",Data, [0,65,0,35],{"lb","lb"})
plotData(2004,"F3D",Data, [0,65,-90,-10],{"lb","lb"})
plotData(2005,"LVDT",Data, [0,65,0,.07],{"in","lb"})
RSquared = calculateUncertainty(Data)'

%Calculate Expected Displacement
P = [0,10,20,30,40,50,40,30,20,10,0]; %Load values to test
[a,nu,F_i] = SingleLoadNonCentered(Case2,P,L,E,I_x,A,c)

figure(2006); hold on; grid on;
scatter(Case2.LoadingCase,Case2.LVDT,'filled');
plot(P,nu,LineStyle='--');
title("Comparison of Measured Data and Model Prediction")
xlabel("Magnitude of Load P")
ylabel("Dispalcement at center of beam")
legend("Measured Deflection", "Modeled Beam Deflection")

%% Task 2
%a -> location of load
%x -> location of interest
function [a,nu,F_i] = SingleLoadNonCentered(data,P,L,E,I_x,A,c)
    
    % Measured reactions
    F_A_meas = data.F0 + data.F1;
    F_B_meas = data.F2;

    % Infer load location
    P_meas = F_A_meas + F_B_meas;
    a=(L .* F_B_meas ./ P_meas);
    a = mean(a(11:end));

    % Theoretical reaction
    R_A = P.*(L-a)/L;

    % Integration constant
    c1 = (P.*(L-a).^3 - R_A.*L.^3)/(6*L);

    % Midspan deflection
    x = L/2;

    if x <= a
        nu = -(R_A.*x.^3/6 + c1.*x)/(E*I_x);
    else
        nu = -(R_A.*x.^3/6 - P.*(x-a).^3/6 + c1.*x)/(E*I_x);
    end
w
    F_i = 0;
end

%% Functions
plotData(2001,"F0",Data, [0,65,0,35],{"N","Kg"})
plotData(2002,"F1",Data, [0,65,0,35],{"N","Kg"})
plotData(2003,"F2",Data, [0,65,0,35],{"N","Kg"})
plotData(2004,"F3D",Data, [0,65,-90,-10],{"N","Kg"})
plotData(2005,"LVDT",Data, [0,65,0,.07],{"m","Kg"})

RSquared = calculateUncertainty(Data);
[Stress_part3,Deflection_part3] = Part3Model(L,E,I,A,c,Case3);


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


function [Stress_Part3, Deflection_Part3] = Part3Model(L,E,I_x,A,c,Case3)
L_in = (L*0.250)*39.37; %Feet to inches here
F_ly = Case3.LoadingCase./2;
%C_y = Case3.F2;

%a = ((F_ry - C_y.*L))./((-F_ly + F_ry)); %This doesn't work and its
%important, because it shows we can't find a just with reaction forces,
%because things will cancel out.
F3D_initial = 86;
Delta_Fi = Case3.F3D+F3D_initial;
a = (Delta_Fi*I_x)./(A*c.*F_ly);
a(~isfinite(a)) = 0;
Deflection_Part3 = ((F_ly .* a)./(6*E*I_x)).*((3*L_in^2)/4 - a.^2);

Stress_Part3 = (F_ly.*a*c)/I_x;
end