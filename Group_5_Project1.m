function Group_5_Project1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project #1
% Payton Downey, Charles Lee, Kade Pizzuto, Nicholas Tate, Jaylin Trice, Ethan Winchester 
% ME 2543--Simulations Methods
% Spring 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; % Clear the variable list and the command window.
options = optimoptions('fsolve','Display','none'); % Removes fsolve notes from command window
% Problem Number: 1
% Define Given Parameters 
L1 = 10; L2 = 10; L3 = 10; % pipe length; ft
D1 = 2/12; D2 = 2/12; D3 = 2/12; % pipe diameter; converting in. to ft
e1 = 0.00085; e2 = 0.00085; e3 = 0.00085; % surface roughness; ft
rho = 1.94*32.174; % desnity; converting slug/ft^3 to lbf*s^2/ft^3
mu = 2.05e-5; % dynamic viscosity; lbf*s/ft^2
A1 = (pi*(D1)^2)/4; A2 = (pi*(D2)^2)/4; A3 = (pi*(D3)^2)/4; % cross-sectional area for all pipes; ft^2
Q = 750*0.002228; % Total Flow Rate; converting gpm to ft^3/s
% Inital Flow Rate and Chnage in Pressure guesses
Q1 = Q/4; Q2 = Q/2; Q3 = Q/4; dP = 50000; 
Q0 = [Q1;Q2;Q3;dP]; 
% Intial Friction Factor Guesses
f1 = 0.25*(log10(((e1/D1)/3.7)+((5.74)/(myRe(myV(Q1,A1),D1)^0.9))))^(-2); %Colebrook Equation
f2 = 0.25*(log10(((e2/D2)/3.7)+((5.74)/(myRe(myV(Q2,A2),D2)^0.9))))^(-2); %Colebrook Equation
f3 = 0.25*(log10(((e3/D3)/3.7)+((5.74)/(myRe(myV(Q3,A3),D3)^0.9))))^(-2); %Colebrook Equation
% Solve for Q1,Q2,Q3, and dP (Change in Pressure) 
z = fsolve(@(x) equations(x,Q),Q0,options); 
fprintf("Test 1: \n" ...
        + "The flow rate for pipe 1 (Q1) is: " + (z(1)/0.002228) + " gpm.\n" ...
        + "The flow rate for pipe 2 (Q2) is: " + (z(2)/0.002228) + " gpm.\n" ...
        + "The flow rate for pipe 3 (Q3) is: " + (z(3)/0.002228) + " gpm.\n" ...
        + "The total flow rate (Q) is: " + ((z(1)/0.002228)+(z(2)/0.002228)+(z(3)/0.002228)) + " gpm.\n" ...
        + "The total change in pressure (dP) is: " + (z(4)/144) + " lbf/in^2 (psi).\n");

% Problem Number: 2
% Define/Update Given Parameters 
L1 = 20; L2 = 10; L3 = 30; % pipe length; ft 
D1 = 2/12; D2 = 2.5/12; D3 = 1.5/12; % pipe diameter; converting in. to ft
A1 = (pi*(D1)^2)/4; A2 = (pi*(D2)^2)/4; A3 = (pi*(D3)^2)/4; % cross-sectional area for all pipes; ft^2
Q = 750*0.002228; % Total Flow Rate; converting gpm to ft^3/s
% Initial Guesses 
Q1 = Q/4; Q2 = Q/2; Q3 = Q/4; dP = 50000; 
Q0 = [Q1;Q2;Q3;dP];
% Intial Friction Factor Guesses
f1 = 0.25*(log10(((e1/D1)/3.7)+((5.74)/(myRe(myV(Q1,A1),D1)^0.9))))^(-2); %Colebrook Equation
f2 = 0.25*(log10(((e2/D2)/3.7)+((5.74)/(myRe(myV(Q2,A2),D2)^0.9))))^(-2); %Colebrook Equation
f3 = 0.25*(log10(((e3/D3)/3.7)+((5.74)/(myRe(myV(Q3,A3),D3)^0.9))))^(-2); %Colebrook Equation
z = fsolve(@(x) equations(x,Q),Q0,options);
fprintf("\nTest 2: \n" ...
        + "The flow rate for pipe 1 (Q1) is: " + (z(1)/0.002228) + " gpm.\n" ...
        + "The flow rate for pipe 2 (Q2) is: " + (z(2)/0.002228) + " gpm.\n" ...
        + "The flow rate for pipe 3 (Q3) is: " + (z(3)/0.002228) + " gpm.\n" ...
        + "The total flow rate (Q) is: " + ((z(1)/0.002228)+(z(2)/0.002228)+(z(3)/0.002228)) + " gpm.\n" ...
        + "The total change in pressure (dP) is: " + (z(4)/144) + " lbf/in^2 (psi).\n");

% Problem Number: 3
% Update Q for a range of values from 100 to 1500 gpm
Q = 100:100:1500;
% Initialize Q1,Q2,Q3 as arrays with the same size as Q. 
Q1 = zeros(size(Q));
Q2 = zeros(size(Q)); 
Q3 = zeros(size(Q)); 
dP = zeros(size(Q)); 
for i=1:length(Q)
    % Initial Guesses for Q1,Q2,Q3,dP
    Q10 = (Q(i)/4)*0.002228; 
    Q20 = (Q(i)/2)*0.002228; 
    Q30 = (Q(i)/4)*0.002228; 
    dP0 = 50000; 
    Q0 = [Q10;Q20;Q30;dP0]; 
    % Intial Friction Factor Guesses
    f1 = 0.25*(log10(((e1/D1)/3.7)+((5.74)/(myRe(myV(Q10,A1),D1)^0.9))))^(-2); %Colebrook Equation
    f2 = 0.25*(log10(((e2/D2)/3.7)+((5.74)/(myRe(myV(Q20,A2),D2)^0.9))))^(-2); %Colebrook Equation
    f3 = 0.25*(log10(((e3/D3)/3.7)+((5.74)/(myRe(myV(Q30,A3),D3)^0.9))))^(-2); %Colebrook Equation
    z = fsolve(@(x) equations(x,(Q(i)*0.002228)),Q0,options); 
    Q1(i) = (z(1)/0.002228);
    Q2(i) = (z(2)/0.002228);
    Q3(i) = (z(3)/0.002228); 
    dP(i) = (z(4)/144); 
end 
fig1 = figure('Name','Pipe Flow Rates vs. Total Flow Rate'); 
ax1 = axes(fig1); 
plot(ax1,Q,Q1,'r-o',MarkerSize=8);
grid on; 
axis([min(Q) max(Q) min(Q3) max(Q2)]);
xlabel("Total Flow Rate (gpm)","FontSize",14);
ylabel("Pipe Flow Rates (gpm)","FontSize",14);
title("Pipe Flow Rates vs. Total Flow Rate","FontSize",14);
hold on; 
plot(ax1,Q,Q2,'b-o',MarkerSize=8);
hold on; 
plot(ax1,Q,Q3,'g-o',MarkerSize=8); 
hold on; 
fig2 = figure('Name','Change in Pressure vs. Total Flow Rate'); 
ax2 = axes(fig2); 
plot(ax2,Q,dP,'g-d',LineWidth=2);
axis([min(Q) max(Q) min(dP) max(dP)])
grid on; 
xlabel("Total Flow Rate (gpm)","FontSize",14);
ylabel("Change in Pressure (psi)","FontSize",14);
title("Change in Pressure vs. Total Flow Rate","FontSize",14);
hold on; 

% Problem Number: 4a) 
% Update (e/D) by an increase of 25% by multiplying the e values by 1.25
e1 = e1*1.25; e2 = e2*1.25; e3 = e3*1.25; 
% Initialize Q1,Q2,Q3 as arrays with the same size as Q. 
Q1 = zeros(size(Q));
Q2 = zeros(size(Q)); 
Q3 = zeros(size(Q)); 
dP = zeros(size(Q)); 
for i = 1:length(Q)
    % Initial Guesses for Q1,Q2,Q3,dP
    Q10 = (Q(i)/4)*0.002228; 
    Q20 = (Q(i)/2)*0.002228; 
    Q30 = (Q(i)/4)*0.002228; 
    dP0 = 50000; 
    Q0 = [Q10;Q20;Q30;dP0];
    % Intial Friction Factor Guesses
    f1 = 0.25*(log10(((e1/D1)/3.7)+((5.74)/(myRe(myV(Q10,A1),D1)^0.9))))^(-2); %Colebrook Equation
    f2 = 0.25*(log10(((e2/D2)/3.7)+((5.74)/(myRe(myV(Q20,A2),D2)^0.9))))^(-2); %Colebrook Equation
    f3 = 0.25*(log10(((e3/D3)/3.7)+((5.74)/(myRe(myV(Q30,A3),D3)^0.9))))^(-2); %Colebrook Equation
    z = fsolve(@(x) equations(x,(Q(i)*0.002228)),Q0,options); 
    Q1(i) = (z(1)/0.002228);
    Q2(i) = (z(2)/0.002228);
    Q3(i) = (z(3)/0.002228); 
    dP(i) = (z(4)/144); 
end
plot(ax1,Q,Q1,'c-^',MarkerSize=8); 
hold on; 
plot(ax1,Q,Q2,'m-^',MarkerSize=8); 
hold on; 
plot(ax1,Q,Q3,'y-^',MarkerSize=8); 
hold on; 
plot(ax2,Q,dP,'r-^',LineWidth=2); 
axis([min(Q) max(Q) min(dP) max(dP)])
hold on; 

% Problem Number: 4b) 
% Update (e/D) by an increase of 35% by multiplying the e values by 1.35
e1 = e1*1.35; e2 = e2*1.35; e3 = e3*1.35; 
% Initialize Q1,Q2,Q3 as arrays with the same size as Q. 
Q1 = zeros(size(Q));
Q2 = zeros(size(Q)); 
Q3 = zeros(size(Q)); 
dP = zeros(size(Q)); 
for i = 1:length(Q)
    % Initial Guesses for Q1,Q2,Q3,dP
    Q10 = (Q(i)/4)*0.002228; 
    Q20 = (Q(i)/2)*0.002228; 
    Q30 = (Q(i)/4)*0.002228; 
    dP0 = 50000;
    Q0 = [Q10;Q20;Q30;dP0]; 
    % Intial Friction Factor Guesses
    f1 = 0.25*(log10(((e1/D1)/3.7)+((5.74)/(myRe(myV(Q10,A1),D1)^0.9))))^(-2); %Colebrook Equation
    f2 = 0.25*(log10(((e2/D2)/3.7)+((5.74)/(myRe(myV(Q20,A2),D2)^0.9))))^(-2); %Colebrook Equation
    f3 = 0.25*(log10(((e3/D3)/3.7)+((5.74)/(myRe(myV(Q30,A3),D3)^0.9))))^(-2); %Colebrook Equation
    z = fsolve(@(x) equations(x,(Q(i)*0.002228)),Q0,options); 
    Q1(i) = (z(1)/0.002228);
    Q2(i) = (z(2)/0.002228);
    Q3(i) = (z(3)/0.002228); 
    dP(i) = (z(4)/144); 
end 
plot(ax1,Q,Q1,'k-*',MarkerSize=8); 
hold on; 
plot(ax1,Q,Q2,'g-*',MarkerSize=8); 
hold on; 
plot(ax1,Q,Q3,'r-*',MarkerSize=8); 
legend(ax1,'Pipe 1 Flow Rate for Test 1','Pipe 2 Flow Rate for Test 1','Pipe 3 Flow Rate for Test 1','Pipe 1 Flow Rate for Test 2','Pipe 2 Flow Rate for Test 2','Pipe 3 Flow Rate for Test 2','Pipe 1 Flow Rate for Test 3','Pipe 2 Flow Rate for Test 3','Pipe 3 Flow Rate for Test 3','Location','best');
plot(ax2,Q,dP,'b-o',LineWidth=2);
axis([min(Q) max(Q) min(dP) max(dP)])
legend(ax2,'Clean Pipes','Fouled Pipes (25% (e/D) Increase)','Fouled Pipes (35% (e/D) Increase)','Location','best','FontSize',14);

% Functions Used %
    function V = myV(Q,A) % Function used to compute the effective flow velocity
        V = (Q/A); % Velocity; ft/s
    end
    function Re = myRe(V,D) % Function used to compute the Reynolds Number
        Re = (rho*V*D)/mu; % Reynolds Number; dimensionless 
    end
    function F = myF(Re,e,D,f) % Function used to compute the friction factor
        if Re <= 2300 % Laminar Flow 
            F = 64/Re; % Friction Factor; Dimensionless 
        else % Turbulent Flow
            F = (inv(-2*log10(((e/D)/3.7)+(2.51/(Re*sqrt(f))))))^2; % Friction Factor; Dimensionless 
        end
    end
    function g = equations(x,Q) % Function used to solve for Q1,Q2,Q3,and dP
        g = [(myF(myRe(myV(x(1),A1),D1),e1,D1,f1)*(L1/D1)*(1/2)*((myV(x(1),A1))^2))-(x(4)/rho);
             (myF(myRe(myV(x(2),A2),D2),e2,D2,f2)*(L2/D2)*(1/2)*((myV(x(2),A2))^2))-(x(4)/rho);
             (myF(myRe(myV(x(3),A3),D3),e3,D3,f3)*(L3/D3)*(1/2)*((myV(x(3),A3))^2))-(x(4)/rho);
             x(1)+x(2)+x(3)-Q]; 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%