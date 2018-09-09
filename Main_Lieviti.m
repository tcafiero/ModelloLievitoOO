% YEAST model solver MAIN
% data scrittura: 28/10/2014
% autore: Fabrizio Carteni

% close all
clear all
clear global all
global t_f t_end mu Initial_glucose

% PARAMETERS TO SET MANUALLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_f=15;                         % Batch end time
t_end=40;                       % End simulation time
mu=0.16;                        % Feeding rate
Initial_glucose=20;             % Initial glucose concentration
t_start=0;
x = inputdlg({'Batch end time','End simulation time','Feeding rate','Initial glucose concentration'},'Set Simulation Parameter', [1 50; 1 50; 1 50; 1 50], {'15','49', '0.16', '20'});
[val status] = str2num(x{1});  
if ~status
    % Handle empty value returned 
    % for unsuccessful conversion
    % ...
end
t_f=val;                         % Batch end time
[val status] = str2num(x{2});  
if ~status
    % Handle empty value returned 
    % for unsuccessful conversion
    % ...
end
t_end=val;                       % End simulation time
[val status] = str2num(x{3});  
if ~status
    % Handle empty value returned 
    % for unsuccessful conversion
    % ...
end
mu=val;                        % Feeding rate
[val status] = str2num(x{4});  
if ~status
    % Handle empty value returned 
    % for unsuccessful conversion
    % ...
end
Initial_glucose=val;             % Initial glucose concentration

% OTHER PARAMETERS (DO NOT CHANGE)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Initial values [V G E P Cm I D R]
x0=[1 Initial_glucose 0 0.00005 0.022 0 0 0];

% Set up simulation and solve equations
YeastObject=YeastClassSpecial;
tspan=t_start:0.1:t_end;
options = odeset('NonNegative',[1 2 3 4 5 6 7 8],'MaxStep',0.1);
[t,x] = ode15s(@(t, y)YeastObject.Lieviti_eqs(YeastObject, t, y), tspan, x0,  options);

% Calculate variables to plot
Mass_conc=(x(:,4)+x(:,5)+x(:,7)+x(:,8))./x(:,1);
Viable_conc=(x(:,4)+x(:,5))./x(:,1);
Eth_conc= x(:,3)./x(:,1);
Glu_conc= x(:,2)./x(:,1);
Dead_conc= x(:,7)./x(:,1);
I_conc= x(:,6)./x(:,1);
Py_conc=x(:,4)./(x(:,4)+x(:,5)+x(:,8))*YeastObject.c;
R_conc=x(:,8)/x(:,1);
Res_Perc=x(:,8)./(x(:,4)+x(:,5)+x(:,8));
NF=YeastObject.sigma_i*I_conc/YeastObject.I_max;
NE=YeastObject.sigma_e*Eth_conc/YeastObject.E_max;
Batch_end_vector=ones(1,size(tspan,2)).*(tspan<t_f)*1000;

% Plots
% Figure Settings [left, bottom, width, height]
Fig=figure('Position',[150 0 1600 1250]);

figure(Fig)
subplot(3,1,1)
area(tspan,Batch_end_vector,'FaceColor',[.9 .9 .9],'EdgeColor','none')
hold on
plot(tspan,Mass_conc,'k','LineWidth',1.2)
hold off
set(gca,'layer','top')
xlim([0 t_end])
ylim([0 max(max(Mass_conc))*1.1])
ylabel('[X] (g/l)');

subplot(3,1,2)
area(tspan,Batch_end_vector,'FaceColor',[.9 .9 .9],'EdgeColor','none')
hold on
plot(tspan,Eth_conc,'k','LineWidth',1.2)
hold off
set(gca,'layer','top')
xlim([0 t_end])
ylim([0 max(max(Eth_conc))*1.1])
ylabel('[Ethanol] (g/l)');

subplot(3,1,3)
area(tspan,Batch_end_vector,'FaceColor',[.9 .9 .9],'EdgeColor','none')
hold on
plot(tspan,Glu_conc,'k','LineWidth',1.2)
hold off
set(gca,'layer','top')
xlim([0 t_end])
ylim([0 max(max(Glu_conc))*1.1])
ylabel('[Glucose] (g/l)');
xlabel('Time (h)')