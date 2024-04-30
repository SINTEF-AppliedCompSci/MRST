%% Create example CO2-brine kr-Pc properties as tables
% Brooks-Corey type properties
% Relative permeability
np = 24;
% Reservoir
sg1 = linspace(0,0.95,np)';
krg1 = sg1.^2;
krw1 = (1-sg1).^5;  krw1(end) = 0;

% Caprock
sg2 = linspace(0,0.8,np-1)';
sg2 = [sg2(1); 1e-5; sg2(2:end)];
krg2 = sg2.^1.5;
krw2 = (1-sg2).^5;  krw2(end) = 0;

% Fault (mix)
sg3 = linspace(0,0.875,np-1)';
sg3 = [sg3(1); 1e-5; sg3(2:end)];
krg3 = sg3.^1.75;
krw3 = (1-sg3).^5;  krw3(end) = 0;


% Capillary pressure
pc1 = zeros(np, 1);

pce2 = 0.06;            % capillary entry pressure seal [bar]
lambda = 2;
pc2 = pce2*(1-sg2).^(-1/lambda); pc2(1) = 0;

pce3 = 0.008;            % capillary entry pressure fault [bar]
lambda = .8;
pc3 = pce3*(1-sg3).^(-1/lambda); pc3(1) = 0;

% Generate tables for input to .DATA file
hdrs = {'SGAS', 'KRG', 'KRW', 'PCOG'};
vartyp  = cellstr(repmat('double', numel(hdrs), 1));
clear t
t.reservoir = table('Size', [np, 4], 'VariableTypes', vartyp, ...
                   'VariableNames', hdrs); 
t.caprock = table('Size', [np, 4], 'VariableTypes', vartyp, ...
                   'VariableNames', hdrs); 
t.fault = table('Size', [np, 4], 'VariableTypes', vartyp, ...
                   'VariableNames', hdrs); 
t.reservoir = [sg1, krg1, krw1, pc1];
t.caprock = [sg2, krg2, krw2, pc2];
t.fault = [sg3, krg3, krw3, pc3];

%% Plot relative permeabilities
figure(1)
subplot(1,3,1)
hold on
plot((1-sg1), krw1, '-b')
plot((1-sg1), krg1, '-r')
grid on
xlabel('S_w [-]')
ylabel('k_r [-]')
xlim([0 1])
ylim([0 1])
grid on
title('Reservoir')
legend('k_{r,w}', 'k_{r,g}')
subplot(1,3,2)
hold on
plot((1-sg2), krw2, '-b')
plot((1-sg2), krg2, '-r')
grid on
xlabel('S_w [-]')
ylabel('k_r [-]')
xlim([0 1])
ylim([0 1])
grid on
title('Seal')
subplot(1,3,3)
hold on
plot((1-sg3), krw3, '-b')
plot((1-sg3), krg3, '-r')
grid on
xlabel('S_w [-]')
ylabel('k_r [-]')
xlim([0 1])
ylim([0 1])
grid on
title('Fault')
%% Plot capillary pressures
figure(2)
hold on
plot((1-sg2), pc2*1e3, '-k')
plot((1-sg3), pc3*1e3, '--k')
grid on
xlabel('S_w [-]')
ylabel('P_c [mbar]')
xlim([0 1])
legend('Seal', 'Fault')
grid on