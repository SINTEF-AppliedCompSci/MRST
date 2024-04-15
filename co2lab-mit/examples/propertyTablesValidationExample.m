%% Validation of pvtBrineWithCO2BlackOil
% Here, we validate our implementation by comparing the output with
% published results. Specifically, we compare our output with:
%   Fig. 4 and Fig. 5 in Hassanzadeh et al, IGGC (2008).
%   Fig. 2 in Spycher and Pruess, GGA (2005)
%    
% NOTE: Although Fig. 5 in Hassanzadeh et al. (2008) states 50C, the plots 
% correspond to T=45C (Hassan Hassanzadeh, personal communication, 2019).
%
% Documentation on the functionality itself is provided within
% pvtBrineWithCO2BlackOil.
%
% Note that the table generation produces warnings due to evaluation of
% functions outside their strict range of validity. These are nothing to
% worry about from a user's perspective, but we include the command for
% disabling warnings in your current MATLAB session if needed:
% warning('off', 'all')              % disable warnings, use with caution


%% Comparison with Fig. 5 in Hassanzadeh et al. IJGGC (2008)
T = {'C', 45};                      % temperature
P = {'MPa','mMn',[0.1, 50, 100]};   % pressure range
S = {'ppm', 'NaCl', 1.5e5};         % salinity
saltVar = false;                    % see pvtBrineWithCO2BlackOil
vapH2O = false;                     % "
figs = true;                        % set to false if you don't want figs

% Generate properties and plot figures (make sure to select correct T, etc)
% Compare Fig. 1 and Fig. 4 with fig 5c and Fig. 5a in Hassanzadeh et al.
% (2008), Fig. 2 with Fig.2 in Spycher and Pruess (2005), Fig. 3 with
% Fig. 4 in Hassanzadeh et al. (2008), and Fig. 5 with Fig. 5d in
% Hassanzadeh et al. (2008) (T=45 C)
pvtBrineWithCO2BlackOil(T, P, S, saltVar, vapH2O, false, figs);



%% Comparison with Figure 2 in Spycher and Pruess, Geochim Cosmochim Ac (2005)
T_val = [30, 60, 90];
s_val = [0 1 2 4];
P = {'MPa', 'vals', linspace(0.1,60,50)};
saltVar = true;
vapH2O = true;
out = cell(numel(T_val),numel(s_val));
for t = 1:numel(T_val)
    for s = 1:numel(s_val) 
        T = {'C', T_val(t)};
        S = {'m', 'NaCl', s_val(s)};
        out{t,s} = pvtBrineWithCO2BlackOil(T, P, S, saltVar, vapH2O, false, false);
    end
end

% Plot
col_co2  = [100, 0, 0; 150, 0, 0; 200, 0, 0; 250, 0, 0]./255;
col_brine = [0, 0, 100; 0, 0, 150; 0, 0, 200; 0, 0, 250]/255;
latx       = {'Interpreter','latex'};
h = figure(10);
tiledlayout(3, 2,'tilespacing','compact','Padding','compact')
nexttile(1)         % --------
hold on
p1 = plot(out{1,1}.P_bar, out{1,1}.aq.m_co2, 'color', col_co2(1,:), ...
          'linewidth', 1, 'DisplayName', '0 m'); 
p2 = plot(out{1,2}.P_bar, out{1,2}.aq.m_co2, 'color', col_co2(2,:), ...
          'LineWidth', 1, 'DisplayName', '1 m'); 
p3 = plot(out{1,3}.P_bar, out{1,3}.aq.m_co2, 'color', col_co2(3,:), ...
          'linewidth', 1, 'DisplayName', '2 m'); 
p4 = plot(out{1,4}.P_bar, out{1,4}.aq.m_co2, 'color', col_co2(4,:), ...
          'LineWidth', 1, 'DisplayName', '4 m'); 
hold off
grid on; ylim([0 2]); xlim([0 600])
yticks([0 0.5 1 1.5 2]); xticks(0:100:600)
ylabel('CO$_{2(aq)}$ [m]','fontsize', 14, latx{:})
xlabel('$p$ [bar]','fontsize', 14, latx{:})
leg = legend([p1 p2 p3 p4], 'fontsize', 12);
leg.BoxFace.ColorType='truecoloralpha';
leg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
nexttile(2)         % -------
hold on
p1 = plot(out{1,1}.P_bar, out{1,1}.gas.Y_h2o*1000, 'color', col_brine(1,:), ...
     'linewidth', 1, 'DisplayName', '0 m');
p2 = plot(out{1,2}.P_bar, out{1,2}.gas.Y_h2o*1000, 'color', col_brine(2,:), ...
     'linewidth', 1, 'DisplayName', '1 m'); 
p3 = plot(out{1,3}.P_bar, out{1,3}.gas.Y_h2o*1000, 'color', col_brine(3,:), ...
     'linewidth', 1, 'DisplayName', '2 m'); 
p4 = plot(out{1,4}.P_bar, out{1,4}.gas.Y_h2o*1000, 'color', col_brine(4,:), ... 
     'linewidth', 1, 'DisplayName', '4 m'); 
hold off
grid on; ylim([0 10]); xlim([0 600])
yticks(0:10); xticks(0:100:600)
ylabel('$y_\mathrm{H_{2}O}\times10^3$ [-]','fontsize', 14, latx{:})
leg = legend([p1 p2 p3 p4], 'fontsize', 12);
leg.BoxFace.ColorType='truecoloralpha';
leg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
nexttile(3)         % --------
hold on
plot(out{2,1}.P_bar, out{2,1}.aq.m_co2, 'color', col_co2(1,:), 'linewidth', 1); 
plot(out{2,2}.P_bar, out{2,2}.aq.m_co2, 'color', col_co2(2,:), 'LineWidth', 1); 
plot(out{2,3}.P_bar, out{2,3}.aq.m_co2, 'color', col_co2(3,:), 'linewidth', 1); 
plot(out{2,4}.P_bar, out{2,4}.aq.m_co2, 'color', col_co2(4,:), 'LineWidth', 1); 
hold off
grid on; ylim([0 2]); xlim([0 600])
yticks([0 0.5 1 1.5 2]); xticks(0:100:600)
nexttile(4)         % -------
hold on
plot(out{2,1}.P_bar, out{2,1}.gas.Y_h2o*1000, 'color', col_brine(1,:), 'linewidth', 1); 
plot(out{2,2}.P_bar, out{2,2}.gas.Y_h2o*1000, 'color', col_brine(2,:), 'linewidth', 1); 
plot(out{2,3}.P_bar, out{2,3}.gas.Y_h2o*1000, 'color', col_brine(3,:), 'linewidth', 1); 
plot(out{2,4}.P_bar, out{2,4}.gas.Y_h2o*1000, 'color', col_brine(4,:), 'linewidth', 1); 
hold off
grid on; ylim([4 14]); xlim([0 600])
yticks(4:14); xticks(0:100:600)
nexttile(5)         % --------
hold on
plot(out{3,1}.P_bar, out{3,1}.aq.m_co2, 'color', col_co2(1,:), 'linewidth', 1); 
plot(out{3,2}.P_bar, out{3,2}.aq.m_co2, 'color', col_co2(2,:), 'LineWidth', 1); 
plot(out{3,3}.P_bar, out{3,3}.aq.m_co2, 'color', col_co2(3,:), 'linewidth', 1); 
plot(out{3,4}.P_bar, out{3,4}.aq.m_co2, 'color', col_co2(4,:), 'LineWidth', 1); 
hold off
grid on; ylim([0 2]); xlim([0 600])
yticks([0 0.5 1 1.5 2]); xticks(0:100:600)
nexttile(6)         % -------
hold on
plot(out{3,1}.P_bar, out{3,1}.gas.Y_h2o*1000, 'color', col_brine(1,:), 'linewidth', 1); 
plot(out{3,2}.P_bar, out{3,2}.gas.Y_h2o*1000, 'color', col_brine(2,:), 'linewidth', 1); 
plot(out{3,3}.P_bar, out{3,3}.gas.Y_h2o*1000, 'color', col_brine(3,:), 'linewidth', 1); 
plot(out{3,4}.P_bar, out{3,4}.gas.Y_h2o*1000, 'color', col_brine(4,:), 'linewidth', 1); 
hold off
grid on; ylim([10 20]); xlim([0 600])
yticks(10:20); xticks(0:100:600)
set(h, 'Position', [600, 600, 800, 600])
%exportgraphics(h,'untitled.pdf','ContentType','vector')


%% Comparison with Fig. 8 and 9 in Spycher et al. Geochim Cosmochim Ac (2003)
T_val = [6.85, 9.85, 16.85, 26.85, 30.98 36.85:10:106.85];
S = {'m', 'NaCl', 1}; 
P = {'MPa', 'vals', [linspace(0.1,25,120) linspace(25.5,60,60)]};
saltVar = false;
vapH2O = false;
out = cell(numel(T_val),1);
for t = 1:numel(T_val)
    T = {'C', T_val(t)};
    out{t} = pvtBrineWithCO2BlackOil(T, P, S, saltVar, vapH2O, false, false);
end

% Plot
col_co2  = repelem([20, 0, 0], numel(T_val), 1);
col_co2(:,1) = [(20:20:240) 255]/255;
h = figure(11); % Compare with Fig. 8 in Spycher et al. (2003)
tiledlayout(1, 2,'tilespacing','compact','Padding','compact')
nexttile(1)         % --------
hold on
plot(out{2}.P_bar,out{2}.gas.Z_co2, 'color', col_co2(2,:),'linewidth', 1, ...
    'displayname', [num2str(T_val(2)) ' ºC'] ); 
plot(out{3}.P_bar,out{3}.gas.Z_co2, 'color', col_co2(3,:),'linewidth', 1, 'displayname', num2str(T_val(3)))
plot(out{4}.P_bar,out{4}.gas.Z_co2, 'color', col_co2(4,:),'linewidth', 1, 'displayname', num2str(T_val(4)))
plot(out{5}.P_bar,out{5}.gas.Z_co2, 'color', col_co2(5,:),'linewidth', 1, 'displayname', num2str(T_val(5)))
plot(out{6}.P_bar,out{6}.gas.Z_co2, 'color', col_co2(6,:),'linewidth', 1, 'displayname', num2str(T_val(6)))
plot(out{7}.P_bar,out{7}.gas.Z_co2, 'color', col_co2(7,:),'linewidth', 1, 'displayname', num2str(T_val(7)))
plot(out{8}.P_bar,out{8}.gas.Z_co2, 'color', col_co2(8,:),'linewidth', 1, 'displayname', num2str(T_val(8)))
plot(out{9}.P_bar,out{9}.gas.Z_co2, 'color', col_co2(9,:),'linewidth', 1, 'displayname', num2str(T_val(9)))
plot(out{10}.P_bar,out{10}.gas.Z_co2, 'color', col_co2(10,:),'linewidth', 1, 'displayname', num2str(T_val(10)))
plot(out{11}.P_bar,out{11}.gas.Z_co2, 'color', col_co2(11,:),'linewidth', 1, 'displayname', num2str(T_val(11)))
plot(out{12}.P_bar,out{12}.gas.Z_co2, 'color', col_co2(12,:),'linewidth', 1, 'displayname', num2str(T_val(12)))
plot(out{13}.P_bar,out{13}.gas.Z_co2, 'color', col_co2(13,:),'linewidth', 1, 'displayname', num2str(T_val(13)))
grid on; 
ylim([0 1.1]); xlim([0 600])
yticks(0:.1:1.1); xticks(0:100:600)
ylabel('$Z$ [-]','fontsize', 14, latx{:})
xlabel('$p$ [bar]','fontsize', 14, latx{:})
leg = legend('fontsize', 12,'location','best');
leg.BoxFace.ColorType='truecoloralpha';
leg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
nexttile(2)         % --------
hold on
plot(out{1}.P_bar,out{1}.gas.phi_co2, 'color', col_co2(1,:), 'linewidth', 1, ...
     'displayname', [num2str(T_val(1)) ' ºC']); 
plot(out{3}.P_bar,out{3}.gas.phi_co2, 'color', col_co2(3,:),'linewidth', 1, 'displayname', num2str(T_val(3)))
plot(out{4}.P_bar,out{4}.gas.phi_co2, 'color', col_co2(4,:),'linewidth', 1, 'displayname', num2str(T_val(4)))
plot(out{6}.P_bar,out{6}.gas.phi_co2, 'color', col_co2(6,:),'linewidth', 1, 'displayname', num2str(T_val(6)))
plot(out{7}.P_bar,out{7}.gas.phi_co2, 'color', col_co2(7,:),'linewidth', 1, 'displayname', num2str(T_val(7)))
plot(out{8}.P_bar,out{8}.gas.phi_co2, 'color', col_co2(8,:),'linewidth', 1, 'displayname', num2str(T_val(8)))
plot(out{9}.P_bar,out{9}.gas.phi_co2, 'color', col_co2(9,:),'linewidth', 1, 'displayname', num2str(T_val(9)))
plot(out{10}.P_bar,out{10}.gas.phi_co2, 'color', col_co2(10,:),'linewidth', 1, 'displayname', num2str(T_val(10)))
plot(out{11}.P_bar,out{11}.gas.phi_co2, 'color', col_co2(11,:),'linewidth', 1, 'displayname', num2str(T_val(11)))
plot(out{12}.P_bar,out{12}.gas.phi_co2, 'color', col_co2(12,:),'linewidth', 1, 'displayname', num2str(T_val(12)))
plot(out{13}.P_bar,out{13}.gas.phi_co2, 'color', col_co2(13,:),'linewidth', 1, 'displayname', num2str(T_val(13)))
hold off
grid on;
ylim([0.1 1]); xlim([0 600])
yticks(0.1:.1:1); xticks(0:100:600)
ylabel('$\phi$ [-]','fontsize', 14, latx{:})
leg = legend('fontsize', 12,'location','best');
leg.BoxFace.ColorType='truecoloralpha';
leg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
set(h, 'Position', [600, 600, 800, 300])


%% Comparison with Fig. 14,15,16 in Yan et al., IJGGC (2011)
T_val = [323.2 373.2] - 273.15;
s_val = [0, 1, 5];
P = {'MPa', 'vals', 5:2.5:40};
saltVar = true;
vapH2O = true;
out = cell(numel(T_val),numel(s_val));
for t = 1:numel(T_val)
    for s = 1:numel(s_val) 
        T = {'C', T_val(t)};
        S = {'m', 'NaCl', s_val(s)};
        out{t,s} = pvtBrineWithCO2BlackOil(T, P, S, saltVar, vapH2O, false, false);
    end
end

%             p    T=323.2   373.2  423.2
data_yan = [5.00 0.99722 0.96370 0.92928        % 0m
            10.00 1.00268 0.96741 0.93367      
            15.00 1.00528 0.97062 0.93760       
            20.00 1.00688 0.97425 0.94108       
            30.00 1.01293 0.97962 0.94700       
            40.00 1.01744 0.98506 0.95282       
            5.00 1.03116 1.00026 0.96883    % 1m NaCl
            10.00 1.03491 1.00321 0.97169
            15.00 1.03968 1.00667 0.97483
            20.00 1.04173 1.00961 0.97778
            30.00 1.04602 1.01448 0.98301
            40.00 1.05024 1.01980 0.98817
            5.00 1.15824 1.12727 1.09559    % 5m NaCl
            10.00 1.16090 1.12902 1.10183
            15.00 1.16290 1.13066 1.10349
            20.00 1.16468 1.13214 1.10499
            30.00 1.16810 1.13566 1.10882
            40.00 1.17118 1.13893 1.11254];

% viscosities
P = {'bar', 'vals', [50 200 600]};
T_v = 20:5:100;
S = {'m', 'NaCl', 1};
S2 = {'m', 'NaCl', 0};
S3 = {'ppm', 'seawater', 35e3};
saltVar = true;
vapH2O = true;
mu_b_cP = zeros(numel(T_v), numel(P{3}));
mu_aq_cP = zeros(numel(T_v), numel(P{3}));
mu_h2o_cP = zeros(numel(T_v), numel(P{3}));
mu_h2o_co2_cP = zeros(numel(T_v), numel(P{3}));
mu_sea_cP = zeros(numel(T_v), numel(P{3}));
for t = 1:numel(T_v)
    T = {'C', T_v(t)};
    data = pvtBrineWithCO2BlackOil(T, P, S, saltVar, vapH2O, false, false);
    mu_b_cP(t, :) = data.aq.mu_b_cP(1:3);
    mu_aq_cP(t, :) = data.aq.mu_aq_cP(1:3);
    data2 = pvtBrineWithCO2BlackOil(T, P, S2, saltVar, vapH2O, false, false);
    mu_h2o_cP(t,:) = data2.aq.mu_b_cP(1:3);
    mu_h2o_co2_cP(t, :) = data2.aq.mu_aq_cP(1:3);
    data3 = pvtBrineWithCO2BlackOil(T, P, S3, saltVar, vapH2O, false, false);
    mu_sea_cP(t,:) = data3.aq.mu_b_cP(1:3);
end

col_b  = zeros(6, 3);
col_b(:,3) = (40:40:250)/255;
latx={'interpreter','latex'};
h = figure(12); % Compare with Fig. 8 in Spycher et al. (2003)
tiledlayout(2, 3,'tilespacing','compact','Padding','compact')
nexttile(1)         % --------
hold on
p1 = plot(data_yan(1:6, 1)*10, data_yan(1:6, 2)*1e3, 'o', 'markeredgecolor', col_b(1,:), ...
          'displayname', 'Yan (2011)');
plot(data_yan(1:6, 1)*10, data_yan(1:6, 3)*1e3, 'o', 'markeredgecolor', col_b(4,:));
p2 = plot(out{1,1}.P_bar,out{1,1}.aq.rho_b_kgm3, '--', 'color', col_b(1,:),'linewidth', 1, ...
    'displayname', [num2str(T_val(1)) ' $^\circ$C, ' 'no CO$_2$'] ); 
p3 = plot(out{1,1}.P_bar,out{1,1}.aq.rho_aq_kgm3, '-', 'color', col_b(1,:),'linewidth', 1, ...
    'displayname', [num2str(T_val(1)) ', CO$_2$ sat.']); 
p4 = plot(out{2,1}.P_bar,out{2,1}.aq.rho_b_kgm3, '--', 'color', col_b(4,:),'linewidth', 1, ...
    'displayname', [num2str(T_val(2))] ); 
p5 = plot(out{2,1}.P_bar,out{2,1}.aq.rho_aq_kgm3, '-', 'color', col_b(4,:),'linewidth', 1, ...
    'displayname', [num2str(T_val(2))]); 
grid on; 
ylim([920 1020]); xlim([50 400])
yticks(920:10:1020); 
xticks(0:50:400)
ylabel('$\rho$ [kg/m$^3$]','fontsize', 14, latx{:})
xlabel('$p$ [bar]','fontsize', 14, latx{:})
leg = legend([p1 p2 p3 p4], 'fontsize', 12,'location','best',latx{:});
leg.BoxFace.ColorType='truecoloralpha';
leg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
title('Density, pure water')
nexttile(2)         % --------
hold on
p1 = plot(data_yan(7:12, 1)*10, data_yan(7:12, 2)*1e3, 'o', 'markeredgecolor', col_b(2,:), ...
          'displayname', 'Yan (2011)');
plot(data_yan(7:12, 1)*10, data_yan(7:12, 3)*1e3, 'o', 'markeredgecolor', col_b(5,:));
p2 = plot(out{1,2}.P_bar,out{1,2}.aq.rho_b_kgm3, '--', 'color', col_b(2,:),'linewidth', 1); 
p3 = plot(out{1,2}.P_bar,out{1,2}.aq.rho_aq_kgm3, '-', 'color', col_b(2,:),'linewidth', 1); 
p4 = plot(out{2,2}.P_bar,out{2,2}.aq.rho_b_kgm3, '--', 'color', col_b(5,:),'linewidth', 1); 
p5 = plot(out{2,2}.P_bar,out{2,2}.aq.rho_aq_kgm3, '-', 'color', col_b(5,:),'linewidth', 1); 
grid on; 
ylim([960 1060]); xlim([50 400])
yticks(960:10:1060); 
xticks(0:50:400)
%ylabel('$\rho$ [kg/m$^3$]','fontsize', 14, latx{:})
%xlabel('$p$ [bar]','fontsize', 14, latx{:})
title('1 m NaCl')
nexttile(3)         % --------
hold on
p1 = plot(data_yan(13:18, 1)*10, data_yan(13:18, 2)*1e3, 'o', 'markeredgecolor', col_b(3,:), ...
          'displayname', 'Yan (2011)');
plot(data_yan(13:18, 1)*10, data_yan(13:18, 3)*1e3, 'o', 'markeredgecolor', col_b(6,:));
p2 = plot(out{1,3}.P_bar,out{1,3}.aq.rho_b_kgm3, '--', 'color', col_b(3,:),'linewidth', 1); 
p3 = plot(out{1,3}.P_bar,out{1,3}.aq.rho_aq_kgm3, '-', 'color', col_b(3,:),'linewidth', 1); 
p4 = plot(out{2,3}.P_bar,out{2,3}.aq.rho_b_kgm3, '--', 'color', col_b(6,:),'linewidth', 1); 
p5 = plot(out{2,3}.P_bar,out{2,3}.aq.rho_aq_kgm3, '-', 'color', col_b(6,:),'linewidth', 1); 
grid on; 
ylim([1080 1180]); xlim([50 400])
yticks(1080:10:1180); 
xticks(0:50:400)
%ylabel('$\rho$ [kg/m$^3$]','fontsize', 14, latx{:})
%xlabel('$p$ [bar]','fontsize', 14, latx{:})
title('5 m NaCl')
nexttile(4)         % --------
hold on
p2 = plot(T_v,mu_h2o_cP(:,1), '-', 'color', col_b(1,:),'linewidth', 1, 'displayName', 'H$_2$O');
p3 = plot(T_v,mu_b_cP(:,1), '--', 'color', col_b(2,:),'linewidth', 1, 'displayName', 'H$_2$O$+$NaCl(1m)');
p4 = plot(T_v,mu_h2o_co2_cP(:,1), ':', 'color', col_b(3,:),'linewidth', 1.5, 'displayName', 'H$_2$O+CO$_2$ sat.');
p5 = plot(T_v,mu_aq_cP(:,1), '-.', 'color', col_b(4,:),'linewidth', 1, 'displayName', 'H$_2$O+NaCl(1m)+CO$_2$ sat.');
p6 = plot(T_v,mu_sea_cP(:,1), '-', 'color', col_b(5,:),'linewidth', 1.5, 'displayName', 'Seawater');
hold off
grid on; 
ylim([0.2 1.4]); xlim([20 100])
yticks(0.2:.2:1.4); 
xticks(20:10:100)
ylabel('$\mu$ [cP]','fontsize', 14, latx{:})
xlabel('$T$ [$^\circ$C]','fontsize', 14, latx{:})
leg = legend([p2 p3 p4 p5 p6], 'fontsize', 12,'location','best',latx{:});
leg.BoxFace.ColorType='truecoloralpha';
leg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
title('Viscosity, 50 bar')
nexttile(5)         % --------
hold on
p2 = plot(T_v,mu_h2o_cP(:,2), '-', 'color', col_b(1,:),'linewidth', 1, 'displayName', 'H$_2$O');
p3 = plot(T_v,mu_b_cP(:,2), '--', 'color', col_b(2,:),'linewidth', 1, 'displayName', 'H$_2$O$+$NaCl');
p4 = plot(T_v,mu_h2o_co2_cP(:,2), ':', 'color', col_b(3,:),'linewidth', 1.5, 'displayName', 'H$_2$O+CO$_2$');
p5 = plot(T_v,mu_aq_cP(:,2), '-.', 'color', col_b(4,:),'linewidth', 1, 'displayName', 'H$_2$O+NaCl+CO$_2$');
p6 = plot(T_v,mu_sea_cP(:,2), '-', 'color', col_b(5,:),'linewidth', 1.5, 'displayName', 'Seawater');
hold off
grid on; 
ylim([0.2 1.4]); xlim([20 100])
yticks(0.2:.2:1.4); 
xticks(20:10:100)
%ylabel('$\mu$ [cP]','fontsize', 14, latx{:})
%xlabel('$T$ [$^\circ$C]','fontsize', 14, latx{:})
title('200 bar')
nexttile(6)         % --------
hold on
p2 = plot(T_v,mu_h2o_cP(:,3), '-', 'color', col_b(1,:),'linewidth', 1, 'displayName', 'H$_2$O');
p3 = plot(T_v,mu_b_cP(:,3), '--', 'color', col_b(2,:),'linewidth', 1, 'displayName', 'H$_2$O$+$NaCl(1m)');
p4 = plot(T_v,mu_h2o_co2_cP(:,3), ':', 'color', col_b(3,:),'linewidth', 1.5, 'displayName', 'H$_2$O+CO$_2$');
p5 = plot(T_v,mu_aq_cP(:,3), '-.', 'color', col_b(4,:),'linewidth', 1, 'displayName', 'H$_2$O+NaCl(1m)+CO$_2$');
p6 = plot(T_v,mu_sea_cP(:,3), '-', 'color', col_b(5,:),'linewidth', 1.5, 'displayName', 'Seawater');
hold off
grid on; 
ylim([0.2 1.4]); xlim([20 100])
yticks(0.2:.2:1.4); 
xticks(20:10:100)
%ylabel('$\mu$ [cP]','fontsize', 14, latx{:})
%xlabel('$T$ [$^\circ$C]','fontsize', 14, latx{:})
title('600 bar')
set(h, 'Position', [600, 600, 800, 500])

%% CO2 density and viscosity

% 1. Density of pure CO2
% Density Data from Holste et al., J Chem Thermodyn (1987)
Mw_co2  = 44.0095/10^3; % kg / mol
data_holste = [300.017 4.82296 2753.1;  % Tab 1 [T_K, P_MPa, 1/V_mol/m3]
               300.014 3.33009 1640.9;  % Tab 1
               300.016 2.16077 978.03;  % Tab 1
               299.924 30.9366 21927;   % Tab 5
               299.943 5.61021 3591.2;  % Tab 5
               300.033 4.26151 2284.6;  % Tab 5
               310.044 39.4082 21924;   % Tab 5
               309.980 21.9169 19839;   % Tab 5
               309.985 8.98531 13946;   % Tab 5
               310.078 8.22716 8871.9;  % Tab 5
               310.055 7.45520 5644.2;  % Tab 5
               310.077 6.06892 3590.7;  % Tab 5
               310.041 4.51573 2284.3;  % Tab 5
               320.019 47.7093 21920;   % Tab 5
               319.956 28.1858 19836;   % Tab 5
               320.061 11.6686 13944;   % Tab 5
               320.033 9.61812 8870.7;  % Tab 5
               320.035 8.24721 5643.4;  % Tab 5
               320.043 6.50816 3590.2;  % Tab 5
               320.121 4.76723 2284.0;  % Tab 5
               319.994 3.31356 1453.1;  % Tab 5
               348.151 6.33575 2760.8;  % Tab 3
               348.154 4.60679 1863.6;  % Tab 3
               348.150 3.27221 1258.0;  % Tab 3
               348.150 2.28670 849.21;  % Tab 3
               373.151 7.07654 2757.2;  % Tab 3 
               373.156 5.07010 1861.3;  % Tab 3
               373.151 3.56755 1256.4;  % Tab 3 
               398.156 7.79652 2753.7;  % Tab 3
               398.151 5.52235 1858.9;  % Tab 3
               398.147 3.85735 1254.8;  % Tab 3
               ]; 
data_holste = [data_holste Mw_co2*data_holste(:,3)]; % add rho_kg/m3
% Density data from Fenghour et al., J Chem Thermodyn (1995)
data_fenghour = [347.54 6.927 3118.41
                 348.27 7.126 3223.15 
                 348.61 11.600 6771.45
                 374.30 7.842 3113.03
                 375.26 14.068 6761.57];
data_fenghour = [data_fenghour Mw_co2*data_fenghour(:,3)]; % add rho_kg/m3
% Density data from Ely et al., J Chem Thermodyn (1989)
data_ely = [310.000 5.43213 2.99001
            310.000 7.12164 4.98887
            310.000 7.86433 6.85824
            310.000 8.06822 7.83874
            310.000 8.20823 8.82938
            310.000 8.31622 9.87864
            310.000 8.44698 11.24963
            310.000 8.62861 12.59945
            310.000 9.02537 14.05004
            310.000 10.37431 15.95375
            310.000 14.10470 17.94770
            310.000 14.09659 18.04778
            310.000 22.01753 19.83996
            320.000 5.77614 2.98536
            320.000 7.79098 4.97857
            320.000 8.87599 6.84063
            320.000 9.26616 7.81680
            320.000 9.59111  8.80289
            320.000 9.89966 9.84698
            320.000 10.31430 11.21009
            320.000 10.82856 12.54955
            320.000 11.67866 13.98434
            320.000 13.81104 15.87138
            320.000 18.74969 17.90012
            320.000 18.73590 17.89656
            320.000 28.10439 19.80854
            330.000 6.11169 2.97716
            330.000 8.43934 4.96743
            330.000 9.85303 6.81963
            330.000 10.42414 7.78929
            330.000 10.94096 8.76818
            330.000 11.46213 9.80457
            330.000 12.18270 11.15849
            330.000 13.05130 12.49189
            330.000 14.37347 13.92777
            330.000 17.32428 15.82941
            330.000 23.43457 17.86949
            330.000 23.41716 17.86591
            330.000 34.18039 19.78174];
data_ely = [data_ely Mw_co2*data_ely(:,3)*1000]; % add rho_kg/m3
% Data from Klimeck et al., J Chem Thermodynamics 2001
% T = 360 K
data_klimeck = [177.676 8.93780
                237.408 10.9306
                273.475 11.9930
                302.223 12.7870
                375.618 14.7057
                479.295 17.5408
                556.172 20.1844
                647.965 24.7640
                686.030 27.4211
                717.383 30.0864];
% Data from Pensado et al., J Supercrit Fluids 2008
% T = 353.15 K
data_pensado = [15 0.4290
                20 0.5946
                25 0.6868
                30 0.7453
                35  0.7885
                40  0.8212
                45 0.8505
                50  0.8748
                55  0.8959
                60 0.9148];
T_d = [300 310 320 330 348.15 353.15 360 373.15 398.15]-273.15;
S = {'m', 'NaCl', 1}; 
P = {'MPa', 'vals', [linspace(0.1,25,120) linspace(25.5,60,60)]};
saltVar = false;
vapH2O = false;
out_d = cell(numel(T_d),1);
for t = 1:numel(T_d)
    T = {'C', T_d(t)};
    out_d{t} = pvtBrineWithCO2BlackOil(T, P, S, saltVar, vapH2O, false, false);
end

% 2. Viscosity of pure CO2
%                P    280K  300   320   340   360   380
data_fengv =    [0.1 14.05 15.02 15.98 16.93 17.87 18.79
                 0.5 14.09 15.06 16.02 16.96 17.90 18.82
                 1.0 14.15 15.11 16.07 17.01 17.94 18.86
                 2.5 14.51 15.41 16.32 17.23 18.13 19.03
                 5.0 90.41 16.72 17.24 17.95 18.73 19.54
                 7.5 96.86 60.47 19.78 19.48 19.85 20.44
                 10.0 102.42 71.13 32.58 22.80 21.8 21.86
                 15.0 111.94 83.74 60.11 40.23 30.29 27.05
                 20.0 120.07 93.06 71.74 54.76 42.46 35.29
                 25.0 127.85 101.08 80.39 64.45 52.36 43.90
                 30.0 134.98 108.29 87.68 72.09 60.18 51.34
                 40.0 148.28 120.65 100.20 84.55 72.56 63.35
                 50.0 160.73 132.55 111.64 95.21 82.77 73.15
                 60.0 172.63 143.55 121.62 104.80 91.94 81.79];
data_fengv(:,2:end) = data_fengv(:,2:end)*1e-3; % convert to cP

%                P   280K  300   320   340   360  380
data_vesovic = [0.1 14.05 15.02 15.98 16.93 17.87 18.79
                0.5 14.08 15.05 16.01 16.96 17.89 18.81
                1.0 14.14 15.11 16.06 17.00 17.93 18.85
                2.5 14.48 15.39 16.30 17.21 18.12 19.01
                5.0 94.82 16.67 17.19 17.91 18.70 19.51
                7.5 101.30 61.00 19.69 19.41 19.80 20.38
                10.0 106.82 72.44 32.40 22.70 21.75 21.79
                12.5 111.72 79.74 50.96 29.95 25.07 23.92
                15.0 116.20 85.55 60.04 40.06 30.14 26.93
                20.0 124.30 94.96 71.88 54.77 42.33 35.14
                25.0 131.64 102.81 80.72 64.71 52.40 43.81
                30.0 138.47 109.76 88.18 72.62 60.46 51.40
                40.0 151.15 122.08 101.03 85.58 73.36 63.86
                50.0 163.02 133.15 112.33 96.63 84.09 74.17
                60.0 174.40 143.47 122.75 106.64 93.69 83.31];
data_vesovic(:,2:end) = data_vesovic(:,2:end)*1e-3;

T_v = (280:20:380)-273.15;
P = {'MPa', 'vals', [linspace(0.1,25,50) linspace(25.5,60,20)]};
saltVar = false;
vapH2O = false;
out_v = cell(numel(T_v),1);
for t = 1:numel(T_v)
    T = {'C', T_v(t)};
    out_v{t} = pvtBrineWithCO2BlackOil(T, P, S, saltVar, vapH2O, false, false);
end

latx = {'interpreter','latex'};
col_co2  = repelem([20, 0, 0], numel(T_d), 1);
col_co2(:,1) = [(20:30:250) 255]/255;
h=figure(23);
tiledlayout(1,4,'TileSpacing','compact','Padding','compact')
nexttile(1)        % ---------
hold on
pd1 = plot(data_holste(1:6,2)*10,data_holste(1:6,end), 'o', 'markerEdgeColor', col_co2(1,:), ...
          'displayname', 'Holste (1987)');
plot(data_holste(7:13,2)*10,data_holste(7:13,end), 'o', 'markerEdgeColor', col_co2(2,:))
plot(data_holste(14:21,2)*10,data_holste(14:21,end), 'o', 'markerEdgeColor', col_co2(3,:))
plot(data_holste(22:25,2)*10,data_holste(22:25,end), 'o', 'markerEdgeColor', col_co2(5,:))
plot(data_holste(26:28,2)*10,data_holste(26:28,end), 'o', 'markerEdgeColor', col_co2(8,:))
plot(data_holste(29:31,2)*10,data_holste(29:31,end), 'o', 'markerEdgeColor', col_co2(9,:))
pd2 = plot(data_fenghour(1:3,2)*10,data_fenghour(1:3,end), 's', 'markerEdgeColor', col_co2(5,:), ...
          'displayname', 'Fenghour (1995)');
plot(data_fenghour(4:5,2)*10,data_fenghour(4:5,end), 's', 'markerEdgeColor', col_co2(8,:));

pd3 = plot(data_ely(1:13,2)*10,data_ely(1:13,end), '*', 'markerEdgeColor', col_co2(2,:), ...
           'displayname', 'Ely (1989)');
plot(data_ely(14:26,2)*10,data_ely(14:26,end), '*', 'markerEdgeColor', col_co2(3,:))
plot(data_ely(27:39,2)*10,data_ely(27:39,end), '*', 'markerEdgeColor', col_co2(4,:))
pd4 = plot(data_klimeck(:,2)*10,data_klimeck(:,1), '>', 'markerEdgeColor', col_co2(7,:), ...
           'displayname', 'Klimeck (2001)');
pd5 = plot(data_pensado(:,1)*10,data_pensado(:,2)*1000, 'p', 'markerEdgeColor', col_co2(6,:), ...
           'displayname', 'Pensado (2008)');
p2 = plot(out_d{1}.P_bar,out_d{1}.gas.rho_g_kgm3, 'color', col_co2(1,:), 'linewidth', 1, ...
          'displayname', [num2str(T_d(1)) ' ºC']); 
p3 = plot(out_d{2}.P_bar,out_d{2}.gas.rho_g_kgm3, 'color', col_co2(2,:),'linewidth', 1, 'displayname', num2str(T_d(2)));
p4 = plot(out_d{3}.P_bar,out_d{3}.gas.rho_g_kgm3, 'color', col_co2(3,:),'linewidth', 1, 'displayname', num2str(T_d(3)));
p5 = plot(out_d{4}.P_bar,out_d{4}.gas.rho_g_kgm3, 'color', col_co2(4,:),'linewidth', 1, 'displayname', num2str(T_d(4)));
p6 = plot(out_d{5}.P_bar,out_d{5}.gas.rho_g_kgm3, 'color', col_co2(5,:),'linewidth', 1, 'displayname', num2str(T_d(5)));
p7 = plot(out_d{6}.P_bar,out_d{6}.gas.rho_g_kgm3, 'color', col_co2(6,:),'linewidth', 1, 'displayname', num2str(T_d(6)));
p8 = plot(out_d{7}.P_bar,out_d{7}.gas.rho_g_kgm3, 'color', col_co2(7, :),'linewidth', 1, 'displayname', num2str(T_d(7)));
p9 = plot(out_d{8}.P_bar,out_d{8}.gas.rho_g_kgm3, 'color', col_co2(8, :),'linewidth', 1, 'displayname', num2str(T_d(8)));
p10 = plot(out_d{9}.P_bar,out_d{9}.gas.rho_g_kgm3, 'color', col_co2(9, :),'linewidth', 1, 'displayname', num2str(T_d(9)));
hold off
grid on;
ylim([0 1100]); xlim([0 600])
yticks(0:100:1100); xticks(0:100:600)
ylabel('$\rho$ [kg/m$^3$]','fontsize', 14, latx{:})
xlabel('$p$ [bar]', 'fontsize', 14, latx{:})
leg = legend([pd1 pd3 pd2 pd4 pd5 p2 p3 p4 p5 p6 p7 p8 p9 p10], 'fontsize', 12,'location','best');
leg.BoxFace.ColorType='truecoloralpha';
leg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
nexttile(2)        % ---------
hold on
pd1 = plot(data_holste(1:6,2)*10,data_holste(1:6,end), 'o', 'markerEdgeColor', col_co2(1,:), ...
          'displayname', 'Holste (1987)');
plot(data_holste(7:13,2)*10,data_holste(7:13,end), 'o', 'markerEdgeColor', col_co2(2,:))
plot(data_holste(14:21,2)*10,data_holste(14:21,end), 'o', 'markerEdgeColor', col_co2(3,:))
plot(data_holste(22:25,2)*10,data_holste(22:25,end), 'o', 'markerEdgeColor', col_co2(5,:))
plot(data_holste(26:28,2)*10,data_holste(26:28,end), 'o', 'markerEdgeColor', col_co2(8,:))
plot(data_holste(29:31,2)*10,data_holste(29:31,end), 'o', 'markerEdgeColor', col_co2(9,:))
pd2 = plot(data_fenghour(1:3,2)*10,data_fenghour(1:3,end), 's', 'markerEdgeColor', col_co2(5,:), ...
          'displayname', 'Fenghour (1995)');
plot(data_fenghour(4:5,2)*10,data_fenghour(4:5,end), 's', 'markerEdgeColor', col_co2(8,:));

pd3 = plot(data_ely(1:13,2)*10,data_ely(1:13,end), '*', 'markerEdgeColor', col_co2(2,:), ...
           'displayname', 'Ely (1989)');
plot(data_ely(14:26,2)*10,data_ely(14:26,end), '*', 'markerEdgeColor', col_co2(3,:))
plot(data_ely(27:39,2)*10,data_ely(27:39,end), '*', 'markerEdgeColor', col_co2(4,:))
pd4 = plot(data_klimeck(:,2)*10,data_klimeck(:,1), '>', 'markerEdgeColor', col_co2(7,:), ...
           'displayname', 'Klimeck (2001)');
pd5 = plot(data_pensado(:,1)*10,data_pensado(:,2)*1000, 'p', 'markerEdgeColor', col_co2(6,:), ...
           'displayname', 'Pensado (2008)');
p2 = plot(out_d{1}.P_bar,out_d{1}.gas.rho_g_kgm3, 'color', col_co2(1,:), 'linewidth', 1, ...
          'displayname', [num2str(T_d(1)) ' ºC']); 
p3 = plot(out_d{2}.P_bar,out_d{2}.gas.rho_g_kgm3, 'color', col_co2(2,:),'linewidth', 1, 'displayname', num2str(T_d(2)));
p4 = plot(out_d{3}.P_bar,out_d{3}.gas.rho_g_kgm3, 'color', col_co2(3,:),'linewidth', 1, 'displayname', num2str(T_d(3)));
p5 = plot(out_d{4}.P_bar,out_d{4}.gas.rho_g_kgm3, 'color', col_co2(4,:),'linewidth', 1, 'displayname', num2str(T_d(4)));
p6 = plot(out_d{5}.P_bar,out_d{5}.gas.rho_g_kgm3, 'color', col_co2(5,:),'linewidth', 1, 'displayname', num2str(T_d(5)));
p7 = plot(out_d{6}.P_bar,out_d{6}.gas.rho_g_kgm3, 'color', col_co2(6,:),'linewidth', 1, 'displayname', num2str(T_d(6)));
p8 = plot(out_d{7}.P_bar,out_d{7}.gas.rho_g_kgm3, 'color', col_co2(7, :),'linewidth', 1, 'displayname', num2str(T_d(7)));
p9 = plot(out_d{8}.P_bar,out_d{8}.gas.rho_g_kgm3, 'color', col_co2(8, :),'linewidth', 1, 'displayname', num2str(T_d(8)));
p10 = plot(out_d{9}.P_bar,out_d{9}.gas.rho_g_kgm3, 'color', col_co2(9, :),'linewidth', 1, 'displayname', num2str(T_d(9)));
hold off
grid on;
ylim([0 350]); xlim([0 125])
yticks(0:50:350); xticks(0:25:125)
nexttile(3)        % ---------
hold on
pv1 = plot(data_fengv(:,1)*10,data_fengv(:,2), 'o', 'color', col_co2(1,:), ...
           'displayname', 'Fenghour (1998)');
plot(data_fengv(:,1)*10,data_fengv(:,3), 'o', 'markerEdgeColor', col_co2(2,:));
plot(data_fengv(:,1)*10,data_fengv(:,4), 'o', 'markerEdgeColor', col_co2(3,:));
plot(data_fengv(:,1)*10,data_fengv(:,5), 'o', 'markerEdgeColor', col_co2(5,:));
plot(data_fengv(:,1)*10,data_fengv(:,6), 'o', 'markerEdgeColor', col_co2(7,:));
plot(data_fengv(:,1)*10,data_fengv(:,7), 'o', 'markerEdgeColor', col_co2(9,:));
pv12 = plot(data_vesovic(:,1)*10,data_vesovic(:,2), '+', 'color', col_co2(1,:), ...
           'displayname', 'Vesovic (1990)');
plot(data_vesovic(:,1)*10,data_vesovic(:,3), '+', 'markerEdgeColor', col_co2(2,:));
plot(data_vesovic(:,1)*10,data_vesovic(:,4), '+', 'markerEdgeColor', col_co2(3,:));
plot(data_vesovic(:,1)*10,data_vesovic(:,5), '+', 'markerEdgeColor', col_co2(5,:));
plot(data_vesovic(:,1)*10,data_vesovic(:,6), '+', 'markerEdgeColor', col_co2(7,:));
plot(data_vesovic(:,1)*10,data_vesovic(:,7), '+', 'markerEdgeColor', col_co2(9,:));
pv2 = plot(out_v{1}.P_bar,out_v{1}.gas.mu_co2_cP, 'color', col_co2(1,:), 'linewidth', 1, ...
          'displayname', [num2str(T_v(1)) ' ºC']); 
pv3 = plot(out_v{2}.P_bar,out_v{2}.gas.mu_co2_cP, 'color', col_co2(2,:), 'linewidth', 1, 'displayname', num2str(T_v(2))); 
pv4 = plot(out_v{3}.P_bar,out_v{3}.gas.mu_co2_cP, 'color', col_co2(3,:), 'linewidth', 1, 'displayname', num2str(T_v(3)));
pv5 = plot(out_v{4}.P_bar,out_v{4}.gas.mu_co2_cP, 'color', col_co2(5,:), 'linewidth', 1, 'displayname', num2str(T_v(4)));
pv6 = plot(out_v{5}.P_bar,out_v{5}.gas.mu_co2_cP, 'color', col_co2(7,:), 'linewidth', 1, 'displayname', num2str(T_v(5)));
pv7 = plot(out_v{6}.P_bar,out_v{6}.gas.mu_co2_cP, 'color', col_co2(9,:), 'linewidth', 1, 'displayname', num2str(T_v(6)));
grid on;
%ylim([0 1100]); yticks(0:100:1100); 
xlim([0 600]); xticks(0:100:600)
ylabel('$\mu$ [cP]','fontsize', 14, latx{:})
%xlabel('$p$ [bar]', 'fontsize', 14, latx{:})
leg = legend([pv12 pv1 pv2 pv3 pv4 pv5 pv6 pv7], 'fontsize', 12,'location','best');
leg.BoxFace.ColorType='truecoloralpha';
leg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
nexttile(4)        % ---------
hold on
plot(data_fengv(:,1)*10,data_fengv(:,2), 'o', 'color', col_co2(1,:));
plot(data_fengv(:,1)*10,data_fengv(:,3), 'o', 'markerEdgeColor', col_co2(2,:));
plot(data_fengv(:,1)*10,data_fengv(:,4), 'o', 'markerEdgeColor', col_co2(3,:));
plot(data_fengv(:,1)*10,data_fengv(:,5), 'o', 'markerEdgeColor', col_co2(5,:));
plot(data_fengv(:,1)*10,data_fengv(:,6), 'o', 'markerEdgeColor', col_co2(7,:));
plot(data_fengv(:,1)*10,data_fengv(:,7), 'o', 'markerEdgeColor', col_co2(9,:));
plot(data_vesovic(:,1)*10,data_vesovic(:,2), '+', 'color', col_co2(1,:));
plot(data_vesovic(:,1)*10,data_vesovic(:,3), '+', 'markerEdgeColor', col_co2(2,:));
plot(data_vesovic(:,1)*10,data_vesovic(:,4), '+', 'markerEdgeColor', col_co2(3,:));
plot(data_vesovic(:,1)*10,data_vesovic(:,5), '+', 'markerEdgeColor', col_co2(5,:));
plot(data_vesovic(:,1)*10,data_vesovic(:,6), '+', 'markerEdgeColor', col_co2(7,:));
plot(data_vesovic(:,1)*10,data_vesovic(:,7), '+', 'markerEdgeColor', col_co2(9,:));
plot(out_v{1}.P_bar,out_v{1}.gas.mu_co2_cP, 'color', col_co2(1,:), 'linewidth', 1); 
plot(out_v{2}.P_bar,out_v{2}.gas.mu_co2_cP, 'color', col_co2(2,:), 'linewidth', 1, 'displayname', num2str(T_v(2))); 
plot(out_v{3}.P_bar,out_v{3}.gas.mu_co2_cP, 'color', col_co2(3,:), 'linewidth', 1, 'displayname', num2str(T_v(3)));
plot(out_v{4}.P_bar,out_v{4}.gas.mu_co2_cP, 'color', col_co2(5,:), 'linewidth', 1, 'displayname', num2str(T_v(4)));
plot(out_v{5}.P_bar,out_v{5}.gas.mu_co2_cP, 'color', col_co2(7,:), 'linewidth', 1, 'displayname', num2str(T_v(5)));
plot(out_v{6}.P_bar,out_v{6}.gas.mu_co2_cP, 'color', col_co2(9,:), 'linewidth', 1, 'displayname', num2str(T_v(6)));
grid on;
ylim([0 0.05]); yticks(0:0.01:0.05); 
xlim([0 125]); xticks(0:25:125)
ylabel('$\mu$ [cP]','fontsize', 14, latx{:})
set(h, 'Position', [600, 600, 1000, 350])
%exportgraphics(h,'untitled.pdf','ContentType','vector')