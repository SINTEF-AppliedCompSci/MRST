function plot_EosPRSW(pres,xliqExp,xliqsw,xliqpr)
% Function which plot the hydrogen or dioxyde carbone solubility 
% from experiment, SW EoS and PR EoS.

patm = 1e5; % Atmospheric pressure in Pa
presbar=pres./patm;
min_xliq=[min(xliqExp),min(xliqsw),min(xliqpr)];
max_xliq=[max(xliqExp),max(xliqsw),max(xliqpr)];

f21=figure('Name','solubility_SwPr','NumberTitle','off');
f21.Position(3:4) = [900 700];

plot(presbar,xliqsw,'k*','MarkerSize',7,'LineWidth',2)
hold on
plot(presbar,xliqpr,'r square','MarkerSize',8,'LineWidth',2)
hold on
plot(presbar,xliqExp,'b o','MarkerSize',8,'LineWidth',2)

title('solubility in pure water PR vs SW','FontSize',16,'FontWeight','normal','Color','k')
xlabel({'pressure (bar)'},'FontWeight','normal','Color','k')
ylabel('molar fraction (-)','FontWeight','normal','Color','k')
ax = gca;
ax.FontSize = 16; 

legend({'SW','PR','Exp'},...
    'FontSize',16,'TextColor','black',...
    'Location','best')
xlim([min(presbar)-10 max(presbar)+10])
ylim([min(min_xliq)-1e-4 max(max_xliq)+1e-4])

end
