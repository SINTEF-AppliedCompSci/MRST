%% Illustrate Riemann solutions for various mobility ratios
% This script shows how the Riemann solution for a case with adsorption and
% full mixing changes with viscosity ratio between the polymer solution and
% the oil phase. This can alternatively be seen as an illustration of the
% effect of permeability reduction caused by polymer retention, since this
% will also increase the viscosity ratio between polymer and oil.
[swr, swor] = deal(0.2,0.8);
[nw,no,ads] = deal(2,4,0.2);
mobRatios   = [2 5 10 20 40];

%% Plot fractional flow and solution profile for waterflooding
s=linspace(0,1,501);
sm = (min(max(s,swr),swor)-swr)/(swor-swr);
fw = sm.^nw./(sm.^nw + (1-sm).^no);
clf, subplot(1,2,1)
plot(s,fw,'k--','LineWidth',1);

subplot(1,2,2),
iw  = find(fw(2:end)./(s(2:end)-swr)-diff(fw)./diff(s)>0,1);
i0 = find(s<swr,1,'last')+1;
xiw = fliplr(diff(fw([i0 iw:end]))./diff(s([i0 iw:end])));
ssw = s(end:-1:iw);
xiw = [xiw xiw(end) max(xiw)+.5];
ssw = [ssw, swr swr];
plot(xiw(ssw<=swor),ssw(ssw<=swor),'--k','LineWidth',1);

%% Plot fractional flow and solution profiles for polymer flooding
for M=mobRatios
    % Compute and plot the fractional flow curve with polymer
    fp = sm.^nw./(sm.^nw + M*(1-sm).^no);
    subplot(1,2,1), hold all
    plot(s,fp,'-','LineWidth',2);

    % Find envelope for polymer curve
    i3 = find(fp(2:end)./(s(2:end)+ads)-diff(fp)./diff(s)>0,1);

    % Find intersection of this envelope and the pure water curve
    i2 = find(min(fp(i3)/(s(i3)+ads).*(s(2:end)+ads),1)<fw(2:end),1)+1;

    % Find envelope for pure water curve
    i0 = find(s<swr,1,'last');
    i1 = find(fw(2:i2)./(s(2:i2)-swr)>diff(fw(1:i2))./diff(s(1:i2)),1);
    if isempty(i1), i1=i2; end

    % Plot self-similar solution
    subplot(1,2,2), hold all
    xi = fliplr(diff(fp(i3:end))./diff(s(i3:end)));
    ss = fliplr(s(i3:end-1));
    df = fliplr(diff(fw([i0 i1:i2]))./diff(s([i0 i1:i2])));
    ds = fliplr(s(i1:i2));
    xi = [xi([1:end end]), df];
    ss = [ss, ds([1 1:end])]; %#ok<*AGROW>
    xi = [xi([1:end end]), max(xi)+.5];
    ss = [ss, s(i0), swr];
    plot(xi(ss<=swor),ss(ss<=swor),'-','LineWidth',2);
    axis tight, set(gca,'Ylim',[0 1]);
end
subplot(1,2,1), set(gca,'FontSize',12)
legend(cellfun(@(x) sprintf('M=%2d', x), num2cell([1, mobRatios]),...
    'UniformOutput', false));
subplot(1,2,2), set(gca,'FontSize',12)

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
