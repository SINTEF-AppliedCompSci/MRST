function fn = grauePcAD(sr, sr_tot, model)
% Twophase capillary pressure function handle for AD style fluids us
%
% INPUT:
%    sr - residual water saturation
%    sr_tot - sum of sr for all phases present
% DESCRIPTION:
% The capillary pressure-saturation curve is approximated as a curve fit to 
% experimental data from Graue SPE56672. The fit is based on the following data:    
%
%   sw_pc	= [0.15	0.2	0.26	0.35	0.41	0.45	0.5	0.58	0.6	0.65	0.7	0.73	0.76];
%   pc_imb_sim = [60	29	19	12	7	5	4	3	2	2	0	-1	-11] * bar_psi;
%   Swe = (sw_pc-Swr)/(1-Swr-Sor); 
%   pc = polyfit(Swe, pc_imb_sim, 3);	% -16.8901098927788	30.4949345867430	-18.0265406816741	3.78772984965509
%   pc_val = polyval(pc, Swe);

    fn = @(s) grauePc(s, sr, sr_tot, model);
end

function pc = grauePc(s, sr, sr_tot, model)
    den = 1 - sr_tot;
    sat = ((s - sr)./den);
    if isa(sat, 'ADI')
        sat.val = max(min(sat.val, 1), 0);
    else
        sat = max(min(sat, 1), 0);
    end
    
    p = [-16.8901098927788, 30.4949345867430, -18.0265406816741, 3.78772984965509];   % in bars
    pc = (p(1)*sat.^3 + p(2)*sat.^2 + p(3)*sat + p(4)) * 1e5;   % in Pa
    
    % Set Pc=0 in the cells, adjacent to the outflow boundary
    %pc(model.G.iout) = 0;
    
end
