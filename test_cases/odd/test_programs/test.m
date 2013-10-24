function test
    gravity on;
    [pc, tc] = CO2CriticalPoint;
    pRefMin = 5.0e5;%0.8e7; %1.0e5;
    pRefMax = 7e7;%pc - (pc - pRefMin) * 0.05;%7e7;%1.0e8;%1.0e8;
    tRefMin = 280;
    tRefMax = 400;%tc - (tc - tRefMin) * 0.05;;%400;

    bounds = [tRefMin tRefMax pRefMin pRefMax];
    
    p_res = 150;
    t_res = 150;

    height = 200;
    sensitivity = 0.01;

    CO2 = CO2props;

    tvals = linspace(tRefMin, tRefMax, t_res);
    pvals = linspace(pRefMin, pRefMax, p_res);

    [T_top, P_top] = meshgrid(tvals, pvals);
    
    top_rho   = reshape(CO2.rho     (P_top(:), T_top(:)), p_res, t_res);
    top_der   = reshape(CO2.rhoder  (P_top(:), T_top(:)), p_res, t_res);
    top_beta  = reshape(CO2.beta    (P_top(:), T_top(:)), p_res, t_res);
    top_bder  = reshape(CO2.betader (P_top(:), T_top(:)), p_res, t_res);
    top_bder2 = reshape(CO2.betader2(P_top(:), T_top(:)), p_res, t_res);
    
    buffered = true;
    if (~buffered)
        int_bot_rho   = zeros(size(top_rho));
        int_bot_press = zeros(size(top_rho));
        i = 0;
        j = 0;
        for t = tvals
            j = j + 1;
            for p = pvals
                i = i + 1;
                % computing depth pressure profile at this temperature and top
                % pressure
                bp = numIntP(0, p, height, @CO2.rho, @(z)t);
                
                int_bot_press(i,j) = bp.y(end);
                int_bot_rho  (i,j) = CO2.rho(bp.y(end), t);
            end
            i=0;
        end
    else
        load('int_150_150_200m');
    end

    lterm = linterm(height, top_rho, top_beta);           % linear term of alpha
    sterm = sqrterm(height, top_rho, top_beta, top_bder); % quadratic term of alpha
    %cterm = cubterm(height, top_rho, top_beta, top_bder, top_bder2); % cubic term
    
    lapx_bot_rho = top_rho .* (1 + lterm);
    sapx_bot_rho = top_rho .* (1 + lterm + sterm);
    %capx_bot_rho = top_rho .* (1 + lterm + sterm + cterm);
    apx_bpress   = apx_bot_press(height, P_top, top_rho, top_beta, top_bder, top_bder2);
    
    % -----------------

    % Computing relative differences
    PRelDiffInt = (int_bot_press - P_top)./P_top;
    PRelDiffApx = (apx_bpress - P_top)./P_top;
    
    RhoRelDiffInt  = ( int_bot_rho - top_rho)./top_rho;
    %RhoRelDiffApxC = (capx_bot_rho - top_rho)./top_rho; % equal to 'lterm + sterm + cterm'
    RhoRelDiffApxS = (sapx_bot_rho - top_rho)./top_rho; % equal to `lterm + sterm`
    RhoRelDiffApxL = (lapx_bot_rho - top_rho)./top_rho; % equal to 'lterm'

    % Determining columns where discontinuities have been crossed
    g2l_int = double(CO2.phaseOf(int_bot_press(:), T_top(:)) == 1 & CO2.phaseOf(P_top(:), T_top(:)) == 2);
    g2l_apx = double(CO2.phaseOf(apx_bpress(:), T_top(:)) == 1 & CO2.phaseOf(P_top(:), T_top(:)) == 2);
    g2l_int = reshape(g2l_int, p_res, t_res);
    g2l_apx = reshape(g2l_apx, p_res, t_res);

    g2l = double(g2l_int | g2l_apx);
    
    RhoRelDiffInt(find(g2l)) = NaN;
    RhoRelDiffApxS(find(g2l)) = NaN;
    RhoRelDiffApxL(find(g2l)) = NaN;
    
    close all;
    
    %figure(1); subPlotAnalyze(RhoRelDiffInt, RhoRelDiffApxS, sensitivity, tvals, pvals, bounds);
    %figure(2); subPlotAnalyze(RhoRelDiffInt, RhoRelDiffApxC, sensitivity, tvals, pvals, bounds);
    %figure(2); subPlotAnalyze(RhoRelDiffInt-RhoRelDiffApxL, sterm, sensitivity, tvals, pvals, bounds);

    
    % figure(3); PVplot(g2l, tvals, pvals, bounds); view(0, 90);

    error0 = abs(int_bot_rho - top_rho)./top_rho;
    error1 = abs(int_bot_rho - lapx_bot_rho)./top_rho;
    error2 = abs(int_bot_rho - sapx_bot_rho)./top_rho;
    
    error_map = zeros(p_res, t_res);
    error_map(error0 > sensitivity) = 1;
    error_map(error1 > sensitivity) = 2;
    error_map(error2 > sensitivity) = 3;

    error_map(find(g2l_apx)) = NaN;

    PVplot(error_map, tvals, pvals, bounds);
    
    %error_map_apx(find((abs(lterm)+abs(sterm)) > sensitivity)) = 1;
    % error_map_apx(find(abs(sterm) > sensitivity)) = 2;
    % error_map_apx(find((abs(sapx_bot_rho - int_bot_rho)./top_rho) > sensitivity)) = 3;
    % error_map_apx(find(g2l_apx)) = NaN;
    

    % figure(4); PVplot(error_map_apx, tvals, pvals, bounds); view(0, 90);
    
    CO2.dispose();       
    keyboard;        

end

function subPlotAnalyze(m1, m2, tol, tvals, pvals, bounds)
    subplot(2, 3, 1); PVplot(m1, tvals, pvals, bounds);    
    subplot(2, 3, 2); PVplot(m2, tvals, pvals, bounds);
    subplot(2, 3, 3); PVplot(double(abs(m1 - m2)>tol), tvals, pvals, bounds);
    subplot(2, 3, 4); PVplot(double(abs(m1) > tol), tvals, pvals, bounds); 
    subplot(2, 3, 5); PVplot(double(abs(m2) > tol), tvals, pvals, bounds);
    subplot(2, 3, 6); PVplot(double((abs(m1) > tol) ~= (abs(m2) > tol)), tvals, pvals, bounds);
end

function PVplot(m, tvals, pvals, bounds)
    surf(tvals, pvals, m, 'faceAlpha', 0.5, 'edgeAlpha', 0.1);
    line_props = {'color', [.3 .3 .3], 'lineWidth', 1.5};
    line(tvals, CO2VaporPressure(tvals), line_props{:});
    axis(bounds);
    %shading interp;
    view(0, 90);
end



function x = linterm(dz, rho, beta)
    pstar = rho.*norm(gravity).*dz;
    x     = beta.*pstar;
end

function x = sqrterm(dz, rho, beta, bder)
    pstar = rho.*norm(gravity).*dz;
    x = 1/2 * (2 * beta.*beta + bder) .* (pstar.^2);
end

% function x = cubterm(dz, rho, beta, bder, bder2)
%     pstar = rho.*norm(gravity).*dz;
%     x = 1/6 * ((6 * beta.^3) + (7 * beta.*bder) + bder2) .* (pstar.^3);
%     end
% 
function P_bot = apx_bot_press(dz, p, rho, beta, bder, bder2)

    pstar = norm(gravity) .* dz .* rho;
    
    P_bot = p + ...
            pstar + ...
            1/2 .* beta .* pstar.^2 + ...
            1/6 * (2*beta.^2 + bder) .* pstar.^3 + ...
            1/24 * (6*beta.^3 + (7 * beta.*bder) + bder2) .* pstar.^4;
end


    % max_rel_err_0 = max(abs(int_bot_rho(:) - top_rho(:))./top_rho(:));
    % max_rel_err_1 = max(abs(int_bot_rho(:) - lapx_bot_rho(:))./top_rho(:));
    % max_rel_err_2 = max(abs(int_bot_rho(:) - sapx_bot_rho(:))./top_rho(:));

    % mean_rel_err_0 = mean(abs(int_bot_rho(:) - top_rho(:))./top_rho(:));
    % mean_rel_err_1 = mean(abs(int_bot_rho(:) - lapx_bot_rho(:))./top_rho(:));
    % mean_rel_err_2 = mean(abs(int_bot_rho(:) - sapx_bot_rho(:))./top_rho(:));
    
    % fprintf('%8.5f | %8.5f | %8.5f\n', max_rel_err_0, max_rel_err_1, max_rel_err_2);
    % fprintf('%8.5f | %8.5f | %8.5f\n', mean_rel_err_0, mean_rel_err_1, mean_rel_err_2);
