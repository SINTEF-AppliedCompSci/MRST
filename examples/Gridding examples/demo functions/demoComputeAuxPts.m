function [pIn, pOut, R] = demoComputeAuxPts(p, bn, m0)
% An excerpt from 'generateVOIGridNodes' to show how to generate the 
% auxiliary points for the Voronoi diagram. The auxiliary points help to
% locate the boundary nodes in p, i.e. p(bn, :).

    pib = p(bn, :);
    pib = [pib; pib(1,:)];
    n   = size(bn,1);
    e2n = [(1:n)', [(2:n)';1]];
    % Compute the radius
    edges = [bn, [bn(2:end); bn(1)]];
    L  = p(edges(:,1),:) - p(edges(:,2),:);
    L  = sqrt(sum(L.^2,2));
    do = true;
    m  = m0;
    while do
        if m > 0.5
            throwError(L)
        end
        try
            m = m + 0.02;
            R = ( [L(1);L(1:n-1)] + [L(n);L(2:n)] ) * m;
            [pIn, pOut] = deal( zeros(size(edges,1), 2) );
            for i = 1 : size(edges,1)
                [x1, y1, x2, y2] = deal(p(edges(i,1),1), p(edges(i,1),2), ...
                    p(edges(i,2),1), p(edges(i,2),2));
                [r1, r2] = deal(R(e2n(i,1)), R(e2n(i,2)));
                pCross = circleCross(x1, y1, r1, x2, y2, r2);
                in = inpolygon(pCross(:,1), pCross(:,2), pib(:,1), pib(:,2));
                pIn(i, :)   = pCross(in, :);
                pOut(i, :)  = pCross(~in, :);
            end
            do = false;
        catch
            do = true;
        end
    end
    
    figure, hold on, axis equal tight off
    demoPlotPoly(p(bn, :), 'ko-', 'y', 5)
    th = linspace(0, 2*pi, 100)';
    for i = 1 : size(edges,1)
        pCri = bsxfun(@plus, [R(i)*cos(th), R(i)*sin(th)], p(bn(i),:));
        demoPlotLine(pCri, 'r-', 'r', 1)
        demoPlotLine([pIn(i, :); pOut(i, :)], 'bs-', 'g', 5)
    end
    legend('WR boundary points','Designed circumcircles of boundary points',....
        'Intersection points of neighboring circles')
end
    

function throwError(L)
    error(['Cannot generate appropriate Voronoi sites, please \n',...
        '   (1) Increase the resolution of well trajectory (add more well points) \n', ...
        'Or (2) Increase the value of ''ly'', ',...
        'the suggested value is %.0f (may require to enlarge the VOI boundary)\n'], 1.1*max(L));
end
