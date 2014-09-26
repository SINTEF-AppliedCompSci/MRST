function [qW, qO, qG, bWw, bOw, bGw, Rw, mobWw, mobOw, mobGw] = simpleWellFlux_BO(fluid, state, W, mob, rho, b, Q, pBHP, p, rs)




    [Tw, dzw, Rw, wc, perf2well, pInx, iInxW, iInxO, iInxG] = getWellStuff(W);

    qWs  = Q{1};     qOs = Q{2};     qGs = Q{3}; 
    bW   = b{1};      bO = b{2};      bG = b{3};
    mobW = mob{1};  mobO = mob{2};  mobG = mob{3};
    
    rhoW = rho{1}; rhoO = rho{2}; rhoG = rho{3}; 
%     rhow = [double(rho{1}(wc)), double(rho{2}(wc)),double(rho{3}(wc))];
    
    verySimple = true;
    % connection mobilities
    bWw = bW(wc);
    bOw = bO(wc);
    bGw = bG(wc);
    
    if verySimple
        grav = norm(gravity(), 2);
        
        pw = p(wc);

        
        mobWw = mobW(wc);
        mobOw = mobO(wc);
        mobGw = mobG(wc);
        
        nPerf = arrayfun(@(x)numel(x.cells), W)';
        inj   = vertcat(W.sign)==1;
        injComp = rldecode(vertcat(W(inj).compi), nPerf(inj));
        iInx  = rldecode(inj, nPerf);
        
        perfInj = wc(inj(perf2well));
         if ~isempty(injComp)
            totMobInj = mobW(perfInj) + mobO(perfInj) + mobG(perfInj);
            
            mobWw(iInx) = injComp(:, 1).*totMobInj;
            mobOw(iInx) = injComp(:, 2).*totMobInj;
            mobGw(iInx) = injComp(:, 3).*totMobInj;
            
         end
        drawdown = pBHP(perf2well) - pw;
        qW = mobWw.*Tw.*(drawdown + grav.*dzw.*rhoW(wc));
        qO = mobOw.*Tw.*(drawdown + grav.*dzw.*rhoO(wc));
        qG = mobGw.*Tw.*(drawdown + grav.*dzw.*rhoG(wc));
        
        
    else
        if(~isempty(W))
            rhow = cellfun(@(x) double(x(wc)), rho, 'UniformOutput', false);
            rhow = horzcat(rhow{:});
            [Hw, alpha]   = computeWellHead(W, state.wellSol, rhow);
            drawdown = pBHP(perf2well) - p(wc) + Hw;
            isInj    = double(drawdown) > 0;
            rsw = rs(wc);
            [mobWw, mobOw, mobGw, crossFlow] = ...
                computeWellMobilities(W, qWs, qOs, qGs, mobW(wc), mobO(wc), mobG(wc), ...
                                      bWw, bOw, bGw, rsw, isInj);

            if any(crossFlow)
               fprintf('Crossflow in %2.0d connections\n', nnz(crossFlow));
            end
        else
            mobWw=0*bOw;
            mobOw=0*bOw;
            mobGw=0*bOw;
        end
        qW = mobWw.*Tw.*drawdown;
        qO = mobOw.*Tw.*drawdown;
        qG = mobGw.*Tw.*drawdown;
%         disp(num2str([double(mobWw), double(mobOw), double(mobGw)]))
%         fprintf('\n')
    end

end


function [mW, mO, mG, crossFlow] = computeWellMobilities(W, qWs, qOs, qGs, mW, mO, mG, bW, bO, bG, rs, isInj)
    % for producer producing connections, m = [mW mO mG] remains unchanged
    % for producer injecting connections the (assumed uniform) flow in the
    % well-bore is calculated by
    %   qW = qWs/bW
    %   qO = qOs/bO
    %   qG = (qGs-rs*qOs)/bG
    % and the mobility as m = [qW qO qG]*(mW+mO+mG)/(qW+qO+qG)
    %
    % for injector injecting connections the mobility is calculated by
    % m = compi*(mW+mO+mG)
    % injector producing connections is treated as above, but supresses a
    % warning

    nPerf = arrayfun(@(x)numel(x.cells), W)';
    inj   = vertcat(W.sign)==1;
    iInx  = rldecode(inj, nPerf);
    % update crossflow flag
    %check for injector producing connections
%     if any(and(iInx, ~isInj))
%         warning('Crossflow detected in injectors, no special treatment for this case')
%     end
    injComp = rldecode(vertcat(W(inj).compi), nPerf(inj));
    if ~isempty(injComp)
        mtInj   = mW(iInx) + mO(iInx) + mG(iInx);
        mW(iInx) = injComp(:,1).*mtInj;
        mO(iInx) = injComp(:,2).*mtInj;
        mG(iInx) = injComp(:,3).*mtInj;
    end
    crossFlow = and(~iInx,isInj);
    if any(crossFlow) && 0
        ix = find(crossFlow); 
        perf2well = rldecode((1:numel(W))', nPerf);

        qW = qWs(perf2well(ix))./bW(ix);
        qO = qOs(perf2well(ix))./bO(ix);
        qG = (qGs(perf2well(ix))-rs(ix).*qOs(perf2well(ix)))./bG(ix);
        qt = qW + qO + qG;

        % don't bother with connections having almost zero flow
    %     zeroInx = abs(qt.val)<1e-6/day;
    %     ix = ix(~zeroInx);
        if ~isempty(ix)
            fprintf('Crossflow detected in %2.0d connections in wells ', numel(ix));
            fprintf('%d ', unique(perf2well(ix)))
            fprintf('\n')
            mt = mW(ix) + mO(ix) + mG(ix);
            mW(ix) = qW.*mt./qt;
            mO(ix) = qO.*mt./qt;
            mG(ix) = qG.*mt./qt;

        end
    end
end

function [Hw, alpha] = computeWellHead(W, wellSol, rho)
    inx = 0;
    for wnr = 1:numel(W)
        if ~isfield(wellSol(wnr), 'flux') || isempty(wellSol(wnr).flux)
            wellSol(wnr).flux = ones(numel(W(wnr).cells), 1)*W(wnr).compi*W(wnr).sign;
        end
    end

    rhoMix = cell(numel(W),1);
    alpha  = cell(numel(W),1);

    % the following could be replaced by e.g. a linear system for multilateral
    % wells
    for wnr = 1:numel(W)
        wbOut = [0 0 0];
        alpha{wnr}  = zeros(numel(W(wnr).cells), 3);
        rhoMix{wnr} = zeros(numel(W(wnr).cells), 1);
        for cnr = numel(W(wnr).cells):-1:1
            fluxOut = wellSol(wnr).flux(cnr,:);
            if cnr == 1 && W(wnr).sign==1
                wbIn    = sum(wbOut + fluxOut)*W(wnr).compi;
            else
                wbIn    = wbOut + fluxOut;
            end
            totIn   = (wbIn    > 0).*wbIn  + ...
                      (wbOut   < 0).*wbOut + ...
                      (fluxOut < 0).*fluxOut;
            if abs(sum(totIn)) < 1e-6/day
                totIn = W(wnr).compi;
            end
            a = totIn/sum(totIn);
            alpha{wnr}(cnr,:)  = a;
            rhoMix{wnr}(cnr) = sum(a.*rho(inx+cnr,:));
            wbOut = wbIn;
        end
        inx = inx + numel(W(wnr).cells);
    end
    rhoMixAvg = cellfun(@(x)[x(1); .5*(x(1:end-1)+x(2:end))], rhoMix, 'UniformOutput', false);

    g   = norm(gravity);
    dzw = arrayfun(@(x)[x.dZ(1); x.dZ(2:end)-x.dZ(1:end-1)], W, 'UniformOutput', false);
    ddp = cellfun(@(x,y)g*x.*y, dzw, rhoMixAvg, 'UniformOutput', false);
    Hw  = cellfun(@cumsum, ddp, 'UniformOutput', false);
    Hw = vertcat(Hw{:});
    alpha = vertcat(alpha{:});
    %dzw = vertcat(W.dZ);
    %dpw = g*dzw.*vertcat(rhoMixAvg{:});
end