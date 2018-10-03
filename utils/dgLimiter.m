function state = dgLimiter(disc, state, bad, type, varargin)

    opt = struct('plot', false);
    opt = merge_options(opt, varargin{:});
    
    nDofMax = disc.basis.nDof;
    G       = disc.G;
    sdof    = state.sdof;
    
    
    if opt.plot
        figure(1); clf; hold on
        plotSaturationDG(disc, state, 'n', 500, 'plot1d', true, 'color', 'b', 'linew', 2);
    end

    switch type
            
         case 'kill'
            % Simple "limiter" that reduces to dG(0) for bad cells
            
            % Fix saturation so that it is inside [0,1]
            ix         = disc.getDofIx(state, 1, bad);
            sdof(ix,:) = min(max(state.s(bad,:), 0), 1);
            % Ensure sum(s,2) = 1
            sdof(ix,:) = sdof(ix,:)./sum(sdof(ix,:),2);
            % Remove dofs for dofNo > 1
            if disc.degree > 0
                ix         = disc.getDofIx(state, 2:nDofMax, bad);
                sdof(ix,:) = [];
                state.degree(bad) = 0;
            end
            % Assign new dofs to state
            state.sdof = sdof;

        case 'tvb'
            % TVB limiter
            
            dofbar = approximatGradient(sdof(:,1), state, disc);
            dofbar = dofbar(bad,:)';
            dofbar = dofbar(:);
            
            ix = disc.getDofIx(state, (1:G.griddim)+1, bad);
            sdof(ix,:) = dofbar.*[1,-1];
            
            if disc.degree > 1
                ix = disc.getDofIx(state, (G.griddim+2):nDofMax, bad);
                sdof(ix,:) = [];
            end
            state.sdof = sdof;
            state.degree(bad) = 1;
            
        case 'scale'
            
            s = state.s;
            [sMin, sMax] = disc.getMinMaxSaturation(state);
            theta = [(s(:,1) - 0)./(s(:,1) - sMin), ...
                     (1 - s(:,1))./(sMax - s(:,1)), ...
                     ones(G.cells.num,1)          ];
            theta(abs(theta) > 1) = 1;
            theta = min(theta, [], 2);
            
            for dofNo = 2:nDofMax%(1:G.griddim)+1
                ix = disc.getDofIx(state, dofNo, (1:G.cells.num)', true);
                sdof(ix(ix>0),:) = sdof(ix(ix>0),:).*theta(ix>0);
            end
            
            satIx = state.degree > 0;
            dofIx = disc.getDofIx(state, 1, satIx);
            sdof(dofIx,:) = (sdof(dofIx,:) - s(satIx,:)).*theta(satIx) + s(satIx,:);
            state.sdof = sdof;
            [sMin, sMax] = disc.getMinMaxSaturation(state);
            if any(sMin < 0 - eps | sMax > 1 + eps)
                warning('Some saturations still outside [0,1]')
            end
            
            ind = theta < 1;
            if disc.degree > 1 && 0
                ix = disc.getDofIx(state, (G.griddim+2):nDofMax, ind);
                sdof(ix,:) = [];
                state.degree(ind & state.degree > 1) = 1;
%                 sdof(ix,:) = 0;
            end
            state.sdof = sdof;
            state.scaled = ind;
            
        case 'orderReduce'
            
            [sMin, sMax] = disc.getMinMaxSaturation(state);
            outside = sMin < 0 - disc.outTolerance | ...
                      sMax > 1 + disc.outTolerance;
                  
            bad = outside & state.degree > 1;
                  
            if any(bad)
                ix = disc.getDofIx(state, (G.griddim + 2):nDofMax, bad);
                sdof(ix,:) = [];
                state.degree(bad) = 1;
%                 sdof(ix,:) = 0;
                state.sdof = sdof;
            end

    end
    
    % Update state dofPos and satuarion
    state   = disc.updateDofPos(state);
    state.s = disc.getCellSaturation(state);
    
    if opt.plot
        plotSaturationDG(disc, state, 'n', 500, 'plot1d', true, 'color', 'r', 'linew', 4, 'LineStyle', '--'); hold off
        legend({['Before ', type], ['After ', type]});
        drawnow
    end

end

function dofbar = approximatGradient(dof, state, disc)

%     ind  = 1:disc.basis.nDof:disc.G.cells.num*disc.basis.nDof;
    dofIx   = disc.getDofIx(state, 1);
    q    = dof(dofIx);
    G    = disc.G;
    nDof = disc.basis.nDof;
    
    sigma = cell(1, disc.dim);
    [sigma{:}] = deal(0);
     
    for d = 1:disc.dim
        b = disc.interp_setup.tri_basis{d};
        for l = 1:disc.dim+1
            
            cells = disc.interp_setup.C(:, l);
%             
%             cells = disc.interp_setup
            
            ccl = disc.interp_setup.tri_cells(cells);
            ds = b(:, l).*q(ccl);
            sigma{d} = sigma{d} + ds;
        end
    end
    
    dofbar = zeros(G.cells.num, disc.G.griddim);
    for cNo = 1:G.cells.num
        
%         gradix = gradpos(cNo):gradpos(cNo+1)-1;
        gradIx = disc.interp_setup.cell_support{cNo};
        for dNo = 1:G.griddim
            dofIx = disc.getDofIx(state, 1+dNo, cNo);
            ss = sigma{dNo}(gradIx);
%             ss = ss(ss~=0);
            dd = dof(dofIx);
            val = [dd; ss];
            db = minmod(val);
            dofbar(cNo, dNo) = db;
        end

    end
    
end

function v = minmod(val)
    v = max(val)*all(val < 0) + min(val)*all(val > 0);
end

%             while any(bad) && it < 4
% 
%                 [sMin, sMax] = disc.getMinMaxSaturation(state);
%                 under = sMin < 0 - disc.outTolerance;
%                 over  = sMax > 1 + disc.outTolerance;
%                 
%                 
%                 
%                 bad   = bad & (over | under) & state.degree > 0;
%                 under = under & bad;
%                 over  = over  & bad;
%                     
%                 if any(under)
%                     badCells = find(under);
%                     ix = disc.getDofIx(state, (1:G.griddim) + 1, badCells);
%                     sd = reshape(sdof(ix), G.griddim, [])';
%                     [~, i] = max(abs(sd), [], 2);
%                     ix0 = disc.getDofIx(state, 1, under);
%                     s0  = state.sdof(ix0,:);
%                     for dofNo = (1:G.griddim) + 1
%                         ii = i == (dofNo - 1);
%                         ix = disc.getDofIx(state, dofNo, badCells(ii));
%                         sdof(ix,1) = -s0(ii,1)/G.griddim;
%                         sdof(ix,2) = -sdof(ix,1);
%                     end
% %                     state.degree(badCells) = 1;
%                 end
% 
%                 if any(over)
%                     badCells = find(over);
%                     ix = disc.getDofIx(state, (1:G.griddim) + 1, badCells);
%                     sd = reshape(sdof(ix), G.griddim, [])';
%                     [~, i] = max(abs(sd), [], 2);
%                     ix0 = disc.getDofIx(state, 1, over);
%                     s0  = state.sdof(ix0,:);
% 
%                     for dofNo = (1:G.griddim) + 1
%                         ii = i == (dofNo - 1);
%                         ix = disc.getDofIx(state, dofNo, badCells(ii));
%                         sdof(ix,1) = (1-s0(ii,1))/G.griddim;
%                         sdof(ix,2) = -sdof(ix,1);
%                     end
% %                     state.degree(badCells) = 1;
%                 end
%                 
%                 if disc.degree > 1
%                     ix = disc.getDofIx(state, (G.griddim+2):nDofMax, badCells);
%                     sdof(ix,:) = [];
%                 end
%                 state.sdof = sdof;
% %                 state = disc.updateDofPos(state);
%                 
%                 [sMin, sMax] = disc.getMinMaxSaturation(state);
%                 it = it+1;
%             end
%             
%             state = disc.getCellSaturation(state);