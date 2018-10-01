function state = dgLimiter(disc, state, bad, type, varargin)

    nDofMax = disc.basis.nDof;
    G       = disc.G;
    sdof    = state.sdof;

    switch type
            
         case 'kill'
            % Simple "limiter" that reduces to dG(0) for bad cells
            
            % Fix saturation so that it is inside [0,1]
            ix         = disc.getDofIx(state, 1, bad);
            sdof(ix,:) = min(max(state.s(bad,:), 0), 1);
            % Ensure sum(s,2) = 1
            sdof(ix,:) = sdof(ix,:)./sum(sdof(ix,:),2);
            % Set saturation
            state.s(bad,:) = sdof(ix,:);
            % Remove dofs for dofNo > 1
            if disc.degree > 0
                ix         = disc.getDofIx(state, 2:nDofMax, bad);
                sdof(ix,:) = [];
            end
%             sdof(ix,:) = 0;
            % Assign new dofs to state
            state.sdof = sdof;
            state.degree(bad) = 0;

        case 'tvb'
            % TVB limiter
            
            dofbar = approximatGradient(sdof(:,1), state, disc);
            dofbar = dofbar(bad,:)';
            dofbar = dofbar(:);
            
            ix = disc.getDofIx(state, (1:G.griddim)+1, bad);
            sdof(ix,:) = dofbar.*[1,-1];
            
%             if disc.degree > 1
%                 ix = disc.getDofIx(state, (G.griddim+2):nDofMax, bad);
%                 sdof(ix,:) = [];
%             end
            state.sdof = sdof;
%             state.degree(bad) = 1;
            
        case 'scale'

            it = 1;
            while any(bad) && it < 4

                [smin, smax] = disc.getMinMaxSaturation(state);
                under = smin < 0 - disc.outTolerance;
                over  = smax > 1 + disc.outTolerance;
                
                bad   = bad & (over | under) & state.degree > 0;
                under = under & bad;
                over  = over  & bad;
                    
                if any(under)
                    badCells = find(under);
                    ix = disc.getDofIx(state, (1:G.griddim) + 1, badCells);
                    sd = reshape(sdof(ix), G.griddim, [])';
                    [~, i] = max(abs(sd), [], 2);
                    ix0 = disc.getDofIx(state, 1, under);
                    s0  = state.sdof(ix0,:);
                    for dofNo = (1:G.griddim) + 1
                        ii = i == (dofNo - 1);
                        ix = disc.getDofIx(state, dofNo, badCells(ii));
                        sdof(ix,1) = -s0(ii,1)/G.griddim;
                        sdof(ix,2) = -sdof(ix,1);
                    end
%                     state.degree(badCells) = 1;
                end

                if any(over)
                    badCells = find(over);
                    ix = disc.getDofIx(state, (1:G.griddim) + 1, badCells);
                    sd = reshape(sdof(ix), G.griddim, [])';
                    [~, i] = max(abs(sd), [], 2);
                    ix0 = disc.getDofIx(state, 1, over);
                    s0  = state.sdof(ix0,:);

                    for dofNo = (1:G.griddim) + 1
                        ii = i == (dofNo - 1);
                        ix = disc.getDofIx(state, dofNo, badCells(ii));
                        sdof(ix,1) = (1-s0(ii,1))/G.griddim;
                        sdof(ix,2) = -sdof(ix,1);
                    end
%                     state.degree(badCells) = 1;
                end
                
                if disc.degree > 1
                    ix = disc.getDofIx(state, (G.griddim+2):nDofMax, badCells);
                    sdof(ix,:) = [];
                end
                state.sdof = sdof;
%                 state = disc.updateDofPos(state);
                
                [smin, smax] = disc.getMinMaxSaturation(state);
                it = it+1;
            end
            
            state = disc.getCellSaturation(state);

    end
    
    state = disc.updateDofPos(state);

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