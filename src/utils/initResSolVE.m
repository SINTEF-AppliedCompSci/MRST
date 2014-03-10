function resSol = initResSolVE(G, p0, s0, varargin)
    % Wrapper for initResSol which adds any extra properties needed by the
    % vertical-equil module solvers. See resSol for details of valid
    % arguments.
    
    % G = varargin{1};

    opt = struct('use_s_form', false, 'use_ADI', false);
    opt = merge_options(opt, varargin{:});
    
    % Create initial reservoir solution
    resSol = initResSol(G, p0, s0, varargin{:});

    % Height of the plume is assumed to be zero at t=0
    resSol.h     = zeros(size(resSol.s));
    resSol.h_max = resSol.h;
    resSol.extSat= repmat(resSol.s, 1, 2);
    
    if ~opt.use_ADI 
        % outside the ADI framework, we need to provide explicit functions
        % for computing the Jacobian.
        if any(strcmp(G.type, 'topSurfaceGrid'))
            resSol.twophaseJacobian = @twophaseJacobianWithVE_s;
        else
            % TODO: figure out region3D option
            resSol.twophaseJacobian = @twophaseJacobianWithVE_coupled;
        end
    end
end
