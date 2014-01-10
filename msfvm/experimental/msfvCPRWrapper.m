function varargout = msfvCPRWrapper(fn, varargin)
    persistent itercount oldsat;
    itercount = 0;

    varargout = cell(nargout, 1);
    [varargout{:}] = fn(varargin{:});


    state = varargout{1};
    if isempty(oldsat);
        oldsat = state.s;
    end

    if varargout{2}.converged
        itercount = itercount + 1;
%         updateTargets = sum(abs(oldsat - state.s) > 0.2, 2) > 0;
%         if any(updateTargets)
%             msfvPressureSolve(updateTargets);
%         end


%         if itercount == 1
          if state.T > 10*day

            % Reset solver
            disp('Converged, resetting fv basis!')
            msfvPressureSolve();
            itercount = 0;
            state.T = 0;
        end
    end
end
