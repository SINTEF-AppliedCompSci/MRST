function jsonstruct = removeJsonStructField(jsonstruct, fdnames, varargin)

    opt = struct('handleMissing', 'warn');
    opt = merge_options(opt, varargin{:});
    fdname = fdnames{1};
    if numel(fdnames) == 1
        if ~isfield(jsonstruct, fdname)
            msgtxt = sprintf('Field %s is missing\n', fdname);
            switch opt.handleMissing
              case 'warn'
                fprintf(msgtxt);
              case 'error'
                error(msgtxt);
              case 'quiet'
                % do nothin
              otherwise
                error('handleMissing case not recognized.');
            end
            return
        end
        jsonstruct = rmfield(jsonstruct, fdname);
        return
    else
        jsonstruct.(fdname) = removeJsonStructField(jsonstruct.(fdname), fdnames(2:end));
    end
    
end

