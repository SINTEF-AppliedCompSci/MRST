function [updata, report] = upRk(block, updata, method, varargin)
% Upscaling of polymer reduction factor Rk
%
% This method applied the methodology for upscaling of polymer as described
% in [1]. First, the absolute permeability and then relative permeabilities
% are upscaled. Then, the polymer reduction factor Rk is upscaled. The
% method uses an important result from [1], which states that at steady
% state, the polymer concentration is constant. Therefore, we can apply a
% two-phase upscaling strategy.
% 
% To save computation, the Rk upscaling is performed together with the
% relative permeability upscaling. Therefore, do not call the two-phase
% upscaling first.
% 
% References:
% [1] Hilden, Lie and Xavier Raynaud - "Steady State Upscaling of Polymer
%     Flooding", ECMOR XIV - 14th European conference on the mathematics of
%     oil recovery, 2014.
% 
%
opt = struct(...
    'npoly', 15 ...
    );
warning('off', 'Option:Unsupported');
[opt, relPermOpt] = merge_options(opt, varargin{:});
warning('on', 'Option:Unsupported');

watopt = []; % Relperm upscaling options
Rk = []; % Upscaled values
nc = opt.npoly;
cvalues = linspace(0, block.fluid.cmax, nc)';


% Run upscaling of relperm and Rk together
[updata, report] = upRelPerm(block, updata, method, ...
    'fun_setup',   @(varargin) fun_setup(varargin{:}), ...
    'fun_upscale', @(varargin) fun_upscale(varargin{:}), ...
    'fun_final',   @(varargin) fun_final(varargin{:}), ...
    relPermOpt{:});


%--------------------------------------------------------------------------
% NESTED FUNCTIONS
%--------------------------------------------------------------------------

function fun_setup(opt)
    % Called at the beginning of the relative permeability upscaling. We
    % use the options to allocate space and setup the Rk upscaling.
    
    watopt = opt;
    ndims = length(opt.dims);
    Rk = cell(1, ndims);
    for d=1:ndims
        Rk{d} = nan(watopt.nsat, nc);
    end
end

function fun_upscale(block, rock_Kkr, krKU, is, id, opt)
    % Called after upscaling of relative permeability for water for each
    % saturation value and each dimension.
    
    d = opt.dims(id); % Current dimension
    
    if krKU==0
        % Water is immobile, and so we cannot run the upscaling. To have
        % some value of Rk for this saturation value, we assume it is equal
        % to the previous saturation.
        if (is > 1)
            Rk{id}(is,:) = Rk{id}(is-1,:);
            return;
        end
    end
    
    % Polymer Rk upscaling
    for ic = 1:nc
        c = cvalues(ic);

        cstate = c.*ones(block.G.cells.num,1);
        Rkval = 1 + ( block.fluid.rrf - 1 ) .* ...
            ( block.fluid.ads(cstate) ./ block.fluid.adsMax );
        block.rock.perm = bsxfun(@rdivide, rock_Kkr.perm, Rkval);

        % Perform one phase upscaling with the altered
        % permeability field
        data = upAbsPerm(block, [], 'dims', d, ...
            'dp', opt.dp, 'method', opt.absmethod);
        krKURk = data.perm;

        RkU = krKU / krKURk;

        Rk{id}(is,ic) = RkU;
    end

end

function updata = fun_final(updata)
    % Called after upscaling is completed. The upscaled values of Rk are
    % added to the updata structure.
    
    updata.Rk.val = Rk;
    updata.Rk.s = cell(1,numel(Rk));
    for d=1:numel(Rk)
        updata.Rk.s{d} = updata.krW{d}(:,1);
    end
    updata.Rk.c = cvalues;
    
end

end



