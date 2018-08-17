function [prop, p] = perturb_rockProps(Gt, prop, interval, varargin)
% Add a perturbation to rock properties of formation

% interval of form: [low high]

    opt.pert_type       = 'gaussian';
    opt.pert_clevel     = 1;        % coarsening level to apply to 
                                    % perturbations. Default is no
                                    % coarsening
    opt.kernal_size     = 3; % default, can pass in [3,3,3], [1,3,1], etc.
    opt.std             = 0.65;
    
    opt = merge_options(opt, varargin{:});

    if numel(opt.kernal_size)==1
       opt.kernal_size = [opt.kernal_size, opt.kernal_size, opt.kernal_size]; 
    end

    % -----------------------------
    % Generate a gaussian field
    % -----------------------------
    % for the rock property with dimensions [m,n] equal to Gt.cartDims
    % (cell dimensions)

    cell_dims = Gt.cartDims;

    % Generate (a possibly coarsened) p field (then use interp to get p of
    % size Gt.cartDims)
    N = opt.pert_clevel;
    coarse_cell_dims = ceil(Gt.cartDims./N);
    if strcmpi(opt.pert_type,'gaussian')
        pc = gaussianField([coarse_cell_dims], interval, opt.kernal_size(1:2), opt.std);
    else
       error('Undefined perturbation type.') 
    end

    [xc,yc] = meshgrid(1:coarse_cell_dims(2), 1:coarse_cell_dims(1));
    r = 1;
    % r becomes smaller if we require a smaller step size to ensure enough
    % resolution for grid
    while 1
        dx = coarse_cell_dims(1)/cell_dims(1) * r;
        dy = coarse_cell_dims(2)/cell_dims(2) * r;
        [x,y] = meshgrid(1:dy:coarse_cell_dims(2), 1:dx:coarse_cell_dims(1));
        if size(x,1) >= cell_dims(1) && size(y,2) >= cell_dims(2)
            break
        else
            r = r*0.95;
        end
    end
    p = interp2(xc,yc,pc,x,y,'cubic');
    assert(all(size(p) >= cell_dims))

    
    % -----------------------------
    % Add the gaussian field values to the cell prop values 
    % -----------------------------
    % This requires finding which index of p(i,j) corresponds to each
    % cell

    % loop through grid cells:
    for c = 1:Gt.cells.num
        
        % ij index of this cell:
        [i,j] = deal(Gt.cells.ij(c,1), Gt.cells.ij(c,2));
        
        % add field value (i.e., p(i,j)) to prop value:
        prop( c ) = prop( c ) + p(i,j);
    end
    
    
    % -----------------------------
    % Handle possible values that fall below zero
    % -----------------------------
    assert(all(prop > 0))

end