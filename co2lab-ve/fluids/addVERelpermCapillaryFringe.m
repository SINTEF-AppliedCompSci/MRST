function fluid = addVERelpermCapillaryFringe(fluid, Gt, rock, invPc3D, kr3D, varargin)

    % type can be 'linear cap.', 'S table', 'P-scaled table' or 'P-K-scaled table'.
    opt = struct('type', 'P-scaled table', 'samples' 2000);
    opt = merge_options(opt, varargin{:});

    % create sampled tables
    table_water = ;
    table_co2 = ;
    
    % create fluid and return
    fluid = tabulatedVERelperm(table_water, table_co2);

end