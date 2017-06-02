function [ P_over ] = computeOverburdenPressure( Gt, rock2D, ...
    seafloor_depth, water_density, varargin )
% Compute overburden pressure acting on top of a formation, which is due to
% the weight of all overlying layers (fluid + solid). Here, we assume
% overlying formation is made up of same rock and fluid type as found in
% storage formation (i.e., rock porosity, water density)


% Inputs:
%   Gt                - 2D grid, as generated using topSurfaceGrid
%   rock2D            - rock structure, containing poro (porosity)
%   seafloor_depth    - depth to seafloor
%   water_density     - density of overlying formation fluid
%
% Optional:
%   seawater_density  - density of seawater
%   rock_density      - dry bulk density of overlying formation rock
%
% Outputs:
%   P_over            - overburden pressure, in Pascals


    opt.surface_pressure = 1*atm;
    opt.seawater_density = 1000 * kilogram/meter^3;
    opt.rock_density = (2 * gram)/(centi * meter)^3;        
    % NB: dry bulk density of sandstone, taken from online reference:
    % Johnson and Olhoeft, Handbook of Physical Properties of Rocks, 1984
    opt = merge_options( opt, varargin{:} );

    gravity on;
    g = norm(gravity());
    
    % Layer 1 properties: sea layer
    H1      = seafloor_depth * meter;

    % Layer 2 properties: rock/fluid between sea bottom to top of formation
    poro2 = mean(rock2D.poro);
    rhoW2 = water_density;
    rhoR2 = opt.rock_density;
    H2    = Gt.cells.z - seafloor_depth;

    % Overburden pressure:
    P_over = opt.surface_pressure + opt.seawater_density * g * H1 + ...
                (poro2 .* rhoW2 + (1 - poro2) .* rhoR2) .* g .* H2;

end

