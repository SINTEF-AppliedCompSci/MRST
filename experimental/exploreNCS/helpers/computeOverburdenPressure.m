function [ P_over, P_limit ] = computeOverburdenPressure( Gt, rock2D, ...
    seafloor_depth, water_density, varargin )
% Compute overburden or fracture pressure at a particular depth in the
% Norwegian Continental Shelf.

    % Rock types of some formations:
    %   Garn      - sandstone, medium to coarse grain, interlayered by
    %               mica-rich zones.
    %   Not       - separates Garn from Ile: claystone
    %   Are/Tilje - alternating layers of sandstone, coal, claystones.
    %   Ror       - seal for Are/Tilje: mudstones, some silty and sandy
    %               sequences

    % Dry bulk density of sandstone ~= 2 grams/cm^3
    % @@ comes from online reference: Johnson and Olhoeft, In: Handbook of
    % Physical Properties of Rocks (1984).
    rhoR_sandstone = (2 * gram)/(centi * meter)^3; 

    P_surf = 1 * atm;
    g      = norm(gravity());

    % sea property:
    rhoW    = 1000 * kilogram/meter^3;
    H1      = seafloor_depth * meter;

    % rock between sea bottom to top of formation: %@@ could be passed in as
    % varargin, otherwise use formation properties
    poro2 = mean(rock2D.poro);
    rhoW2 = water_density;
    rhoR2 = rhoR_sandstone;
    H2    = Gt.cells.z - seafloor_depth;


    % rock between top of formation to center of cell (depth of injection):
    poro3 = rock2D.poro;
    rhoW3 = water_density;
    rhoR3 = rhoR_sandstone;
    H3    = Gt.cells.z + Gt.cells.H./2;


    %% Compute over pressure: forces exerted by overlying rock and water:
    P_over = P_surf + ...
                rhoW*g*H1 + ...
                (poro2 .* rhoW2 + (1 - poro2) .* rhoR2) .* g .* H2 + ...
                (poro3 .* rhoW3 + (1 - poro3) .* rhoR3) .* g .* H3;

    %% Compute a limit that the formation pressure should not exceed:

    % take 90 percent of the minimum overpressure to be conservative 
    P_limit = 0.9*(min(P_over)); % Pascals


end

