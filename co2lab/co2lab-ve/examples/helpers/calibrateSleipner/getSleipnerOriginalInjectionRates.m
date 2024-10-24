function [ inj_year, inj_rates ] = getSleipnerOriginalInjectionRates()

    % The following data corresponds to 32 years of CO2 entry into Sleipner
    % Layer 9. Data is cummulative. Retrieved from the Preliminary Sleipner
    % benchmark model (in-house from Statoil).
    
    % Surface density was reported to be 1.87 kg/m3
    % Reservoir density was reported to be 695.89 kg/m3

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


CumSurfVol_m3 = [0;0;9726800.00000003;29532000.0000000;50765000.0000001; ...
    77537600.0000000;113961600.000000;164148800.000000;232211000; ...
    322260000.000000;438407600.000000;584765600.000000;765445800.000000; ...
    984560000.000000;1246220000.00000;1554537600.00000;1913624600.00000; ...
    2327592800.00000;2800554000.00000;3336620000.00000;3939902600.00000; ...
    4614513600.00000;5364564800.00000;6194168000.00000;7107435000.00000; ...
    8108477600.00000;9201407600.00000;10390336800.0000;11679377000.0000; ...
    13072640000.0000;14574237600.0000;16188281600.0000;17918883800.0000; ...
    19770156000.0000];

CumMass_kg = [0;0;18189116.0000001;55224840.0000001;94930550.0000001; ...
    144995312.000000;213108192.000000;306958256.000000;434234570; ...
    602626200.000000;819822212.000000;1093511672.00000;1431383646.00000; ...
    1841127200.00000;2330431400.00000;2906985312.00000;3578478002.00000; ...
    4352598536.00000;5237035980.00000;6239479400.00000;7367617862.00000; ...
    8629140432.00000;10031736176.0000;11583094160.0000;13290903450.0000; ...
    15162853112.0000;17206632212.0000;19429929816.0000;21840434990.0000; ...
    24445836800.0000;27253824312.0000;30272086592.0000;33508312706.0000; ...
    36970191720.0000];

CumReservoirVol_m3 = [0;0;26137.9183491645;79358.5767865612;136416.028395293; ...
    208359.528086336;306238.330770668;441101.691359267;623998.864763109; ...
    865979.105893173;1178091.66966044;1571385.81097587;2056910.78475046; ...
    2645715.84589519;3348850.24932101;4177363.24993893;5142304.10265990; ...
    6254722.06239492;7525666.38405495;8966186.32255098;10587331.1327940; ...
    12400150.0696949;14415692.3881648;16645007.3431146;19099144.1894552; ...
    21789152.1820977;24726080.5759531;27920978.6259323;31384895.5869462; ...
    35128880.7139059;39163983.2617224;43501252.4853066;48151737.6395695; ...
    53126487.9794220];

year_Jan1st = [1998;1999;2000;2001;2002;2003;2004;2005;2006;2007;2008; ...
    2009;2010;2011;2012;2013;2014;2015;2016;2017;2018;2019;2020;2021; ...
    2022;2023;2024;2025;2026;2027;2028;2029;2030;2031];

    % NOTES:
    %
    % 1998 is first year reported. Cummulative CO2 injected by Jan 1st,
    % 1998 was 0. Cummulative CO2 injected by Jan 1st 1999 was 0. Thus
    % there was no injection taking place in 1998.

    % Cummulative CO2 injected by Jan 1st, 2000 was 9e6 m3. Thus
    % injection rate for 1999 was (9e6 - 0) m3/year.

    % Cummulative CO2 injected by Jan 1st, 2001 was 29.5e6 m3. Thus
    % injection rate for 2000 was (29.5e6 - 9e6) m3/year.

    % Cummulative CO2 injected by Jan 1st, 2030 was 17918e6 m3, and
    % cummulative CO2 injected by Jan 1st, 2031 was 19770e6 m3. Thus
    % injection rate for 2030 was (19770e6 - 17918e6) m3/year.

    % No cummulative CO2 volume is reported for Jan 1st, 2032, thus we
    % cannot determine the injection rate for the year 2031.

    % Thus the injection years (that contain non-zero rates) are 1999 to
    % 2030, a total of 32 years.
    
    
    % Injection years (with non-zero rate):
    inj_year = year_Jan1st(2:end-1);
    
    
    % Injection rates, computed using cummulative reservoir volume:
    inj_rates          = zeros(numel(inj_year),1);
    CumReservoirVol_m3 = CumReservoirVol_m3(2:end);
    
    for i = 1:numel(CumReservoirVol_m3)-1
        inj_rates(i,1) = ( CumReservoirVol_m3(i+1) - CumReservoirVol_m3(i) ); % annual rate
    end
    
    inj_rates = inj_rates * meter^3/year; % in meter^3/second


end

