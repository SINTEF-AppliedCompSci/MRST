function [w, ih] = getRestartWellInfo(rstrt, repStep)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

ih = getINTEHEAD(rstrt.INTEHEAD{repStep});

if ih.nwell > 0
    w = getZWEL([], rstrt.ZWEL{repStep}, ih);
    w = getIWEL(w, rstrt.IWEL{repStep}, ih);
    if isfield(rstrt, 'SWEL')
        w = getSWEL(w, rstrt.SWEL{repStep}, ih);
    end
    w = getXWEL(w, rstrt.XWEL{repStep}, ih);
    w = getICON(w, rstrt.ICON{repStep}, ih);
    if isfield(rstrt, 'SCON')
        w = getSCON(w, rstrt.SCON{repStep}, ih);
    end
    w = getXCON(w, rstrt.XCON{repStep}, ih);
else
    w = [];
end
end

    
function w = getZWEL(w, zwel, ih)
d = @(k)(k-1)*ih.nzwel;
for k = 1:ih.nwell
    w(k).name = zwel{1 + d(k)};
end
end

function w = getIWEL(w, iwel, ih)
d = @(k)(k-1)*ih.niwel;
for k = 1:ih.nwell
    w(k).ijk  = iwel((1:3) + d(k))';
    w(k).ncon = iwel( 5 + d(k));
    w(k).ginx = iwel( 6 + d(k));
    w(k).type = iwel( 7 + d(k));
    w(k).cntr = iwel( 9 + d(k));
    w(k).stat = iwel(11 + d(k))>0;
end
end

function w = getSWEL(w, swel, ih)
d = @(k)(k-1)*ih.nswel;
for k = 1:ih.nwell
    w(k).depth = swel(10 + d(k));
end
end

function w = getXWEL(w, xwel, ih)
d = @(k)(k-1)*ih.nxwel;
for k = 1:ih.nwell
    w(k).qOs  = - xwel( 1 + d(k));
    w(k).qWs  = - xwel( 2 + d(k));
    w(k).qGs  = - xwel( 3 + d(k));
    w(k).lrat = - xwel( 4 + d(k));
    w(k).qr   = - xwel( 5 + d(k));
    w(k).bhp  =   xwel( 7 + d(k));
end
end

function w = getICON(w, icon, ih) 
d = @(k)(k-1)*ih.nicon*ih.ncwma;
for k = 1:ih.nwell
    dc = (0:(w(k).ncon-1))*ih.nicon;
    cijk = icon(bsxfun(@plus, dc(:), 2:4) + d(k));
    w(k).cijk  = reshape(cijk, [w(k).ncon, 3]);
    w(k).cstat = icon( 6 + dc + d(k));
    w(k).cdir  = icon( 14 + dc + d(k));
end
end

function w = getSCON(w, scon, ih) 
d = @(k)(k-1)*ih.nscon*ih.ncwma;
for k = 1:ih.nwell
    dc = (0:(w(k).ncon-1))*ih.nscon;
    w(k).cwi    = scon( 1 + dc + d(k));
    w(k).cdepth = scon( 2 + dc + d(k));
    w(k).cdiam  = scon( 3 + dc + d(k));
    w(k).ckh    = scon( 4 + dc + d(k));
end
end

function w = getXCON(w, xcon, ih) 
d = @(k)(k-1)*ih.nxcon*ih.ncwma;
for k = 1:ih.nwell
    dc = (0:(w(k).ncon-1))*ih.nxcon;
    w(k).cqs  = - reshape(xcon(bsxfun(@plus, dc(:), 1:3) + d(k)), [], 3);
    w(k).cqr  = - xcon( 50 + dc + d(k));
    w(k).press = xcon(35 + dc + d(k));
end
end

function ih = getINTEHEAD(intehead)
ih.	unit	=	intehead(	3	);
ih.	nx      =	intehead(	9	);
ih.	ny      =	intehead(	10	);
ih.	nz      =	intehead(	11	);
ih.	nactive	=	intehead(	12	);
ih.	iphs	=	intehead(	15	);
ih.	nwell	=	intehead(	17	);
ih.	ncwma	=	intehead(	18	);
ih.	nwgmax	=	intehead(	20	);
ih.	ngmaxz	=	intehead(	21	);
ih.	niwel	=	intehead(	25	);
ih.	nswel	=	intehead(	26	);
ih.	nxwel	=	intehead(	27	);
ih.	nzwel	=	intehead(	28	);
ih.	nicon	=	intehead(	33	);
ih.	nscon	=	intehead(	34	);
ih.	nxcon	=	intehead(	35	);
ih.	nigrpz	=	intehead(	37	);
ih.	iday	=	intehead(	65	);
ih.	imon	=	intehead(	66	);
ih.	iyear	=	intehead(	67	);
ih.	iprog	=	intehead(	95	);
end

