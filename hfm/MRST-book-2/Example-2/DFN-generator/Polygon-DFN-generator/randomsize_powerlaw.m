function diam=randomsize_powerlaw(mindiam,maxdiam,exponent)
% This function generates random sizes based on power law scaling. The
% inputs are the minimum and maximum fracture diams. Exponent is the power
% law scaling exponent. All sizes are in diameters.

mindiam_a=mindiam^(-exponent+1);
maxdiam_a=maxdiam^(-exponent+1);
diam=(maxdiam_a + rand*(mindiam_a - maxdiam_a))^(1/(1-exponent));


end