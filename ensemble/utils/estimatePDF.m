function pdf = estimatePDF(mean, median, std, scale)
    gamma3 = 3*(mean-median)/std;
    gamma3 = sign(gamma3)*min(0.99, abs(gamma3));
%     gamma3 = -0.8
    a     = abs(gamma3)^(2/3);
    delta = sign(gamma3).*sqrt(pi/2*a/(a+((4-pi)/2)^(2/3)));
    alpha = delta/sqrt(1-delta^2);
    
    std = std/sqrt((1 - 2*delta^2/pi));
    mean = mean - std*delta*sqrt(2/pi);
    
    xi = @(x) (x-mean)./std;
    
    pdf = @(x) 2/(std*sqrt(2*pi))*exp(-0.5*xi(x).^2).*0.5.*(1 + erf(alpha*xi(x)/sqrt(2))).*scale;
    
end