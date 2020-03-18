function plotConvergence(ep,ef,x,mylegend,varargin)
%Plot convergence plots for the given L2 norm of pressure and flux errors
%   ep - L2 norm of pressure error matrix, each column corresponding to one
%   discretization method. Number of rows equal to the length of h;
%   ef - L2 norm of flux error matrix, similar to ep
%   h - mesh size;
%   nc - number of cells
%   mylegend - my legend


if(all(x>1))
    figure, loglog(x,ep,'o-','linewidth',2);
    legend(mylegend,'location','best','fontsize',12);
    xlabel('\itn\rm_c');ylabel('\ite_p\rm');grid on
    
    figure, loglog(x,ef,'o-','linewidth',2);
    legend(mylegend,'location','best','fontsize',12);
    xlabel('\itn\rm_c');ylabel('\ite_f\rm');grid on
else
    x=log2(1./x);
    set(figure,'color','w');
    if(isempty(varargin))
        plot(x,log2(ep),'o-','linewidth',2);
    else
        plot(x,log2(ep),varargin{:});
    end
    legend(mylegend,'location','best','fontsize',12);
    xlabel('log_{2}(1/\ith\rm)');
    ylabel('log_{2}(\ite_p\rm)');grid on
    
    set(figure,'color','w');
    if(isempty(varargin))
        plot(x,log2(ef),'o-','linewidth',2);
    else
        plot(x,log2(ef),varargin{:});
    end
    legend(mylegend,'location','best','fontsize',12);
    xlabel('log_{2}(1/\ith\rm)');
    ylabel('log_{2}(\ite_f\rm)');grid on
end
end

