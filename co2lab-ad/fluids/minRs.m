function min_rs= minRs(p,sG,sGmax,f, G)
    
    % computing |g|.(rho_w - rho_g)
    drho=norm(gravity)*(f.rhoWS.*f.bW(p)-f.rhoGS.*f.bG(p));
    
    % Computing position of the h_max interface, consistent with the current
    % value of sGmax.
    pcmax=f.pcWG(sGmax, p,'sGmax',double(sGmax));       
    h_max=pcmax./drho;
    assert(all(double(h_max)>=0));

    % If the computed value for h_max ever exceeds H, set it to H (while
    % keeping ADI structure if applicable).
    ind=double(h_max)>G.cells.H;
    h_max(ind)=G.cells.H(ind)+h_max(ind)*0;
    
    % Minimal dissolved quantity equals the 'dis_max' value multiplied by the
    % content of brine in the zone above 'h_max' (i.e. residual brine within the
    % CO2 plume, as well as the brine in the residual CO2 zone).  The content of
    % brine in this zone equals the total content of brine in the column (1-sG),
    % minus the content of brine below 'h_max'.
    min_rs=((1-sG)-(G.cells.H-h_max)./G.cells.H).*f.dis_max;   

end

% -------------------------------- DEPRECATED --------------------------------
%     %{
%     pcWG=f.pcWG(sG,p,'sGmax',double(sGmax));
%     h=pcWG./drho;
%     assert(all(double(h)>=0));
%     h_max(h<=0)=sG(h<=0).*G.cells.H(h<=0)/(f.res_gas);
%     h(h<=0)=0*h(h<=0);
%     min_rs1=(1-f.res_gas).*(h_max-h).*f.dis_max...% rs in the oil/water sone
%                 +(f.res_water).*h.*f.dis_max;
%     min_rs1=min_rs1/G.cells.H;        
%     %}

%     %assume all water above h_max is saturated
%     %min_rs=((1-sG)-(G.cells.H-h_max)./G.cells.H).*f.dis_max;   
%     %min_rs=min_rs./G.cells.H;
% end