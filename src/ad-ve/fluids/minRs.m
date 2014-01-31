function min_rs= minRs(p,sG,sGmax,f, G)
    drho=norm(gravity)*(f.rhoOS.*f.bO(p)-f.rhoGS.*f.bG(p));
    pcmax=f.pcOG(sGmax, p,'sGmax',double(sGmax));       
    h_max=pcmax./drho;
    assert(all(double(h_max)>=0));
    %h_max=min(h_max,G.cells.H);
    ind=double(h_max)>G.cells.H;
    h_max(ind)=G.cells.H(ind)+h_max(ind)*0;
    min_rs=((1-sG)-(G.cells.H-h_max)./G.cells.H).*f.dis_max;   
    %{
    pcOG=f.pcOG(sG,p,'sGmax',double(sGmax));
    h=pcOG./drho;
    assert(all(double(h)>=0));
    h_max(h<=0)=sG(h<=0).*G.cells.H(h<=0)/(f.res_gas);
    h(h<=0)=0*h(h<=0);
    min_rs1=(1-f.res_gas).*(h_max-h).*f.dis_max...% rs in the oil/water sone
                +(f.res_oil).*h.*f.dis_max;
    min_rs1=min_rs1/G.cells.H;        
    %}
    %assume all water above h_max is saturated
    %min_rs=((1-sG)-(G.cells.H-h_max)./G.cells.H).*f.dis_max;   
    %min_rs=min_rs./G.cells.H;
end