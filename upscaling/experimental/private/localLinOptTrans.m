function upscaled = localLinOptTrans(upscaled, state, CG, W_cg, fluid, opt)
   %T_cg=upscale.trans;% transmissibilities of coarse grid
   %W_cg=W_cg_w;% well optimization not done

   f_press  = state.pressure;% fine cale pressure
   cg_press = upscaled.pressure;% original cg_pressur
   cg_flux  = upscaled.flux;% cg_flux

   % prepare to make condition
   internal = all(CG.faces.neighbors > 0, 2);
   cg_idp   = cg_press(CG.faces.neighbors(internal,2)) - ...
              cg_press(CG.faces.neighbors(internal,1));

   tok = (cg_idp .* cg_flux(internal)) <= 0;

   % find max and min in a coarse cell
   p_max = accumarray(CG.partition, f_press, size(cg_press), @max, NaN);
   p_min = accumarray(CG.partition, f_press, size(cg_press), @min, NaN);

   % optimization is only done on internal faces
   % make maping from faces to internal faces
   nfi     = sum(internal);
   ncg     = CG.cells.num;
   tmp     = 1 : CG.faces.num;
   if_to_f = tmp(internal);

   assert(all(cg_press-p_min>=-eps*abs(p_min) & cg_press-p_max<=eps*abs(p_max)))

   ff = cg_flux(internal);

   NN = -(sparse(1:nfi,CG.faces.neighbors(internal,2),sign(ff),nfi,ncg)-...
      sparse(1:nfi,CG.faces.neighbors(internal,1),sign(ff),nfi,ncg));

   NN_f=-(sparse(1:nfi,CG.faces.neighbors(internal,2),1,nfi,ncg)-...
      sparse(1:nfi,CG.faces.neighbors(internal,1),1,nfi,ncg));

   p_scale=max(abs(cg_press));

   if mrstVerbose,
      cvx_quiet(0);
   else
      cvx_quiet(1);
   end

   well_cells=[W_cg.cells];
   switch opt.opt_trans_method
      case 'linear_simple'
         c = ones(ncg,1);
         cc = ones(nfi,1);
         cvx_begin
         variable p(ncg)
         variable y(ncg)
         variable fs(nfi)
         %dual variables pmin pmax ptrans yeq yp yn
         dual pmax pmin yep yen yp ptrans pfs pwell
         minimize( c'*(y) + 5e6*cc'*fs )
         subject to
         yep : y >= p -  cg_press/p_scale ;%#ok
         yen : y >= cg_press/p_scale - p ;%#ok
         %p <= p_max
         ptrans : NN*p+fs >= 0;%#ok
         pfs  : fs >=0;%#ok
         pmin : p <= p_max/p_scale;%#ok
         pmax : p >= p_min/p_scale;%#ok
         pwell  : p(well_cells) == cg_press(well_cells)/p_scale;
         cvx_end

         tmp = cg_flux(internal)./(NN_f*(p*p_scale));
         tmp(~isfinite(tmp))=max(upscaled.trans(internal));
         %tmp(NN_f*(p*p_scale)==0)=max(upscaled.trans(internal));
         %tmp(find(abs(NN_f*(p*p_scale))<1e-10*max(abs((NN_f*(p*p_scale))))))=max(upscaled.trans(internal));
         upscaled.trans(internal)=fluid.properties()*tmp;
         upscaled.pressure=p*p_scale;

      case 'convecs_simple'
         %cvx_solver
         %cvx_perecision low
         scale_f=max(abs(cg_flux));
         cvx_begin
         variable p(ncg)
         variable Tinv(nfi)
         variable cgflux(nfi)
         dual ptrans pcon ptinv incomp pwell
         minimize( 1e5*norm(cgflux) + norm(p-cg_press/p_scale) )
         subject to
         pcon : p_min/p_scale <= p <= p_max/p_scale;%#ok
         %ptrans : NN*p == Tinv.*abs(cg_flux(internal))/scale_f
         ptrans : NN_f*p - Tinv.*cg_flux(internal)/scale_f == cgflux;%#ok
         ptinv : Tinv >=0;%#ok
         pwell  : p(well_cells) == cg_press(well_cells)/p_scale;
         cvx_end
         upscaled.trans(internal)= fluid.properties()./(Tinv*p_scale/scale_f);
         upscaled.pressure=p*p_scale;

      case 'convecs_new'
         % cvx_quiet(0);
         %cvx_solver
         %cvx_perecision low
         scale_f=max(abs(cg_flux));

         cvx_begin
         variable p(ncg)
         variable Tinv(nfi)
         variable cgflux(nfi)
         dual ptrans pcon ptinv incomp pwell
         %minimize(max(Tinv)+1e5*norm(cgflux))
         minimize(-min(Tinv)+1e5*norm(cgflux))
         %minimize(norm(p-cg_press/p_scale))
         %maximize(min(NN*p))
         %minimize(max(1./Tinv)+1e5*norm(cgflux)/scale_f )
         %maximize(min(Tinv*p_scale./abs(cg_flux(internal))))
         subject to
         pcon : p_min/p_scale <= p <= p_max/p_scale;%#ok
         %ptrans : NN*p == Tinv.*abs(cg_flux(internal))/scale_f
         ptrans : NN_f*p - Tinv.*cg_flux(internal)/scale_f == cgflux;%#ok
         ptinv : Tinv >=0;%#ok
         pwell  : p(well_cells) == cg_press(well_cells)/p_scale;
         cvx_end
         upscaled.trans(internal)= fluid.properties()./(Tinv*p_scale/scale_f);
         upscaled.pressure=p*p_scale;

      case 'nonlinear'
         error('Not properly implemented')
         cellNo  = rldecode(1:CG.cells.num, diff(CG.cells.facePos), 2) .';
         fcg=CG.faces.num;
         CC= sparse(cellNo,CG.cells.faces(:,1),2*(cellNo==CG.faces.neighbors(CG.cells.faces(:,1),1))-1,ncg,fcg);
         CC=CC(:,internal);
         % inteded to minimize flux error with incompressible
         % constraint
         %incomp : CC*cgflux==0;
         %ptrans : NN*p == Tinv.*abs(cg_flux(internal))/scale_f

      otherwise
         error('No such implementation done');
   end

   cg_press_new=p*p_scale;
   %
   cg_idp_new=cg_press_new(CG.faces.neighbors(internal,2))-cg_press_new(CG.faces.neighbors(internal,1));
   tok_new=cg_idp_new.*cg_flux(internal)<=0;
   %
   disp(['Internal faces ', num2str(numel(tok))])
   disp(['Internal not not ok ', num2str(sum(NN*cg_press<0))])
   disp(['Internal faces new ',  num2str(numel(tok_new))])
   disp(['Internal not not ok new ', num2str(sum(NN*cg_press_new<0))])

   % Prepare to analyse optimisation
   if opt.plot_opt_grid,
      df_new=cg_press_new(CG.partition)-state.pressure;
      df_old=cg_press(CG.partition)-state.pressure;
      cc_new=nan(ncg,1);
      cc_old=nan(ncg,1);
      df_cc_new=inf(ncg,1);
      df_cc_old=inf(ncg,1);
      for i=1:CG.parent.cells.num
         if(abs(df_new(i))<df_cc_new(CG.partition(i)))
            df_cc_new(CG.partition(i))=abs(df_new(i));
            cc_new(CG.partition(i))=i;
         end
         if(abs(df_old(i))<df_cc_old(CG.partition(i)))
            df_cc_old(CG.partition(i))=abs(df_old(i));
            cc_old(CG.partition(i))=i;
         end
      end

      cf_old=if_to_f(find(NN*cg_press<0));%#ok
      cf_new=if_to_f(find(NN*cg_press_new<0));%#ok
      %mcolon(CG.faces.connPos(cf_old),CG.faces.connPos(cf_old+1)-1);
      ff_old=CG.faces.fconn(mcolon(CG.faces.connPos(cf_old),CG.faces.connPos(cf_old+1)-1));
      fc_old=CG.parent.faces.neighbors(ff_old,:);
      ff_new=CG.faces.fconn(mcolon(CG.faces.connPos(cf_new),CG.faces.connPos(cf_new+1)-1));
      fc_new=CG.parent.faces.neighbors(ff_new,:);

      figure(33); cla
      %outlineCoarseGrid(CG.parent,CG.partition,'LineWidth',1,'Color','k'); box on; axis equal tight;

      %
      cc_changed=find(~(cc_new==cc_old));
      plotGrid(CG.parent,cc_new(cc_changed),'FaceColor','r')
      plotGrid(CG.parent,cc_old(cc_changed),'FaceColor','b')
      plotGrid(CG.parent,cc_old,'FaceColor','g')
      plotGrid(CG.parent,fc_old(:))
      plotGrid(CG.parent,fc_new(:),'FaceColor','r')
      drawnow;
      %{
            gg={'uniform','uniform_ct', 'uniform_fa', ...
    'tof_ka','tof_ka_h1','tof_ka_h2', ...
    'flow_tof','flow_tof_cs','flow_tof_facies'};
      %}
      disp('hei')
   end
end
