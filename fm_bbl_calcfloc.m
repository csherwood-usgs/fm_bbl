%% fm_bbl_calcfloc
f_dt=dt;
Gval = G(iz);
dtmin=dt;
dttemp=0.0;

cvtotmud=sum(cv_tmp);
NNin=cv_tmp./mass;
%fprintf(1,'z, G, cvtotmud: %f %f %f\n',zc(iz), Gval, cvtotmud)
if( any (NNin<0.0) )
   fprintf(1,'***************************************\n')
   fprintf(1,'CAUTION, negative mass at t = %f\n', t(n))
   fprintf(1,'***************************************\n') 
end

if (cvtotmud > f_clim)
   while (dttemp <= dt)
      fm_bbl_comp_fsd % NNin -> NNout
      % fm_mass_control
      ineg = find(NNout<0.0);
      mneg = sum( -NNout(ineg).*mass(ineg) );
      
      %     fprintf(1, 'mneg',mneg
      if (mneg > f_mneg_param)
         while (mneg > f_mneg_param)
            f_dt=min(f_dt/2.0,dt-dttemp);
            fm_bbl_comp_fsd % NNin -> NNout
            % fm_mass_control
            ineg = find(NNout<0.0);
            mneg = sum( -NNout(ineg).*mass(ineg) );
         end
      else       
         if (f_dt<dt)
            while (mneg <=f_mneg_param)
               
               if (dttemp+f_dt == dt)
                  fm_bbl_comp_fsd % NNin -> NNout
                  break
               else
                  dt1=f_dt;
                  f_dt=min(2.0*f_dt,dt-dttemp);
                  fm_bbl_comp_fsd % NNin -> NNout
                  % fm_mass_control
                  ineg = find(NNout<0.0);
                  mneg = sum( -NNout(ineg).*mass(ineg) );
                  if (mneg > f_mneg_param)
                     f_dt=dt1;
                     fm_bbl_comp_fsd % NNin -> NNout
                     break
                  end
               end
            end
         end
      end
      
  %    if (f_dt<10); fprintf(1,'fdt=%f3.0\n',f_dt); end
      
      dtmin = min(dtmin,f_dt);
      dttemp = dttemp+f_dt;
      NNin = NNout; % update new Floc size distribution
      
     while any(NNin<0.0); fm_bbl_mass_distribute; end % redistribute negative masses in NNin (if any) over positive classes,
      % depends on f_mneg_param
      
      if (abs(sum(NNin.*mass)-cvtotmud)>epsilon*10000.0)
         fprintf(1, 'CAUTION flocculation routine not conservative!\n')
         fprintf(1, 'time = %f\n',t);
         fprintf(1, 'f_dt= %f\n',f_dt);
         fprintf(1, 'before : cvtotmud= %f\n',cvtotmud);
         fprintf(1, 'after  : cvtotmud= %f\n',sum(NNin.*mass));
         fprintf(1, 'absolute difference  : cvtotmud= %f\n',abs(cvtotmud-sum(NNin.*mass)));
         fprintf(1, 'absolute difference reference  : epsilon= %f\n',epsilon);
         fprintf(1, 'before redistribution %f\n', sum(NNout.*mass));
         fprintf(1, 'after redistribution %f\n', sum(NNin.*mass));
         error('Simultation stopped')
      end
      
      if (dttemp == dt); break; end
   end % loop on full dt
end
if (abs( sum( NNin.*mass )-cvtotmud) > epsilon*10000.0)
   fprintf(1, 'CAUTION flocculation routine not conservative!\n');
   fprintf(1, 'time = %g\n',t);
   fprintf(1, 'before : cvtotmud= %f\n',cvtotmud)
   fprintf(1, 'after  : cvtotmud= %f\n',sum( NNin.*mass ))
   fprintf(1, 'absolute difference  : cvtotmud= %f\n',...
      abs(cvtotmud-sum( NNin.*mass )))
   fprintf(1, 'absolute difference reference  : espilon= %f\n',epsilon);
   error('Simultation stopped')
end

% update mass concentration for all mud classes
cv_wat  = NNin.*mass;

if(0)
   % compute floc distribution statistics before output
   f_csum=0.0;
   f_ld50=1;
   f_ld10=1;
   f_ld90=1;
   
   f_davg = sum(NNin.*mass.*Df)./(sum(NNin.*mass)+eps);
   f_dtmin = dtmin;
   
   for iv1=1:np
      f_csum=f_csum + NNin(iv1)*mass(iv1)/((sum(NNin.*mass))+eps);
      if (f_csum > 0.1 && f_ld10)
         f_d10 =Df(iv1);
         f_ld10 = 0;
      end
      
      if (f_csum > 0.5 && f_ld50)
         f_d50 = Df(iv1);
         f_ld50=0;
      end
      
      if (f_csum > 0.9 && f_ld90)
         f_d90=Df(iv1);
         f_ld90=0;
      end
   end
   f_d_area_weighted=(NNin.*f_area)'*Df/(sum(NNin.*f_area));
   fprintf(fid,'%f %f %f %f %f\n',...
      t, f_dt, Gval, f_d50*1e6, f_d_area_weighted*1e6);
end