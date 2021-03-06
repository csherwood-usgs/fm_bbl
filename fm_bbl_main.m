% fm_bbl_main - No user input needed below here
% Floc model of R. Verney
% Implemented by R. Verney and C.R. Sherwood
%
% Reference: Verney, R., R. Lafite, J C. Brun-Cottan, and P. Le Hir (2011)
%   Behaviour of a floc population during a tidal cycle: Laboratory
%   experiments and numerical modelling.
%   Continental Shelf Research 31:S64-S83
%
%% Initialize floc parameters and coefficients
fm_bbl_initfloc

%% Initalize arrays

totmass_m = sum(sum(massconc)*dzz)
Cm = zeros( nzc, npmud, nt );
Gsave=zeros(nzc, nt);
Cm(:,:,1) = massconc;
totmass_m = sum( sum( Cm(:,:,1) )*dzz)
Cvs = zeros( nzc, nps, nt );
erate = zeros(nt,1);
if (do_wcbbl_model)
   ustc = zeros(nt,1);
   ustw = zeros(nt,1);
   ustcw = zeros(nt,1);
   zoa = zeros(nt,1);
   taub = zeros(nt,1);
end

for n=2:nt
   totmass_m = sum( sum( Cm(:,:,n-1) )*dzz);
   totmass_s = sum( sum( rhoss*Cvs(:,:,n-1) )*dzz);
   fprintf(1,'t = %f, tot. mass sand= %f and tot. mass mud =%f\n',t(n),totmass_s,totmass_m)
   if(do_wcbbl_model)
      m=m94(uw(n),omegaw,uc(n),1,0.,kN,0);
      ustc(n) = [m.ustrc];
      ustw(n) = [m.ustrwm];
      ustcw(n) = [m.ustrr];
      zoa(n) = [m.zoa];
   end
   erate(n) = max(0., E0*Cb*(taub(n)-tauc)/tauc ); 
   
   % Calc bbl mixing and turbulence
   K = vk*ustc(n)*zf(2:end-1);              % K is on faces
   e = 1.0 ./(vk*h) * (1-zc./h)./(zc./h);    % e, and G are on centers
   e = abs(ustc(n)).^3 .* e;
   G = sqrt(e./nu );                         % 1/s shear rate parameter
   
   if(do_wibergK)
      % Pat Wiberg's u*cw(z) profile
      lc= min(h/2, max(ustc(n)./omegac,zoc));
      lw= max(ustw(n)./omegaw,zoc);
      uscfzf2=ustc(n)^2 * exp(-2*zf./lc);
      uswfzf2=ustw(n)^2 * exp(-2*zf./lw);
      usfzf=(uscfzf2.^2 + uswfzf2.^2 + 2.*uscfzf2.*uswfzf2).^0.25;
      % Wiberg K profile
      K = nu + vk*usfzf.*zf;
   end
   if(do_altG)
      % alternative w-c G profile on cell centers
      lc= min(h/2, max(ustc(n)./omegac,zoc));
      lw= max(ustw(n)./omegaw,zoc);
      uscfzc2=ustc(n)^2 *exp(-2*zc./lc);
      uswfzc2=ustw(n)^2 *exp(-2*zc./lw);
      usfzc=(uscfzc2.^2 + uswfzc2.^2 + 2.*uscfzc2.*uswfzc2).^0.25;
      e1 = usfzc.^3 ./(vk*zc);
      G = sqrt(e./nu ); % 1/s shear rate parameter
   end
   Gsave(:,n)=G;
   
   if(do_sand)
      tint = (0:t(n)-t(n-1));
      % Settle sand
      for k=1:nps
         % right now, wss and erate are for only one sand class
         w = wss;
         Erate = erate(n); % passed as global
         Drate = 1;
         C0 = squeeze( Cvs(:,k,n-1) );
         [t2,C2]=ode45('resus_settle_leer',tint,C0); % general ode equation solver
         %[t2,C2]=ode45('resus_settle_cv',tint,C0); % general ode equation solver
         %[t2,C2]=ode23('resus_settle_cv',tint,C0); % stiff equation solver
         Cvs(:,k,n)=C2(end,:);
      end
   else
      Cvs(:,:,n)=Cvs(:,:,n-1);
   end
   
   if(do_floc)
      if(do_settle)
         tint = (0:t(n)-t(n-1));
         % Settle all particles
         for k=1:npmud
            w = wsf(k);
            Drate = 0; Erate = 0; % no flux bottom bc
            C0 = squeeze( Cm(:,k,n-1) );
            [t2,C2]=ode45('resus_settle_leer',tint,C0); % general ode equation solver
            Cm(:,k,n)=C2(end,:);
         end
      end
      % Floc all depths
      for iz = 1:nzc
         cv_tmp = Cm(iz,:,n)';
         fm_bbl_calcfloc
         Cm(iz,:,n) = cv_wat;
      end
   else
      Cm(:,:,n)=Cm(:,:,n-1);
   end
end
