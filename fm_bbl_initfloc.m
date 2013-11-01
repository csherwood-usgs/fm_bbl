% fm_bbl_init_floc - Initialize floc parameters and coefficients

% fm_print_init - Print out initial values
  fprintf(1,'\n');
  fprintf(1,'***********************\n')
  fprintf(1,'    FLOCMOD\n')
  fprintf(1,'***********************\n')

  for iv=1:npmud
     fprintf(1,'% 3d % 7.1f % 14g % 6.1f % 14g % 12f\n',...
        iv,Df(iv)*1e6,volf(iv),rhof(iv),mass(iv),wsf(iv))
  end
  fprintf(1,'\n')
  fprintf(1,' *** PARAMETERS ***\n')
  fprintf(1,'Primary particle size (Dp)                                   : %f\n',Dp)
  fprintf(1,'Fractal dimension (nf)                                       : %f\n',nf)
  fprintf(1,'Flocculation efficiency (alpha)                              : %f\n',alpha)
  fprintf(1,'Floc break up parameter (beta)                               : %f\n',beta)
  fprintf(1,'Nb of fragments (f_nb_frag)                                  : %f\n',f_nb_frag)
  fprintf(1,'Ternary fragmentation (f_ater)                               : %f\n',f_ater)
  fprintf(1,'Floc erosion (pct of mass) (f_ero_frac)                      : %f\n',f_ero_frac)
  fprintf(1,'Nb of fragments by erosion (f_ero_nbfrag)                    : %f\n',f_ero_nbfrag)
  fprintf(1,'fragment class (f_ero_iv)                                    : %f\n',f_ero_iv)
  fprintf(1,'negative mass tolerated before redistribution (f_mneg_param) : %f\n',f_mneg_param)
  fprintf(1,'Boolean for differential settling aggregation (L_ADS)        : %d\n',l_ADS)
  fprintf(1,'Boolean for shear aggregation (L_ASH)                        : %d\n',l_ASH)
  fprintf(1,'Boolean for collision fragmenation (L_COLLFRAG)              : %d\n',l_COLLFRAG)
  fprintf(1,'Collision fragmentation parameter (f_collfragparam)          : %f\n',f_collfragparam)
  fprintf(1,'\n')
  fprintf(1,'*** END FLOCMOD INIT *** \n')    

  if ~(l_ADS+l_ASH)
     fprintf(1,'CAUTION : incompatible flocculation kernel options : \n')
     fprintf(1,'*****************************************************\n')
     fprintf(1,'l_ADS=%d\n',l_ADS)
     fprintf(1,'l_ASH=%d\n',l_ASH)
     error('simulation stopped')
  end

% dimension arrays
f_coll_prob_sh=zeros(npmud,npmud);
f_coll_prob_ds=zeros(npmud,npmud);
f_g1_sh = zeros(npmud,npmud,npmud);
f_g1_ds = zeros(npmud,npmud,npmud);
f_g3 = zeros(npmud,npmud);
f_l3 = zeros(npmud);
f_g4=zeros(npmud,npmud,npmud);
f_l4=zeros(npmud,npmud);

%% floc kernals
mass = [mass; mass(npmud)*2.0+eps];
%% fm_flocmod_aggregation_statistics
% flocmod_agregation_statistics
for iv1=1:npmud
   for iv2=1:npmud
      % Shear (eqn 9)
      f_coll_prob_sh(iv1,iv2)=1.0/6.0*(Df(iv1)+Df(iv2))^3.0;
      % Differential settling
      f_coll_prob_ds(iv1,iv2)=0.250*pi*(Df(iv1)+Df(iv2))^2.0 ...
         *g/mu*abs((rhof(iv1)-rhow)*Df(iv1)^2.0 ...
         -(rhof(iv2)-rhow)*Df(iv2)^2.0);      
   end
end
%% fm_aggregation_gain
diffmass = diff([mass; 0]);
for iv1=1:npmud
   for iv2=1:npmud
      for iv3=iv2:npmud
         if(iv1==1)
            masslo = 0.0;
         else
            masslo = mass(iv1-1);
         end
         if((mass(iv2)+mass(iv3)) > masslo ...
               && ((mass(iv2)+mass(iv3)) <= mass(iv1)))
            
            %f_weight=(mass(iv2)+mass(iv3)-mass(iv1-1))/(mass(iv1)-mass(iv1-1));
            f_weight=(mass(iv2)+mass(iv3)-masslo)/(diffmass(iv1));
            
         elseif ((mass(iv2)+mass(iv3)) > mass(iv1) ...
               && ((mass(iv2)+mass(iv3)) < mass(iv1+1)))
            
            if (iv1 == npmud)
               f_weight=1.0;
            else
               %f_weight=1.0-(mass(iv2)+mass(iv3)-mass(iv1))/(mass(iv1+1)-mass(iv1));
               f_weight=1.0-(mass(iv2)+mass(iv3)-mass(iv1))/(diffmass(iv1+1));
            end
            
         else
            f_weight=0.0;
         end
         
         f_g1_sh(iv2,iv3,iv1)=f_weight*alpha*f_coll_prob_sh(iv2,iv3)*(mass(iv2)+mass(iv3))/mass(iv1);
         f_g1_ds(iv2,iv3,iv1)=f_weight*alpha*f_coll_prob_ds(iv2,iv3)*(mass(iv2)+mass(iv3))/mass(iv1);
         
      end
   end
end
clear f_weight
%% fm_shear_frag_gain
%****************************************
% Shear fragmentation : GAIN : f_g3
%****************************************
for iv1=1:npmud
   for iv2=iv1:npmud
      
      if(iv1==1)
         masslo = 0.0;
      else
         masslo = mass(iv1-1);
      end
      
      if (Df(iv2)>dfragmax)
         % binary fragmentation   
         if (mass(iv2)/f_nb_frag > masslo ...
               && mass(iv2)/f_nb_frag <= mass(iv1))
            
            if (iv1 == 1)
               f_weight=1.0;
            else
               f_weight=(mass(iv2)/f_nb_frag-masslo)/(mass(iv1)-masslo);
            end
         elseif (mass(iv2)/f_nb_frag > mass(iv1) ...
               && mass(iv2)/f_nb_frag < mass(iv1+1))
            f_weight=1.0-(mass(iv2)/f_nb_frag-mass(iv1))/(mass(iv1+1)-mass(iv1));
         else            
            f_weight=0.0;
         end
      else
         f_weight=0.0;
      end
      
      f_g3(iv2,iv1)=f_g3(iv2,iv1)+(1.0-f_ero_frac)*(1.0-f_ater)*f_weight*beta ...
         *Df(iv2)*((Df(iv2)-Dp)/Dp)^(3-nf)           ...
         *mass(iv2)/mass(iv1);
      
      % ternary fragmentation
      if (Df(iv2)>dfragmax)
         if (mass(iv2)/(2.0*f_nb_frag) > masslo ...
               && mass(iv2)/(2.0*f_nb_frag) <= mass(iv1))
            if (iv1 == 1)
               f_weight=1.0;
            else
               f_weight=(mass(iv2)/(2.0*f_nb_frag)-mass(iv1-1))/(mass(iv1)-masslo);
            end
         elseif (mass(iv2)/(2.0*f_nb_frag) > mass(iv1) ...
               && mass(iv2)/(2.0*f_nb_frag) < mass(iv1+1))
            f_weight=1.0-(mass(iv2)/(2.0*f_nb_frag)-mass(iv1))/(mass(iv1+1)-mass(iv1));            
         else
            f_weight=0.0;
            
         end
         % update for ternary fragments
         f_g3(iv2,iv1)=f_g3(iv2,iv1)+(1.0-f_ero_frac)*(f_ater)*f_weight*beta ...
            *Df(iv2)*((Df(iv2)-Dp)/Dp)^(3-nf)           ...
            *mass(iv2)/mass(iv1);
         
         % Floc erosion
         if ((mass(iv2)-mass(f_ero_iv)*f_ero_nbfrag) > mass(f_ero_iv))
            if (((mass(iv2)-mass(f_ero_iv)*f_ero_nbfrag) >masslo) ...
                  && (mass(iv2)-mass(f_ero_iv)*f_ero_nbfrag) <= mass(iv1))              
               if (iv1 == 1)
                  f_weight=1.0;
               else
                  f_weight=(mass(iv2)-mass(f_ero_iv)*f_ero_nbfrag-masslo)/(mass(iv1)-masslo);
               end
            elseif ((mass(iv2)-mass(f_ero_iv)*f_ero_nbfrag) > mass(iv1) ...
                  && (mass(iv2)-mass(f_ero_iv)*f_ero_nbfrag) < mass(iv1+1))
               f_weight=1.0-(mass(iv2)-mass(f_ero_iv)*f_ero_nbfrag-mass(iv1))/(mass(iv1+1)-mass(iv1));       
            else
               f_weight=0.0;
            end
            
            % update for eroded floc masses            
            f_g3(iv2,iv1)=f_g3(iv2,iv1)+f_ero_frac*f_weight*beta                    ...
               *Df(iv2)*((Df(iv2)-Dp)/Dp)^(3-nf)           ...
               *(mass(iv2)-mass(f_ero_iv)*f_ero_nbfrag)/mass(iv1);
            
            if (iv1 == f_ero_iv)               
               f_g3(iv2,iv1)=f_g3(iv2,iv1)+f_ero_frac*beta                           ...
                  *Df(iv2)*((Df(iv2)-Dp)/Dp)^(3-nf)           ...                                      ...
                  *f_ero_nbfrag*mass(f_ero_iv)/mass(iv1);
            end
         end
      end % condition on dfragmax
   end
end
clear f_weight
%% fm_aggregation_loss
%**************************************************************************
%   Shear agregation : LOSS : f_l1
%**************************************************************************
for iv1=1:npmud
   for iv2=1:npmud     
      if(iv2 == iv1)
         mult=2.0;
      else
         mult=1.0;
      end
      f_l1_sh(iv2,iv1)=mult*alpha*f_coll_prob_sh(iv2,iv1);
      f_l1_ds(iv2,iv1)=mult*alpha*f_coll_prob_ds(iv2,iv1);      
   end
end
clear mult
%% f m_shear_frag_loss
%**************************************************************************
%  Shear fragmentation : LOSS : f_l2 (stet...f13)
%**************************************************************************
% TODO - this is easy to vectorize
for iv1=1:npmud
   if (Df(iv1)>dfragmax)
      % shear fragmentation
      f_l3(iv1)=f_l3(iv1)+(1.0-f_ero_frac)*beta*Df(iv1)*((Df(iv1)-Dp)/Dp)^(3-nf);
      % shear erosion
      if ((mass(iv1)-mass(f_ero_iv)*f_ero_nbfrag) > mass(f_ero_iv))
         f_l3(iv1)=f_l3(iv1)+f_ero_frac*beta*Df(iv1)*((Df(iv1)-Dp)/Dp)^(3-nf);
      end
   end
end

% I think we can delete fmass(np+1) now
mass = mass(1:npmud,1);
%% fm_kernal_stats
knames = {'f_coll_prob_sh',...
'f_coll_prob_ds',...
'f_g1_sh',...
'f_g1_ds',...
'f_g3',...
'f_l3'};

for i=1:length(knames);
fprintf('%s ',knames{i})
%eval(['size( ',char(knames{i}),')'])
eval(['val(1)=min( ',char(knames{i}),'(:));'])
eval(['val(2)=max( ',char(knames{i}),'(:));'])
eval(['val(3)=sum( ',char(knames{i}),'(:));'])
eval(['val(4)=sum(abs( ',char(knames{i}),'(:)));'])
%fprintf(1,'Min: %g Max: %g Sum: %g SumAbs: %g\n',val)
fprintf('%g\n',val(3))
end


