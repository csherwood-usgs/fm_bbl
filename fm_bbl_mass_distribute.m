% fm_bbl_mass_redistribute -  based on a tolerated negative mass parameter, negative masses
% are redistributed linearly towards remaining postive masses

% TODO - Vectorize this
mneg=0.0;
npos=0.0;
NNtmp=NNin;

% for iv=1:np
%    if (NNin(iv) < 0.0) 
%       mneg=mneg-NNin(iv)*mass(iv);
%       NNtmp(iv)=0.0;
%    else
%       npos=npos+1.0;
%    end
% end

mneg=-sum(NNin((NNin<0)).*mass((NNin<0)));
npos=sum(any(NNin>0));



if (mneg > 0.0) 
   if (npos == 0.0) 
      error('All floc sizes have negative mass!')
   else
      for iv=1:np
         if (NNin(iv) > 0.00) 
            NNin(iv)=NNin(iv)-mneg/sum(NNtmp)*NNin(iv)/mass(iv);
         else
            NNin(iv)=0.0;
         end         
      end
   end
end % and negative masses are set to 0



