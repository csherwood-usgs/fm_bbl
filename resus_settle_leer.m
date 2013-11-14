function dcdt = resus_settle_leer(t,c)
% kernal for resuspension/settling/diffusion model
% global dzc dzf w K Erate
% This uses control volume formulation with flux limiters
% ws is a uniform scalar, K is on interfaces
% Erate is flux into cell 1
% dzc are distances between cell centers
% dzf are distances between cell faces

% TODO - Either simplify this for the case at had (w<0, uniform grid) to
%        speed up, or test varying grid to make more general.

% csherwood@usgs.gov
% Thanks to Gretar Tryggvason's lecture notes:
% http://www.nd.edu/~grtyggva/CFD-Course
global dzc dzf w K Erate Drate
nzc = length(c);
dcdt = zeros(nzc,1);
% Diffusion
dcdt(1) =       K(1)*(c(2)-c(1))/dzc(1);
dcdt(2:nzc-1) = K(2:nzc-1).*(c(3:nzc)-c(2:nzc-1))./dzc(2:nzc-1)...
               -K(1:nzc-2).*(c(2:nzc-1)-c(1:nzc-2))./dzc(1:nzc-2);
dcdt(nzc) =     K(nzc-1)*(c(nzc-1)-c(nzc))/dzc(nzc-1);

% Advection
% r = (c(j+1)-c(j))./(c(j)-c(j-1)) but is calc'd as below to prevent NaNs
% and Infs

r  = zeros(nzc,1); % default r=0 is first-order upwind
f = zeros(nzc+1,1);
if(w>eps)
   top = w*(c(3:nzc)-c(2:nzc-1))./dzc(2:nzc-1);
   bot = w*(c(2:nzc-1)-c(1:nzc-2))./dzc(1:nzc-2);
   r(2:nzc-1,1) = (top.*bot)./(bot.^2+eps);
   pr = (r+abs(r))./(1+r);                   % PHI(r) van Leer
   % pr = max(0,min(2*r,(r+1)./2.2))          % PHI(r) MUSCL
   % pr = max(0,max( min(2*r,1),min(r, 2)));    % PHI(r) Superbee
   f(3:nzc) =  w*( c(2:nzc-1)+0.5*pr(2:nzc-1).*(c(2:nzc-1)-c(1:nzc-2)) );
   f(2) = w*c(1);
elseif(w<-eps)
   % TODO - Change indices for w<0 to avoid flipping arrays
   c=flipud(c);
   dzc=flipud(dzc);
   top = -w*(c(3:nzc)-c(2:nzc-1));
   bot = -w*(c(2:nzc-1)-c(1:nzc-2));
   r(2:nzc-1,1) = (top.*bot)./(bot.^2+eps);
   pr = (r+abs(r))./(1+r) ;              % PHI(r) van Leer
   % pr = max(0,min(2*r,(r+1)./2.2));    % PHI(r) MUSCL
   % pr = max(0,max( min(2*r,1),min(r, 2)));    % PHI(r) Superbee
   f(3:nzc) =  w*( c(2:nzc-1)+0.5*pr(2:nzc-1).*(c(2:nzc-1)-c(1:nzc-2)) );
   f=flipud(f);
   c=flipud(c);
   dzc=flipud(dzc);
   f(nzc)=w*c(nzc);
end
% flux difference 
dcdt = dcdt - diff(f);

% Settling and resuspension
dcdt(1) = dcdt(1)-Drate*(w<0)*abs(w)*c(1) + Erate;
dcdt = dcdt./dzf;