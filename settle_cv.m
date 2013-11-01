function dcdt = settle_cv(t,c)
% kernal for settling/diffusion model
% This uses control volume formulation
% ws is a uniform scalar, K is on interfaces
% dz should be distances between cell centers
global dz w K
nz = length(c);
dcdt = zeros(nz,1);
%
dcdt(1) =      K(1)*(c(2)-c(1))/dz(1);
dcdt(2:nz-1) = K(2:nz-1).*(c(3:nz)-c(2:nz-1))./dz(2:nz-1)...
              -K(1:nz-2).*(c(2:nz-1)-c(1:nz-2))./dz(1:nz-2);
dcdt(nz) =     K(nz-1)*(c(nz-1)-c(nz))/dz(nz-1);

% upwind flux 
dcdt(2:nz) =   dcdt(2:nz)   +(w>0)*abs(w)*c(1:nz-1) - (w<0)*abs(w)*c(2:nz);
dcdt(1:nz-1) = dcdt(1:nz-1) -(w>0)*abs(w)*c(1:nz-1) + (w<0)*abs(w)*c(2:nz);

