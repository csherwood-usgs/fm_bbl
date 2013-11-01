% fm_bbl_comp_fsd
% This processes NNin and returns NNout
NNout = zeros(size(NNin));
tmp_g1=0.0;
tmp_g3=0.0;
tmp_g4=0.0;
tmp_l1=0.0;
tmp_l3=0.0;
tmp_l4=0.0;
f_g1_tmp=zeros(npmud,npmud,npmud);
f_l1_tmp=zeros(npmud,npmud);

if (l_COLLFRAG)
   fm_collfrag
end

for iv1=1:npmud
   f_g1_tmp(:,:,iv1)=f_g1_tmp(:,:,iv1)+l_ASH*f_g1_sh(:,:,iv1)*Gval;
   f_g1_tmp(:,:,iv1)=f_g1_tmp(:,:,iv1)+l_ADS*f_g1_ds(:,:,iv1);
   tmp_g1=tmp_g1+(NNin'*(f_g1_tmp(:,:,iv1))*NNin);
   tmp_g4=tmp_g4+l_COLLFRAG*(NNin'*(f_g4(:,:,iv1)*Gval)*NNin);
   tmp_g3=tmp_g3+f_g3(:,iv1)'*NNin*Gval^1.5;
   
   f_l1_tmp(:,iv1)=f_l1_tmp(:,iv1)+l_ASH*f_l1_sh(:,iv1)*Gval;
   f_l1_tmp(:,iv1)=f_l1_tmp(:,iv1)+l_ADS*f_l1_ds(:,iv1)*Gval;
   tmp_l1=tmp_l1+(f_l1_tmp(:,iv1))'*NNin;
   tmp_l4=tmp_l4+l_COLLFRAG*(f_l4(:,iv1)*Gval)'*NNin;
   tmp_l1=tmp_l1*NNin(iv1);
   tmp_l4=tmp_l4*NNin(iv1); 
   tmp_l3=f_l3(iv1)*Gval^1.50*NNin(iv1);
   
   NNout(iv1)=NNin(iv1)+f_dt*(tmp_g1+tmp_g3+tmp_g4-(tmp_l1+tmp_l3+tmp_l4));
   
   tmp_g1=0.0;
   tmp_g3=0.0;
   tmp_g4=0.0;
   tmp_l1=0.0;
   tmp_l3=0.0;
   tmp_l4=0.0;
end
