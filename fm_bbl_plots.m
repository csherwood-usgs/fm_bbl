%% fm_bbl_plots
%    size(Df)
%    size(C(:,:,n))
%    size(sum(C(:,:,n),2))
figure(2); clf;
subplot(121)
hold on
for k=1:np
   plot(squeeze(Cm(:,k,n)),zc ,'Color',cmap(k,:),'LineWidth',2 )
end

plot(sum(squeeze(Cm(:,:,n)),2),zc,'k' ,'LineWidth',2)
legend([num2str(Df*1e6,'%3.0f');'SUMC'])
title('FLOCMOD model - Concentration','FontSize',16)
axis([.001 1 .1 12])
grid on
set(gca,'FontSize',14,'XScale','log','YScale','log','box','on')
xlabel('Concentration (g/l)','FontSize',16);   ylabel('distance from the bed (m)','FontSize',16);

subplot(122)
hold on
plot(squeeze(Cm(:,:,n)*1e6*Df./sum(Cm(:,:,n),2)),zc,'b' ,'LineWidth',2)
plot(1e6*(sum(squeeze(C(:,:,n)),2).*ka./kb./sqrt(G)+Dp),zc,'r','LineWidth',2)
xlabel('mean (mass weighted) floc diameter (um)','FontSize',16);   ylabel('distance from the bed (m)','FontSize',16);
set(gca,'FontSize',14,'XScale','log','YScale','log','box','on')
title('FLOCMOD - Mass averaged diameter','FontSize',16)
legend ({'FLOCMOD';'EQUI WIN'},'Location','south')
axis([20 1500 .1 12])
grid on
drawnow
shg
pause(.1)

%%
%load('..\fm_bbl_win\win_settling_ustar_01.mat','Gsave')

cmap=jet(np);
t2=t/86400;
tmin=min(t2);
tmax=max(t2);
dts=13;
fs=16;
figure

subplot(411)
imagesc(t2,zc,Gsave)
title('FLOCMOD Settling mid depth -- FSD','FontSize',fs)
axis([tmin tmax 0 12])
colormap(cmap)
set(gca,'YDir','Normal','XTick',(tmin:(tmax-tmin)/12:tmax),'XTickLabel',datestr((tmin:(tmax-tmin)/12:tmax),dts),'FontSize',12);
%colorbar('Location','East','YTick',(1:3:np),'YTickLabel',num2str(1e6*Df(1:3:np),'%3.0f'))

subplot(412)
bar(t2,squeeze(Cm(nzc,:,:))','Stack')
title('FLOCMOD Settling surface -- FSD','FontSize',fs)
axis([tmin tmax 0 0.3])
colormap(cmap)
set(gca,'YDir','Normal','XTick',(tmin:(tmax-tmin)/12:tmax),'XTickLabel',datestr((tmin:(tmax-tmin)/12:tmax),dts),'FontSize',12);
colorbar('Location','East','YTick',(1:3:np),'YTickLabel',num2str(1e6*Df(1:3:np),'%3.0f'))

subplot(413)
bar(t2,squeeze(Cm(floor(nzc/2),:,:))','Stack')
title('FLOCMOD Settling mid depth -- FSD','FontSize',fs)
axis([tmin tmax 0 0.3])
colormap(cmap)
set(gca,'YDir','Normal','XTick',(tmin:(tmax-tmin)/12:tmax),'XTickLabel',datestr((tmin:(tmax-tmin)/12:tmax),dts),'FontSize',12);
colorbar('Location','East','YTick',(1:3:np),'YTickLabel',num2str(1e6*Df(1:3:np),'%3.0f'))

subplot(414)
bar(t2,squeeze(Cm(1,:,:))','Stack')
title('FLOCMOD Settling Bottom -- FSD','FontSize',fs)
axis([tmin tmax 0 1])
colormap(cmap)
set(gca,'YDir','Normal','XTick',(tmin:(tmax-tmin)/12:tmax),'XTickLabel',datestr((tmin:(tmax-tmin)/12:tmax),dts),'FontSize',12);
colorbar('Location','East','YTick',(1:3:np),'YTickLabel',num2str(1e6*Df(1:3:np),'%3.0f'))

cmap=jet(20);
figure
subplot(311)
imagesc(t2,zc,Gsave)
title('FLOCMOD Settling -- G (/s) values','FontSize',fs)
axis([tmin tmax 0 12])
colormap(cmap)
set(gca,'YDir','Normal','XTick',(tmin:(tmax-tmin)/12:tmax),'XTickLabel',datestr((tmin:(tmax-tmin)/12:tmax),dts),'FontSize',12);
colorbar('Location','East')
%colorbar('Location','East','YTick',(1:3:np),'YTickLabel',num2str(1e6*Df(1:3:np),'%3.0f'))
caxis([0 10])


for i=1:length(t)
   Dfsave(:,i)=squeeze((Cm(:,:,i)*1e6*Df)./sum(Cm(:,:,i),2));
end

subplot(312)
imagesc(t2,zc,log10(Dfsave))
title('FLOCMOD Settling surface -- Davg (um)','FontSize',fs)
axis([tmin tmax 0 12])
colormap(cmap)
set(gca,'YDir','Normal','XTick',(tmin:(tmax-tmin)/12:tmax),'XTickLabel',datestr((tmin:(tmax-tmin)/12:tmax),dts),'FontSize',12);
caxis([1 3])
colorbar('Location','East','YTick',(1:1:3),'YTickLabel',{'10';'100';'1000'})


subplot(313)
imagesc(t2,zc,log10(squeeze(sum(C,2))))
title('FLOCMOD Settling surface -- SSC (g/l)','FontSize',fs)
axis([tmin tmax 0 12])
colormap(cmap)
set(gca,'YDir','Normal','XTick',(tmin:(tmax-tmin)/12:tmax),'XTickLabel',datestr((tmin:(tmax-tmin)/12:tmax),dts),'FontSize',12);
caxis([-2 0])
colorbar('Location','East','YTick',(-3:1:0),'YTickLabel',{'0.001';'0.01';'0.1';'1'})

ik=find(Df>d50s,1,'first');
Cvs_Df=zeros(size(Cmud));
Cvs_Df(:,ik-1,:)=Cvs(:,1,:)*rhos*(Df(ik)-d50s)/(Df(ik)-Df(ik-1));
Cvs_Df(:,ik,:)=Cvs(:,1,:)*rhos*(1-(Df(ik)-d50s)/(Df(ik)-Df(ik-1)));


figure(3); clf; subplot(311)
plot(t2,uw,'-b',t2,uc,'-r')
axis([tmin tmax 0 0.4])
set(gca,'YDir','Normal','XTick',(tmin:(tmax-tmin)/12:tmax),'XTickLabel',datestr((tmin:(tmax-tmin)/12:tmax),dts),'FontSize',12);
legend('Uw','Uc')
   
   
figure
subplot(311)
imagesc(t2,(1:npmud),squeeze(Cmud(6,:,:)+Cvs_Df(6,:,:))) % 0.5mab
title('FLOCMOD sand/mud 0.55mab -- FSD','FontSize',fs)
axis([tmin tmax 0 npmud])
colormap(cmap)
caxis([0 0.05])
set(gca,'YDir','Normal','XTick',(tmin:(tmax-tmin)/12:tmax),'YTick',(1:3:npmud),'YTickLabel',num2str(1e6*Df(1:3:npmud),'%3.0f'),'XTickLabel',datestr((tmin:(tmax-tmin)/12:tmax),dts),'FontSize',12);
% colorbar('Location','East','YTick',(1:3:np),'YTickLabel',num2str(1e6*Df(1:3:np),'%3.0f'))

subplot(312)
%bar(t2,squeeze(C(floor(nzc/2),:,:)+Cvs_Df(floor(nzc/2),:,:))','Stack')
imagesc(t2,(1:npmud),squeeze(Cmud(3,:,:)+Cvs_Df(3,:,:))) % 0.25 mab
title('FLOCMOD sand/mud 0.25mab -- FSD','FontSize',fs)
axis([tmin tmax 0 npmud])
colormap(cmap)
caxis([0 0.05])
set(gca,'YDir','Normal','XTick',(tmin:(tmax-tmin)/12:tmax),'YTick',(1:3:npmud),'YTickLabel',num2str(1e6*Df(1:3:npmud),'%3.0f'),'XTickLabel',datestr((tmin:(tmax-tmin)/12:tmax),dts),'FontSize',12);
% colorbar('Location','East','YTick',(1:3:np),'YTickLabel',num2str(1e6*Df(1:3:np),'%3.0f'))

subplot(313)
imagesc(t2,(1:npmud),squeeze(Cmud(1,:,:)+Cvs_Df(1,:,:))) % 0.05mab
title('FLOCMOD Settling 0.05mab -- FSD','FontSize',fs)
axis([tmin tmax 0 npmud])
colormap(cmap)
caxis([0 1])
set(gca,'YDir','Normal','XTick',(tmin:(tmax-tmin)/12:tmax),'YTick',(1:3:npmud),'YTickLabel',num2str(1e6*Df(1:3:npmud),'%3.0f'),'XTickLabel',datestr((tmin:(tmax-tmin)/12:tmax),dts),'FontSize',12);
% colorbar('Location','East','YTick',(1:3:np),'YTickLabel',num2str(1e6*Df(1:3:np),'%3.0f'))




%% plot specific events

specevent=[3600/86400 3*3600/86400 4*3600/86400];
nspec=length(specevent);
fs=14;
size(Cvs_Df)
figure
for i=1:nspec
    tspec=find(t2>specevent(i),1,'first');
    
subplot(3,nspec,i)
hold on
plot(1e6*Df,squeeze(Cmud(6,:,tspec)),'r','LineWidth',2)
plot(1e6*Df,squeeze(Cvs_Df(6,:,tspec)),'b','LineWidth',2)
plot(1e6*Df,squeeze(Cmud(6,:,tspec)+Cvs_Df(6,:,tspec)),'k','LineWidth',2)
axis([4 1500 0 0.2])
xlabel('Particle size (um)','FontSize',fs); ylabel('Concentration (g/l)','FontSize',fs);
if i==1; legend('mud flocs','sand','total'); ylabel({'0.5mab';'Concentration (g/l)'},'FontSize',fs);end
title(['date=' datestr(t2(tspec),13)],'FontSize',fs)
set(gca,'XScale','log','FontSize',12)

subplot(3,nspec,i+nspec)
hold on
plot(1e6*Df,squeeze(Cmud(3,:,tspec)),'r','LineWidth',2)
plot(1e6*Df,squeeze(Cvs_Df(3,:,tspec)),'b','LineWidth',2)
plot(1e6*Df,squeeze(Cmud(3,:,tspec)+Cvs_Df(3,:,tspec)),'k','LineWidth',2)
axis([4 1500 0 0.2])
set(gca,'XScale','log','FontSize',12)
xlabel('Particle size (um)','FontSize',fs); ylabel('Concentration (g/l)','FontSize',fs);
if i==1; ylabel({'0.25mab';'Concentration (g/l)'},'FontSize',fs);end

subplot(3,nspec,i+2*nspec)
hold on
plot(1e6*Df,squeeze(Cmud(1,:,tspec)),'r','LineWidth',2)
plot(1e6*Df,squeeze(Cvs_Df(1,:,tspec)),'b','LineWidth',2)
plot(1e6*Df,squeeze(Cmud(1,:,tspec)+Cvs_Df(1,:,tspec)),'k','LineWidth',2)
axis([4 1500 0 0.2])
set(gca,'XScale','log','FontSize',12)
xlabel('Particle size (um)','FontSize',fs); ylabel('Concentration (g/l)','FontSize',fs);
if i==1; ylabel({'0.05mab';'Concentration (g/l)'},'FontSize',fs);end
end

toc



% profile viewer