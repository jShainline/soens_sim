%%
clc;
clear all;
close all;
[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;

p = f_physicalConstants;

%% 
gamma = 0.5772;
NVec = linspace(1e1,1e4,10000);
APL_vec = [2 3 4 5];%vector of average path lengths to consider
diameter_vec = [2 3 4 5];%vector of network diameters to consider

%% Erdos-Renyi APL
kExpectation_mat = zeros(length(NVec),length(APL_vec));

for ii = 1:length(APL_vec)

    kExpectation_mat(:,ii) = exp( ( (log(NVec(:))-gamma) / APL_vec(ii) ) + 0.5 );
    
end

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
loglog(NVec(:),kExpectation_mat(:,1),'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
hold on
loglog(NVec(:),kExpectation_mat(:,2),'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
loglog(NVec(:),kExpectation_mat(:,3),'Color',bRGY(13,:),'LineStyle','-','LineWidth',3)
loglog(NVec(:),kExpectation_mat(:,4),'Color',bRGY(18,:),'LineStyle','-','LineWidth',3)
xlabel('N_{tot}' ,'FontSize',fontSize,'FontName','Times')
ylabel('<k>','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
legend(sprintf('APL = %g',APL_vec(1)),sprintf('APL = %g',APL_vec(2)),sprintf('APL = %g',APL_vec(3)),sprintf('APL = %g',APL_vec(4)))
title('Average node degree versus network size for a given APL, random','FontSize',fontSize,'FontName',fontName)
% xlim([1e3 1e8])
% ylim([1e-6 1e2])
grid on

%% Erdos-Renyi diameter
% kExpectation_mat2 = zeros(length(NVec),length(diameter_vec));
% 
% for ii = 1:length(diameter_vec)
% 
%     kExpectation_mat2(:,ii) = exp( log(NVec(:)) / diameter_vec(ii) );
%     
% end
% 
% figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
% loglog(NVec(:),kExpectation_mat2(:,1),'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
% hold on
% loglog(NVec(:),kExpectation_mat2(:,2),'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
% loglog(NVec(:),kExpectation_mat2(:,3),'Color',bRGY(13,:),'LineStyle','-','LineWidth',3)
% loglog(NVec(:),kExpectation_mat2(:,4),'Color',bRGY(18,:),'LineStyle','-','LineWidth',3)
% xlabel('N_{tot}' ,'FontSize',fontSize,'FontName','Times')
% ylabel('<k>','FontSize',fontSize,'FontName','Times')
% set(gca,'FontSize',fontSize,'FontName',fontName)
% legend(sprintf('d = %g',diameter_vec(1)),sprintf('d = %g',diameter_vec(2)),sprintf('d = %g',diameter_vec(3)),sprintf('d = %g',diameter_vec(4)))
% % title(sprintf('kMax vs number of nodes in network, gamma = %g',gamma),'FontSize',fontSize,'FontName',fontName)
% % xlim([1e3 1e8])
% % ylim([1e-6 1e2])
% grid on

%% Scale-free APL
m = 2;
kExpectation_mat = zeros(length(NVec),length(APL_vec));

for ii = 1:length(APL_vec)

    kExpectation_mat(:,ii) = exp( (log(NVec(:))+(1.5-APL_vec(ii))*log(log(NVec(:)))-1-gamma)/(APL_vec(ii)-0.5) );
    
end

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
loglog(NVec(:),kExpectation_mat(:,1),'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
hold on
loglog(NVec(:),kExpectation_mat(:,2),'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
loglog(NVec(:),kExpectation_mat(:,3),'Color',bRGY(13,:),'LineStyle','-','LineWidth',3)
loglog(NVec(:),kExpectation_mat(:,4),'Color',bRGY(18,:),'LineStyle','-','LineWidth',3)
xlabel('N_{tot}' ,'FontSize',fontSize,'FontName','Times')
ylabel('<k>','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)
legend(sprintf('APL = %g',APL_vec(1)),sprintf('APL = %g',APL_vec(2)),sprintf('APL = %g',APL_vec(3)),sprintf('APL = %g',APL_vec(4)))
title('Average node degree versus network size for a given APL, scale-free','FontSize',fontSize,'FontName',fontName)
% xlim([1e3 1e8])
% ylim([1e-6 1e2])
grid on
