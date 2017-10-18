% Non-unitarity project
% (C) Timo Karkkainen 2017 (last modified 10.10.17)
clear;
fs=20; % Font size for the plots
Gf = 1.1664e-5*1.0e-18; % Fermi constant in eV^-2
% Values from 0907.0097 (Blennow-Fernandez-Martinez) Model-independent NSI
% (Experimental limits)
eps_ee = 4.2;
eps_emu = 0.3;
eps_etau = 3.0;
eps_mutau = 0.04;
eps_tautau = 0.15;
eps_all = [eps_ee eps_emu eps_etau eps_mutau eps_tautau];

Ymax = 3.5449; % Raja ei-perturbatiivisuudelle: sqrt(4*pi) = 3.5449...
Yrange = 0.01:0.01:4;
M = zeros(1,length(Yrange));
figure; hold on;
for eps=eps_all
    j=1;
    for Y=Yrange
        M(j) = Y/(8*Gf^2*eps^2)^(1/4);
        if(eps == eps_ee || eps == eps_tautau)
           M(j) = M(j)*sqrt(2); 
        end
        j=j+1;
    end
    plot(Yrange,M*1e-12);
    fprintf('Upper limit for triplet mass: %.2f TeV\n', M(354)*1e-12);
    %plot(Yrange,M*1e-12,'-r','LineWidth',2);
end
set(gca,'FontSize',fs);

y1 = ones(1,length(Yrange))*0.75;
yl = ylim; % Minimum and maximum value of y-axis
xl = xlim; % and x
xm = (xl(1)+xl(2))/2; % middle value of x-axis
xm = 0.5*xm;
yupper = ones(1,length(Yrange))*0;
Yrange = [(2*Yrange(1) - Yrange(2)) Yrange];
X=[Yrange,fliplr(Yrange)];
Y=[y1,y1(1),yupper(1),fliplr(yupper)];
h=fill(X,Y,'b');
set(h,'facealpha',.2)


legend('{\fontsize{15} ee-\mu\mu}','{\fontsize{15} e\mu}','{\fontsize{15} e\tau}','{\fontsize{15} \mu\tau}','{\fontsize{15} \tau\tau}','Location','NorthEast');
xlabel('Y_{\alpha\beta}','FontSize',fs);
ylabel('M_{\Delta} (TeV)','FontSize',fs);
title('Triplet mass as a function of Yukawa coupling','FontSize',fs);

%ylim([1e9 1e10])
%xlim([0 1])