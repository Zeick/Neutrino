% Non-unitarity project
% (C) Timo Karkkainen 2017 (last modified 10.10.17)
clear;
fs=20; % Font size for the plots
Gf = 1.1664e-5*1.0e-18; % Fermi constant in eV^-2
% Values from 0907.0097 (Blennow-Fernandez-Martinez) Model-independent NSI
% (Experimental limits)
eps_1211 = 3.5e-7;
eps_1311 = 1.4e-4;
eps_2322 = 1.2e-4;
eps_1312 = 1.0e-4;
eps_2321 = 1.0e-4;
eps_1322 = 1.0e-4;
eps_1321 = 9.9e-5;
eps_all = [eps_1211 eps_1311 eps_2322 eps_1312 eps_2321 eps_1322 eps_1321];

Ymax = 3.54; % Raja ei-perturbatiivisuudelle: sqrt(4*pi) = 3.5449...
Yrange = 0.01:0.01:4;
M = zeros(1,length(Yrange));
figure; hold on;
eps=eps_1211;
for eps=eps_all
    j=1;
    for Y=Yrange
        M(j) = Y/(8*Gf^2*eps^2)^(1/4);
        j=j+1;
    end
    plot(Yrange,M*1e-12);
    fprintf('Upper limit for triplet mass: %.2f TeV\n', M(354)*1e-12);
    %plot(Yrange,M*1e-12,'-r','LineWidth',2);
end
set(gca,'FontSize',fs);
xlabel('Y_{\alpha\beta}','FontSize',fs);
ylabel('M_{\Delta} (TeV)','FontSize',fs);
title(['Triplet mass as a function of Yukwawa coupling'],'FontSize',fs);

%ylim([1e9 1e10])
%xlim([0 1])