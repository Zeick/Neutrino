% Non-unitarity project (related to PMNS neutrino mixing matrix)
% (C) Timo Karkkainen 2017 (last modified 13.9.17)

%%%%%%%%%%%%%%%%%%
% INITIALIZATION %
%%%%%%%%%%%%%%%%%%
clear               % Destroy all variables in the beginning
t = cputime;        % See the time elapsed (DON'T TOUCH!)
counter=1;          % See the progress     (DON'T TOUCH!)
stages=15;          % 
y_upper = 14;       % Upper limit of y-axis
y_lower = 10;       % Lower limit of y-axis
fs = 20;            % Font size of title and axes
nh = true;          % Normal hierarchy (false -> inverse hierarchy)
s23Fix = false;     % If true, set best-fit value for s23
deltaFix = false;   % likewise for deltaCP

% Connect lepton flavour to integers
e = 1;              % Use exp(1) for Neper's number if necessary
mu = 2;
tau = 3;

%%%%%%%%%%%%%%%%%%
% SET PARAMETERS %
%%%%%%%%%%%%%%%%%%
mDeltaRange = (500:100:10000)*10^9;
m1 = 0.2; % Lightest neutrino mass

% Below in all cases
% eps_ee_max is actually (eps_ee - eps_mumu)_max and same for
% eps_tautau_max = (eps_tautau - eps_mumu)_max

% Values from 0907.0097 (Blennow-Fernandez-Martinez) Model-independent NSI
% (Experimental limits)
eps_ee_max = 4.2;      eps_emu_max = 0.3;   eps_etau_max = 3.0;
eps_mumu_max = 0;  eps_mutau_max = 0.04; eps_tautau_max = 0.15;
eps_exp = [[eps_ee_max         eps_emu_max         eps_etau_max];
           [conj(eps_emu_max)  eps_mumu_max        eps_mutau_max];
           [conj(eps_etau_max) conj(eps_mutau_max) eps_tautau_max]];

% Values from 1609.08637 (Blennow-Coloma) Table I nonunitarity column
% (Nonunitarity limits)
% Note: Original limits are 2sigma limits -> Transformaxg to 90% CL
eps_ee_max = 1.3e-3;    eps_emu_max = 0.5*6.8e-4;     eps_etau_max = 0.5*2.7e-3;
eps_mumu_max = 2.0e-4;  eps_mutau_max = 0.5*1.2e-3;   eps_tautau_max = 2.8e-3;
eps_nonunit = [[eps_ee_max         eps_emu_max         eps_etau_max];
              [conj(eps_emu_max)  eps_mumu_max        eps_mutau_max];
              [conj(eps_etau_max) conj(eps_mutau_max) eps_tautau_max]];
eps_nonunit = 1.64485/2.0*eps_nonunit; % mathworld.wolfram.com/ConfidenceInterval.html

% Values from 1606.08851 (Blennow-Ohlsson) Table I all NSI 
% (DUNE limits)
eps_ee_max = 0.9;      eps_emu_max = 0.074;     eps_etau_max = 0.19;
eps_mumu_max = 0;   eps_mutau_max = 0.038;   eps_tautau_max = 0.08;
eps = [[eps_ee_max         eps_emu_max         eps_etau_max];
       [conj(eps_emu_max)  eps_mumu_max        eps_mutau_max];
       [conj(eps_etau_max) conj(eps_mutau_max) eps_tautau_max]];

% Check s23 and deltaCP fixing. If fixed, best-fit. Hierarchy-dependent.
if s23Fix
   if nh
       s23range = sqrt(0.441);
   else
       s23range = sqrt(0.587);
   end
else
    s23range = sqrt(0.385:0.01:0.635);
end
if deltaFix
    deltarange = -pi/2;
else
    deltarange = 0:0.04:6.28;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FOR LBNO/DUNE DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('###############\n     DUNE\n###############\n');
figure;
subplot(2,2,1);
ee_dune = RajatJaPlotti_ML(e,e,m1,mDeltaRange,s23range,deltarange,eps,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1; hold on;
emu_dune = RajatJaPlotti_ML(e,mu,m1,mDeltaRange,s23range,deltarange,eps,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
etau_dune = RajatJaPlotti_ML(e,tau,m1,mDeltaRange,s23range,deltarange,eps,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
mutau_dune = RajatJaPlotti_ML(mu,tau,m1,mDeltaRange,s23range,deltarange,eps,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
tautau_dune = RajatJaPlotti_ML(tau,tau,m1,mDeltaRange,s23range,deltarange,eps,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
% Beautifying the plots (showLegend, ymax, yMax, fontSize)
set(gca,'FontSize',fs);
xlabel('M_{\Delta} (eV)','FontSize',fs);
ylabel('log_{10}|\lambda_{\phi}/eV|','FontSize',fs);
title('DUNE (m1=0.2 eV)','FontSize',fs);
hold off;

[maxIndexValues, maxValuesDune] = FindMaxIndex(ee_dune, emu_dune, etau_dune, mutau_dune, tautau_dune);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FOR MAXIMAL NON-UNITARY DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('#################\n# NON-UNITARITY #\n#################\n');

%figure;
subplot(2,2,2);
ee_nonunit = RajatJaPlotti_ML(e,e,m1,mDeltaRange,s23range,deltarange,eps_nonunit,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1; hold on;
emu_nonunit = RajatJaPlotti_ML(e,mu,m1,mDeltaRange,s23range,deltarange,eps_nonunit,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
etau_nonunit = RajatJaPlotti_ML(e,tau,m1,mDeltaRange,s23range,deltarange,eps_nonunit,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
mutau_nonunit = RajatJaPlotti_ML(mu,tau,m1,mDeltaRange,s23range,deltarange,eps_nonunit,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
tautau_nonunit = RajatJaPlotti_ML(tau,tau,m1,mDeltaRange,s23range,deltarange,eps_nonunit,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;

set(gca,'FontSize',fs);
xlabel('M_{\Delta} (eV)','FontSize',fs);
ylabel('log_{10}|\lambda_{\phi}/eV|','FontSize',fs);
title('Nonunitarity (m1=0.2 eV)','FontSize',fs);
hold off;

maxValuesNonunit = GetValuesByIndex(maxIndexValues, ee_nonunit, emu_nonunit, etau_nonunit, mutau_nonunit, tautau_nonunit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FOR CURRENT EXPERIMENTAL DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('##################\n# CURRENT LIMITS #\n##################\n');
%figure;
subplot(2,2,3);
ee_exp = RajatJaPlotti_ML(e,e,m1,mDeltaRange,s23range,deltarange,eps_exp,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1; hold on;
emu_exp = RajatJaPlotti_ML(e,mu,m1,mDeltaRange,s23range,deltarange,eps_exp,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
etau_exp = RajatJaPlotti_ML(e,tau,m1,mDeltaRange,s23range,deltarange,eps_exp,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
mutau_exp = RajatJaPlotti_ML(mu,tau,m1,mDeltaRange,s23range,deltarange,eps_exp,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
tautau_exp = RajatJaPlotti_ML(tau,tau,m1,mDeltaRange,s23range,deltarange,eps_exp,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;

set(gca,'FontSize',fs);
xlabel('M_{\Delta} (eV)','FontSize',fs);
ylabel('log_{10}|\lambda_{\phi}/eV|','FontSize',fs);
title('Experimental limit (m1=0.2 eV)','FontSize',fs);
hold off;

maxValuesExp = GetValuesByIndex(maxIndexValues, ee_exp, emu_exp, etau_exp, mutau_exp, tautau_exp);

%%%%%%%%%%%%%%%%%%%
% COMBINE RESULTS %
%%%%%%%%%%%%%%%%%%%
%figure;
subplot(2,2,4);
%area(m1range,log10(maxValuesExp),'FaceColor','g');
plot(mDeltaRange,maxValuesExp,'Color',[0 0.5 0]);
hold on;
%area(m1range,log10(maxValuesTot),'FaceColor', [0.5 0.5 0.5]); % Grey = [0.5 0.5 0.5]
%area(m1range,log10(maxValuesTot),'FaceColor', 'c');
%area(m1range,log10(maxValuesNonUnit),'FaceColor','r');
plot(mDeltaRange,maxValuesDune,'Color', 'b');
plot(mDeltaRange,maxValuesNonunit,'Color','r');
%plot(m1range,emu_ee,'Color','k');
set(gca,'FontSize',fs);
xlabel('M_{\Delta} (eV)','FontSize',fs);
ylabel('log_{10}|\lambda_{\phi}/eV|','FontSize',fs);
title('Combined plot','FontSize',fs);
%text(0.01,15.0,'Excluded','FontSize',fs)
%text(0.01,13.1,'DUNE coverage','FontSize',fs)
%text(0.01,11,'Nonunitarity','FontSize',fs)
legend('{\fontsize{15}Experimental limit}', '{\fontsize{15}DUNE coverage}','{\fontsize{15}Non-unitary limit}','Location','NorthEast');
hold off;