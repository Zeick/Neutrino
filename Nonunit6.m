% Non-unitarity project (related to PMNS neutrino mixing matrix)
% (C) Timo Karkkainen 2017 (last modified 22.6.17)

% IMPORTANT NOTE: tautau_dune limit is incomputable, since
% eps_tautau - eps_mumu = 0 exactly in eps-matrix. The values
% are from paper 1601.07730. Otherwise everything works.

%%%%%%%%%%%%%%%%%%
% INITIALIZATION %
%%%%%%%%%%%%%%%%%%
clear               % Destroy all variables in the beginning
t = cputime;        % See the time elapsed (DON'T TOUCH!)
counter=1;          % See the progress     (DON'T TOUCH!)
stages=15;          % 
y_upper = 16;       % Upper limit of y-axis
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
m1range = 0.0:0.002:0.20;   % Lightest neutrino mass range
initialMin = 16*ones(1,length(m1range)); % Choose the highest of 5 curves
%initialMin = zeros(1,length(m1range)); % Choose the lowest of 5 curves

% Values from 0907.0097 (Blennow-Fernandez-Martinez) Model-independent NSI
eps_ee_max = 4.2;      eps_emu_max = 0.5*0.33;   eps_etau_max = 0.5*3.0;
eps_mumu_max = 0.068;  eps_mutau_max = 0.5*0.33; eps_tautau_max = 21;
eps_exp = [[eps_ee_max         eps_emu_max         eps_etau_max];
           [conj(eps_emu_max)  eps_mumu_max        eps_mutau_max];
           [conj(eps_etau_max) conj(eps_mutau_max) eps_tautau_max]];

% Values from 1609.08637 (Blennow-Coloma) Table I nonunitarity column
eps_ee_max = 1.3e-3;    eps_emu_max = 0.5*6.8e-4;     eps_etau_max = 0.5*2.7e-3;
eps_mumu_max = 2.0e-4;  eps_mutau_max = 0.5*1.2e-3;   eps_tautau_max = 2.8e-3;
eps_nonunit = [[eps_ee_max         eps_emu_max         eps_etau_max];
              [conj(eps_emu_max)  eps_mumu_max        eps_mutau_max];
              [conj(eps_etau_max) conj(eps_mutau_max) eps_tautau_max]];

% Values from our paper 1601.07730 table II case for LBNO 50 kt
eps_ee_max = 0.18;      eps_emu_max = 0.5*0.025;     eps_etau_max = 0.5*0.021;
eps_mumu_max = 0.063;   eps_mutau_max = 0.5*0.008;   eps_tautau_max = 0.063;
eps = [[eps_ee_max         eps_emu_max         eps_etau_max];
       [conj(eps_emu_max)  eps_mumu_max        eps_mutau_max];
       [conj(eps_etau_max) conj(eps_mutau_max) eps_tautau_max]];

%%%%
% Another possibility for defining the matrix, given limits for
% non-diagnonal entries and differencees of diagonal entries and eps_mumu
% Given max-values of  ee-mumu, tautau-mumu, emu, etau and mutau
% eps_new = [[eps_ee_max          eps_emu_max           eps_etau_max
%            [conj(eps_emu_max)   0                     eps_mutau_max
%            [conj(eps_etau_max)  conj(eps_mutau_max)   eps_tautau_max]];
% where eps_ee_max is actually (eps_ee - eps_mumu)_max and same for
% eps_tautau_max
%%%%%%
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
minValues=initialMin;
figure;
subplot(2,2,1);
ee_dune = RajatJaPlotti(e,e,m1range,s23range,deltarange,eps,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1; hold on;
emu_dune = RajatJaPlotti(e,mu,m1range,s23range,deltarange,eps,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
etau_dune = RajatJaPlotti(e,tau,m1range,s23range,deltarange,eps,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
mutau_dune = RajatJaPlotti(mu,tau,m1range,s23range,deltarange,eps,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
% MUST CHANGE LIMITS FOR EPS
tautau_dune = RajatJaPlotti(tau,tau,m1range,s23range,deltarange,eps,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;

% Beautifying the plots (showLegend, yMin, yMax, fontSize)
Beautify(true,y_lower,y_upper,fs);
text(0.01,15,'Our limit (Preliminary)','FontSize',fs)
hold off;

% Find minimum values of all datasets
minValues = min(minValues,ee_dune);
minValues = min(minValues,emu_dune);
minValues = min(minValues,etau_dune);
minValues = min(minValues,mutau_dune);
minValuesDune = min(minValues,tautau_dune);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FOR MAXIMAL NON-UNITARY DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('#################\n# NON-UNITARITY #\n#################\n');
minValues=initialMin;

% figure;
subplot(2,2,2);
ee_nonunit = RajatJaPlotti(e,e,m1range,s23range,deltarange,eps_nonunit,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1; hold on;
emu_nonunit = RajatJaPlotti(e,mu,m1range,s23range,deltarange,eps_nonunit,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
etau_nonunit = RajatJaPlotti(e,tau,m1range,s23range,deltarange,eps_nonunit,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
mutau_nonunit = RajatJaPlotti(mu,tau,m1range,s23range,deltarange,eps_nonunit,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
tautau_nonunit = RajatJaPlotti(tau,tau,m1range,s23range,deltarange,eps_nonunit,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;

Beautify(true,y_lower,y_upper,fs);
text(0.01,15,'Non-unitarity limit','FontSize',fs)
hold off;

minValues = min(minValues,ee_nonunit);
minValues = min(minValues,emu_nonunit);
minValues = min(minValues,etau_nonunit);
minValues = min(minValues,mutau_nonunit);
minValuesNonunit = min(minValues,tautau_nonunit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FOR CURRENT EXPERIMENTAL DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('##################\n# CURRENT LIMITS #\n##################\n');

minValues=initialMin;
subplot(2,2,3);
ee_exp = RajatJaPlotti(e,e,m1range,s23range,deltarange,eps_exp,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1; hold on;
emu_exp = RajatJaPlotti(e,mu,m1range,s23range,deltarange,eps_exp,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
etau_exp = RajatJaPlotti(e,tau,m1range,s23range,deltarange,eps_exp,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
mutau_exp = RajatJaPlotti(mu,tau,m1range,s23range,deltarange,eps_exp,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
tautau_exp = RajatJaPlotti(tau,tau,m1range,s23range,deltarange,eps_exp,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;

Beautify(true,y_lower,y_upper,fs);
text(0.01,15,'Current experimental limit','FontSize',fs)
hold off;

minValues = min(minValues,ee_exp);
minValues = min(minValues,emu_exp);
minValues = min(minValues,etau_exp);
minValues = min(minValues,mutau_exp);
minValuesExp = min(minValues,tautau_exp);

%%%%%%%%%%%%%%%%%%%
% COMBINE RESULTS %
%%%%%%%%%%%%%%%%%%%
%figure;
subplot(2,2,4);
%area(m1range,log10(minValuesExp),'FaceColor','g');
plot(m1range,minValuesExp,'Color','g');
hold on;
%area(m1range,log10(minValuesTot),'FaceColor', [0.5 0.5 0.5]); % Grey = [0.5 0.5 0.5]
%area(m1range,log10(minValuesTot),'FaceColor', 'c');
%area(m1range,log10(minValuesNonUnit),'FaceColor','r');
plot(m1range,minValuesDune,'Color', 'b');
plot(m1range,minValuesNonunit,'Color','r');
Beautify(false,y_lower,y_upper,20);
text(0.01,15.0,'Excluded','FontSize',fs)
text(0.01,13.1,'DUNE coverage','FontSize',fs)
text(0.01,11,'Nonunitarity','FontSize',fs)
%legend('{\fontsize{10}Experimental limit}', '{\fontsize{10}DUNE coverage}','{\fontsize{10}Non-unitary limit}','Location','NorthEast');
hold off;