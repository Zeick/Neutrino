% Non-unitarity project
% (C) Timo Karkkainen 2017 (last modified 26.4.17)
%%%%%%%%%%%%%%%%%%
% INITIALIZATION %
%%%%%%%%%%%%%%%%%%
clear               % Destroy all variables in the beginning
t = cputime;        % See the time elapsed (DON'T TOUCH!)
counter=1;          % See the progress     (DON'T TOUCH!)
y_upper = 16;       % Upper limit of y-axis
y_lower = 10;       % Lower limit of y-axis
fs = 20;            % Font size of title and axes
nh = true;          % Normal hierarchy
ih = not(nh);       % Inverse hierarchy    (DON'T TOUCH!)
s23Fix = false;     % If true, set best-fit value for s23
deltaFix = false;   % likewise for deltaCP

% Connect lepton flavour to integers
e = 1;              % Use exp(1) for Neper's number if necessary
mu = 2;
tau = 3;

%%%%%%%%%%%%%%%%%%
% SET PARAMETERS %
%%%%%%%%%%%%%%%%%%

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

m1range = 0.0:0.002:0.20;   % Lightest neutrino mass range

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

D=zeros(length(m1range),2);
D(:,1) = m1range;
minValues=[m1range' zeros(length(m1range),1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FOR LBNO/DUNE DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,2,1);
[p_ee, minValues] = RajatJaPlotti(e,e,m1range,s23range,deltarange,eps,D,minValues,ih);
fprintf('%d/18 (%.2f seconds elapsed)\n',counter, cputime-t); counter = counter+1;
hold on;
[p_emu, minValues] = RajatJaPlotti(e,mu,m1range,s23range,deltarange,eps,D,minValues,ih);
fprintf('%d/18 (%.2f seconds elapsed)\n',counter, cputime-t); counter = counter+1;
[p_etau, minValues] = RajatJaPlotti(e,tau,m1range,s23range,deltarange,eps,D,minValues,ih);
fprintf('%d/18 (%.2f seconds elapsed)\n',counter, cputime-t); counter = counter+1;
[p_mumu, minValues] = RajatJaPlotti(mu,mu,m1range,s23range,deltarange,eps,D,minValues,ih);
fprintf('%d/18 (%.2f seconds elapsed)\n',counter, cputime-t); counter = counter+1;
[p_mutau, minValues] = RajatJaPlotti(mu,tau,m1range,s23range,deltarange,eps,D,minValues,ih);
fprintf('%d/18 (%.2f seconds elapsed)\n',counter, cputime-t); counter = counter+1;
[p_tautau, minValues] = RajatJaPlotti(tau,tau,m1range,s23range,deltarange,eps,D,minValues,ih);
fprintf('%d/18 (%.2f seconds elapsed)\n',counter, cputime-t); counter = counter+1;

% Beautifying the plots (showLegend, yMin, yMax, fontSize)
Beautify(true,y_lower,y_upper,20);
ylabel('log_{10}(M_{\Delta}/|\lambda_{\phi}|)','FontSize',fs);
text(0.01,15,'Our limit (Preliminary)','FontSize',fs)
%hold off;

minValuesTot = minValues(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FOR MAXIMAL NON-UNITARY DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minValues=[m1range' zeros(length(m1range),1)];

% figure;
subplot(2,2,2);
[p_ee, minValues] = RajatJaPlotti(e,e,m1range,s23range,deltarange,eps_nonunit,D,minValues,ih);
hold on;
fprintf('%d/18 (%.2f seconds elapsed)\n',counter, cputime-t); counter = counter+1;
[p_emu, minValues] = RajatJaPlotti(e,mu,m1range,s23range,deltarange,eps_nonunit,D,minValues,ih);
fprintf('%d/18 (%.2f seconds elapsed)\n',counter, cputime-t); counter = counter+1;
[p_etau, minValues] = RajatJaPlotti(e,tau,m1range,s23range,deltarange,eps_nonunit,D,minValues,ih);
fprintf('%d/18 (%.2f seconds elapsed)\n',counter, cputime-t); counter = counter+1;
[p_mumu, minValues] = RajatJaPlotti(mu,mu,m1range,s23range,deltarange,eps_nonunit,D,minValues,ih);
fprintf('%d/18 (%.2f seconds elapsed)\n',counter, cputime-t); counter = counter+1;
[p_mutau, minValues] = RajatJaPlotti(mu,tau,m1range,s23range,deltarange,eps_nonunit,D,minValues,ih);
fprintf('%d/18 (%.2f seconds elapsed)\n',counter, cputime-t); counter = counter+1;
[p_tautau, minValues] = RajatJaPlotti(tau,tau,m1range,s23range,deltarange,eps_nonunit,D,minValues,ih);
fprintf('%d/18 (%.2f seconds elapsed)\n',counter, cputime-t); counter = counter+1;

Beautify(true,y_lower,y_upper,fs);
ylabel('log_{10}(M_{\Delta}/|\lambda_{\phi}|)','FontSize',fs);
text(0.01,15,'Non-unitarity limit','FontSize',fs)
%hold off;

minValuesNonUnit = minValues(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FOR CURRENT EXPERIMENTAL DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minValues=[m1range' zeros(length(m1range),1)];

%figure;
subplot(2,2,3);
[p_ee, minValues] = RajatJaPlotti(e,e,m1range,s23range,deltarange,eps_exp,D,minValues,ih);
hold on;
fprintf('%d/18 (%.2f seconds elapsed)\n',counter, cputime-t); counter = counter+1;
[p_emu, minValues] = RajatJaPlotti(e,mu,m1range,s23range,deltarange,eps_exp,D,minValues,ih);
fprintf('%d/18 (%.2f seconds elapsed)\n',counter, cputime-t); counter = counter+1;
[p_etau, minValues] = RajatJaPlotti(e,tau,m1range,s23range,deltarange,eps_exp,D,minValues,ih);
fprintf('%d/18 (%.2f seconds elapsed)\n',counter, cputime-t); counter = counter+1;
[p_mumu, minValues] = RajatJaPlotti(mu,mu,m1range,s23range,deltarange,eps_exp,D,minValues,ih);
fprintf('%d/18 (%.2f seconds elapsed)\n',counter, cputime-t); counter = counter+1;
[p_mutau, minValues] = RajatJaPlotti(mu,tau,m1range,s23range,deltarange,eps_exp,D,minValues,ih);
fprintf('%d/18 (%.2f seconds elapsed)\n',counter, cputime-t); counter = counter+1;
[p_tautau, minValues] = RajatJaPlotti(tau,tau,m1range,s23range,deltarange,eps_exp,D,minValues,ih);
fprintf('%d/18 (%.2f seconds elapsed)\n',counter, cputime-t); counter = counter+1;

Beautify(true,y_lower,y_upper,fs);
ylabel('log_{10}(M_{\Delta}/|\lambda_{\phi}|)','FontSize',fs);
text(0.01,15,'Experimental NSI limit','FontSize',fs)
%hold off;

minValuesExp = minValues(:,2);


%%%%%%%%%%%%%%%%%%%
% COMBINE RESULTS %
%%%%%%%%%%%%%%%%%%%
figure;
%subplot(2,2,4);
area(m1range,log10(minValuesExp),'FaceColor','g');
hold on;
%area(m1range,log10(minValuesTot),'FaceColor', [0.5 0.5 0.5]); % Grey = [0.5 0.5 0.5]
area(m1range,log10(minValuesTot),'FaceColor', 'c');
area(m1range,log10(minValuesNonUnit),'FaceColor','r');
Beautify(false,y_lower,y_upper,20);
ylabel('log_{10}(M_{\Delta}/|\lambda_{\phi}|)','FontSize',fs);
text(0.08,15.5,'Excluded','FontSize',fs)
text(0.08,14.5,'Accessible by DUNE','FontSize',fs)
text(0.08,13.6,'NSI not from nonunitarity','FontSize',fs-5)
text(0.08,11,'Nonunitarity','FontSize',fs)
hold off;