% Non-unitarity project (related to PMNS neutrino mixing matrix and CLFV)
% (C) Timo Karkkainen 2017 (last modified 28.6.17)

%%%%%%%%%%%%%%%%%%
% INITIALIZATION %
%%%%%%%%%%%%%%%%%%
clear               % Destroy all variables in the beginning
t = cputime;        % See the time elapsed (DON'T TOUCH!)
counter=1;          % See the progress     (DON'T TOUCH!)
stages=8;          % 
y_upper = 14;       % Upper limit of y-axis
y_lower = 8;       % Lower limit of y-axis
fs = 20;            % Font size of title and axes
nh = false;          % Normal hierarchy (false -> inverse hierarchy)
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
initialMin = 16*ones(1,length(m1range)); % Choose the lowest of 5 curves
%initialMin = zeros(1,length(m1range)); % Choose the highest of 5 curves

% Below in all cases
% eps_ee_max is actually (eps_ee - eps_mumu)_max and same for
% eps_tautau_max = (eps_tautau - eps_mumu)_max


% Values from XXXX.XXXXX 
% (LFV limits)
eps_emu_ee = 3.5e-7;
eps_etau_ee = 1.6e-4;
eps_mutau_mumu = 1.5e-4;
eps_etau_emu = 1.2e-4;
eps_mutau_mue = 1.3e-4;
eps_etau_mumu = 1.2e-4;
eps_etau_mue = 9.9e-5;
eps_mue_mue = 3e-3;

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
% PLOT FOR LFV DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('###############\n     LFV\n###############\n');
minValues=initialMin;
figure;
emu_ee = RajatJaPlotti_FLV(e,mu,e,e,m1range,s23range,deltarange,eps_emu_ee,nh);
hold on;
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1; hold on;
etau_ee = RajatJaPlotti_FLV(e,tau,e,e,m1range,s23range,deltarange,eps_etau_ee,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
mutau_mumu = RajatJaPlotti_FLV(mu,tau,mu,mu,m1range,s23range,deltarange,eps_mutau_mumu,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
etau_emu = RajatJaPlotti_FLV(e,tau,e,mu,m1range,s23range,deltarange,eps_etau_emu,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
mutau_mue = RajatJaPlotti_FLV(mu,tau,mu,e,m1range,s23range,deltarange,eps_mutau_mue,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
etau_mumu = RajatJaPlotti_FLV(e,tau,mu,mu,m1range,s23range,deltarange,eps_etau_mumu,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
etau_mue = RajatJaPlotti_FLV(e,tau,mu,e,m1range,s23range,deltarange,eps_etau_mue,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;
mue_mue = RajatJaPlotti_FLV(mu,e,mu,e,m1range,s23range,deltarange,eps_mue_mue,nh);
fprintf('%d/%d (%.2f seconds elapsed)\n',counter, stages, cputime-t); counter = counter+1;

% Beautifying the plots (showLegend, yMin, yMax, fontSize)
Beautify(false,y_lower,y_upper,fs);
legend('{\fontsize{15}e\mu,ee}','{\fontsize{15}e\tau,ee}','{\fontsize{15}\mu\tau,\mu\mu}','{\fontsize{15}e\tau,e\mu}','{\fontsize{15}\mu\tau,\mu e}','{\fontsize{15}e\tau,\mu\mu}','{\fontsize{15}e\tau,\mu e}','{\fontsize{15}\mu e,\mu e}','Location','NorthEast');
hold off;