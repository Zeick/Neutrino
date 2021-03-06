% Non-unitarity project (related to PMNS neutrino mixing matrix)
% (C) Timo Karkkainen 2017-18 (last modified 28.2.18)

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

% Values from 1612.07377 (Valle et al.) Table I
eps_ee_max = 0.02;      eps_emu_max = 1.0e-2;     eps_etau_max = 4.2e-2;
eps_mumu_max = 0.01;   eps_mutau_max = 9.8e-3;   eps_tautau_max = 0.07;
eps_nonunit = [[eps_ee_max         eps_emu_max         eps_etau_max];
       [conj(eps_emu_max)  eps_mumu_max        eps_mutau_max];
       [conj(eps_etau_max) conj(eps_mutau_max) eps_tautau_max]];

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
else % 90%CL = 1.645*1sigma CL
    if nh
        s23range = sqrt(0.424:0.01:0.592);
    else
        s23range = sqrt(0.500:0.01:0.592);
    end
end
if deltaFix
    deltarange = -pi/2;
else
    deltarange = 0:0.04:6.28;
end

for k=0:20
    m1 = 0.01*k;
    mDeltaRange = (0:100:10000)*10^9;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT FOR LBNO/DUNE DATA %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('###############\n     DUNE\n###############\n');
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

    [minIndexValues, minValuesDune] = FindMaxIndex(ee_dune, emu_dune, etau_dune, mutau_dune, tautau_dune);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT FOR MAXIMAL NON-UNITARY DATA %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('#################\n# NON-UNITARITY #\n#################\n');
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

    minValuesNonunit = GetValuesByIndex(minIndexValues, ee_nonunit, emu_nonunit, etau_nonunit, mutau_nonunit, tautau_nonunit);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT FOR CURRENT EXPERIMENTAL DATA %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('##################\n# CURRENT LIMITS #\n##################\n');
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

    minValuesExp = GetValuesByIndex(minIndexValues, ee_exp, emu_exp, etau_exp, mutau_exp, tautau_exp);

    %%%%%%%%%%%%%%%%%%%
    % COMBINE RESULTS %
    %%%%%%%%%%%%%%%%%%%

    %area(m1range,log10(minValuesExp),'FaceColor','g');
    mDeltaRange = mDeltaRange*1e-12;
    g = figure('visible', 'off');
    plot(mDeltaRange,minValuesExp,'Color',[0 0.5 0]);
    hold on;
    plot(mDeltaRange,minValuesDune,'Color', 'k');
    plot(mDeltaRange,minValuesNonunit,'Color','k');

    set(gca,'Color','c'); % Background color of the plot
    set(gcf,'color','w'); % Background color of the plot window
    area(mDeltaRange,minValuesNonunit,'FaceColor','g');
    area(mDeltaRange,minValuesDune,'FaceColor', 'y');
    area(mDeltaRange,minValuesExp,'FaceColor', 'w'); % Grey = [0.5 0.5 0.5]

    %text(0.01,15.0,'Excluded','FontSize',fs)
    %text(0.01,13.1,'DUNE coverage','FontSize',fs)
    %text(0.01,11,'Nonunitarity','FontSize',fs)

    set(gca,'FontSize',fs);
    xlim([0 10]);
    ylim([0 3]);
    xlabel('M_{\Delta} (TeV)','FontSize',fs);
    ylabel('\lambda_{\phi} (eV)','FontSize',fs);
    if nh
        text(8,0.5,'NH','FontSize',fs+20);
    else
        text(8,0.5,'IH','FontSize',fs+20);
    end
    title(['m_1 = ' num2str(m1) ' eV'],'FontSize',fs);
    %legend('{\fontsize{15}Experimental limit}', '{\fontsize{15}DUNE coverage}','{\fontsize{15}Non-unitary limit}','Location','NorthEast');
    hold off;
    if k < 10
        filename = sprintf('Fig3_nh0%d.png',k);
    else
        filename = sprintf('Fig3_nh%d.png',k);
    end
    g.InvertHardcopy = 'off';
    saveas(gcf,filename);
    %saveas(gcf,filename,'epsc');
    close(g);
end