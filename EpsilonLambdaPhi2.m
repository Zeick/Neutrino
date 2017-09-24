clear;
nh=true; m1=0.2; fs = 20;

% Linear scale
lambdaRange = 0.02:0.02:5;

alpha = 1; beta = 1;    % We look for the best limit, eps_m(e,e)
ne = 7.645e-18*1.0e27;  % In eV^3
Gf = 1.1664e-5*1.0e-18; % Fermi constant in eV^-2
Enu = 2.0e9; 		    % Neutrino energy in eV
A = 2*sqrt(2)*Gf*Enu*ne;% Matter density potential of electrons
e = 1; mu = 2;
m21 = 7.5e-5;           
s12 = sqrt(0.306); s23 = sqrt(0.441); delta = 0;

% Values from 0907.0097 (Blennow-Fernandez-Martinez) Model-independent NSI
% (Experimental limits)
eps_ee_max = 4.2;      eps_emu_max = 0.3;   eps_etau_max = 3.0;
eps_mumu_max = 0;  eps_mutau_max = 0.04; eps_tautau_max = 0.15;
eps = [[eps_ee_max         eps_emu_max         eps_etau_max];
           [conj(eps_emu_max)  eps_mumu_max        eps_mutau_max];
           [conj(eps_etau_max) conj(eps_mutau_max) eps_tautau_max]];

% Values from 1609.08637 (Blennow-Coloma) Table I nonunitarity column
% (Nonunitarity limits)
% Note: Original limits are 2sigma limits -> Transforming to 90% CL
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
eps_dune = [[eps_ee_max         eps_emu_max         eps_etau_max];
       [conj(eps_emu_max)  eps_mumu_max        eps_mutau_max];
       [conj(eps_etau_max) conj(eps_mutau_max) eps_tautau_max]];

if(nh) % Normal hierarchy
    m31 = 2.524e-3;     % Squared mass differences in eV^2
    s13 = sqrt(0.02166);
else % Inverse hierarchy
    m31 = 2.514e-3;     % Squared mass differences in eV^2
    s13 = sqrt(0.02179);
end
Gf = 1.1664e-5*1.0e-18; % Fermi constant in eV^-2
v = 174.0e9;            % SM Higgs field VEV in eV
if(nh)
    m2 = sqrt(m21 + m1^2);    m3 = sqrt(m31 + m1^2); % NORMAL HIERARCHY
else
    m3 = sqrt(m31 + m1^2);    m2 = sqrt(m3^2 - m21); % INVERSE HIERARCHY
end
mD2 = diag([m1^2 m2^2 m3^2]); % Diagonal mass matrix squared
U = GenerateMixingMatrix(s12,s13,s23,delta);
mDeltaRange = [1 2 5 10 20 50]*1.0e12;
%mDeltaRange = (0.5:0.5:10)*1.0e12;
%mDeltaRange = 1.0e12;
values = zeros(length(mDeltaRange),length(lambdaRange));
k=0;
for mDelta = mDeltaRange
    k=k+1;
    j=1;
    for lambda = lambdaRange
        mnu2 = GenerateMassMatrix(U,m1,m21,m31,nh);
        if (alpha == e)
            mnu2(alpha, beta) = mnu2(alpha, beta) + A*eps(alpha, beta);
            if (alpha ~= beta) % For non-diagonal elements both (alpha,beta) and (beta,alpha)
                mnu2(beta, alpha) = mnu2(beta, alpha) + A*conj(eps(alpha,beta));
            end
        end
        mnu2dag = ctranspose(mnu2);
        if(alpha == beta) % Case ee -> ee - mumu; Case tautau -> tautau - mumu
            result = -mDelta^2*(sqrt(mnu2(e,alpha)*mnu2dag(alpha,e)) - sqrt(mnu2(e,mu)*mnu2dag(mu,e)))/(8*sqrt(2)*Gf*v^4*lambda^2);
        else
            result = -mDelta^2*sqrt(mnu2(e,alpha)*mnu2dag(alpha,e))/(8*sqrt(2)*Gf*v^4*lambda^2);
        end
        values(k,j) = log10(sqrt(abs(result))); % lambda_phi
        j=j+1;
    end
    plot(lambdaRange, values(k,:));
    hold on;
end
set(gca,'FontSize',fs);
xlabel('\lambda_\phi (eV)','FontSize',fs);
ylabel(['log_{10}|' char(949) '|'],'FontSize',fs);
legend('{\fontsize{15} 1 TeV}','{\fontsize{15} 2 TeV}','{\fontsize{15} 5 TeV}','{\fontsize{15} 10 TeV}','{\fontsize{15} 20 TeV}','{\fontsize{15} 50 TeV}','Location','NorthEast');
title(['m_1 = ' num2str(m1) ' eV']);

y1 = ones(1,length(lambdaRange))*log10(eps(alpha,beta));
y2 = ones(1,length(lambdaRange))*log10(eps_dune(alpha,beta));
y3 = ones(1,length(lambdaRange))*log10(eps_nonunit(alpha,beta));
%plot(lambdaRange, y1,'-k','LineWidth',3);
%plot(lambdaRange, y2,'-k','LineWidth',3);
%plot(lambdaRange, y3,'-k','LineWidth',3);

yl = ylim; % Minimum and maximum value of y-axis
xl = xlim; % and x
xm = (xl(1)+xl(2))/2; % middle value of x-axis
xm = 0.5*xm;
yupper = ones(1,length(lambdaRange))*6;
lambdaRange = [0 lambdaRange];
X=[lambdaRange,fliplr(lambdaRange)];
Y=[y1,y1(1),yupper(1),fliplr(yupper)];
h=fill(X,Y,'b');
set(h,'facealpha',.2)

Y=[y2,y2(1),y1(1),fliplr(y1)];
h=fill(X,Y,'g');
set(h,'facealpha',.2)

Y=[y3,y3(1),y2(1),fliplr(y2)];
h=fill(X,Y,'r');
set(h,'facealpha',.2)
text(xm,2,'Excluded','FontSize',fs)
text(xm,0.3,'DUNE','FontSize',fs)
text(xm,-2,'New interactions (NSI)','FontSize',fs)
text(xm,-3.4,'Possible misinterpretation of NSI','FontSize',fs)
