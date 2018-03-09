% (C) Timo Karkkainen 2017-18
% Last modified: 9.3.2018
% Helper function. Returns maximum of M/Lambda acquired from minimum of 
%                  neutrino mass matrix elements
% 
% NOTE: s23range and deltarange can be single-valued too.
function values = RajatJaPlotti(alpha,beta,m1range,s23range,deltarange,eps,nh,anti)
% Electron number density of Earth: 10^-30 m^-3 in natural units
% Note: 1/GeV = 0.197e-15 m and 1 m = 5076142.13198 eV
ne = 7.645e-18*1.0e27;  % In eV^3
if anti % Potential has different sign for antineutrinos
    ne = -ne;
end
Gf = 1.1664e-5*1.0e-18; % Fermi constant in eV^-2
Enu = 2.0e9; 		    % Neutrino energy in eV
A = 2*sqrt(2)*Gf*Enu*ne;% Matter density potential of electrons

e = 1; mu = 2;
m21 = 7.5e-5;           
s12 = sqrt(0.306);
if(alpha == beta)
    fprintf('NSI parameter: %d%d - 22. ',alpha,alpha);
else
    fprintf('NSI parameter: %d%d. ', alpha, beta);
end
if(nh) % Normal hierarchy
    %fprintf('Using normal mass ordering.\n');
    m31 = 2.524e-3;     % Squared mass differences in eV^2
    s13 = sqrt(0.02166);
else % Inverse hierarchy
    %fprintf('Using inverted mass ordering.\n');
    m31 = 2.514e-3;     % Squared mass differences in eV^2
    s13 = sqrt(0.02179);
end
Gf = 1.1664e-5*1.0e-18; % Fermi constant in eV^-2
v = 174.0e9;            % SM Higgs field VEV in eV
values = zeros(1,length(m1range));
j=1;
for m1 = m1range
%    fprintf('%d\t%d\t%d/100 ready\n',alpha,beta,j-1);
    if(nh)
        m2 = sqrt(m21 + m1^2);    m3 = sqrt(m31 + m1^2); % NORMAL HIERARCHY
    else
        m3 = sqrt(m31 + m1^2);    m2 = sqrt(m3^2 - m21); % INVERSE HIERARCHY
    end
    mD2 = diag([m1^2 m2^2 m3^2]); % Diafgonal mass matrix squared
    for s23 = s23range
        for delta = deltarange
            U = GenerateMixingMatrix(s12,s13,s23,delta);
            Udag = ctranspose(U);
            mnu2 = U*mD2*Udag + diag([A 0 0]);
            if (alpha == e)
               mnu2(alpha, beta) = mnu2(alpha, beta) + A*eps(alpha, beta);
               if (alpha ~= beta) % For non-diagonal elements both (alpha,beta) and (beta,alpha)
                   mnu2(beta, alpha) = mnu2(beta, alpha) + A*conj(eps(alpha,beta));
               end
            end
            mnu2dag = ctranspose(mnu2);
            if(alpha == beta) % Case ee -> ee - mumu; Case tautau -> tautau - mumu
                result = eps(alpha,alpha)*8*sqrt(2)*Gf*v^4/abs(sqrt(mnu2(e,alpha)*mnu2dag(alpha,e)) - sqrt(mnu2(e,mu)*mnu2dag(mu,e)));
            else
                result = eps(alpha,beta)*8*sqrt(2)*Gf*v^4/abs(sqrt(mnu2(e,beta))*sqrt(mnu2dag(alpha,e)));
            end
%            fprintf('%d\t%g\n',j,result);
            values(j) = sqrt(result);
%            if(result > minValues(j,2)) % Use '<' for maxValues. Then, also must change initial values.
%                minValues(j,2) = result; 
%            end
        end
    end
    j=j+1;
end
values = log10(values);
plot(m1range, values);