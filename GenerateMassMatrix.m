% (C) Timo Karkkainen 2017 (last chage 20.9.2017)
% Helper-function. Creates Neutrino mass matrix squared from PMNS
% matrix, lightest neutrino mass and mass squared differences.
% Requires mass ordering choice (normal true/false)
function MassMatrix = GenerateMassMatrix(U,m1,m21,m31,normalOrdering)
    ne = 7.645e-18*1.0e27;  % In eV^3
    Gf = 1.1664e-5*1.0e-18; % Fermi constant in eV^-2
    Enu = 2.0e9; 		    % Neutrino energy in eV
    A = 2*sqrt(2)*Gf*Enu*ne;% Matter density potential of electrons
    Udag = ctranspose(U);
    if(normalOrdering)
        m2 = sqrt(m21 + m1^2);    m3 = sqrt(m31 + m1^2); % NORMAL HIERARCHY
    else
        m3 = sqrt(m31 + m1^2);    m2 = sqrt(m3^2 - m21); % INVERSE HIERARCHY
    end
    mD2 = diag([m1^2 m2^2 m3^2]); % Diagonal mass matrix squared
    MassMatrix = U*mD2*Udag + diag([A 0 0]);
end