% (C) Timo Karkkainen 2017
% Helper-function. Creates PMNS Neutrino mixing matrix from oscillation
% parameters. No Majorana phases acknowledged (yet).
%
% TODO: If 6 input parmeters, include Majorana phases, if 4, don't.
% Unnecessary for oscillation studies, however.
function MixingMatrix = GenerateMixingMatrix(s12,s13,s23,delta)
c12 = sqrt(1-s12^2);	c13 = sqrt(1-s13^2);	c23 = sqrt(1-s23^2);
MixingMatrix = [[c12*c13 s12*c13 s13*exp(-1i*delta)]; % PMNS matrix
	 [-s12*c23-c12*s23*s13*exp(1i*delta) c12*c23-s12*s23*s13*exp(1i*delta) s23*c13];
	 [s12*s23-c12*c23*s13*exp(1i*delta) -c12*s23-s12*c23*s13*exp(1i*delta) c23*c13]];
end