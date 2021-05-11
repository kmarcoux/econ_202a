function VU = BellmanU(par, be, bu, A, Y, Ap )
% VE = BellmanE( Par, be, bu, A, Y, Ap )
%   Evaluate the RHS of the BellmanE equation
%
% Inputs
% Par   Parameter structure
% be    6 x 1 coefficients in polynomial for E[ V(A',Y') | Ye ]
% bu    6 x 1 coefficients in polynomial for E[ V(A',Y') | Yu ]
% A     n x 1 array of current wealth
% Y     n x 1 array of current income
% Ap    n x 1 array of savings
%
% Output
% VE     n x 1 array of value employed function
%

C = par.b+A-(1/(1+par.r))*Ap;
if C<0
    u = -inf;
else
    u = C.^(1-par.gamma) / (1-par.gamma);
end
VU = u + par.beta * ((1-par.q)*PolyBasis(Ap,Y) * bu+(par.q)*PolyBasis(Ap,Y) * be);

end