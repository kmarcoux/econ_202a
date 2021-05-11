function b = PolyGetCoef(A,Y,Z)
% b = PolyGetCoef(Grid,Z)
%   Fits the polynomial from PolyBasis to the function(s) in column(s) of
%   Z.
%
% inputs
% A    n x 1   points for A
% Y    n x 1   points for Y
% Z    n x 1   values for function at (A,Y)
%
% outputs
% b    6 x 1   basis coefficients

b = PolyBasis(A,Y) \ Z;

end