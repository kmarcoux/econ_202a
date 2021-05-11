function B  = PolyBasis(A,Y)
% B  = PolyBasis(A,Y)
% Polynomial basis functions.  Using 2nd order polynomial
%
% inputs
% A    n x 1   points for A
% Y    n x 1   points for Y
%     or scalar for Y
%
% outputs
% B    n x 6   array of basis functions: 1, A, Y, A^2, A*Y, Y^2

%B = [ones(size(A)) A Y A.^2 A.*Y Y.^2 ];

Yb = Y.*ones(size(A));
B = [ones(size(A)) A Yb A.^2 A.*Yb Yb.^2 A.^3];

end