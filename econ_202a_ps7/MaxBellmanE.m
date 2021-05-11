function [VE, Ap] = MaxBellmanE(par,be,bu,Grid)
% [VE, Ap] = MaxBellmanE(Par,b,Grid)
%   Maximizes the RHS of the Bellman equation using golden section search
%
% Inputs
% Par       Parameter structure
% b     6 x 1 coefficients in polynomial for E[ V(A',Y') | Y ]
% Grid      Grid structure


p = (sqrt(5)-1)/2;

K = Grid.A(1) * ones(size(Grid.AA));
D = min(f(par,Grid.AA,Grid.YY) - 1e-3, Grid.A(end)); % -1e-3 so we always have positve consumption.

B = p*K+(1-p)*D;
C = (1-p)*K + p * D;


fB = BellmanE(par,be,bu,Grid.AA,Grid.YY,B);
fC = BellmanE(par,be,bu,Grid.AA,Grid.YY,C);


MAXIT = 1000;
for it_inner = 1:MAXIT

    if all(D-K < 1e-6)
        break
    end

    I = fB > fC;

    D(I) = C(I);
    C(I) = B(I);
    fC(I) = fB(I);
    B(I) = p*C(I) + (1-p)*K(I);
    fB(I) = BellmanE(par,be,bu,Grid.AA(I),Grid.YY(I),B(I));

    K(~I) = B(~I);
    B(~I) = C(~I);
    fB(~I) = fC(~I);
    C(~I) = p*B(~I) + (1-p)*D(~I);
    fC(~I) = BellmanE(par,be,bu,Grid.AA(~I),Grid.YY(~I),C(~I));

end

% At this stage, K, B, C, and D are all within a small epsilon of one
% another.  We will use the average of B and C as the optimal level of
% savings.
Ap = (B+C)/2;

% evaluate the Bellman equation at the optimal policy to find the new
% value function.
VE = BellmanE(par,be,bu,Grid.AA,Grid.YY,Ap);

end