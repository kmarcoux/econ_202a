function y  = f(par,A,Y)
% y  = f(A,Y)
%   Function for next period wealth

y = (1+par.r)*(Y+A);
end