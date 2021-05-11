function Sim = Simulate(par,bApE,bApU,T)
% Sim = Simulate(par,bApE,bApU,T)
% 
% Simulates the model.
%
% Inputs:
% par       Parameter structure
% bAp       Polynomial coefficients for polynomial for Ap policy rule



Sim.A = zeros(T,1);
Sim.Y = zeros(T,1);

Sim.A(1) = 1;


Sim.Y_tilde(1) = par.mu;
Eps = par.sigma * randn(T,1);

x=[0,1];
Sim.status(1)=x(randi(length(x)));
if Sim.status(1)==0
    Sim.Y=par.b;
else 
    Sim.Y=Sim.Y_tilde(1);
end

for t = 2:T
    if Sim.status(t-1)==0
        w=[(1-par.q),par.q];
    else 
        w=[par.p,(1-par.p)];
    end
    Sim.status(t)=randsample(x,1,true,w);
    Sim.Y_tilde(t) = (1-par.rho)*par.mu+par.rho*Sim.Y_tilde(t-1) + Eps(t);
    if Sim.status(t)==0
        Sim.Y(t) = par.b;
        Sim.A(t) = PolyBasis(Sim.A(t-1),Sim.Y(t)) * bApU;
        Sim.C(t) = Sim.A(t-1)+par.b - (Sim.A(t)/(1+par.r));
    else 
        Sim.Y(t) = Sim.Y_tilde(t);
        Sim.A(t) = PolyBasis(Sim.A(t-1),Sim.Y(t)) * bApE;
        Sim.C(t) = Sim.A(t-1)+Sim.Y(t) - (Sim.A(t)/(1+par.r));
    end
%     if Sim.A(t)<0
%         Sim.A(t)=0;
%     else 
%         Sim.A(t)=Sim.A(t);
%     end

end


% Compute quantities from state variables
Ti = 2:T-1;
Ap = Sim.A(Ti+1);
Sim.A = Sim.A(Ti);
Sim.C = Sim.C(Ti);
Sim.Y = Sim.Y(Ti)';
Sim.status = Sim.status';
Sim.Y_tilde = Sim.Y_tilde';
end