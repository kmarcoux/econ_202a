% Store the parameter values in a structure 
par.gamma=2;
par.beta=.98;
par.mu=5;
par.b=1;
par.rho=.95;
par.sigma=sqrt(.01);
par.p=.05;
par.q=.55;
par.r=.01;


% Create a grid for Y
Grid.nY=7;
numStdY=4;
[Grid.Y, Grid.PY]=tauchen(Grid.nY,par.mu,par.rho,par.sigma,numStdY);
Grid.PY=Grid.PY';

% Create grid for A
Grid.nA=100;
Grid.A=linspace(-5,20,Grid.nA)';

% Create a product of the two grids
[YY,AA]=meshgrid(Grid.Y,Grid.A);
Grid.AA=AA(:);
Grid.YY=YY(:);

be=zeros(7,1);
bu=zeros(7,1);

ApE0 = zeros(size(Grid.AA));
ApU0 = zeros(size(Grid.AA));
MAXIT = 2000;
for it = 1:MAXIT

    [VE, ApE] = MaxBellmanE(par,be,bu,Grid);
    
    % take the expectation of the value function from the perspective of
    % the previous Y
    EVE = reshape(VE,Grid.nA,Grid.nY) * Grid.PY;

    % update our polynomial coefficients
    be = PolyGetCoef(Grid.AA,Grid.YY,EVE(:));

    [VU, ApU] = MaxBellmanU(par,be,bu,Grid);
    % take the expectation of the value function from the perspective of
    % the previous Y
    EVU = reshape(VU,Grid.nA,Grid.nY) * Grid.PY;

    % update our polynomial coefficients
    bu = PolyGetCoef(Grid.AA,Grid.YY,EVU(:));

    % see how much our policy rule has changed
    testE = max(abs(ApE0 - ApE));
    ApE0 = ApE; 

    testU = max(abs(ApU0 - ApU));
    ApU0 = ApU;

    disp(['iteration ' num2str(it) ', test = ' num2str(testE) ', ' num2str(testU)])
    if (testE < 1e-7 && testU < 1e-7)
        break
    end
end

bApE =  PolyGetCoef(Grid.AA,Grid.YY,ApE);
bApU =  PolyGetCoef(Grid.AA,Grid.YY,ApU);

%% Simulate

T = 10000;
Sim = Simulate(par,bApE,bApU,T);

%% Histograms
figure;
subplot(1, 3, 1)
histogram(Sim.A);
xlabel('A');
ylabel('Frequency');
title('Histogram of A');
subplot(1, 3, 2)
histogram(Sim.Y);
xlabel('Y');
ylabel('Frequency');
title('Histogram of Y');
subplot(1, 3, 3)
histogram(Sim.C);
xlabel('C');
ylabel('Frequency');
title('Histogram of C');

% Summary stats for table
stddevA=std(Sim.A);
meanA=mean(Sim.A);
stddevY=std(Sim.Y);
meanY=mean(Sim.Y);
stddevC=std(Sim.C);
meanC=mean(Sim.C);


%% Time series

figure;
subplot(1, 3, 1)
plot(Sim.A(1000:1100));
xlabel('Time');
ylabel('A');
title('Time Series of A');
xlim([0 100]);
subplot(1, 3, 2)
plot(Sim.Y(1000:1100));
xlabel('Time');
ylabel('Y');
title('Time Series of Y');
xlim([0 100]);
subplot(1, 3, 3)
plot(Sim.C(1000:1100));
xlabel('Time');
ylabel('C');
title('Time Series of C');
xlim([0 100]);

%% Part L

AL=linspace(0,meanA+.5*stddevA,Grid.nA)';

% Create a product of the two grids
[YYL,AAL]=meshgrid(Grid.Y,AL);

AAL=AAL(:);
YYL=YYL(:);

SavE= 100.*(PolyBasis(AAL,YYL)*bApE-AAL)./AAL;
SavE= reshape(SavE,Grid.nA, Grid.nY);

SavU= 100.*(PolyBasis(AAL,YYL)*bApU-AAL)./AAL;
SavU= reshape(SavU,Grid.nA, Grid.nY);

figure;
subplot(1,2,1)
plot(AL,SavE);
hold on
plot(AL,SavU(:,4), "black");
xlim([meanA-.5*stddevA,meanA+.5*stddevA]);
xlabel('[Mean A -.5SD, Mean A +.5SD]');
ylabel('Savings');
hold off
subplot(1,2,2)
plot(AL,SavE);
hold on
plot(AL,SavU(:,4), "black");
xlim([0,meanA+.5*stddevA]);
xlabel('[0, Mean A +.5SD]');
ylabel('Savings');

legend('Employed 1','Employed 2','Employed 3','Employed 4','Employed 5','Employed 6','Employed 7','Unemployed');

%% Part M

CE = YYL+AAL-(PolyBasis(AAL,YYL)*bApE)/(1+par.r);
CU = par.b+AAL-(PolyBasis(AAL,YYL)*bApU)/(1+par.r);
bCE = PolyGetCoef(AAL,YYL,CE);
bCU = PolyGetCoef(AAL,YYL,CU);

CE1=PolyBasis(AAL,YYL)*bCE;
CU1=PolyBasis(AAL,YYL)*bCU;

CE2=PolyBasis(AAL+.01*meanA,YYL)*bCE;
CU2=PolyBasis(AAL+.01*meanA,YYL)*bCU;

mpcE = (CE2-CE1)./(.01*meanA);
mpcU = (CU2-CU1)./(.01*meanA);

mpcE = reshape(mpcE,Grid.nA,Grid.nY);
mpcU = reshape(mpcU,Grid.nA,Grid.nY);

figure;
subplot(1,2,1)
plot(AL,mpcE);
hold on
plot(AL,mpcU(:,4), "black");
xlim([meanA-.5*stddevA,meanA+.5*stddevA]);
ylabel('Marginal Propensity to Consume');
xlabel('[Mean A -.5SD, Mean A +.5SD]');
hold off
subplot(1,2,2)
plot(AL,mpcE);
hold on
plot(AL,mpcU(:,4), "black");
xlim([0,meanA+.5*stddevA]);
xlabel('[0, Mean A +.5SD]');
ylabel('Marginal Propensity to Consume');

legend('Employed 1','Employed 2','Employed 3','Employed 4','Employed 5','Employed 6','Employed 7','Unemployed');


