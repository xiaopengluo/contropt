% Contractility of continuous optimization
% by Xiaopeng Luo and Xin Xu 
% 
% Fig. 6
% 
% Comparison of CO to BO and GD on 2D quadratic function plus a Gaussian
clear; clc;

global a;
a = 20; 


plus = @(x)10*exp(-4*(x(:,1)-0.5).^2-4*(x(:,2)+2).^2);
obj = @(x)a*x(:,1).^2+x(:,2).^2 - plus(x);

xstar = fminunc(obj,[0.5 -2]);
fstar = obj(xstar);

%{%
% Contraction Optimization
N = 2; maxN = 500; minL = 5; K = 8; 
lb = [-5 -5]; ub = [5 5];
type={'median','medium',0.2};
[best,mdl,thr] = contropt(obj,type,K,minL,N,maxN,lb,ub);

figure('units','normalized','position',[0.1,0.2,0.3,0.4])
scatter(mdl{1}.X(:,1),mdl{1}.X(:,2),25,'k.')
hold on
for i=2:sum(thr(:,1)~=0)
    scatter(mdl{i}.X(:,1),mdl{i}.X(:,2),25,'k.')
end
hold off
box on
axis([-5 5 -5 5])
title('x trace of CO')
saveas(gcf,'Q2P-COx','epsc')

figure('units','normalized','position',[0.1,0.2,0.3,0.4])
semilogy([0;best.N],[max(mdl{1}.Y);best.y]-fstar,'k-');
hold on
xlabel('Function evaluations');
ylabel('Absolute error');
axis([0 150 1e-8 1e2])
xticks([0 50 100 150])
yticks([1e-8 1e-6 1e-4 1e-2 1e0 1e2])
title('Comparison')
saveas(gcf,'Q2P-CR','epsc')

% Bayesian Optimization (BO)
% BO-LCB
fun = @obj_bo; 
x1 = optimizableVariable('x1',[-5,5]);
x2 = optimizableVariable('x2',[-5,5]);
vars = [x1,x2];
results = bayesopt(fun,vars,'Verbose',0,'PlotFcn',[],...
    'AcquisitionFunctionName','lower-confidence-bound',...
    'MaxObjectiveEvaluations',150);

REC_BO=(1:results.NumObjectiveEvaluations)';
REC_BO=[REC_BO results.ObjectiveTrace];
for i=1:results.NumObjectiveEvaluations
    REC_BO(i,3)=min(REC_BO(1:i,2));
end
semilogy(REC_BO(:,1),REC_BO(:,3)-fstar,'k:','LineWidth',0.75);

% BO-EI
vars = [x1,x2];
results = bayesopt(fun,vars,'Verbose',0,'PlotFcn',[],...
    'AcquisitionFunctionName','expected-improvement',...
    'MaxObjectiveEvaluations',150);

REC_BO=(1:results.NumObjectiveEvaluations)';
REC_BO=[REC_BO results.ObjectiveTrace];
for i=1:results.NumObjectiveEvaluations
    REC_BO(i,3)=min(REC_BO(1:i,2));
end
semilogy(REC_BO(:,1),REC_BO(:,3)-fstar,'k--');

% BO-PI
vars = [x1,x2];
results = bayesopt(fun,vars,'Verbose',0,'PlotFcn',[],...
    'AcquisitionFunctionName','probability-of-improvement',...
    'MaxObjectiveEvaluations',150);

REC_BO=(1:results.NumObjectiveEvaluations)';
REC_BO=[REC_BO results.ObjectiveTrace];
for i=1:results.NumObjectiveEvaluations
    REC_BO(i,3)=min(REC_BO(1:i,2));
end
semilogy(REC_BO(:,1),REC_BO(:,3)-fstar,'k-.');
legend('CO','BO-LCB','BO-EI','BO-PI','Location','southwest')
hold off
saveas(gcf,'Q2P-C','epsc')

figure('units','normalized','position',[0.1,0.2,0.3,0.4])
scatter(results.XTrace.x1,results.XTrace.x2,25,'k.')
axis([-5 5 -5 5])
box on
title('x trace of BO-LCB')
saveas(gcf,'Q2P-BOx','epsc')

function f = obj_bo(x)
global a;
x = [x.x1 x.x2];
f = a*x(:,1).^2+x(:,2).^2-10*exp(-4*(x(:,1)-0.5).^2-4*(x(:,2)+2).^2);
end