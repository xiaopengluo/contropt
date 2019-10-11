% Contractility of continuous optimization
% by Xiaopeng Luo and Xin Xu 
% 
% Fig. 7
% 
% Comparison of CO to BO on 2D Rosenbrock function
clear; clc;

obj = @(x)100*(x(:,2)-x(:,1).^2).^2+(x(:,1)-1).^2;
fstar = 0;

N = 2; maxN = 500; minL = 5; K = 7; 
lb = [-5 -5]; ub = [5 5];
type={'median','low',0.4};
[best,mdl,thr] = contropt(obj,type,K,minL,N,maxN,lb,ub);

figure('units','normalized','position',[0.1,0.2,0.3,0.4])
semilogy([0;best.N],[max(mdl{1}.Y);best.y]-fstar,'k-');
hold on
xlabel('Function evaluations');
ylabel('Absolute error');
axis([0 600 1e-9 1e3])
xticks([0 200 400 600])
yticks([1e-9 1e-6 1e-3 1e0 1e3])
title('Comparison')

% Bayesian Optimization (BO)
% BO-LCB
fun = @obj_bo; 
x1 = optimizableVariable('x1',[-5,5]);
x2 = optimizableVariable('x2',[-5,5]);
vars = [x1,x2];
results1 = bayesopt(fun,vars,'Verbose',0,'PlotFcn',[],...
    'AcquisitionFunctionName','lower-confidence-bound',...
    'MaxObjectiveEvaluations',600);

REC_BO1=(1:results1.NumObjectiveEvaluations)';
REC_BO1=[REC_BO1 results1.ObjectiveTrace];
for i=1:results1.NumObjectiveEvaluations
    REC_BO1(i,3)=min(REC_BO1(1:i,2));
end
semilogy(REC_BO1(:,1),REC_BO1(:,3)-fstar,'k:','LineWidth',0.75);

% BO-EI
vars = [x1,x2];
results2 = bayesopt(fun,vars,'Verbose',0,'PlotFcn',[],...
    'AcquisitionFunctionName','expected-improvement',...
    'MaxObjectiveEvaluations',600);

REC_BO2=(1:results2.NumObjectiveEvaluations)';
REC_BO2=[REC_BO2 results2.ObjectiveTrace];
for i=1:results2.NumObjectiveEvaluations
    REC_BO2(i,3)=min(REC_BO2(1:i,2));
end
semilogy(REC_BO2(:,1),REC_BO2(:,3)-fstar,'k--');

% BO-PI
vars = [x1,x2];
results3 = bayesopt(fun,vars,'Verbose',0,'PlotFcn',[],...
    'AcquisitionFunctionName','probability-of-improvement',...
    'MaxObjectiveEvaluations',600);

REC_BO3=(1:results3.NumObjectiveEvaluations)';
REC_BO3=[REC_BO3 results3.ObjectiveTrace];
for i=1:results3.NumObjectiveEvaluations
    REC_BO3(i,3)=min(REC_BO3(1:i,2));
end
semilogy(REC_BO3(:,1),REC_BO3(:,3)-fstar,'k-.');
legend('CO','BO-LCB','BO-EI','BO-PI','Location','southwest')
hold off
saveas(gcf,'PTC-CR','epsc')

figure('units','normalized','position',[0.1,0.2,0.3,0.4])
scatter(mdl{1}.X(:,1),mdl{1}.X(:,2),25,'k.')
Nk = length(mdl{1}.Y);
hold on
for i=2:sum(thr(:,1)~=0)
    scatter(mdl{i}.X(:,1),mdl{i}.X(:,2),25,'k.')
    Nk = [Nk; length(mdl{i}.Y)];
end
hold off
box on
axis([-5 5 -5 5])
title('x trace of CO')
saveas(gcf,'PTC-COx','epsc')

figure('units','normalized','position',[0.1,0.2,0.3,0.4])
plot((0:length(best.N)-1)',Nk);
hold on
plot((0:length(best.N)-1)',Nk,'b.');
hold off
box on
xlabel('k-th model');
ylabel('Sample size');
axis([0 length(best.N)-1 0 300])
saveas(gcf,'PTC-size','epsc')

function f = obj_bo(x)
x = [x.x1 x.x2];
f = 100*(x(:,2)-x(:,1).^2).^2+(x(:,1)-1).^2;
end