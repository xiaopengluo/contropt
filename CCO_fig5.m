% Contractility of continuous optimization
% by Xiaopeng Luo and Xin Xu 
% 
% Fig. 5
% 
% Comparison of CO to BO and GD on 2D quadratic function
clear; clc;

global a;
a = 20; 

%{%
% Contraction Optimization
obj = @(x)a*x(:,1).^2+x(:,2).^2;
N = 2; maxN = 500; minL = 5; K = 6; 
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
title('x trace of Contraction Optimization')
saveas(gcf,'Q2-COx','epsc')

figure('units','normalized','position',[0.1,0.2,0.3,0.4])
fstar = obj([0 0]);
semilogy([0;best.N],[max(mdl{1}.Y);best.y]-fstar,'k-');
hold on
plot([0;best.N],[max(mdl{1}.Y);best.y]-fstar,'b.');
hold off
xlabel('Function evaluations');
ylabel('Absolute error');
axis([0 100 1e-12 1e4])
yticks([1e-12 1e-8 1e-4 1e0 1e4])
title('Contraction Optimization')
saveas(gcf,'Q2-CO','epsc')

% Bayesian Optimization (BO)
% BO-LCB
fun = @obj_bo; 
x1 = optimizableVariable('x1',[-5,5]);
x2 = optimizableVariable('x2',[-5,5]);
vars = [x1,x2];
results = bayesopt(fun,vars,'Verbose',0,'PlotFcn',[],...
    'AcquisitionFunctionName','lower-confidence-bound',...
    'MaxObjectiveEvaluations',100);

figure('units','normalized','position',[0.1,0.2,0.3,0.4])
scatter(results.XTrace.x1,results.XTrace.x2,25,'k.')
axis([-5 5 -5 5])
box on
title('x trace of BO-LCB')
saveas(gcf,'Q2-BOx','epsc')

REC_BO=(1:results.NumObjectiveEvaluations)';
REC_BO=[REC_BO results.ObjectiveTrace];
for i=1:results.NumObjectiveEvaluations
    REC_BO(i,3)=min(REC_BO(1:i,2));
end
figure('units','normalized','position',[0.1,0.2,0.3,0.4])
semilogy(REC_BO(:,1),REC_BO(:,3)-fstar,'k-');
xlabel('Function evaluations');
ylabel('Absolute error');
axis([0 100 1e-12 1e4])
yticks([1e-12 1e-8 1e-4 1e0 1e4])
hold on

% BO-EI
vars = [x1,x2];
results = bayesopt(fun,vars,'Verbose',0,'PlotFcn',[],...
    'AcquisitionFunctionName','expected-improvement',...
    'MaxObjectiveEvaluations',100);

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
    'MaxObjectiveEvaluations',100);

REC_BO=(1:results.NumObjectiveEvaluations)';
REC_BO=[REC_BO results.ObjectiveTrace];
for i=1:results.NumObjectiveEvaluations
    REC_BO(i,3)=min(REC_BO(1:i,2));
end
semilogy(REC_BO(:,1),REC_BO(:,3)-fstar,'k-.');
legend('LCB','EI','PI')
hold off
title('Bayesian Optimization')
saveas(gcf,'Q2-BO','epsc')
%}

% Gradient Descent
x = [5 5];
e = [.047 .045 .0475];
for i=1:99
    t = x(i,:);
    x(i+1,:) = t - e(1) * [2*a 2] .* t;
end

figure('units','normalized','position',[0.1,0.2,0.3,0.4])
scatter(x(:,1),x(:,2),25,'k.')
axis([-5 5 -5 5])
box on
title('x trace of Gradient Descent')
saveas(gcf,'Q2-GDx','epsc')

figure('units','normalized','position',[0.1,0.2,0.3,0.4])
semilogy((1:100),x.^2*[a;1],'k-');
xlabel('Function evaluations');
ylabel('Absolute error');
hold on

x = [5 5];
for i=1:99
    t = x(i,:);
    x(i+1,:) = t - e(2) * [2*a 2] .* t;
end
semilogy((1:100),x.^2*[a;1],'k--');
x = [5 5];
for i=1:99
    t = x(i,:);
    x(i+1,:) = t - e(3) * [2*a 2] .* t;
end
semilogy((1:100),x.^2*[a;1],'k-.');
hold off
axis([0 100 1e-12 1e4])
yticks([1e-12 1e-8 1e-4 1e0 1e4])
legend(['\eta=' num2str(e(1))],['\eta=' num2str(e(2))],['\eta=' num2str(e(3))])
title('Gradient Descent')
saveas(gcf,'Q2-GD','epsc')

function f = obj_bo(x)
global a;
x = [x.x1 x.x2];
f = a*x(:,1).^2+x(:,2).^2;
end