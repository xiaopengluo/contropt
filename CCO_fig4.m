% Contractility of continuous optimization
% by Xiaopeng Luo and Xin Xu 
% 
% Fig. 4
% 
% A typical logarithmic time contractile example

clear; clc;

obj = @(x)(x-pi).^2 + sin(50*(x-pi)).^2;

t = (0:0.01:2*pi)';
figure('units','normalized','position',[0.1,0.2,0.3,0.4])
plot(t,obj(t),'k')
xticks([0 0.5*pi pi 1.5*pi 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
title('f')
saveas(gcf,'HLFDF-f','epsc')

N = 2; maxN = 500; minL = 5; K = 3; 
lb = 0; ub = 7;
type={'median','low',0.4};
[best,mdl,thr] = contropt(obj,type,K,minL,N,maxN,lb,ub);

figure('units','normalized','position',[0.3,0.2,0.3,0.4])
scatter(mdl{1}.X,0*mdl{1}.Y,25,'k.')
hold on
for i=2:sum(thr(:,1)~=0)
    scatter(mdl{i}.X,0*mdl{i}.Y,25,'k.')
end
hold off
box on
axis([0 2*pi -1 1])
xticks([0 0.5*pi pi 1.5*pi 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
title('x trace')
saveas(gcf,'HLFDF-x','epsc')

figure('units','normalized','position',[0.5,0.2,0.3,0.4])
semilogy([0;best.N],[max(mdl{1}.Y);best.y],'k-');
hold on
plot([0;best.N],[max(mdl{1}.Y);best.y],'b.');
hold off
xlabel('Function evaluations');
ylabel('Absolute error');
axis([0 100 1e-8 1e2])
yticks([0 20 40 60 80 100])
yticks([1e-8 1e-6 1e-4 1e-2 1e0 1e2])
title('Convergence')
saveas(gcf,'HLFDF-cr','epsc')