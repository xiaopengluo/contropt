% Contractility of continuous optimization
% by Xiaopeng Luo and Xin Xu 
% 
% Fig. 1
% 
% Two-dimensional illustration of the performance of 
% the recursive algorithm

clear; clc;

n      = 2;
p      = haltonset(2);
x      = net(p,50);
t      = norminv(x(2:21,:));

r      = 0.25;
theta  = 0:2*pi/3600:2*pi;
xs     = x(sqrt(sum(bsxfun(@minus,x,[0.5 0.5]).^2,2))<=r,:);

N      = size(xs,1);
s      = zeros(N,1);
S      = [];
for i=1:N
    D    = sort(sqrt(sum(bsxfun(@minus,xs,xs(i,:)).^2,2)));
    s(i) = D(2*n+1)/2;
    S    = [S; xs(i,:)+s(i)*t];
end

S      = S(sqrt(sum(bsxfun(@minus,S,[0.5 0.5]).^2,2))<=r,:);
SN     = size(S,1);
D  = sum(bsxfun(@minus,S,xs(1,:)).^2,2);
PT = ones(SN,1);
for i = 2:N
    [D,C] = min([D sum(bsxfun(@minus,S,xs(i,:)).^2,2)],[],2);
    PT(C==2) = i;
end
I  = [];
for i = 1:N
    [~,id] = max(D);
    [D,C] = min([D sum(bsxfun(@minus,S,S(id,:)).^2,2)],[],2);
    [PT(C==2),I] = deal(N+1,[I;id]);
end

figure('units','normalized','position',[0.1,0.2,0.3,0.4])
scatter(x(:,1),x(:,2),'.')
hold on
scatter(xs(:,1),xs(:,2),10,'k*')
plot(0.5+r*cos(theta),0.5+r*sin(theta),'k','Linewidth',0.5);
hold off
box on
axis([0 1 0 1])
% text(-0.15,1,'A','FontWeight','bold')
saveas(gcf,['S1-' num2str(1)],'epsc')

figure('units','normalized','position',[0.3,0.2,0.3,0.4])
scatter(xs(:,1),xs(:,2),10,'k*')
hold on
plot(0.5+r*cos(theta),0.5+r*sin(theta),'k','Linewidth',0.5);
for i=1:N
    scatter(xs(i,1)+s(i)*t(:,1),xs(i,2)+s(i)*t(:,2),10,'o')
end
hold off
box on
axis([0 1 0 1])
% text(-0.15,1,'B','FontWeight','bold')
saveas(gcf,['S1-' num2str(2)],'epsc')

figure('units','normalized','position',[0.5,0.2,0.3,0.4])
scatter(xs(:,1),xs(:,2),10,'k*')
hold on
plot(0.5+r*cos(theta),0.5+r*sin(theta),'k','Linewidth',0.5);
scatter(S(I,1),S(I,2),10,'bo')
hold off
box on
axis([0 1 0 1])
% text(-0.15,1,'C','FontWeight','bold')
saveas(gcf,['S1-' num2str(3)],'epsc')