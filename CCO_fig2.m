% Contractility of continuous optimization
% by Xiaopeng Luo and Xin Xu 
% 
% Fig. 2
% 
% 2D illustration of the performance of the median-type algorithm

clear;clc;

obj = @(x) (x(:,2)-5.1*x(:,1).^2/(4*pi^2) + ...
    5*x(:,1)/pi-6).^2+10*(1-1/(8*pi))*cos(x(:,1))+10;
lb = [-5 0]; ub = [10 15];

N = 2; maxN = 500; minL = 5; K = 7; 
type={'median','medium',0.25};
[best,mdl,thr] = contropt(obj,type,K,minL,N,maxN,lb,ub);

alphbet = ['A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N'];
figure('units','normalized','position',[0.1+6/20,0.3,0.3,0.4])
Nk = length(mdl{1}.Y);
scatter(mdl{1}.X(:,1),mdl{1}.X(:,2),25,'k.')
hold on
for i=2:sum(thr(:,1)~=0)
    scatter(mdl{i}.X(:,1),mdl{i}.X(:,2),25,'k.')
    Nk = [Nk; length(mdl{i}.Y)];
end
hold off
box on
saveas(gcf,['2-' num2str(K+1)],'epsc')

figure('units','normalized','position',[0.1+7/20,0.3,0.3,0.4])
fstar = obj([pi 2.275]);
semilogy([0;best.N],[max(mdl{1}.Y);best.y]-fstar);
hold on
semilogy([0;best.N],[max(mdl{1}.Y);thr(:,5)]-fstar,'r');
plot([0;best.N],[max(mdl{1}.Y);best.y]-fstar,'b.');
plot([0;best.N],[max(mdl{1}.Y);thr(:,5)]-fstar,'r.');
hold off
legend('ERR(N)','UB(N)');
xlabel('number of function evaluations (N)');
axis([0 350 1e-12 1e4])
yticks([1e-12 1e-8 1e-4 1e0 1e4])
saveas(gcf,['2-' num2str(K+2)],'epsc')

figure('units','normalized','position',[0.1+8/20,0.3,0.3,0.4])
plot((0:K-1)',Nk);
hold on
plot((0:K-1)',Nk,'k.');
hold off
ylabel('N^{(k)}');
xlabel('k');
axis([0 K-1 20 100])
saveas(gcf,['2-' num2str(K+3)],'epsc')

function [best,mdl,thr] = contropt(obj,type,K,minL,N,maxN,lb,ub,cons)
%CONTROPT finds a global minima of a continuous obj on a compact set 
%   by using a contracting algorithm. 
%   
%   CONTROPT attempts to solve problems of the form:
%    min obj(X) subject to: lb <= X <= ub (bounds)
%     X                     cons(X) == 1  (nonlinear constraints)
%   
%   type  -  cell(3,1): type{1} - three contracting types: 'median', 
%            'mean', or 'middle'; type{2} - three levels of the residual 
%            bound: 'low','medium' or 'high'; type{3} - bound of the 
%            sampling ratio deviation
%   level -  levels of the residual bound: 'low','medium' or 'high'
%   K     -  maximum number of contractions
%   minL  -  minimum number of inner loops (i.e., explorations)
%   N     -  number of samples per update
%   maxN  -  maximum number for reducing the contracting condition


%   References
%   X Luo, X Xu, "Contractility of Continuous Optimization", 
%   https://arxiv.org/abs/1909.01123.

%   Xiaopeng Luo and Xin Xu, 7/24/2019

%   global variable
global cr MUL;
%   initialization
if nargin > 8
    px = net(haltonset(length(lb)),tN);
    px = bsxfun(@plus,bsxfun(@times,px,(ub-lb)),lb);
    cr = tN / sum(cons(px));
else
    cr = 1; cons = [];
end
[TN,X,Y,sig] = deal(0,[],[],1);
[cn,mdl,thr] = deal(0,cell(1,1),[0 1 0 0 0]);
best         = struct('x',[],'y',[],'N',[]);
%   start the contractions until the tolerance is reached
k  = 0; 
MUL = 100;
CX  = [];
t0 = norminv(net(scramble(haltonset(length(lb),'Skip',1),'RR2'),1e5));
while k < K
    iter = 0;
    cut  = -Inf;
    while 1
        iter         = iter + 1;
        %   sampling in D^(k)
        [x,cn,cx]    = samplingDk(cn,k,N,X,lb,ub,mdl,thr,cons);
        if isempty(cx), return; end
        [X,Y,TN]     = updateData(obj,x,X,Y,TN);
        CX           = [CX; cx];
        
        % normalization, use the factor 20 here
        if iter<=minL
            [E,STD]      = deal(mean(Y),std(Y)/1);
        end
        YY           = (Y-E) / STD;
        
        %   modeling in D^(k)
        [Mdl,ek]     = modelingDk(X,YY,sig,type);
        if ~isempty(Mdl)
            %   calculate the judgment parameters
            x        = CX;
            y        = predict(Mdl,x);
            
            [~,id]   = min(y);
            if k==0
                Yn       = obj(x(id,:));
                [X,Y,TN] = deal([X; x(id,:)],[Y; Yn],TN+1);
                YY       = [YY; (Yn-E)/STD];
            else
                sD       = sort(sqrt(sum(bsxfun(@minus,X,x(id,:)).^2,2)));
                sc       = sD(min(2*length(lb)+1,length(sD)))/3;
                s0       = bsxfun(@plus,sc*t0,x(id,:));
                s0       = inDk(s0,k,lb,ub,mdl,thr,cons);
                sy       = predict(Mdl,s0);
                [~,id]   = min(sy);
                Yn       = obj(s0(id,:));
                [X,Y,TN] = deal([X; s0(id,:)],[Y; Yn],TN+1);
                YY       = [YY; (Yn-E)/STD];
            end
            
            %   get uk
            [uk,cut,rc]  = get_uk(iter,minL,maxN,cut,YY,ek,type);
            
            %   issuing a judgment
            Iy       = (y(length(y)-size(cx,1)+1:end)<=cut);
            sr       = sum(Iy) / size(cx,1);
            rd       = abs(sr-rc);
            if ek <= uk && rd < type{3} && iter>=minL
                I        = (YY<=cut); 
                [X,Y]    = deal(X(I,:),Y(I,:));
                mdl{k+1} = Mdl;
                sig      = mdl{k+1}.KernelInformation.KernelParameters(1);
                [cb,id]  = min(Y);
                best.x(k+1,:) = X(id,:);
                best.y(k+1,:) = cb;
                best.N(k+1,:) = TN;
                thr(k+1,:)    = [cut sr E STD STD*cut+E];
                CX            = [];
                fprintf(' - iteration %d with the ub = %d;\n',k,thr(k+1,5));
                stepshow(k,obj,lb,ub,mdl,thr,x,y);
                k             = k + 1;
                break;
            end
        end
    end
end
end  %end main function
%   subfunction: 
%   sampling in Dk
function [dx,cn,T] = samplingDk(cn,k,N,X,lb,ub,mdl,thr,cons)
if k==0
    [dx,cn,T] = sampling0(cn,k,N,lb,ub,mdl,thr,cons);
else
    [dx,cn,T] = samplingX(cn,k,N,X,lb,ub,mdl,thr,cons);
end
end
function [dx,cn,T] = sampling0(cn,k,N,lb,ub,mdl,thr,cons)
global cr MUL;
[dx,SN] = deal([],ceil(N/cr));
while SN > 0
    p   = haltonset(length(lb),'Skip',cn);
    tx  = net(scramble(p,'RR2'),SN); 
    tx  = bsxfun(@plus,bsxfun(@times,tx,(ub-lb)),lb);
    dx  = [dx; inDk(tx,k,lb,ub,mdl,thr,cons)];
    cn  = cn + SN;
    SN  = ceil((N-size(dx,1))/cr);
end
T      = [];
cn1    = cn;
SN     = ceil(MUL*N/cr);
while SN > 0
    p   = haltonset(length(lb),'Skip',1e3+cn1,'Leap',1e2);
    tx  = net(scramble(p,'RR2'),SN); 
    tx  = bsxfun(@plus,bsxfun(@times,tx,(ub-lb)),lb);
    T   = [T; inDk(tx,k,lb,ub,mdl,thr,cons)];
    cn1 = cn1 + SN;
    SN  = ceil((MUL*N-size(T,1))/cr);
end
end
function [dx,cn,T] = samplingX(cn,k,N,X,lb,ub,mdl,thr,cons)
global MUL;
[XN,n]  = size(X);
sn      = min(XN,2*n+1);
T       = [];
needN   = max(10,MUL*N);
SN      = min(1e6,5*ceil(needN/XN));
X_NBR   = zeros(size(X,1),sn);
parfor i=1:XN
    [~,I]      = sort(sqrt(sum(bsxfun(@minus,X,X(i,:)).^2,2)));
    X_NBR(i,:) = I(1:sn);
end
iter = 0;
while size(T,1)<MUL*N && iter<1000
    iter   = iter + 1;
    if log2(cn)>53, cn = 1; end
    p      = haltonset(sn,'Skip',cn);
    W      = net(scramble(p,'RR2'),SN);
    W      = W.^2/2;
    W      = [-W; W] + 1/2;
    W      = W./sum(W,2);
    cn     = cn + SN;
    tT     = cell(XN,1);
    parfor i=1:XN
        t     = W*X(X_NBR(i,:),:);
        tm    = mean(X(X_NBR(i,:),:));
        dX    = sqrt(mean(sum(bsxfun(@minus,X(X_NBR(i,:),:),tm).^2,2)));
        dt    = sqrt(mean(sum(bsxfun(@minus,t,tm).^2,2)));
        t2    = 1.5*dX/dt*(t-tm)+tm;
        tT{i} = inDk(t2,k,lb,ub,mdl,thr,cons);
    end
    tT     = cell2mat(tT);
    T      = [T; tT];
    needN  = max(10,MUL*N-size(T,1));
    SN0    = 5*ceil(needN/XN);
    SN     = min(1e6,SN0);
end
T      = unique(T,'rows');
if isempty(T), return; end
D      = sum(bsxfun(@minus,T,X(1,:)).^2,2);
PT     = ones(size(T,1),1);
for i=2:XN
    [D,C]    = min([D sum(bsxfun(@minus,T,X(i,:)).^2,2)],[],2);
    PT(C==2) = i;
end
I      = [];
for i=1:N
    [~,id] = max(D);
    [D,C]  = min([D sum(bsxfun(@minus,T,T(id,:)).^2,2)],[],2);
    [PT(C==2),I] = deal(XN+1,[I;id]);
end
dx     = T(I,:);
dx     = dx(1:N,:);
end
%   select samples in Dk
function [dx,I] = inDk(dx,k,lb,ub,mdl,thr,cons)
I = (sum(bsxfun(@lt,dx,lb),2)==0) & (sum(bsxfun(@gt,dx,ub),2)==0);
if ~isempty(cons)
    I = I & (cons(dx)==1);
end
ks = max(k-1,1);
for i=ks:k
    I = I & (predict(mdl{i},dx)<thr(i,1));
end
dx = dx(I,:);
end
%   create a model and estimate the error bound ek on Dk
function [mdl,ek] = modelingDk(X,Y,t,type)
warning('off');
switch type{2}
    case 'low',     Ce = 2;
    case 'medium',	Ce = 3;
    case 'high',	Ce = 5;
end
ek = 0;
[mdl,s]   = deal([],std(Y));
[XN,n]    = size(X);
sn        = min(XN,2*n+1);
cc        = zeros(XN,1);
for i=1:XN
    tc    = sort(sqrt(sum(bsxfun(@minus,X,X(i,:)).^2,2)));
    cc(i) = tc(sn)/2;
end
mc        = mean(cc);
[m,CV,L]  = deal(cell(2,1),cell(2,1),[Inf;Inf]);
try
    m{1}  = fitrgp(X,Y,'KernelFunction','squaredexponential',...
        'KernelParameters',[mc;s],'Sigma',s/2,'Standardize',1);
    CV{1} = crossval(m{1});
    L(1)  = sqrt(kfoldLoss(CV{1}));
end
try
    m{2}  = fitrgp(X,Y,'KernelFunction','squaredexponential',...
        'KernelParameters',[t;s],'Sigma',s/2,'Standardize',1);
    CV{2} = crossval(m{2});
    L(2)  = sqrt(kfoldLoss(CV{2}));
end
[L,id]    = min(L);
if L<Inf
    eY    = kfoldPredict(CV{id});
    eI    = ~isnan(eY);
    mu    = abs(mean(eY(eI)-Y(eI)));
    ek    = mu + Ce*L;
    mdl   = m{id};
end
end
function [uk,thr,ratio] = get_uk(iter,minL,maxN,thr,Y,ek,type)
NY      = length(Y);
Y0      = Y((Y-mean(Y))<2*std(Y));
switch type{1}
    case 'median'
        if iter<=minL
            thr = median(unique(Y0));
        end
        J1 = NY>=1.0*maxN && min(Y)+ek<prctile(Y0,75);
        J2 = NY>=1.5*maxN && min(Y)+ek<prctile(Y0,90);
        J3 = NY>=2.0*maxN && min(Y)+ek<prctile(Y0,95);
        if J1 || J2 || J3
            thr = min(Y) + ek;
        end
        uk    = thr - min(Y);
        ratio = sum(Y<=thr)/NY;
    case 'mean'
        if NY>=2*maxN && min(Y)+ek<prctile(Y0,95)
            thr = min(Y) + ek;
        end
        if NY>=1.5*maxN && min(Y)+ek<prctile(Y0,90)
            thr = min(Y) + ek;
        end
        if iter<=minL
            thr = mean(unique(Y0));
        end
        uk    = thr - min(Y);
        ratio = sum(Y<=thr)/NY;
    case 'middle'
        uk    = ( max(Y) - min(Y) )/2;
        thr   = min(Y) + ek;
        ratio = sum(Y<=thr)/NY;
end
end
function [X,Y,TN] = updateData(obj,x,X,Y,TN)
[y,Np] = deal(obj(x),size(X,1));
[x,y]  = deal(x(~isnan(y),:),y(~isnan(y)));
[X,Y]  = deal([X; x],[Y; y]);
[X,ix] = unique(X,'rows');
nN     = size(x,1) - (size(x,1)+Np-size(X,1));
[Y,TN] = deal(Y(ix),TN+nN);
end

function stepshow(k,obj,lb,ub,mdl,thr,x,y)
if length(lb)==2
alphbet = ['A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N'];
figure('units','normalized','position',[0.1,0.2,0.3,0.4])
scatter(x(y<=thr(k+1),1),x(y<=thr(k+1),2),2,'o')
hold on
scatter(mdl{k+1}.X(:,1),mdl{k+1}.X(:,2),25,'k.')
t = [(ub(1)-lb(1))/200 (ub(2)-lb(2))/200];
[mx,my] = meshgrid(lb(1):t(1):ub(1),lb(2):t(2):ub(2));
z = reshape(obj([mx(:) my(:)]),size(mx));
level = [1 2 3 10 20 30 60 120 180];
TL    = [3 10 30 60 120 180];
contour(mx,my,z,level,'ShowText','on','TextList',TL);
box on
hold off
saveas(gcf,['2-' num2str(k+1)],'epsc')
end
end