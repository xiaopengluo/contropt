% Contractility of continuous optimization
% by Xiaopeng Luo and Xin Xu 
% 
% Fig. 8
% 
% Three examples for noncontractile functions
clear; clc

f1 = @(x) 2*rand(size(x))-1;
f2 = @(x) sin(50*x);
f3 = @(x) abs(x).^2 - exp(-10000*(x+.7).^2);

n = 1000;
p = haltonset(1);
x = 2*net(p,n)-1;
x = sort(x);

y1 = f1(x);
I = ~isnan(y1);
x1 = x(I);
y1 = y1(I);
F1 = dct1d(y1);

y2 = f2(x);
I = ~isnan(y2);
x2 = x(I);
y2 = y2(I);
F2 = dct1d(y2);

y3 = f3(x);
I = ~isnan(y3);
x3 = x(I);
y3 = y3(I);
F3 = dct1d(y3);

figure('units','normalized','position',[0.1,0.1,0.3,0.4])
plot(x1,y1)
axis([-1 1 -1.5 1.5])
saveas(gcf,['non-' num2str(1)],'epsc')

figure('units','normalized','position',[0.3,0.1,0.3,0.4])
histogram(y1)
axis([-1 1 0 150])
saveas(gcf,['non-' num2str(2)],'epsc')

figure('units','normalized','position',[0.5,0.1,0.3,0.4])
plot(F1)
saveas(gcf,['non-' num2str(3)],'epsc')

figure('units','normalized','position',[0.1,0.2,0.3,0.4])
plot(x2,y2)
axis([-1 1 -1.5 1.5])
saveas(gcf,['non-' num2str(4)],'epsc')

figure('units','normalized','position',[0.3,0.2,0.3,0.4])
histogram(y2)
axis([-1.2 1.2 0 200])
saveas(gcf,['non-' num2str(5)],'epsc')

figure('units','normalized','position',[0.5,0.2,0.3,0.4])
plot(F2)
saveas(gcf,['non-' num2str(6)],'epsc')

figure('units','normalized','position',[0.1,0.3,0.3,0.4])
plot(x3,y3)
saveas(gcf,['non-' num2str(7)],'epsc')

figure('units','normalized','position',[0.3,0.3,0.3,0.4])
histogram(y3)
axis([-0.6 1 0 350])
saveas(gcf,['non-' num2str(8)],'epsc')

figure('units','normalized','position',[0.5,0.3,0.3,0.4])
plot(F3)
saveas(gcf,['non-' num2str(9)],'epsc')

function data=dct1ds(data)
% computes the discrete cosine transform of the column vector data
[nrows,~]= size(data);
% Compute weights to multiply DFT coefficients
weight = [1;2*(exp(-1i*(1:nrows-1)*pi/(2*nrows))).'];
% Re-order the elements of the columns of x
data = [ data(1:2:end,:); data(end:-2:2,:) ];
% Multiply FFT by weights:
data= real(weight.* fft(data));
end