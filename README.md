clc;
close all;
clear all;
T = 1;
t=-4*T:0.01*T:4*T;
x=sinc(pi*t/T); 
Fs=1;
subplot(3,1,1);
title('Sinc function');
[N,X, fvals] = FFTAnalog(x, Fs, 'single-sided");
subplot(3,1,2);
plot(fvals,abs(X(1:N/2)));
title('Single-Sided Amplitude Spectrum of x(t))
xlabel('f (Hz)');
ylabel('|Px)|');
grid on;
[N.X, fvals] = FFTAnalog(x, Fs,'double-sided");
subplot(3,1,3);
plot(fvals,fftshift(abs(X))):
title('Double-Sided Amplitude Spectrum of x(t)');
xlabel('f (Hz)');
ylabel('|Px)|');
grid on;
function [N.X, fvals] FFTAnalog(x, Fs,type)
L = length(x):
N = 1;
for i=0:100
if L>=2^i
continue 
elseif L<2^i
N = 2^i;
break
end
end
X = fft(x,N);
if strcmp(type,'single-sided") fvals Fs (0:(N/2)-1)/N;
else
fvals Fs (1/N) (-1 (N/2):(N/2)-1);
end
end





clear all,
close all;
T=1;
L=100;
alpha=[0,0,4,0,6,0,8,1);
t=-4*T:1/L:4*T;
N=2^nextpow2(length(t));
RV=[];
for i=1:length(alpha)
figure();
[time,k1]=RCfunc(t,Talpha(i).L.N);
plot(t,k1);
title('roll of alpha(i));
grid on
RV(i,:)k1; 
figure();
subplot(211);
[f2,fftr2]=FFT_Analog(RCspec(time,T,alpha(i)),1/L,'double-sided");
plot(f2,fftr2) 
grid on
subplot(212)
[F1,fftr]=FFT_Analog(RV(i,:),1/L,'double-sided');
plot(f1,fftr)
grid on
end
function specvals = RCspec(f,T,alp)
specvals=[]:
specvals(f<=(1-alp)/(2*T))=T;
specvals((f>=(1-alp)/(2*T)) & (f<=(1-alp)/(2*T)))=T/2*(1+cos(pi*T/alp (f((f>=(1-alp)/(2*T)) & (f<=(1-alp)/(2*T)))-(1-alp)/(2*T))});
specvals((1+alp)/(2*T))=0;
end
function [tvals,rvals]=RCfunc(t,T,alp,L,N) 
rvals (double(sin(t*pi/T)).*cos(pi*alp*V/T))./((t*pi/T).*(1-(2*alp*t/T)));
tvals=-N/2:N/2-1; 
rvals(t==0)=1;
rvals (t==abs(T/(2*alp)))=alp*sin(pi/(2*alp))/2;
end
function [freq,fvals]=FFT_Analog(x,Fs,type)
N=2^nextpow2(length(x));
X=abs(fftshift(fft(x,N))):
X=X/sum(X);
if type="single-sided"
fvals=X(1:N/2);
freq=Fs*(0:N/2-1)/(N); 
elseif type="double-sided"
fvals=X;
end
freq=Fs*(-N/2:(N/2)-1)/(N/2);
end
