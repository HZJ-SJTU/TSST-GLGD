clc; clear; close all; 
%Parameters
num_dB = 7;
N = 2000;
fs = 200;
t = (0:N-1)/fs;
f = (0:N/2)*fs/N;
%Mode1
A_f1 = exp(0.008*f);
Phi_f1 = sin(6*pi*f.*f/10000)+2.5*(f);
GD_t1 = 12*pi*f/10000.*cos(6*pi*f.*f/10000)+2.5;
X1 = A_f1.*exp(-1i*2*pi*Phi_f1);
X1(end) = -A_f1(end);
Y1 = [X1  conj(fliplr(X1(2:end-1)))];    
y1 = ifft(Y1);
%Mode2
A_f2 = exp(0.005*f);
Phi_f2 = 4*f+0.07/2*f.^2-0.0007/3*f.^3;
GD_t2 = 4+0.07*f-0.0007*f.^2;
X2 = A_f1.*exp(-1i*2*pi*Phi_f2);
X2(end) = -A_f1(end);
Y2 = [X2  conj(fliplr(X2(2:end-1)))];    
y2 = ifft(Y2);
%Test Signal
ideal_signal = y1+y2;
noise_data = zeros(num_dB,N);
noise_dB = zeros(num_dB,1);
for k =1:num_dB
    noise_dB(k) = (k-1)*5;
    y = ideal_signal;
    y = awgn(y,noise_dB(k),'measured'); 
    noise_data(k,:) = y-ideal_signal;
end
save('noise.mat','noise_data','noise_dB')