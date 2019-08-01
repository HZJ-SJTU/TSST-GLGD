    clc; clear; close all; 
%% Open File
    folder1 = 'ridge';
    addpath(genpath(folder1));
    folder2 = 'other';
    addpath(genpath(folder2));
    folder3 = 'method';
    addpath(genpath(folder3));
    folder4 = 'indicator';
    addpath(genpath(folder4));
    folder5 = 'mfiles';
    addpath(genpath(folder5));
    folder6 = 'ST';
    addpath(genpath(folder6));
    load('noise.mat');
%% Test Signal
    %Parameters
    N = 2000;%2000
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
    y = y1 + y2;
    GD_t(1,:)=GD_t1;
    GD_t(2,:)=GD_t2;
    ideal_y(1,:)=y1;
    ideal_y(2,:)=y2;
%% Parameters
    exband = 20;%25
    map = jet ;
    %Window Parameter
    WindowOpt = struct('type','gauss','s',0.1);
    %Frequency axis Parameter
    Parameter = struct('L',N/2+1,'fmin',0,'fmax',fs/2);
%%  TSST 
    SNRoutput1 = zeros(7,1);
    r1 = zeros(7,1);
    for m = 1:7
        y = y1+y2+noise_data(m,:);
        [Ts1,t,f,xMean,~] = tsst(y , fs,  WindowOpt, Parameter, '1Ord');
        r1(m) = renyi(abs(Ts1),t,f',3);
        [L,N] = size(Ts1);
        for k = 1:2
            Ee = zeros(L,N);
            for i =1:L
                Ee(i,round(GD_t(k,i)*fs)-exband:round(GD_t(k,i)*fs)+exband) = 1;
            end
            Tx = Ts1.*Ee;
            b = itsst(Tx,fs, xMean);
            if k == 1
                SNRoutput1(m) = SNR(ideal_y(k,:),b);
            end
        end
    end    
%% TSST2
    SNRoutput2 = zeros(7,1);
    r2 = zeros(7,1);
    for m = 1:7
        y = y1+y2+noise_data(m,:);
        [Ts2,t,f,xMean,~] = tsst(y , fs,  WindowOpt, Parameter, '2Ord');
        r2(m) = renyi(abs(Ts2),t,f',3);
        [L,N] = size(Ts2);
        for k = 1:2
            Ee = zeros(L,N);
            for i =1:L
                Ee(i,round(GD_t(k,i)*fs)-exband:round(GD_t(k,i)*fs)+exband) = 1;
            end
            Tx = Ts2.*Ee;
            b = itsst(Tx,fs, xMean);
            if k == 1
                SNRoutput2(m) = SNR(ideal_y(k,:),b);
            end
        end
    end
%% S-method
	G=tftb_window(9,'hanning'); 
    h=tftb_window(61,'hanning'); 
    r3 = zeros(7,1);
    for m = 1:7
        y = y1+y2+noise_data(m,:);
        t =1:5:2000;
        f_num = 512;
        [~,rtfr,~] = tfrrstan(y',t,f_num,G,h,1);
        t = (t-1)/fs;
        f = (0:f_num/2-1)/f_num*fs;
        rtfr = rtfr(1:f_num/2,:); 
        r3(m) = renyi(rtfr,t,f',3);
    end
%% the discrete-time Levenberg-Marquard synchrosqueezed Gabor Transform      
    r4 = zeros(7,1);
    for m = 1:7
        y = y1+y2+noise_data(m,:);
        f_num = 512;
        t = (0:N-1)/fs;
        [~, stfr, ~] = tfrvsgab(y, f_num, 10);
        stfr = stfr(1:f_num/2,:);
        f = (0:f_num/2-1)/f_num*fs;
        r4(m) = renyi(abs(stfr),t,f',3);
    end
%% plot
    figure
    plot(noise_dB,r3,'-d')
    hold on
    plot(noise_dB,r4,'-^')
    plot(noise_dB,r1,'-o')
    plot(noise_dB,r2,'-p')  
    
    xlabel({'SNR (dB)','(a)'},'FontSize',20);set(gca,'XTick',0:5:30);
    ylabel('RE (bits)','FontSize',20); set(gca,'YTick',0:1:7);
    legend('RSM','SST2','TSST','TSST2');
    set(gca,'FontSize',20);axis([0 30 0 7])
    grid on

    figure
    plot(noise_dB,SNRoutput1,'-o')
    hold on
    plot(noise_dB,SNRoutput2,'-p')
    xlabel({'SNR (dB)','(b)'},'FontSize',20);set(gca,'XTick',0:5:30);
    ylabel('RQF (dB)','FontSize',20); set(gca,'YTick',8:2:20);
    legend('TSST','TSST2','Location','SouthEast');
    set(gca,'FontSize',20);axis([0 30 8 20])
    grid on
%% Close File
rmpath(genpath(folder1));
rmpath(genpath(folder2));
rmpath(genpath(folder3));
rmpath(genpath(folder4));
rmpath(genpath(folder5));
rmpath(genpath(folder6));