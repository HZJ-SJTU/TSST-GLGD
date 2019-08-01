    clc; clear; close all; 
    figure_num = 5;%3 or 5
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
    if figure_num == 5
        load('noise.mat');
        y = y + noise_data(1,:);
    end
    GD_t(1,:)=GD_t1;
    GD_t(2,:)=GD_t2;
    ideal_y(1,:)=y1;
    ideal_y(2,:)=y2;
%% Parameters
    map = jet ;
    %Window Parameter
    WindowOpt = struct('type','gauss','s',0.10);
    %Frequency axis Parameter
    Parameter = struct('L',N/2+1,'fmin',0,'fmax',fs/2);
%% S-method
    t =1:5:2000;
    f_num = 512;
	G=tftb_window(9,'hanning'); 
    h=tftb_window(61,'hanning'); 
    [~,rtfr,~] = tfrrstan(y',t,f_num,G,h,1);
    t = (t-1)/fs;
    f = (0:f_num/2-1)/f_num*fs;
    rtfr = rtfr(1:f_num/2,:);
    figure;
    imagesc(t,f,rtfr);axis xy;colormap(map);axis tight; 
    xlabel({'Time (s)','(a)'},'FontSize',20);set(gca,'XTick',0:2:10);
    ylabel('Frequency (Hz)','FontSize',20); set(gca,'YTick',0:20:100);
    set(gca,'FontSize',20);axis([0 10 0 100])
    r1 = renyi(rtfr,t,f',3);
    text(5.25,10,['RE:',num2str(r1,3),'(bits)'],'Color','black','FontWeight','bold','FontSize',26) 
%%  STFT 
%     [Wx,t1,f1,~] = stft(y, fs, WindowOpt, Parameter, 'modify');
%     figure;
%     imagesc(t1,f1,abs(Wx));axis xy;
%     colormap(map)
%     axis tight; xlabel({'Time (s)','(a)'},'FontSize',20);set(gca,'FontSize',20);
%     ylabel('Fre (Hz)','FontSize',20);  
%     r1 = renyi(abs(Wx),t1,f1',3);
%     text(5.5,10,['R¨¦nyi:',num2str(r1,3)],'Color','black','FontWeight','bold','FontSize',26)  
%% the discrete-time Levenberg-Marquard synchrosqueezed Gabor Transform    
    f_num = 512;
    t = (0:N-1)/fs;
    [~, stfr, ~] = tfrvsgab(y, f_num, 10);
    stfr = stfr(1:f_num/2,:);
    f = (0:f_num/2-1)/f_num*fs;
    figure;
    imagesc(t,f,abs(stfr));axis xy;colormap(map);axis tight; 
    xlabel({'Time (s)','(b)'},'FontSize',20);set(gca,'XTick',0:2:10);
    ylabel('Frequency (Hz)','FontSize',20); set(gca,'YTick',0:20:100);
    set(gca,'FontSize',20);axis([0 10 0 100])
    r2 = renyi(abs(stfr),t,f',3);
    text(5.25,10,['RE:',num2str(r2,3),'(bits)'],'Color','black','FontWeight','bold','FontSize',26) 
%%  TSST
    [Ts,t,f,xMean,~] = tsst(y , fs,  WindowOpt, Parameter, '1Ord');
    figure;
    imagesc(t,f,abs(Ts));axis xy;colormap(map);axis tight; 
    xlabel({'Time (s)','(c)'},'FontSize',20);set(gca,'XTick',0:2:10);
    ylabel('Frequency (Hz)','FontSize',20); set(gca,'YTick',0:20:100);
    set(gca,'FontSize',20);axis([0 10 0 100])
    r3 = renyi(abs(Ts),t,f',3);
    text(5.25,10,['RE:',num2str(r3,3),'(bits)'],'Color','black','FontWeight','bold','FontSize',26) 
%% TSST2
    [Ts,t,f,xMean,~] = tsst(y , fs,  WindowOpt, Parameter, '2Ord');
    figure;
    imagesc(t,f,abs(Ts));axis xy;colormap(map);axis tight; 
    xlabel({'Time (s)','(d)'},'FontSize',20);set(gca,'XTick',0:2:10);
    ylabel('Frequency (Hz)','FontSize',20); set(gca,'YTick',0:20:100);
    set(gca,'FontSize',20);axis([0 10 0 100])
    r4 = renyi(abs(Ts),t,f',3);
    text(5.25,10,['RE:',num2str(r4,3),'(bits)'],'Color','black','FontWeight','bold','FontSize',26) 
%% Close File
rmpath(genpath(folder1));
rmpath(genpath(folder2));
rmpath(genpath(folder3));
rmpath(genpath(folder4));
rmpath(genpath(folder5));
rmpath(genpath(folder6));