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
    exband = 25;
    map = jet ;
    %Window Parameter
    WindowOpt = struct('type','gauss','s',0.1);
    %Frequency axis Parameter
    Parameter = struct('L',N/2+1,'fmin',0,'fmax',fs/2);
%%  TSST
    [Ts,t,f,xMean,~] = tsst(y , fs,  WindowOpt, Parameter, '1Ord');
%% ITSST
    [L,N] = size(Ts);
    for k = 1:2
        Ee = zeros(L,N);
        for i =1:L
            Ee(i,round(GD_t(k,i)*fs)-exband:round(GD_t(k,i)*fs)+exband) = 1;
        end
        Tx = Ts.*Ee;
%         figure;
%         imagesc(t,f,abs(Tx));axis xy;
        b = itsst(Tx,fs, xMean);
        figure
        plot(t,b,'LineWidth',2);
        hold on
        plot(t,ideal_y(k,:)-b,'color','black','LineWidth',2);
        if k == 1
            xlabel({'Time (s)','(a)'},'FontSize',20);set(gca,'XTick',0:2:10);
        else
            xlabel({'Time (s)','(b)'},'FontSize',20);set(gca,'XTick',0:2:10);
        end
        ylabel('AMP','FontSize',20); set(gca,'YTick',-0.4:0.2:0.4);
        legend('Reconstruction','Error');
        set(gca,'FontSize',20);axis([0 10 -0.4 0.4])

        SNRoutput1(k) = SNR(ideal_y(k,:),b);
        text(5.5,-0.3,['RQF:',num2str(SNRoutput1(k),3),'dB'],'Color','black','FontWeight','bold','FontSize',26)
    end
%% TSST2
    [Ts,t,f,xMean,~] = tsst(y , fs,  WindowOpt, Parameter, '2Ord');
%% ITSST2
    [L,N] = size(Ts);
    for k = 1:2
        Ee = zeros(L,N);
        for i =1:L
            Ee(i,round(GD_t(k,i)*fs)-exband:round(GD_t(k,i)*fs)+exband) = 1;
        end
        Tx = Ts.*Ee;
%         figure;
%         imagesc(t,f,abs(Tx));axis xy;
        b = itsst(Tx,fs, xMean);
        figure
        plot(t,b,'LineWidth',2);
        hold on
        plot(t,ideal_y(k,:)-b,'color','black','LineWidth',2);
        if k == 1
            xlabel({'Time (s)','(c)'},'FontSize',20);set(gca,'XTick',0:2:10);
        else
            xlabel({'Time (s)','(d)'},'FontSize',20);set(gca,'XTick',0:2:10);
        end
        ylabel('AMP','FontSize',20); set(gca,'YTick',-0.4:0.2:0.4);
        legend('Reconstruction','Error');
        set(gca,'FontSize',20);axis([0 10 -0.4 0.4])

        SNRoutput2(k) = SNR(ideal_y(k,:),b);
        text(5.5,-0.3,['RQF:',num2str(SNRoutput2(k),3),'dB'],'Color','black','FontWeight','bold','FontSize',26)
    end
%% Close File
rmpath(genpath(folder1));
rmpath(genpath(folder2));
rmpath(genpath(folder3));
rmpath(genpath(folder4));
rmpath(genpath(folder5));
rmpath(genpath(folder6));