    clc; clear all; close all;
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
%% signal
    load('batdata.mat');
    x = data';
    N = length(data);           
    fs = 1000000/7;          
%% Parameters
    map = jet;
    exband = 20;
    WindowOpt = struct('type','gauss','s',0.00015);
    Parameter = struct('L',N/2+1,'fmin',0,'fmax',fs/2);
%% figure 7
    figure('color',[1 1 1]);
    t = (0:N-1)/fs;
    plot(t*1000,x);
    axis tight; xlabel({'Time (ms)','(a)'},'FontSize',20);set(gca,'XTick',0:1:3);
    ylabel('Amp','FontSize',20);set(gca,'YTick',-0.25:0.1:0.25);
    set(gca,'FontSize',20);axis([0 2.8 -0.25 0.15])
    
    figure('color',[1 1 1]);
    ffty = abs(fft(x))/N*2;
    f = (0:N/2)/N*fs;
    plot(f/1000,ffty(1:N/2+1))
    axis tight; xlabel({'Frequency (kHz)','(b)'},'FontSize',20);set(gca,'XTick',0:10:70);
    ylabel('Amp','FontSize',20);set(gca,'FontSize',20);set(gca,'YTick',0:0.005:0.02);
    set(gca,'yTickLabel',num2str(get(gca,'yTick')','%.3f'));axis([0 70 0 0.02])
%% S-method
    t =1:1:400;
    f_num = 256;
	G=tftb_window(9,'hanning'); 
    h=tftb_window(61,'hanning'); 
    [~,rtfr,~] = tfrrstan(x',t,f_num,G,h,1);
    t = (t-1)/fs;
    f = (0:f_num/2-1)/f_num*fs;
    rtfr = rtfr(1:f_num/2,:);
    figure;
    imagesc(t*1000,f/1000,rtfr);axis xy;colormap(map);axis tight; 
    xlabel({'Time (ms)','(a)'},'FontSize',20);set(gca,'XTick',0:0.5:3);
    ylabel('Frequency (kHz)','FontSize',20); set(gca,'YTick',0:20:70);
    set(gca,'FontSize',20);axis([0 2.79 0 70])
    r1 = abs(renyi(rtfr,t,f',3));
    text(0,5,['RE:',num2str(abs(r1),3),'(bits)'],'Color','black','FontWeight','bold','FontSize',26) 
%% the discrete-time Levenberg-Marquard synchrosqueezed Gabor Transform    
    f_num = 256;
    t = (0:N-1)/fs;
    [~, stfr, ~] = tfrvsgab(x, f_num, 10);
    stfr = stfr(1:f_num/2,:);
    f = (0:f_num/2-1)/f_num*fs;
    figure;
    imagesc(t*1000,f/1000,abs(stfr));axis xy;colormap(map);axis tight; 
    xlabel({'Time (ms)','(b)'},'FontSize',20);set(gca,'XTick',0:0.5:3);
    ylabel('Frequency (kHz)','FontSize',20); set(gca,'YTick',0:20:70);
    set(gca,'FontSize',20);axis([0 2.79 0 70])
    r2 = abs(renyi(abs(stfr),t,f',3));
    text(0,5,['RE:',num2str(abs(r2),3),'(bits)'],'Color','black','FontWeight','bold','FontSize',26) 
%%  TSST
    [Ts1,t,f,xMean,~] = tsst(x , fs,  WindowOpt, Parameter, '1Ord');
    figure;
    imagesc(t*1000,f/1000,abs(Ts1));axis xy;colormap(map);axis tight; 
    xlabel({'Time (ms)','(c)'},'FontSize',20);set(gca,'XTick',0:0.5:3);
    ylabel('Frequency (kHz)','FontSize',20); set(gca,'YTick',0:20:70);
    set(gca,'FontSize',20);axis([0 2.79 0 70])
    r3 = abs(renyi(abs(Ts1),t,f',3));
    text(0,5,['RE:',num2str(abs(r3),3),'(bits)'],'Color','black','FontWeight','bold','FontSize',26) 
%% TSST2
    [Ts2,t,f,xMean,~] = tsst(x , fs,  WindowOpt, Parameter, '2Ord');
    figure;
    imagesc(t*1000,f/1000,abs(Ts2));axis xy;colormap(map);axis tight; 
    xlabel({'Time (ms)','(d)'},'FontSize',20);set(gca,'XTick',0:0.5:3);
    ylabel('Frequency (kHz)','FontSize',20); set(gca,'YTick',0:20:70);
    set(gca,'FontSize',20);axis([0 2.79 0 70])
    r4 = abs(renyi(abs(Ts2),t,f',3));
    text(0,5,['RE:',num2str(abs(r4),3),'(bits)'],'Color','black','FontWeight','bold','FontSize',26) 
%% ITSST2
    [Cs, Es] = brevridge_mult(abs(Ts2'), 0:N-1, 4, 0.009, 10);
%     figure
%     plot(Cs/fs*1000,f/1000)
%     axis([0 2.79 0 70])
    Tx = zeros(Parameter.L,N);
    for k = 1:4
        for i =1:Parameter.L
            min = round(Cs(k,i))-exband;
            min(min<=1) = 1;
            max = round(Cs(k,i))+exband;
            max(max>=N) = N;
            Tx(i,min:max) = 1;
        end
        Tx = Ts2.*Tx;
%         figure;
%         imagesc(t,f,abs(Tx));axis xy;
        b(k,:) = itsst(Tx,fs, 0);
        figure('color',[1 1 1]);
        plot(t*1000,x);hold on
        plot(t*1000,b(k,:));
        if k == 1
            axis tight; xlabel({'Time (ms)','(b)'},'FontSize',20);set(gca,'XTick',0:1:3);
            legend('Original','Mode2','Location','SouthEast');
        elseif k == 2
            axis tight; xlabel({'Time (ms)','(a)'},'FontSize',20);set(gca,'XTick',0:1:3);
            legend('Original','Mode1','Location','SouthEast');
        elseif k == 3
            axis tight; xlabel({'Time (ms)','(c)'},'FontSize',20);set(gca,'XTick',0:1:3);
            legend('Original','Mode3','Location','SouthEast');
        else
            axis tight; xlabel({'Time (ms)','(d)'},'FontSize',20);set(gca,'XTick',0:1:3);
            legend('Original','Mode4','Location','SouthEast');
        end
        ylabel('Amp','FontSize',20);set(gca,'YTick',-0.25:0.1:0.25);
        set(gca,'FontSize',20);axis([0 2.8 -0.25 0.15])
    end
    sum_b = b(1,:)+b(2,:)+b(3,:)+b(4,:);
    figure('color',[1 1 1]);
    plot(t*1000,x);hold on
    plot(t*1000,sum_b);
    axis tight; xlabel({'Time (ms)','(e)'},'FontSize',20);set(gca,'XTick',0:1:3);
    ylabel('Amp','FontSize',20);set(gca,'YTick',-0.25:0.1:0.25);
    set(gca,'FontSize',20);axis([0 2.8 -0.25 0.15])
    legend('Original','Reconstruction','Location','SouthEast');
    
    figure('color',[1 1 1]);
    plot(t*1000,x);hold on
    plot(t*1000,x-sum_b,'Color','black','LineWidth',1);
    axis tight; xlabel({'Time (ms)','(f)'},'FontSize',20);set(gca,'XTick',0:1:3);
    ylabel('Amp','FontSize',20);set(gca,'YTick',-0.25:0.1:0.25);
    set(gca,'FontSize',20);axis([0 2.8 -0.25 0.15])
    legend('Original','Error','Location','SouthEast');
%% Close File
    rmpath(genpath(folder1));
    rmpath(genpath(folder2));
    rmpath(genpath(folder3));
    rmpath(genpath(folder4));
    rmpath(genpath(folder5));
    rmpath(genpath(folder6));