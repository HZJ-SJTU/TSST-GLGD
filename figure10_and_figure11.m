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
    load('vib_data')
    fs = 12000; N = 1200; 
    x=data(N+1:N*2);
%% Parameters
    map = jet;
    x1 = 0.032;x2 = 0.042;
    y1 = 2.55;y2 = 3.9;

    WindowOpt = struct('type','gauss','s',0.0015);

    Parameter = struct('L',N/2+1,'fmin',0,'fmax',fs/2);
%% figure 11
    figure('color',[1 1 1]);
    t = (0:N-1)/fs;
    plot(t,x);
    axis tight; xlabel({'Time (s)','(a)'},'FontSize',20);set(gca,'XTick',0:0.02:0.1);
    ylabel('Amp (g)','FontSize',20);set(gca,'YTick',-3:1:3);
    set(gca,'FontSize',20);axis([0 0.1 -3 3])
    
    figure('color',[1 1 1]);
    ffty = abs(fft(x))/N*2;
    f = (0:N/2)/N*fs;
    plot(f/1000,ffty(1:N/2+1))
    axis tight; xlabel({'Frequency (kHz)','(b)'},'FontSize',20);set(gca,'XTick',0:1:6);
    ylabel('Amp (g)','FontSize',20);set(gca,'FontSize',20);set(gca,'YTick',0:0.05:0.3);
    axis([0 6 0 0.3])
%% S-method
    t =1:5:N;
    f_num = 512;
	G=tftb_window(9,'hanning'); 
    h=tftb_window(61,'hanning'); 
    [~,rtfr,~] = tfrrstan(x,t,f_num,G,h,1);
    t = (t-1)/fs;
    f = (0:f_num/2-1)/f_num*fs;
    rtfr = rtfr(1:f_num/2,:);
    
    figure;
    imagesc(t,f/1000,rtfr);axis xy;colormap(map);axis tight; 
    xlabel({'Time (s)'},'FontSize',26);set(gca,'XTick',0:0.01:0.1);
    ylabel('Frequency (kHz)','FontSize',26); set(gca,'YTick',2.5:0.5:4); 
    set(gca,'FontSize',26);axis([0 0.1 2.5 4])
    r1 = abs(renyi(rtfr,t,f',3));
    text(0,3.8,['RE:',num2str(abs(r1),3),'(bits)'],'Color','black','FontWeight','bold','FontSize',40) 
    set(gcf,'Units','centimeter','Position',[5 5 40 12]);
    rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','r','Linewidth',3);
    
    figure('color',[1 1 1]);
    colormap(map);
    imagesc(t,f/1000,rtfr);axis xy;
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
    set(gcf,'Units','centimeter','Position',[5 5 5 12]);
    xlim([x1 x2]);ylim([y1 y2]);  
%% the discrete-time Levenberg-Marquard synchrosqueezed Gabor Transform    
    f_num = 512;
    t = (0:N-1)/fs;
    [~, stfr, ~] = tfrvsgab(x, f_num, 10);
    stfr = stfr(1:f_num/2,:);
    f = (0:f_num/2-1)/f_num*fs;
    
    figure;
    imagesc(t,f/1000,abs(stfr));axis xy;colormap(map);axis tight; 
    xlabel({'Time (s)'},'FontSize',26);set(gca,'XTick',0:0.01:0.1);
    ylabel('Frequency (kHz)','FontSize',26); set(gca,'YTick',2.5:0.5:4); 
    set(gca,'FontSize',26);axis([0 0.1 2.5 4])
    r2 = abs(renyi(abs(stfr),t,f',3));
    text(0,3.8,['RE:',num2str(abs(r2),3),'(bits)'],'Color','black','FontWeight','bold','FontSize',40) 
    set(gcf,'Units','centimeter','Position',[5 5 40 12]);
    rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','r','Linewidth',3);
    
    figure('color',[1 1 1]);
    colormap(map);
    imagesc(t,f/1000,abs(stfr));axis xy;
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
    set(gcf,'Units','centimeter','Position',[5 5 5 12]);
    xlim([x1 x2]);ylim([y1 y2]);  
%%  TSST
    [Ts1,t,f,xMean,~] = tsst(x , fs,  WindowOpt, Parameter, '1Ord');
    
    figure;
    imagesc(t,f/1000,abs(Ts1));axis xy;colormap(map);axis tight; 
    xlabel({'Time (s)'},'FontSize',26);set(gca,'XTick',0:0.01:0.1);
    ylabel('Frequency (kHz)','FontSize',26); set(gca,'YTick',2.5:0.5:4); 
    set(gca,'FontSize',26);axis([0 0.1 2.5 4])
    r3 = abs(renyi(abs(Ts1),t,f',3));
    text(0,3.8,['RE:',num2str(abs(r3),3),'(bits)'],'Color','black','FontWeight','bold','FontSize',40) 
    set(gcf,'Units','centimeter','Position',[5 5 40 12]);
    rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','r','Linewidth',3);
    
    figure('color',[1 1 1]);
    colormap(map);
    imagesc(t,f/1000,abs(Ts1));axis xy;
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
    set(gcf,'Units','centimeter','Position',[5 5 5 12]);
    xlim([x1 x2]);ylim([y1 y2]); 
%% TSST2
    [Ts2,t,f,xMean,~] = tsst(x , fs,  WindowOpt, Parameter, '2Ord');
    
    figure;
    imagesc(t,f/1000,abs(Ts2));axis xy;colormap(map);axis tight; 
    xlabel({'Time (s)'},'FontSize',26);set(gca,'XTick',0:0.01:0.1);
    ylabel('Frequency (kHz)','FontSize',26); set(gca,'YTick',2.5:0.5:4); 
    set(gca,'FontSize',26);axis([0 0.1 2.5 4])
    r4 = abs(renyi(abs(Ts2),t,f',3));
    text(0,3.8,['RE:',num2str(abs(r4),3),'(bits)'],'Color','black','FontWeight','bold','FontSize',40) 
    set(gcf,'Units','centimeter','Position',[5 5 40 12]);
    rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','r','Linewidth',3);
    
    figure('color',[1 1 1]);
    colormap(map);
    imagesc(t,f/1000,abs(Ts2));axis xy;
    set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
    set(gcf,'Units','centimeter','Position',[5 5 5 12]);
    xlim([x1 x2]);ylim([y1 y2]); 
%% Close File
    rmpath(genpath(folder1));
    rmpath(genpath(folder2));
    rmpath(genpath(folder3));
    rmpath(genpath(folder4));
    rmpath(genpath(folder5));
    rmpath(genpath(folder6));