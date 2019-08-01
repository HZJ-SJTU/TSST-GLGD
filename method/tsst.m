function [Tx,t,f,xMean,GroupDelay] = tsst(x , fs,  WindowOpt, Parameter, Mode)
%% 
    N = length(x);
%% 
    s = WindowOpt.s; type = WindowOpt.type;
    L = Parameter.L; fmin = Parameter.fmin; fmax = Parameter.fmax;
    gamma = sqrt(eps); 
%% 
 
    [Wx,t,f,xMean] = stft(x, fs, WindowOpt, Parameter, 'normal');

    if strcmp(Mode, '1Ord')
        WindowOpt.type = '1ord(w)_gauss';
        [dWx,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
        GroupDelay = imag(dWx./Wx);
        for ptr = 1:N
            GroupDelay(:,ptr) = GroupDelay(:,ptr) + t(ptr);
        end
        GroupDelay( abs(Wx) < gamma ) = Inf;
    elseif strcmp(Mode, '2Ord')
        WindowOpt.type = '1ord(w)_gauss';
        [dWx,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
        WindowOpt.type = '2ord(w)_gauss';
        [ddWx,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
        WindowOpt.type = 'w*gauss';
        [wWx,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
        WindowOpt.type = 'w*1ord(w)_gauss';
        [wdWx,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
        Denominator = wdWx.*Wx-dWx.*wWx;
        Numerator = ddWx.*wWx-dWx.*wdWx;
        p = Numerator./Denominator;
        for ptr = 1:N
            p(:,ptr) = p(:,ptr) - 1i*t(ptr);
        end
        GroupDelay = -imag(p);
        GroupDelay( abs(Denominator) < gamma ) = Inf;
    else
        error('Unknown SST Mode: %s', Mode);
    end

    dt = 1/fs;

    [gf,~] = windowf(s,type);

    g0 = gf(0);
    g0 = conj(g0);
    if(g0 == 0)
        error('window must be non-zero and continuous at 0 !');
    end

    Wx(isinf(GroupDelay)) = 0;
    Tx = zeros(L,N);
    for prt=1:L
        for b=1:N
            m = min(max(1 + round((GroupDelay(prt,b)-0)/dt),1),N);
            Tx(prt, m) = Tx(prt, m) + Wx(prt, b)*dt;
        end
    end
    Tx = Tx / g0;
end
