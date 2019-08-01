function [Wx,t,f,xMean] = stft(x, fs, WindowOpt, Parameter, Mode)
%% 
    N = length(x);          
    t = (0:N-1)/fs;         
    [XPad, NUp, NL, ~] = padsignal(x, 'symmetric');    
    xMean = mean(XPad);     
    XPad = XPad-xMean;      
    XPad = hilbert(XPad);     
    xh = fft(XPad);          
%% 
    if nargin<5, Mode = 'modify'; end               

    if nargin<4, Parameter = struct(); end            
    if ~isfield(Parameter, 'L'), Parameter.L = round(N/2); end
    if ~isfield(Parameter, 'fmin'), Parameter.fmin = 0; end
    if ~isfield(Parameter, 'fmax'), Parameter.fmax = fs/2; end

    if nargin<3, WindowOpt = struct(); end            
    if ~isfield(WindowOpt, 's'), WindowOpt.s = 0.01; end
    if ~isfield(WindowOpt, 'type'), WindowOpt.type = 'gauss'; end

    s = WindowOpt.s; type = WindowOpt.type;
    L = Parameter.L; fmin = Parameter.fmin; fmax = Parameter.fmax;
%% 
    f = linspace(fmin, fmax, L);
    
    [gf,~] = windowf(s,type);
    
    wi = zeros(1, NUp);
    wi(1:NUp/2+1) = 2*pi*(0:NUp/2)*fs/NUp;
    wi(NUp/2+2:end) = 2*pi*(-NUp/2+1:-1)*fs/NUp;
    
    Wx = zeros(L, NUp);
    for ptr = 1:L
        gh = gf(wi-2*pi*f(ptr));
        gh = conj(gh);
        xcpsi = ifft(gh .* xh);
        Wx(ptr, :) = xcpsi;
    end
    
    Wx = Wx(:, NL+1:NL+length(x));
    
    if strcmp(Mode, 'normal')
        for i = 1:L
            for j = 1:N
                Wx(i,j) = Wx(i,j)*exp(-1i*2*pi*f(i)*t(j));
            end
        end
    end
end