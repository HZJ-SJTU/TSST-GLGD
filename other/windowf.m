function [gf,gt] = windowf(s,type)

    switch type
        case 'gauss'
            gt = @(t) s^(-1/2)*pi^(-1/4).*exp(-t.^2/s^2/2);
            gf = @(w) sqrt(2 * s)*pi^(1/4)*exp(-(s * w).^2/2);
        case '1ord(t)_gauss'
            gt = @(t) s^(-5/2)*pi^(-1/4).*(-t).*exp(-t.^2/s^2/2);
            gf = @(w) 1i.*sqrt(2 * s)*pi^(1/4)*exp(-(s * w).^2/2).*w;
        case '2ord(t)_gauss'
            gt = @(t) s^(-9/2)*pi^(-1/4).*(t^2-s^2).*exp(-t.^2/s^2/2);
            gf = @(w) -sqrt(2 * s)*pi^(1/4)*exp(-(s * w).^2/2).*w.^2;
        case 't*gauss'
            gt = @(t) s^(-1/2)*pi^(-1/4).*(t).*exp(-t.^2/s^2/2);
            gf = @(w) -1i.*sqrt(2)*s^(5/2)*pi^(1/4)*exp(-(s * w).^2/2).*w;
        case 't*1ord(t)_gauss'
            gt = @(t) s^(-5/2)*pi^(-1/4).*(-t.^2).*exp(-t.^2/s^2/2);
            gf = @(w) sqrt(2 * s)*pi^(1/4)*exp(-(s * w).^2/2).*((s * w).^2-1);
        case '1ord(w)_gauss'
            gt = @(t) 1i*sqrt(2)*s^(-1/2)*pi^(1/4).*(t).*exp(-t.^2/s^2/2);
            gf = @(w) -sqrt(2)*s^(5/2)*pi^(1/4)*exp(-(s * w).^2/2).*w;
        case '2ord(w)_gauss'
            gt = @(t) s^(-9/2)*pi^(-1/4).*(t^2-s^2).*exp(-t.^2/s^2/2);  
            gf = @(w) -sqrt(2)*s^(5/2)*pi^(1/4)*exp(-(s * w).^2/2).*(1-(s * w).^2);
        case 'w*gauss'
            gt = @(t) s^(-1/2)*pi^(-1/4).*(t).*exp(-t.^2/s^2/2);       
            gf = @(w) sqrt(2 * s)*pi^(1/4)*exp(-(s * w).^2/2).*w;
        case 'w*1ord(w)_gauss'
            gt = @(t) s^(-5/2)*pi^(-1/4).*(-t.^2).*exp(-t.^2/s^2/2);   
            gf = @(w) -sqrt(2)*s^(5/2)*pi^(1/4)*exp(-(s * w).^2/2).*w.^2;
        otherwise
            error('Unknown window type: %s', type);
    end 
end