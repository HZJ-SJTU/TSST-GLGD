function [XPad, NUp, NL, NR] = padsignal(X, PadType)
    N = length(X);
    if iscolumn(X)
        X = X';
    end
    [NUp, NL, NR] = p2up(N);
    if strcmpi(PadType,'symmetric')
        Temp=[X flip(X)];               
        XL = Temp(mod((0:NL-1),2*N)+1);
        XL =flip(XL);
        Temp=[flip(X) X];
        XR = Temp(mod((0:NR-1),2*N)+1);
    elseif strcmpi(PadType,'replicate')
        XL = X(1)*ones(NL,1);
        XR = X(end)*ones(NR,1);
    end
    XPad = [XL X XR];
end