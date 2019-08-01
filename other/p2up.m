function [NUp,NL,NR] = p2up(N)

    NUp = 2^(ceil(log2(N+eps))+1);
    NL = floor((NUp-N)/2);
    NR = NL;
    if (mod(2*NL+N,2)==1)
        NR = NL + 1; 
    end
end