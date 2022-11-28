% bcinit.m : initialize b_change and B_S for FP4 and spectral_dtb, respectively  -CF
function [b_change, B_S] = bcinit(B,t,dt,type,uInf,tau05)

% bound change over time (or not)

switch type
    case 'linear' % linear collapsing
        b_change = zeros(length(t),2);
        bslope = 0.01;
        b_change(:,1) = bslope*t;
        b_change(:,2) = -bslope*t;
        B_S = (B+b_change(:,2))./sqrt(1000);
        B_S(end+1) = (B-bslope*(t(end)+dt))/sqrt(1000);
    case 'hyperbolic' % hyperbolic collapsing (Churchland et al. 2008, Hanks et al. 2011)
        b_change = zeros(length(t),2);
%         uInf = 30;
%         tau05 = 350;
        b_change(:,2) = -uInf*t./(t+tau05);
        b_change(:,1) = -b_change(:,2);
        B_S = (B+b_change(:,2))./sqrt(1000);
        B_S(end+1) = (B-uInf*(t(end)+dt)/((t(end)+dt)+tau05)) / sqrt(1000);
    case 'flat'
        b_change = zeros(length(t),2);
        B_S = zeros(length(t),1);
        B_S(:) = B/sqrt(1000);
end
