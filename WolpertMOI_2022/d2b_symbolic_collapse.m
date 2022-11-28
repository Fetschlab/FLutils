%%
% this calculates the new location (x1,y1) and log-amplitude of  new sources to
% cancel out in x==0 and y==0 respectively given an initial location
% (x0,y0),  with process covariance  rho and a drift

clear all
syms x0 y0 rho drift real  % source location, covariance and drift
syms x1 y1  logA  real  % image location and log-amplitude
syms x y t  real  % generic location and time
syms k_urg real
assume(x0<0); % iniial x position as bound is at 0
assume(y0<0); % intial y position as bound is at 0
assume(rho<0); % covariance must be negative
assumeAlso(rho>-1)
assume(k_urg>=0); % urgency is positive

urgency= k_urg*t;

%mean for decision process starting at location (x0,y0) after time t
mx0 = x0 + drift*t + urgency;
my0 = y0 - drift*t + urgency;

%mean for an image process starting at location (x1,y1)
mx1 = x1 + drift*t + urgency;
my1 = y1 - drift*t + urgency;

%calculate log of pdf  of bivariate gaussian from each ignoring scaling
%constant  https://mathworld.wolfram.com/BivariateNormalDistribution.html
z0=(x-mx0)^2 - 2*rho*(x-mx0)*(y-my0) + (y-my0)^2;
z1=(x-mx1)^2 - 2*rho*(x-mx1)*(y-my1) + (y-my1)^2;

%the t is here as z has a denominator of t for variances so we multiply by t
p0 = log(1)*t +  (-z0/(2*(1-rho^2))); %unit amplitude source
p1 =   logA*t +  (-z1/(2*(1-rho^2))) ;

%
% %  logA*t+ (x1+f(t))^2=  (-B+f(t))^2
%
% logA*t+x1^2+2*x1*f(t)=B^2-2*f(t)*B
%
% x1^2=B^2

%the source and image must  cancel on the x=0 line
s1=subs(p0-p1,'x',0);

%s1 must be zero for all y and all t so find coeff to each solve to make zero
[d,b]=coeffs(s1,[y,t]);

% solve for the location and strength of the image (make sure not at same
% place as source also
Sy=solve(d(1)==0,d(2)==0,d(3)==0,x1~=x0,x1,y1,logA);
Sy.amp=-exp(Sy.logA);

%the source and image must  cancel on the y=0 line
s1=subs(p0-p1,'y',0);

%this must be zero for all x and all t so find coeff to each solve to make zero
[d,b]=coeffs(s1,[x,t]);
Sx=solve(d(1)==0,d(2)==0,d(3)==0,y1~=y0,x1,y1,logA);
Sx.amp=-exp(Sx.logA); %the actual amplitude has to be negative of the original source

Sx.x1
Sx.y1
Sx.amp


save reflections Sx Sy

%% this finds rho values which work for different number of reflections
clear all
syms B k_urg rho drift x0 y0  k_urg real
assume(rho<0); % covariance must be negative
assumeAlso(rho>-1)
assume(k_urg>=0); % urgency is positive

load reflections

N=6:2:8 ; % must be even - can extend at least to 20 - not all work properly for some reason
rho_used=[];

count=0;
for j=1:length(N)
    clear cent amp
    N(j)
    %initial starting location
    cent(1,:)=[-B -B];
    amp(1)=sym(1);
    
    for k=1:N(j) % do reflections alternating x and y
        %reflections
        x0=cent(k,1);
        y0=cent(k,2);
        
        if rem(k,2)==1
            cent(k+1,1)=subs(Sx.x1);
            cent(k+1,2)=subs(Sx.y1);
            amp(k+1)=subs(Sx.amp)*amp(k);
        else
            cent(k+1,1)=subs(Sy.x1);
            cent(k+1,2)=subs(Sy.y1);
            amp(k+1)=subs(Sy.amp)*amp(k);
        end
    end
    
    cent=simplify(cent);
    amp=simplify(amp)';
    
    d1=simplify(solve(cent(1,1)==cent(end,1),cent(1,2)==cent(end,2),amp(end)==1,rho));
    
    drho=unique(d1(double(d1)<0))
    u=vpa(drho)
    [C,IA] =setdiff(u,rho_used);
    rho_used=[rho_used ;C]
    for i=1:length(C)
        count=count+1;
        N_ref(count,1)=N(j);
        R(count,1)=u(IA(i));
        Rs(count,1)=drho(IA(i))
    end
end

[R,i]=sort(R)
N_ref=N_ref(i);
Rs=Rs(i)

save rho_n  R N  Rs



%% calculates centres and weights for a particular Nref reflections

clear all
load reflections
load rho_n
syms B rho drift x0 y0 real
syms xc yc mag  real
syms k_urg t x  y  drift  rho real

assume(t>0)
assume(rho<0)
assumeAlso(rho>-1)
assume(k_urg>=0)

mn=[xc+drift*t yc-drift*t]' + k_urg*t;
cov=[1 rho; rho 1]*t;

X=[x y]';
z=-((X-mn)'*inv(cov)*(X-mn))/2;
p=  mag*exp(z)/(2*pi*sqrt((det(cov))))  ;

px=diff(p,x);
pxflux=-0.5*subs(px,x,0);
pxint=int(pxflux,y);
pxpdf=subs(pxint,y,0)-subs(pxint,y,-1000);

py=diff(p,y);
pyflux=-0.5*subs(py,y,0);
pyint=int(pyflux,x);
pypdf=subs(pyint,x,0)-subs(pyint,x,-100);

%vpa(subs(pxpdf,{y,mx,my,k_urg,mag,x0,y0,rho,t,drift},{0    ,3, 2,3,10,1,2,-0.7,1,1}))

j=1; %chose tje jth rho value
Nref=N(j);

%initial starting location  and amplitude
cent(1,:)=[-B -B];
amp(1)=sym(1); %source amplitude

for k=1:Nref % do Nref reflections alternating x and y
    %reflections
    x0=cent(k,1);
    y0=cent(k,2);
    
    if rem(k,2)==1
        cent(k+1,1)=subs(Sx.x1);
        cent(k+1,2)=subs(Sx.y1);
        amp(k+1)=subs(Sx.amp)*amp(k);
    else
        cent(k+1,1)=subs(Sy.x1);
        cent(k+1,2)=subs(Sy.y1);
        amp(k+1)=subs(Sy.amp)*amp(k);
    end
end
cent(end,:)=[];
amp(end)=[];

cent=simplify(cent);
amp=simplify(amp)';

P=0;
Pxflux=0;
Pyflux=0;
Pxpdf=0;
Pypdf=0;

for k=1:size(cent,1)
    P=P+subs(p,               {mag,xc,yc,},{amp(k),cent(k,1),cent(k,2)});
    Pxflux=Pxflux+subs(pxflux,{mag,xc,yc},{amp(k),cent(k,1),cent(k,2)});
    Pxpdf=Pxpdf+subs(pxpdf,   {mag,xc,yc},{amp(k),cent(k,1),cent(k,2)});
    Pyflux=Pyflux+subs(pyflux,{mag,xc,yc},{amp(k),cent(k,1),cent(k,2)});
    Pypdf=Pypdf+subs(pypdf,   {mag,xc,yc},{amp(k),cent(k,1),cent(k,2)});
end

RHO=Rs(j);

P=subs(P,rho,RHO);
Pxflux=subs(Pxflux,rho,RHO);
Pyflux=subs(Pyflux,rho,RHO);
Pxpdf=subs(Pxpdf,rho,RHO);
Pypdf=subs(Pypdf,rho,RHO);


str=sprintf('collapse_flux_%i',7);

matlabFunction(Pxpdf,Pypdf,Pxflux,Pyflux,'vars',[B k_urg drift t x y],'file',str)


str=sprintf('collapse_notabs_%i',7);

matlabFunction(P,'vars',[B k_urg drift t x y],'file',str)

