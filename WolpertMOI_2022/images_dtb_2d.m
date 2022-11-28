function P =  images_dtb_2d(R)
% P =  images_dtb_2d(R)
% Method of images solution to bounded drift diffusion of two races with
% correlated (-1/sqrt(0.5)) noise  i.e.  2D drift to bound.
% The winner of the race is the one that reaches its upper bound first
% There is  no lower bound and the bounds are flat
% the grid size (ngrid) for simulation is set at 500 in DV space but can be changed
%
% ~~~~~~~~~~~~~~~
% Inputs is a structure R (all in SI units)
% drift     vector of drift rates [length ndrift]
% t         vector time series to simulate [length nt]
% Bup       bound height [scalar]
% k_urg     limnear urgency sclaing
% grid       grid of x and y points
% lose_flag  whether we want the losing densities as well (slower)
% ~~~~~~~~~~~~~~~
% Outputs:  a structure P (which includes the inputs as well)
% P.up and P.lo (for two bounds) contains
%       pdf_t:  pdf of the winner [ndrift x nt]
%       cdf_t:  cdf of the winner [ndrift x nt]
%       p:  total prob up/lo  [ndrift]
%       mean_t:  mean time of up/lo [ndrift]
%       notabs  containe mean and variance for unabsorbed
%       distr_loser: pdf of the losing race [ndrift x nt x ngrid ]
% y:    dv values of grid (bounds are at zero and intial dv is -Bup) [ngrid]
% dy:   grid spacing
%
if ~isfield(R,'lose_flag')
    R.lose_flag=0;
end
P=R; % copy inputs to outputs

ndrift=length(R.drift); %number of drifts
nt=length(R.t); %number of time points to simulate
dt=R.t(2)-R.t(1); %sampling interval

ngrid=length(R.grid); %grid size can change
g=R.grid;  %dv values can change lower
dg=g(2)-g(1); %grid spacing

%  pdf of the losing races if needed
if R.lose_flag
    P.up.distr_loser = zeros(ndrift,nt,ngrid);
    P.lo.distr_loser = zeros(ndrift,nt,ngrid);
end
%  P.notabs = zeros(ndrift,nt,ngrid,ngrid);

for id=1:ndrift  %loop over drifts
    for k=2:nt  %loop over time starting at t(2)
        if R.lose_flag
            [P.up.pdf_t(id,k) P.lo.pdf_t(id,k) P.up.distr_loser(id,k,:) P.lo.distr_loser(id,k,:)]...
                = collapse_flux_7(R.Bup,R.k_urg,R.drift(id),R.t(k),g,g);
        else
            [P.up.pdf_t(id,k) P.lo.pdf_t(id,k) ]=collapse_flux_7(R.Bup,R.k_urg,R.drift(id),R.t(k),g,g);
        end

        [G1, G2]=meshgrid(g,g);
        notabs = collapse_notabs_7(R.Bup,R.k_urg,R.drift(id),R.t(k),G1,G2);

        %lower bound issues
        if 0
            low_indx=find(g>R.low_th,1);

            notabs(low_indx,:)=notabs(low_indx,:)+sum(notabs.*(G1<R.low_th),1);
            notabs(G1<R.low_th)=0;
            notabs(:,low_indx)=notabs(:,low_indx)+sum(notabs.*(G2<R.low_th),2);
            notabs(G2<R.low_th)=0;
        end
        %lower bound issue done

        [GT1, GT2]=meshgrid(g,g);

        GT1(GT1<R.low_th)=R.low_th;
        GT2(GT2<R.low_th)=R.low_th;

        sP=sum(notabs,'all');
        P.up.notabs.mean(id,k)=sum(notabs.*GT1,'all')/sP;
        m2=sum(notabs.*(GT1.^2),'all')/sP;
        P.up.notabs.var(id,k)=m2- P.up.notabs.mean(id,k)^2;

        P.lo.notabs.mean(id,k)=sum(notabs.*GT2,'all')/sP;
        m2=sum(notabs.*(GT2.^2),'all')/sP;
        P.lo.notabs.var(id,k)=m2- P.lo.notabs.mean(id,k)^2;
    end
end
P.up.pdf_t= P.up.pdf_t*dt;
P.lo.pdf_t= P.lo.pdf_t*dt;

if R.lose_flag %put any missing density into the most negative bin
    P.up.distr_loser = P.up.distr_loser*dg*dt;
    P.lo.distr_loser = P.lo.distr_loser*dg*dt;

    P.up.distr_loser(:,:,1)= P.up.pdf_t- nansum(P.up.distr_loser(:,:,2:end),3);
    P.lo.distr_loser(:,:,1)= P.lo.pdf_t- nansum(P.lo.distr_loser(:,:,2:end),3);
end

P.up.cdf_t = cumsum(P.up.pdf_t,2,'omitnan');
P.lo.cdf_t = cumsum(P.lo.pdf_t,2, 'omitnan');

P.up.p = P.up.cdf_t(:,end);
P.lo.p = P.lo.cdf_t(:,end);
P.up.pdf_t(isnan(P.up.pdf_t))=0;
P.lo.pdf_t(isnan(P.lo.pdf_t))=0;

P.up.mean_t = P.up.pdf_t * R.t  ./P.up.p;
P.lo.mean_t = P.lo.pdf_t * R.t  ./P.lo.p;

P.y=g;
P.dy=dg;





