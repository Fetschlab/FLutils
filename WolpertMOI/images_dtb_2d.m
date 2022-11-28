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
% drift:      vector of drift rates [length ndrift]
% t:          vector time series to simulate [length nt]
% Bup:        bound height [scalar]
% lose_flag:  whether we want the losing densities as well (slower)
%
% drift_freq: frequencies of each drift rate (stimulus strength) in the data, for calculation of marginals - CF
% plotflag:   make some plots (or not) [1 = plot, 2 = plot and save] - CF

% ~~~~~~~~~~~~~~~
% Outputs:  a structure P (which includes the inputs as well)
% P.up and P.lo (for two bounds) contains
%       pdf_t:  pdf of the winner [ndrift x nt]
%       cdf_t:  cdf of the winner [ndrift x nt]
%       p:  total prob up/lo  [ndrift]
%       mean_t:  mean time of up/lo [ndrift]
%       distr_loser: pdf of the losing race [ndrift x nt x ngrid ]
% y:    dv values of grid (bounds are at zero and intial dv is -Bup) [ngrid]
% dy:   grid spacing
% 
% logOddsCorrMap: log posterior odds of a correct response,
% as a function of state of losing accumulator (Kiani 2014) -- CF


if ~isfield(R,'lose_flag')
    R.lose_flag=0;
end
if ~isfield(R,'plotflag')
    R.plotflag=0;
end
P=R; % copy inputs to outputs


ndrift=length(R.drift); %number of drifts
nt=length(R.t); %number of time points to simulate
ngrid=500; %grid size can change

dt=R.t(2)-R.t(1); %sampling interval

g=linspace(-4*R.Bup,0,ngrid); %plausible range of dv values

dg=g(2)-g(1); %grid spacing

%  pdf of the losing races if needed
if R.lose_flag
    P.up.distr_loser = zeros(ndrift,nt,ngrid);
    P.lo.distr_loser = zeros(ndrift,nt,ngrid);
end

for id=1:ndrift %loop over drifts
    for k=2:nt  %loop over time starting at t(2)
        if R.lose_flag
            [P.up.pdf_t(id,k) P.lo.pdf_t(id,k) P.up.distr_loser(id,k,:) P.lo.distr_loser(id,k,:)]...
                = flux_img7(R.Bup,R.drift(id),R.t(k),g,g);
        else
            [P.up.pdf_t(id,k) P.lo.pdf_t(id,k) ]=flux_img7(R.Bup,R.drift(id),R.t(k),g,g); %#ok<*NCOMMA>
        end
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

P.up.cdf_t = cumsum(P.up.pdf_t,2);
P.lo.cdf_t = cumsum(P.lo.pdf_t,2);

P.up.p = P.up.cdf_t(:,end);
P.lo.p = P.lo.cdf_t(:,end);

try
    P.up.mean_t = P.up.pdf_t * R.t  ./P.up.p;
    P.lo.mean_t = P.lo.pdf_t * R.t  ./P.lo.p;
catch
    P.up.mean_t = P.up.pdf_t * R.t'  ./P.up.p;
    P.lo.mean_t = P.lo.pdf_t * R.t'  ./P.lo.p;
end

P.y=g;
P.dy=dg;

% CF: log posterior odds of a correct response (Eq. 3 in Kiani et al. 2014)
I = R.drift>=0; % In case drift is signed, calculate only for positives,
                % then the kluge with separate marginals (Pxt's) can be done elsewhere
% unlike 1D code, here we'll marginalize over drift in one step               
odds = (squeeze(sum(P.up.distr_loser(I,:,:),1)) / length(R.drift(I))) ./ ...
       (squeeze(sum(P.lo.distr_loser(I,:,:),1)) / length(R.drift(I)));
odds(odds<1) = 1; % fix some stray negatives/zeros (what about infs?)
odds(isinf(odds)) = nan;
P.logOddsCorrMap = log(odds);
P.logOddsCorrMap = P.logOddsCorrMap';

if R.plotflag

    % (1) plot choice and RT
    figure(111); set(gcf,'Color',[1 1 1],'Position',[600 600 450 700],'PaperPositionMode','auto'); clf;
    subplot(3,1,1); plot(R.drift,P.up.p,'o-'); title('Prob correct bound crossed before tmax') % meaningless w signed drift
    subplot(3,1,2); plot(R.drift,P.up.p./(P.up.p+P.lo.p),'o-'); title('Relative prob of corr vs. incorr bound crossed (Pcorr)'); % actually pRight with signed drift
    subplot(3,1,3); plot(R.drift,P.up.mean_t,'o-'); title('mean RT'); xlabel('drift rate');

    % remake Fig. 5c-d of Kiani et al 2014 (steps analogous to makeLogOddsCorrMap_*)
    n = 100; % set n to 100+ for smooth plots, lower for faster plotting
    
    % (2) first an example PDF
    c = round(length(R.drift(I))/2) - 1 + sum(I==0); % pick an intermediate drift rate, or make a loop to see all of them
    q = 30; % exponent for log cutoff (redefine zero as 10^-q, for better plots)
    Pmap = squeeze(P.up.distr_loser(c,:,:))';
    Pmap(Pmap<10^-q) = 10^-q;
    logPmap = (log10(Pmap)+q)/q;

    figure(c*100); set(gcf, 'Color', [1 1 1], 'Position', [200 300 500 850/R.plotflag], 'PaperPositionMode', 'auto');
    if R.plotflag==1; subplot(2,1,1); end
    [~,h] = contourf(R.t,P.y,logPmap,n); colormap(jet);
    caxis([0 1]); % because we log transformed such that 10^-q is 0 and 1 is 1
    colorbar('YTick',0:.2:1,'YTickLabel',{['10^-^{' num2str(q) '}']; ['10^-^{' num2str(q*.8) '}']; ['10^-^{' num2str(q*.6) '}']; ['10^-^{' num2str(q*.4) '}']; ['10^-^{' num2str(q*.2) '}']; '1'}); 
    set(gca,'XLim',[0 R.t(end)],'XTick',0:0.5:floor(R.t(end)*2)/2,'TickDir','out');
    set(gca,'YLim',[-4*R.Bup 0],'YTick',floor(-4*R.Bup):round(R.Bup*2)/2:0);
    set(h,'LineColor','none');
    xlabel('Time (s)'); ylabel('Accumulated evidence of losing accumulator');
    title('Probability density of losing accumulator');
    changeAxesFontSize(gca,15,15);
    if R.plotflag==2
        ylim([-2.5 0]);
        xlabel([]);ylabel([]);title([]);
        changeAxesFontSize(gca,24,24);
        export_fig('2dacc_PDF','-eps');
    end
        
    % (3) then the log odds corr map
    if R.plotflag==1
        subplot(2,1,2);
    else
        figure(c*1000); set(gcf, 'Color', [1 1 1], 'Position', [700 400 450 700/R.plotflag], 'PaperPositionMode', 'auto');
    end
    [~,h] = contourf(R.t,P.y,P.logOddsCorrMap,n); % colormap(parula);
    caxis([0 3]);
    colorbar('YTick',0:0.5:3); 
    set(gca,'XLim',[0 R.t(end)],'XTick',0:0.5:floor(R.t(end)*2)/2,'TickDir','out');
    set(gca,'YLim',[-4*R.Bup 0],'YTick',floor(-4*R.Bup):round(R.Bup*2)/2:0); 
    set(h,'LineColor','none');
    xlabel('Time (s)'); ylabel('Accumulated evidence of losing accumulator');
    title('Log odds correct vs. state of losing accumulator');
    changeAxesFontSize(gca,15,15);
    
    if R.plotflag==2
        ylim([-2.5 0])
        xlabel([]);ylabel([]);title([]);
        changeAxesFontSize(gca,24,24);
        export_fig('2dacc_logoddscorr','-eps')
    end

end



