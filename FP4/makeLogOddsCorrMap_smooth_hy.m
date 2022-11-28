function [logOddsMapR, logOddsMapL, logOddsCorrMap, PMap, tAxis, vAxis] = makeLogOddsCorrMap_smooth(k,B,sigma,theta,tAxis,plotflag)

% % to run as script, uncomment:
% % clear;
% % close all;
% % plotflag=2;
% % k=0.4; % sensitivity parameter (mean drift rate = k*Coherence)
% % B=40; % bound height
% % sigma=1; % standard deviation of momentary evidence
% % theta=0.8; % threshold, in log odds correct, for making a direction choice as opposed to sure-bet

% keyboard

cohs = [-0.512 -0.256 -0.128 -0.064 -0.032 -eps eps 0.032 0.064 0.128 0.256 0.512];
uvect = k*cohs; % vector of mean drift rates
buf = 4*sigma; % buffer for bound crossings
if isempty(tAxis)
    dT = 1; % time step (ms)
    tAxis = (0:dT:800)'; % time axis
else
    dT = mode(diff(tAxis));
end

y0 = 0;
Bup = ones(1,length(tAxis))*(B+buf); % upper bound
Blo = ones(1,length(tAxis))*-(B+buf); % lower bound

% use Fokker-Planck equation to propagate the probability density of DV;
% xtDist is the probability as a function of x (DV) and t (time)
[~,~,~,~,~,~,~,~,xtDist] = runFPonSeries(uvect,tAxis,Bup,Blo,y0,sigma);


% remake Fig. 4a-b of Kiani and Shadlen 2009
vAxis = linspace(-(B+buf), B+buf, size(xtDist{1},1));
mapLeft = 0; mapRight = 0;
% for each coherence, the "map" is the probability density conditioned on
% the stimulus being either "right" or "left" (Fig. 4A)

% this marginalizes over motion strength for a given direction (see text below Equation in paper)
for c = 1:length(cohs)
    if cohs(c)>0
        mapRight = mapRight + xtDist{c};
    else  % careful if using coh=0!
        mapLeft = mapLeft + xtDist{c};
    end
end
mapRight = mapRight / sum(cohs>0); 
mapLeft = mapLeft / sum(cohs<0);

% clip PDFs at some minimum, for plotting purposes
mapRightClipped = mapRight; mapRightClipped(mapRightClipped<1e-6) = 1e-6;
% mapLeftClipped = mapLeft; mapLeftClipped(mapLeftClipped<1e-9) = 1e-9;

n = 60; % set n to 100+ for Roozbeh's smooth plots, lower for faster plotting
if plotflag==2
    % log-transform P to enable log scale for contourf
    logPmap = (log10(mapRightClipped)+6)/6;
    figure; set(gcf, 'Color', [1 1 1], 'Position', [100 100 450 350], 'PaperPositionMode', 'auto'); hold on;
    [~,h] = contourf(tAxis,vAxis,logPmap,n); caxis([0 1]); % log transformed such that 1e-6 is 0 and 1 is 1.
    colorbar('YTick',[0 1/6 2/6 3/6 4/6 5/6 1],'YTickLabel',{'1e-6';'1e-5';'1e-4';'1e-3';'1e-2';'1e-1';'1'}); 
    set(h,'LineColor','none'); xlabel('time (ms)'); ylabel('DV'); title('Prob of DV for rightward motion');
    colormap(jet);

%     logPmap = (log10(mapLeftClipped)+6)/6;
%     figure; set(gcf, 'Color', [1 1 1], 'Position', [100 100 450 350], 'PaperPositionMode', 'auto'); hold on;
%     [~,h] = contourf(tAxis,dvAxis,logPmap,n); caxis([0 1]); % log transformed such that 1e-6 is 0 and 1 is 1.
%     colorbar('YTick',[0 1/6 2/6 3/6 4/6 5/6 1],'YTickLabel',{'1e-6';'1e-5';'1e-4';'1e-3';'1e-2';'1e-1';'1'}); 
%     set(h,'LineColor','none'); xlabel('time (ms)'); ylabel('DV');  title('PDF of DV for leftward motion');
end


% log odds is log(P/(1-P)), which in this case is Right/Left or vice versa
logOddsMapR = log(mapRight./mapLeft); 
logOddsMapL = log(mapLeft./mapRight);
if mod(size(logOddsMapL,1),2)==0
    logOddsCorrMap = [logOddsMapL(1:(size(vAxis,2))/2,:) ; logOddsMapR((size(vAxis,2))/2+1:end,:)];
    PMap = [mapLeft(1:(size(vAxis,2))/2,:) ; mapRight((size(vAxis,2))/2+1:end,:)];
else
    logOddsCorrMap = [logOddsMapL(1:(size(vAxis,2)+1)/2,:) ; logOddsMapR((size(vAxis,2)+1)/2+1:end,:)];
    PMap = [mapLeft(1:(size(vAxis,2)+1)/2,:) ; mapRight((size(vAxis,2)+1)/2+1:end,:)];
end
if plotflag
    figure; set(gcf, 'Color', [1 1 1], 'Position', [100 100 450 350], 'PaperPositionMode', 'auto'); hold on;
    [~,h] = contourf(tAxis,vAxis,logOddsCorrMap,n);
    colorbar; set(h,'LineColor','none');
%     figure; [c,h] = contourf(tAxis,vAxis,logOddsMapL,n);
%     figure; [c,h] = contourf(tAxis,vAxis,logOddsMapR,n);
    caxis([0 round(B*k/2.2)]);
    [~,hh] = contour(tAxis,vAxis,logOddsCorrMap,[theta -theta]);
    set(hh,'LineColor','k','LineWidth',2);
    colorbar('YTick',linspace(0,round(B*k/2.2),6)); 
    xlabel('Time (ms)'); ylabel('DV');
    set(gca,'XTick',0:200:1000,'YTick',-B:round(B/2):B,'XTickLabel',0:200:1000,'TickDir','out');
%     changeAxesFontSize(gca,24,24); % for paper fig
    title('Log odds correct');
    colormap(jet);
end

