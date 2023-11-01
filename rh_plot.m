function [fr,edges,hhist,hrast]=rh_plot(spktimes,align,varargin)
% generate (and plot if desired) raster and PSTH of neural firing
% plotting splits current axis into two subplots, upper one for raster and lower one for histogram
%
% standard boxcar will look like an old-school histogram
% to plot alignment to multiple events, create a wrapper that loops over
% multiple alignment events, with corresponding tmin,tmax, otherevents, 
%
% to split rasters and PSTHs by e.g. condition, specify 'split', a grouping
% variable matched to the number of trials

%=============
%INPUT PARAMETERS

% REQUIRED 
%
% spktimes      times of spikes (seconds)
%
% align         event times for alignment (seconds) e.g. motionOnset times (should be 1 x ntrials)
%
% OPTIONAL
%
% tmin, tmax    limits for x-axis (seconds) relative to alignment event, default -2s to +3s.
%
% otherevents   other events to plot on raster, cell array of arrays,
%               should each match the length of align (i.e. 1 x ntrials)
%
% eventsnames   names of events, the first one should be the name of alignemnt
%               variable, so length will be length of otherevents+1
%
% ratprc        proportion of current axis which will be devoted to raster plot, default 0.5
%
% bin           bin size for histograms (seconds)
%
% title         title which will appear on top of raster plot
%
% sort          an array, length equal to number of trials, contains a variable
%               according to which perform sorting of trials (in descending
%               order)
%
% split         split trials into groups (grouped in colored blocks on raster, and traces of matching color on PSTH)
%               like sort, should be 1 x ntrials, with each entry
%               indicating the group for each trial e.g. coherence
%
% if sort and split are both provided, trials will first by split, and then
% sorted within each group.
%
% alltrials     include all good trials with valid alignment event, or only
% include trials within spktimes
%
% hist_method   boxcar only for now, TODO add gaussian, causal methods
%
% hist_smbins   how many bins to run moving average over, default is 0,
% which means no moving average
% 
% rejecteventsoutliers      option for trial rejection (falling 1.5x
% outside IQR) e.g. RT. should be a vector matching the length of align
%
% symbol and color          strings defining symbol and color of otherevents
%
% alignlinecolor            color of line for alignment event (time 0), default black
%
% plotsingletrialevents     if 1, show timing of other events on each trial  on raster 
%                           (if 0, then just show median time. PSTH always shows median time)

%% parse inputs 

p = inputParser;
p.addRequired('spktimes');
p.addRequired('align');

var_names={'tmin','tmax','spikecolor','symbols','markersize','eventscolor','otherevents','eventsnames','title','ratprc','bin',...
    'sort','legend','plotsingletrialevents',...
    'rejecteventsoutliers','split','splitcolors','splitlabels','plotsplitevents',...
    'histmethod','alignlinecolor','alltrials',...
    'hist_smbins','plot'};
defaults={-2,3,'k','o*x+^vso*x+^v',4,'cmgbrycmgbry',{},{},'',0.5,0.05,...
    '',1,0,...
    [],[],'kmrcbgykmrcbgy',{},0,...
    'boxcar','k',0,...
    0,1};

for v=1:length(var_names)
    p.addParameter(var_names{v},defaults{v});
end

parse(p,spktimes,align,varargin{:});
p = p.Results;


%% convert colors from char to rgb values, if necessary

if p.plot
    if ischar(p.alignlinecolor), p.alignlinecolor=char2rgb(p.alignlinecolor); end
    if ischar(p.spikecolor), p.spikecolor=char2rgb(p.spikecolor); end
    Nov=length(p.otherevents);
    
    if ischar(p.eventscolor)
        eventscolor=p.eventscolor;
        p.eventscolor=NaN*ones(Nov,3);
        for ie=1:Nov, p.eventscolor(ie,:)=char2rgb(eventscolor(ie));end
        
    end
    
    %% set axes
    
    cla
    pos=get(gca,'Position');
    hhist=subplot('position',[pos(1)*0.9 pos(2) pos(3)*1.1 pos(4)*(1-p.ratprc)*0.95]); %axis off;
    set(hhist,'ylim',[0 1e-10]);
    hrast=subplot('position',[pos(1)*0.9 pos(2)+pos(4)*(1-p.ratprc) pos(3)*1.1 pos(4)*p.ratprc*0.95]); %axis off;
    set(hrast,'xtick',[]); hrast.XAxis.Visible = 'off';
%     ylabel('trial #','fontsize',14)
end

fr=[];
edges=[];

%% discard NaNs in align (i.e. 'bad trials')

nans=isnan(align);
align(nans)=[];
for ii=1:length(p.otherevents)
    p.otherevents{ii}(nans)=[];
end

if ~isempty(p.sort), p.sort(nans)=[]; end
if ~isempty(p.split), p.split(nans)=[]; end
Ntr=length(align);
if (Ntr==0), return; end

%% plot event lines
% what is this for?
% if p.plot
%     Nov=length(p.otherevents);
%     line(0,0,'color',p.alignlinecolor,'linewidth',0.5);
%     for ie=1:Nov, line(0,0,'color',p.eventscolor(ie,:),'marker',p.symbols(ie),'markersize',p.markersize); end
% end

%% check some inputs
% reject trials outliers based on specified variable
% below Q1-1.5*IQR, or above Q3+1.5*IQR

if ~isempty(p.rejecteventsoutliers)
    w=1.5;
    [Q]=prctile(p.rejecteventsoutliers,[25 75]);
    q1=Q(1); q3=Q(2);

    rej=find((p.rejecteventsoutliers<q1-w*(q3-q1))|(p.rejecteventsoutliers>q3+w*(q3-q1)));
    if ~isempty(rej)
        fprintf('%d trials were rejected as outliers\n',length(rej));
        align(rej)=[];
        Ntr=length(align);
        
        for ie=1:Nov
            if isempty(p.otherevents{ie}), continue; end
            p.otherevents{ie}(rej)=[];
        end
        if ~isempty(p.split), p.split(rej)=[]; end
        if ~isempty(p.sort), p.sort(rej)=[]; end
    end
end

% find first and last trials with at least one spike, if alltrials=0;
% if p.alltrials 
%     tr_start=1;
%     tr_end=Ntr;
% else
%     t=find((spktimes(1)-(align+p.tmin))<0);
%     if ~isempty(t), tr_start=t(1); else, tr_start=1; end
%     t=find((spktimes(end)-(align+p.tmax))>0);
%     if ~isempty(t), tr_end=t(end); else, tr_end=Ntr; end
% end

tr_start=1;
tr_end=Ntr;
if ~p.alltrials
    [~,t] = min(abs(align+p.tmin-spktimes(1)));
    if ~isempty(t), tr_start=t; end
    [~,t] = min(abs(align+p.tmax-spktimes(end)));
    if ~isempty(t), tr_end=t; end
end
    
    
% if p.plot && ~isempty(p.sort) && ~isempty(p.split)
%     fprintf('Can''t have both sorting and splitting\n');
%     return;
% end

if p.tmax<p.tmin
    fprintf('tmin (%1.2f) should be smaller than tmax (%1.2f)\n',p.tmin,p.tmax);
    return;
end

if p.plot && ~isempty(p.split)
    if isempty(p.sort)
        p.sort = p.split;
    end
    cas = unique(p.split);
    cas = cas(~isnan(cas));
    clear t
    if ischar(p.splitcolors), for ic=1:length(cas), t(ic,:)=char2rgb(p.splitcolors(ic)); end
        p.splitcolors=t;
    end
end

%% sort trials if needed (for raster viz)

% can't sort and split at the same time, this still needs work
% to be flexible

if ~isempty(p.sort) 
   p.sort(1:tr_start-1) = NaN;
   p.sort(tr_end+1:end) = NaN;

%    temp = [p.split p.sort];
%    [ttt,ind] = sortrows(temp);

   [ttt,ind] = sort(p.sort);

% % why is this necessary??
%    tr_start = 1;
%    t = find(isnan(ttt(:,1)));
%    if isempty(t)
%        tr_end = length(ttt,1);
%    else
%        tr_end = t(1)-1;
%    end

   if ~isempty(p.split)
       p.split=p.split(ind);
   end
   p.sort=p.sort(ind);
   for ie=1:Nov
       if isempty(p.otherevents{ie}), continue; end
       p.otherevents{ie}=p.otherevents{ie}(ind);
   end
   align=align(ind);
end

%%  define the x-axis 'edges' for PSTH, and pre-allocate fr

switch p.histmethod
    case 'boxcar'
        
        edges=0:-p.bin:p.tmin; % backwards from zero
        edges1=0:p.bin:p.tmax; % forwards from zero
        
        % deal with the ends
        if edges(end)~=p.tmin, edges(end+1)=edges(end)-p.bin; end
        if edges1(end)~=p.tmax, edges1(end+1)=edges1(end)+p.bin; end
        
        if p.hist_smbins>1 %
            if p.hist_smbins/2==floor(p.hist_smbins/2), p.hist_smbins=p.hist_smbins+1; end
            smbins2 = floor(p.hist_smbins/2);

            % extend the ends here to avoid edge effects when smoothing
            edges(end+1:end+smbins2)=edges(end)-1*p.bin:-p.bin:edges(end)-smbins2*p.bin;
            edges1(end+1:end+smbins2)=edges1(end)+1*p.bin:+p.bin:edges1(end)+smbins2*p.bin;
            
        end
        
        edges=[edges(end:-1:1) edges1(2:end)];
        ne=length(edges);
        fr=nan(Ntr,ne-1);

    case 'causal'
        
        
    case 'gaussian'
        
end

%% begin plotting raster 

%plot other event lines, or markers
ot=cell(1,Nov);
for ie=1:Nov
   if isempty(p.otherevents{ie}), continue; end
   ot{ie}=nan(1,Ntr);
   for itr=tr_start:tr_end
       ot{ie}(itr)=p.otherevents{ie}(itr)-align(itr);
      if p.plotsingletrialevents
          line(ot{ie}(itr),itr,'color',p.eventscolor(ie,:),'marker',p.symbols(ie),'linestyle','none','markersize',p.markersize);
      end
   end
   if ~p.plotsingletrialevents
       line(ones(1,2)*nanmedian(ot{ie}),[tr_start tr_end],'color',p.eventscolor(ie,:),'linewidth',0.5);
   end
end
% plot alignment event line
line([0 0],[tr_start-0.5 tr_end+0.5],'color',p.alignlinecolor,'linewidth',1);

% plot actual rasters - spike times
for itr=tr_start:tr_end
   if isnan(align(itr)), continue; end
   be=align(itr)+p.tmin;
   en=align(itr)+p.tmax;
   
   if ~isempty(p.split)
       c=find(cas==p.split(itr));
       c=c(1);
       p.spikecolor(1,1:3)=p.splitcolors(c,1:3);
   end
   
   inds=find((spktimes>be-p.hist_smbins*p.bin)&(spktimes<en+p.hist_smbins*p.bin));%spike indices for the trial
   nspk=length(inds);
   
   if isempty(inds), continue; end
   
   index1=spktimes(inds)-align(itr);

   % create the x,y vectors to plot raster 'ticks'
   t=spktimes(inds);
   tt=(ones(2,1)*(t'-be+p.tmin));
   tt(3,:)=NaN;
   ttx=tt(:)';
   tty=ttx;
   tty(1:3:end)=itr-0.3;
   tty(2:3:end)=itr+0.3;
   
   line(ttx,tty,'color',p.spikecolor,'linewidth',0.5); % PLOT THE TRIAL RASTER
   
  % populate the firing rates for histogram
   switch p.histmethod
       case 'boxcar'
           fr(itr,:) = histcounts(index1,edges);
           if p.hist_smbins>1
               fr(itr,:) = smooth(fr(itr,:),p.hist_smbins);
           end
       case 'causal'
       case 'gaussian'
   end
      
end

set(gca,'ydir','reverse');
xlim([p.tmin p.tmax]);
ylim([tr_start-0.5 tr_end+0.5]);
% if ~isempty(p.eventsnames)&&(p.legend)&&(~p.plotsingletrialevents)
%    legend(p.eventsnames,'Location','NorthEastOutside');
%    if length(p.eventsnames)~=length(p.otherevents)+1 
%       warning('plot_raster_hist:EventsNames','Check EventsNames');
%    end
% end

set(gca,'xtick',[]);
tpos=get(gca,'Position');
set(gca,'Box','off');
if ~isempty(p.title)
    ht=title(p.title,'fontsize',14);
    set(ht,'position',[mean([p.tmin p.tmax]) ht.Position(2) ht.Position(3)])
end


%% plot histograms

if p.plot
hhist=subplot('position',[tpos(1) pos(2) tpos(3) pos(4)*(1-p.ratprc)*0.98]);cla; hold on;
histpos=get(gca,'Position');
% ylabel('spikes s^{-1}','fontsize',14)   
end

% if we're splitting trials, need to loop over splits

if ~isempty(p.split)
    for ic=1:length(cas)
        plot(0,0,'color',p.splitcolors(ic,:),'linewidth',0.5);
    end
    
    mx=0;
    for ic=1:length(cas)
        t=nanmean(fr(p.split==cas(ic),:),1); % get fr just for t
        switch p.histmethod
            case 'boxcar'
                t=t/p.bin; % normalize by bin width to get firing rate
                if p.hist_smbins>1
                    x=edges+diff(edges(1:2))/2;x(end)=[]; % if smoothing, set x as the middle of the bins
                else % if true boxcar, double up the x and y vals for histogram
                    t=[t;t];
                    t=t(:);
                    x=[edges;edges];
                    x=x(:);
                    x(1)=[];x(end)=[];
                end
        end
        mx=max(mx,max(t)*1.1);
        if mx==0, mx=1; end

        if p.plot
        if ~all(isnan(t)) % if fr is not nan
            line(x,t,'color',p.splitcolors(ic,:),'linewidth',1.5,'Tag','hist');
            set(gca,'ylim',[0 mx],'xlim',[p.tmin p.tmax]);
            % plot the 'other events'
            for ie=1:Nov
                if isempty(p.otherevents{ie}), continue; end
                if p.plotsplitevents
                    line(ones(1,2)*nanmedian(ot{ie}(p.split==cas(ic))),[0 max(mx,200)],'color',p.splitcolors(ic,:),'linewidth',0.5,'Tag','eventline');
                end
            end
        end
        
        % plot other event lines, IQR range, labels
        if ~p.plotsplitevents
            for ie=1:Nov
                line(ones(1,2)*nanmedian(ot{ie}),[0 max(200,mx)],'color',p.eventscolor(ie,:),'linewidth',0.5,'Tag','eventline');
                %             prcev=prctile(ot{ie},[25 75]);
                %             patch([prcev fliplr(prcev)],[0 0 max(200,mx) max(200,mx)],p.eventscolor(ie,:),'facealpha',0.2,'edgealpha',0)
%                 if nanmedian(ot{ie})>p.tmin && nanmedian(ot{ie})<p.tmax
%                     text(nanmedian(ot{ie}),5+max(200,mx),p.eventsnames{ie+1},'color',p.eventscolor(ie,:),'fontsize',8,'horizontalalignment','left','rotation',45)
%                 end
            end
        end
        line([0 0],[0 max(200,mx)],'color',p.alignlinecolor,'linewidth',1,'Tag','eventline');
%         text(0,5+max(200,mx),p.eventsnames{1},'color',p.alignlinecolor,'fontsize',8,'horizontalalignment','left','rotation',45)
        end
    end
    if p.plot && ~isempty(p.splitlabels)&&(p.legend)
        legend(p.splitlabels,'Location','NorthEast');
        set(gca,'Position',histpos);
    end
    
else % no splitting, this is simpler
    t=nanmean(fr,1);
    switch p.histmethod
        case {'boxcar'}
            t=t/p.bin;
            
            if p.hist_smbins>1
                x=edges+p.bin/2;x(end)=[];
            else
                t=[t;t];
                t=t(:);
                x=[edges;edges];
                x=x(:);
                x(1)=[];x(end)=[];
            end
            
    end
    
    if p.plot && ~all(isnan(t))
        mx=max(t)*1.1;
        plot(x,t,'color','k','linewidth',1.5,'Tag','hist');
        set(gca,'ylim',[0 mx+1e-10],'xlim',[p.tmin p.tmax]);
        for ie=1:Nov
            if isempty(p.otherevents{ie}), continue; end
            line(ones(1,2)*nanmedian(ot{ie}),[0 max(200,mx)],'color',p.eventscolor(ie,:),'linewidth',0.5,'Tag','eventline');
        end
        line([0 0],[0 max(mx,200)],'color',p.alignlinecolor,'linewidth',1,'Tag','eventline');
%         text(0,0-(mx*0.02),p.eventsnames{1},'color',p.alignlinecolor,'fontsize',8,'horizontalalignment','right','rotation',45)
    end
end
if p.plot, set(gca,'Box','off'); end

end

function rgbc=char2rgb(cc)
%convert colors from char to rgb values
rgbc=zeros(length(cc),3);
for irgbc=1:length(cc)
    rgbc(irgbc,1:3)=rem(floor((strfind('kbgcrmyw', cc(irgbc)) - 1) * [0.25 0.5 1]), 2);
end

end
