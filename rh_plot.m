function [fr,edges,hhist,hrast] = rh_plot(spktimes,align,varargin)
% plot raster and PSTH of neural activity for a single array of spiketimes
% and array of corresponding 
%
% the current axis when this function is called gets split in two
% - upper one for raster and lower one for histogram
%
% standard boxcar will look like an old-school histogram
%
% to plot alignment to multiple events, create a wrapper that loops over
% multiple alignment events, with corresponding tmin,tmax, otherevents, etc
%
% to split rasters and PSTHs by e.g. condition, specify 'split', a grouping
% variable matched to the number of trials
%
% note that defaults for tmin, tmax and binsize assume that spike counts
% are specified in seconds. However, it should be possible to use other
% units e.g. milliseconds, as long as one is consistent across all these
% parameters!
%
%
%
%=============
%INPUT PARAMETERS

% REQUIRED 
%
% spktimes      times of spikes
%
% align         event times for alignment (seconds) e.g. motionOnset times (should be 1 x ntrials)
%
% OPTIONAL
%
% tmin, tmax    limits for x-axis (seconds) relative to alignment event, default -2s to +2s.
%
% ratprc        proportion of current axis which will be devoted to raster plot, default 0.5
%
% bin           bin size for spike counts (in seconds), default is 0.01s (10ms)
%
% sm_method     smooth method (boxcar or none only for now) TODO add gaussian, causal
%
% sm_width      width of smoothing filter (default 0.25s/250ms)
%
% sm_sigma      sigma of smoothing kernel (only applicable for gaussian, default 0.05s/50ms)
%
% split         split trials into groups (grouped in colored blocks on raster, and traces of matching color on PSTH)
%               like sort, should be 1 x ntrials, with each entry
%               indicating the group for each trial e.g. coherence
%
% alltrials     true - include all good trials with valid alignment event,
%               false - only include trials within spktimes,
%               default is false
%
% otherevents   times of other events to plot on raster relative to alignment event, 
%               cell array of arrays, should each match the length of align (i.e. 1 x ntrials)
%               default is empty
% othersyms     symbols for otherevents
% othercols     colors for otherevents
%
% event_names   names of events, the first one should be the name of alignemnt
%               variable, so length will be length of otherevents+1, default is empty
%
% plotsingletrialevents     true (default): show timing of other events on each trial on raster 
%                           false: then just show median time. PSTH always shows median time)
%  

%% parse inputs 

p = inputParser;
p.addRequired('spktimes');
p.addRequired('align');

var_names={'tmin', 'tmax', 'ratprc','bin', 'sm_method', 'sm_width', 'sm_sigma', ...
    'hue_style', 'hue_colors', 'style', 'style_lines', 'splitlabels', 'plotsplitevents',...
    'alltrials', ...
    'otherevents', 'event_syms','event_colors','event_names', 'plotsingletrialevents', ...
};

defaults = {-2, 3, 0.5, 0.05, 'boxcar', 0.25, 0.05, ...
    [], 'kmrcbgykmrcbgy', [], {'-', '--', ':', '.-'}, {}, 0,...
    false, ... 
    {}, 'o*x+^vso*x+^v', 'cmgbrycmgbry', {}, true, ...
};

for v=1:length(var_names)
    p.addParameter(var_names{v},defaults{v});
end

parse(p,spktimes,align,varargin{:});
p = p.Results;

assert(p.tmin < p.tmax, 'tmin (%1.2f) must be less than tmax (%1.2f)', p.tmin, p.tmax);

if ~isempty(p.hue_style)
    assert(size(p.hue_style, 2)==1 || size(p.hue_style, 2)==2, 'hue_style must be a 1 or 2-column vector')
    assert(size(p.hue_style, 1)==length(align), 'hue_style must have the same number of rows as align')
end

% placeholders in case of break, return empties
fr = [];
edges = [];

%% discard NaNs in align (i.e. 'bad trials'), if any

nans=isnan(align);
align(nans)=[];
for ii=1:length(p.otherevents)
    p.otherevents{ii}(nans)=[];
end

if ~isempty(p.hue_style), p.hue_style(nans, :)=[]; end

Ntr = length(align);
if Ntr == 0
    disp('Data does not seem to contain any valid trials')
    return; 
end

%% convert colors from char to rgb values, if necessary

Nov=length(p.otherevents);
if ischar(p.event_colors)
    event_colors = p.event_colors;
    p.event_colors = nan(Nov,3);
    for ie=1:Nov
        p.event_colors(ie,:)=char2rgb(event_colors(ie));
    end
end

assert(Nov <= size(p.event_colors, 1), 'Insufficient colors specified for number of other events provided');

if ~isempty(p.hue_style)
    if ischar(p.hue_colors)
        Nhue = length(p.hue_colors);
        hue_colors = p.hue_colors;
        p.hue_colors = nan(Nhue, 3);
        for ic=1:Nhue
            p.hue_colors(ic,:) = char2rgb(hue_colors(ic));
        end
    end
    assert(length(unique(p.hue_style(:,1))) <= size(p.hue_colors, 1), 'Insufficient colors specified for number of hues provided');

    % get the unique styles. if no styles, then just set all to the same
    if size(p.hue_style,2)==2
        unq_styles = unique(p.hue_style(:,2));
    else
        unq_styles = 1;
        p.hue_style(:, 2) = 1;
    end
    unq_hues = unique(p.hue_style(:,1));

    assert(length(unique(p.hue_style(:,2))) <= length(p.style_lines), 'Insufficient line types specified for number of styles provided');
end


%% find first and last trials with at least one spike, if alltrials is false
% this stops us from plotting a bunch of empty raster rows if the unit
% wasn't present in a bunch of trials at the start or the end of the behavior

tr_start=1;
tr_end=Ntr;
if ~p.alltrials
    % this [~, t] = min(abs(... syntax is the classic way of finding a 'nearest' index
    % in this case, t will be the trial index closest to the first spike time
    [~,t] = min(abs(align+p.tmin-spktimes(1)));
    if ~isempty(t), tr_start=t; end

    % repeat for the last trial, overwrite tr_end if valid
    [~,t] = min(abs(align+p.tmax-spktimes(end)));
    if ~isempty(t), tr_end=t; end
end
    

%% for plotting the rasters and conditional PSTHs, re-order the trials here according to hue/style

if ~isempty(p.hue_style)
    [sorted_hs, ind] = sortrows(p.hue_style, [2, 1]); % sort by style first, then hue

    % in case any NaNs were passed as hue, skip these trials
    % do this by updating tr_end, since NaNs will have been sorted to the end
    t = find(isnan(sorted_hs(:,1)));
    if ~isempty(t)
        tr_end = t(1)-1;
    end

    % now re-order the trials in align and otherevents, so that later on we can just plot them in order
    align = align(ind);

    for ie=1:Nov
        if isempty(p.otherevents{ie}), continue; end
        p.otherevents{ie} = p.otherevents{ie}(ind);
    end
end

%%  define the x-axis bin 'edges' for PSTH

% if we are going to smooth the firing rates, extend x so that the first
% and last bins smoothed values get to use the real data before it's cut

x_tmin = p.tmin;
x_tmax = p.tmax;
sm_nbins = floor(p.sm_width / p.bin);

if ~strcmp(p.sm_method, 'none')
    x_tmin = p.tmin - p.bin*sm_nbins/2;
    x_tmax = p.tmax + p.bin*sm_nbins/2;
end

if x_tmin < 0 && x_tmax > 0
    edges  = 0:-p.bin:x_tmin; % backwards from zero
    edges1 = 0:p.bin:x_tmax; % forwards from zero

    % deal with the ends, in case they go too far
    if edges(end)~=x_tmin, edges(end+1)=edges(end)-p.bin; end
    if edges1(end)~=x_tmax, edges1(end+1)=edges1(end)+p.bin; end

    edges  = [fliplr(edges) edges1(2:end)];
else
    edges = x_tmin:p.bin:x_tmax;
end
        
% pre-allocate fr
fr = nan(Ntr, length(edges)-1);

%% set up axes

cla
pos=get(gca,'Position');

hhist = subplot('position',[pos(1)*0.9 pos(2) pos(3)*1.1 pos(4)*(1-p.ratprc)*0.95]); 
hhist.YLim = [0 1e-10];

hrast = subplot('position',[pos(1)*0.9 pos(2)+pos(4)*(1-p.ratprc) pos(3)*1.1 pos(4)*p.ratprc*0.95]); 
hrast.XTick = []; hrast.XAxis.Visible = 'off';

%% begin actual plotting, starting with raster plot

% fix the axes ranges, and flip the raster y-axis - 'trial 1' will then be at the top
set(hrast,'ydir','reverse');
xlim([p.tmin p.tmax]);
ylim([tr_start-0.5 tr_end+0.5]);

% loop over trials to plot raster ticks
for itr = tr_start:tr_end
%    if isnan(align(itr)), continue; end % not really necessary, nans removed earlier
   
   tmin_tr = align(itr) + x_tmin;
   tmax_tr = align(itr) + x_tmax;
   
   % set the tick color 
   if ~isempty(p.hue_style)
       c = find(unq_hues(:,1)==sorted_hs(itr, 1), 1);
       spike_color(1, 1:3) = p.hue_colors(c,1:3);
   else
       % default raster tick color is black
       spike_color = 'k';
   end
   
   inds = (spktimes > tmin_tr) & (spktimes<tmax_tr); %spike indices for the trial
   
   if sum(inds)==0, continue; end % no spikes on this trial
   
   % spiketimes relative to alignment event on this trial
   index1 = spktimes(inds) - align(itr);

   % now create the x,y vectors to plot raster 'ticks'

   % x position is dictated by relative spiketimes
   ttx = repmat(index1', 2, 1); % duplicate to have a top and bottom for each tick
   ttx(3,:) = NaN; % add a nan row so that each tick is separated from the next when we plot
   ttx = ttx(:)';

   % y position array should be the same size as x position
   tty = nan(size(ttx));
   tty(1:3:end) = itr-0.3; % y position is trial number, with top and bottom for each tick slightly separated
   tty(2:3:end) = itr+0.3;

   % now plot all the ticks for this trial, in their color
   line(ttx, tty, 'color', spike_color, 'linewidth',0.5);

   % populate the firing rates for histogram
   fr(itr,:) = histcounts(index1,edges);

   % TODO switch this to a convolution with a kernel, then we can define the
   % kernel outside the loop for each type of smoothing
   switch p.sm_method
       case 'boxcar'
           fr(itr,:) = smooth(fr(itr,:), sm_nbins);
       case 'causal'
       case 'gaussian'
   end
end

%% plot other event lines (or markers) on the rasters

ot = cell(1,Nov); % set a handle for event lines/markers

for ie=1:Nov
    if isempty(p.otherevents{ie}), continue; end
    ot{ie} = nan(1, Ntr);
    
    % for each trial, get the relative time of the 'other' event from the
    % alignment event, (and plot it if we want to show single trial events)
    for itr = tr_start:tr_end
        ot{ie}(itr) = p.otherevents{ie}(itr)-align(itr);
        if p.plotsingletrialevents
            line(ot{ie}(itr), itr, 'color', p.event_colors(ie,:),'marker',p.event_syms(ie),'linestyle', 'none', 'markersize', 4);
        end
    end

    % if only showing average, plot the median relative event time 
    if ~p.plotsingletrialevents
        line(ones(1,2) * nanmedian(ot{ie}), [tr_start tr_end], 'color', p.event_colors(ie,:), 'linewidth', 0.5);
    end
end

% plot alignment event line, at x=0 by definition
line([0 0],[tr_start-0.5 tr_end+0.5], 'color', 'k','linewidth',1);


% if ~isempty(p.event_names)&&(p.legend)&&(~p.plotsingletrialevents)
%    legend(p.event_names,'Location','NorthEastOutside');
%    if length(p.event_names)~=length(p.otherevents)+1 
%       warning('plot_raster_hist:EventsNames','Check EventsNames');
%    end
% end
%         text(0,5+max(200,mx),p.event_names{1},'color',p.alignlinecolor,'fontsize',8,'horizontalalignment','left','rotation',45)


set(gca,'xtick',[], 'box', 'off');
rastpos = get(gca, 'Position'); % we will use this below to set the histogram's position


%% plot PSTHs below

hhist = subplot('position',[rastpos(1) pos(2) rastpos(3) pos(4)*(1-p.ratprc)*0.98]); 
cla; hold on;
set(hhist,'Box','off');
histpos = get(gca,'Position');

mx = 1; % default max firing rate, for setting y-axis lims (will be overwritten below)

if ~strcmp(p.sm_method, 'none')
    x = edges + diff(edges(1:2))/2; 
    x(end) = []; % if smoothing, set x as the middle of the bins
else % if true boxcar, double up the x and y vals for a 'blocky' histogram
    x= [edges; edges];
    x=x(:);
    x(1)=[]; x(end)=[];
end


if ~isempty(p.hue_style)

    % loop over unique hues and styles
    for ih = 1:length(unq_hues)

        for ist=1:length(unq_styles)
            
            hs_inds = (sorted_hs(:,1) == unq_hues(ih)) & (sorted_hs(:,2) == unq_styles(ist));

            % get the average firing rate across trials of this condition
            % (normalize by bin width to get firing rate)
            t = nanmean(fr(hs_inds,:), 1) / p.bin; 

%             if strcmp(p.sm_method, 'none')
%                 t = [t; t]; % TODO use repmat
%                 t=t(:);
%             end

            if ~all(isnan(t)) % if t is not NaN, plot it!
                line(x, t, 'color', p.hue_colors(ih,:), 'linestyle', p.style_lines{ist}, 'linewidth', 1.5, 'Tag', 'hist');
                mx=max(mx, max(t)*1.1); % update max firing rate limit for y-axis lims
            end
        end
    end
    
else 
    % no splitting, just plot one black line

    t = nanmean(fr,1) / p.bin;
            
    if ~all(isnan(t))
        line(x, t, 'color','k','linewidth',1.5,'Tag','hist');
        mx=max(mx, max(t)*1.1); % update max firing rate limit for y-axis lims
    end
end

% plot other event lines
for ie=1:Nov
    line(ones(1,2)*nanmedian(ot{ie}),[0 max(200,mx)], 'color', p.event_colors(ie,:), 'linewidth', 0.5,'Tag','eventline');
end

% plot alignment event line
line([0 0], [0 max(200,mx)], 'color', 'k','linewidth',1,'Tag','eventline');



set(hhist,'ylim',[0 mx+1e-5], 'xlim', [p.tmin p.tmax]);


end


function rgbc=char2rgb(cc)
    % convert colors from char to rgb values
rgbc=zeros(length(cc),3);
for irgbc=1:length(cc)
    rgbc(irgbc,1:3)=rem(floor((strfind('kbgcrmyw', cc(irgbc)) - 1) * [0.25 0.5 1]), 2);
end

end
