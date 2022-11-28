function [Yt, x, Y, durs,whichEventAlign] = trial_psth(spktimes,alignEvent,varargin)
% [Yt, x, Y, durs, whichEventAlign] = TRIAL_PSTH(spktimes,alignEvent,varargin)
%
% This code generates peri-event time histograms  for each trial using times in spktimes, aligned to a single event (alignEvent). 
% It also returns the total activity in the same interval, and duration of
% these intervals. 
% Currently, PSTHs are only computed using boxcar, and should be smoothed
% afterwards for visualization if desired.
%
% INPUTS: 
% required arguments
%   spktimes   - vector of spike times [nspk x 1]
%   alignEvent - times of task events for aligning spikes to
%         [ntr x 1, or ntr x 2] 
%  IMPORTANT - spktimes will always be aligned to the *FIRST* column in
%  alignEvent, regardless of the actual temporal order of the two events in
%  alignEvent
%   
% optional arguments
%       tStart, tEnd  - time relative to alignEvent to start and end psth
%                       if alignEvent is Ntrx1, tStart and tEnd will both be relative to that same alignEvent
%                       if alignEvent is Ntrx2, tStart will be relative to 1st event in alignEvent, tEnd to 2nd event
%       binSize     - the size of each spike count bin, default 0.05s (50ms)
%       alltrials   - default 0, only include trials after the first spike and before the final spike (inclusive)
%                   - 1 will keep all trials. if unit was not recorded for all trials (e.g. due to drifting signal), 
%                       this could artificially lower firing rate when later averaging across trials!
%       countORrate   - default 1, divide by binsize to give firing rate over time, and divide overall 'y' count by durs to get average
%                     firing rate. 0 returns raw counts
%
% OUTPUTS:
%       Yt   - firing rate or spike counts over time [Ntrs x nbins]
%       x    - time of Yt bins (center of bin) relative to alignment event [nbins x 1]
%       Y    - mean firing rate in interval [Ntrs x 1]
%       durs - duration of intervals [Ntrs x 1]
%       whichEventAlign - which event did we end up aligning to (i.e. was
%       the first column in alignEvent the earlier event or the later one)
% ----------
%
% A short essay on the usage of alignEvent...
%
% Note 1. The code has been tested with spike times in seconds, as the defaults for tStart, tEnd, and binSize are also all in seconds. 
% % However, explicitly specifying all necessary variables in different units should work (as long as one is consistent!). This hasn't been tested yet though.
%
% Note 2. In all cases, and for greater flexibility, the alignEvent argument should contain the event times themselves, not
% some rig/lab specific label or variable name (i.e. parse your data structure in a wrapper function in order to generate your alignEvent time values, then pass this in).
%
% Note 3. What this code does not do:
% a. align activity to multiple events
% % if you wish to do this, create a wrapper using this function, and loop over
% different alignEvent vectors (and tStart and tEnd values) to create a
% cell array of firing rate matrices, one matrix per alignment. Later
% splicing or stitching of these matrices is possible if desired (but note
% that each entry in x will be relative to its own 'time zero', and so the
% timing of each event relative to a global time zero would be needed for
% stitching.
%
% (Of course, if an event is a fixed interval from the alignEvent on every
% trial, implicit alignment to this event will also take place, at some
% non-zero time.)
%
% e.g.
% > alignEvent = {stimOn,stimOff,reward};
% > tStart = [-1,-0.5,-0.5];
% > tEnd   = [+2,+1,+0.5];
% > for iae=1:length(alignEvent)
%       [Yt{iae},x{iae}] = trial_psth(spktimes,alignEvent{iae},'tStart',tStart(iae),'tEnd',tEnd(iae));
%   end
%
% b. As of 09/2022, it does not perform any smoothing of psths or averaging across trials, and psth
% calculation is performed by simple boxcar. Future versions will include
% the option for more advanced windowing methods. 
%
% Note 4. Finally, some example use cases of different alignEvents
%
% 1. an Ntrx1 vector in alignEvent to align spikes to a fixed duration around one event 
% e.g.  if stimOn is a vector containing stimulus onset times
%   >> [Yt, x] = psth(spktimes, [stimOn],'tStart',-1,'tEnd',2,'binSize',0.1);
% will produce a matrix of firing rates aligned to stimOn, extending backwards by 1s and forwards by 2s, in 100ms bins.
%
% 2. an Ntrx2 array in alignEvent will align spikes to a fixed time (tStart) before one event, and extend to a fixed time (tEnd) after a second event. 
% % Time '0' will be the time of the first column in alignEvent.
% In this case, the psths will have variable durations on each trial (so there will be NaNs on one 'side' of the matrix for trials shorter than the max dur).
% These NaNs should make it straightforward to later average across trials with attrition (or perform other computations).
% 
% e.g. (1)
%  > [Yt, x] = psth(spktimes, [stimOn stimOff],'tStart',-1,'tEnd',1,'binSize',0.1);
% firing rates aligned to stimulus Onset in 100ms bins. PSTHs will extend back 1s from time zero (stimOn), and run until 1s after stimOff on
% each trial. Trials with shorter stimOn-stimOff intervals will have more NaNs at the end of their psth.
%
% e.g (2) 
%  > [Yt, x] = psth(spktimes, [stimOff stimOn],'tStart',-1,'tEnd',1,'binSize',0.1);
% firing rates aligned to stimulus Onset in 100ms bins. PSTHs will extend forward 1s from time zero (stimOff), and run back until 1s before stimOn on
% each trial. Trials with shorter stimOn-stimOff intervals will have more NaNs at the start of their psth.

% 3. an Ntrx2 array in alignEvent can also be used to calculate average firing rate between two events
% e.g. (1) from start to end of stimulus period
%  > [~,~,Y,durs] = psth(spktimes, [stimOn stimOff],'tStart',0,'tEnd',0);
% 
% e.g. (3) during the middle 1s of the stimulus
% > midStim = stimOn + (stimOff - stimOn)/2;
% > [~,~,Y,durs] = psth(spktimes, midStim,'tStart',-0.5,'tEnd',0.5);
% 
%
% -----------
%
% SJ 08-2022 created
% SJ 09-2022 added functionality for variable duration psths, and aligning
% to first column in alignEvent, regardless of order of events
% SJ 09-2022 cleaned up documentation
%
% TODOs: 
%
% - add alternatives to boxcar method of psth computation
% - plotAsHist option - if we want to plot psth as true 'boxcar', we need
% to duplicate every value in x and y
% - use some flags to skip binning sections of code if user only requests average rate in interval?
% - spin-off a function which calculates a fixed number of bins for trials
% of variable duration, by varying the bin size on each trial as needed
% (i.e. a form of time-warping)

%% parse the inputs

p = inputParser;
p.addRequired('spktimes');
p.addRequired('alignEvent');

assert(isvector(spktimes),'User must provide spike times as a vector')
assert(( isvector(alignEvent) || ismatrix(alignEvent)) && size(alignEvent,2)<=2,'alignEvent must be an Nx1 or Nx2 array')

var_names = {'tStart','tEnd','binSize','alltrials','countORrate'};
defaults = {-1,1,0.05,false,true};

for v=1:length(var_names)
    p.addParameter(var_names{v},defaults{v});
end

parse(p,spktimes,alignEvent,varargin{:});
p = p.Results;

%% deal with the alignEvent var

Ntr=length(alignEvent);
if (Ntr==0), disp('No good trials...'); return; end

% TODO: allow user to specify explicitly which column to align to?
if size(alignEvent,2)==2
    [alignEvent,EvOrder] = sort(alignEvent,2);
    whichEventAlign = EvOrder(1); % align to whichever event was originally first column in alignEvent
    % add a check in here to make sure that alignEvent temporal order is
    % consistent on each trial
else     % if alignEvent is a vector, tStart and tEnd are both relative to same event              
    alignEvent(:,2) = alignEvent(:,1);
    whichEventAlign = 1;
end

% get the start and end times from alignment event matrix
% AFTER sorting, so tStart is always relative to earlier event, tEnd to later event
tr_starts = alignEvent(:,1) + p.tStart;
tr_ends   = alignEvent(:,2) + p.tEnd;

durs   = tr_ends - tr_starts;

%% get indices of first and last trials that we actually want
itr_start = 1; itr_end   = Ntr;

if ~p.alltrials && ~isempty(spktimes)
    [~,t] = min(abs(tr_starts-spktimes(1)));
    if ~isempty(t), itr_start=t; end
    [~,t] = min(abs(tr_ends-spktimes(end)));
    if ~isempty(t), itr_end=t; end
end

%%  define the x-axis 'edges' for PSTH, and pre-allocate fr

% to get the overall extent of bins we need, we need to re-calculate the start OR end relative to the desired alignment event, since tStart and tEnd are
% initially defined relative to separate events in alignEvent. (if alignEvent has just one column, this is redundant, but still fine)
% e.g. if we set tStart and tEnd as zero, relative to stimOn and stimOff respectively, and align to stimOn, tStart_new will remain 0, 
% whereas tEnd_new will become the maximum time of stimOff(+tEnd)

if whichEventAlign==2
    tStarts_new = tr_starts - alignEvent(:,2); % tStarts_new will be <0
    tStart_new  = min(tStarts_new);
    tEnd_new    = p.tEnd;
elseif whichEventAlign==1
    tEnds_new   = tr_ends - alignEvent(:,1);
    tEnd_new    = max(tEnds_new);
    tStart_new  = p.tStart;
end

if p.tStart<0

    % this way, the event you are aligning to has a bin on either side, it
    % does not fall within a bin 
    x = 0:-p.binSize:tStart_new; % backwards from zero
    x1= 0:p.binSize:tEnd_new; % forwards from zero

    % deal with the ends
    if abs(x(end)-tStart_new)>0.01, x(end+1)=x(end)-p.binSize; end
    if abs(x1(end)-tEnd_new)>0.01, x1(end+1)=x1(end)+p.binSize; end

    % for some reason, tEnd_new-x1 has some small non-zero component
    % if x(end)~=p.tStart, x(end+1)=x(end)-p.binSize; end
    % if x1(end)~=tEnd_new, x1(end+1)=x1(end)+p.binSize; end

    % stitch them together now
    x=[fliplr(x) x1(2:end)];

else
    % if both tStart_new and tEnd_new are >0, then just do it the simple way
    x = tStart_new:p.binSize:tEnd_new;
end


nBins = length(x)-1;
Yt = nan(Ntr,nBins);
Y  = nan(Ntr,1);

% now compute the psth within dur limits on valid trials
for itr=itr_start:itr_end
   if isnan(alignEvent(itr,whichEventAlign)), continue; end
   
   spkInd = spktimes >= tr_starts(itr) & spktimes <= tr_ends(itr);
   
   Y(itr) = sum(spkInd);    
   if Y(itr)==0, Yt(itr,:) = 0; continue; end
   
   index1=spktimes(spkInd)-alignEvent(itr,whichEventAlign);
   
   Yt(itr,:) = histcounts(index1,x);

   % if we have two events (with variable delay between them) in alignEvent, discard (i.e. set to NaN) the
   % bins on trials which extend beyond the dur on that trial 
   % if aligning to earlier event,  bins will be at the back end of the trials. 
   % if aligning to the later event, these bins will be at the start.
   if whichEventAlign==1
       [~,endPos] = min(abs(x-tEnds_new(itr)));
       Yt(itr,endPos:end) = NaN;
   elseif whichEventAlign==2
       [~,startPos] = min(abs(x-tStarts_new(itr)));
       Yt(itr,1:startPos) = NaN;
   end
   
end

% shift 'x' from edges of bins to middle of bins, helpful for later
% plotting of Yt at the center of each bin
% TODO, make this flexible, we might want to duplicate bin x and vals to
% plot old-school histogram style
% and if we have a causal or gaussian filter for calculating spike counts,
% this could be different

x=x+diff(x(1:2))/2;
x(end)=[]; % cut off the end, so that x is nBins long
% note that this means there will be no zero in x - the closest values will be -binSize/2 and +binsize/2. 


if p.countORrate==1  
    Yt = Yt / p.binSize; % normalize spike count by bin width to get firing rate!
    Y = Y ./ durs;       % while we're at it, normalize total spike count by total interval time
end