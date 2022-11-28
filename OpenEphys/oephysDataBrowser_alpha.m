%% OpenEphys data browser

% in progress, begun 10-19-20 CF
% must set localDir and currentFileList before running (getDataFromServer does the latter)


keyboard

% 2021-02-16_16-21-35_map_Ch16_Sorting.psort


close all

channelMapping = [fliplr(9:24) 25:32];
channelMapping = [channelMapping channelMapping+32]; % for 2-probe expts

for m = 1:length(currentFileList)
    files = dir([localDir currentFileList{m}]);
    disp(['beginning folder /' currentFileList{m} ' -- press a key to continue']); pause; 
    
%     % continuous snapshot:
%     figure(m); set(gcf,'Color',[1 1 1],'Position',[300+20*m 300 1000 1100],'PaperPositionMode','auto'); clf;
%     
%     % tuning curves
%     figure(m*10); set(gcf,'Color',[1 1 1],'Position',[1000+20*m 200 1000 1100],'PaperPositionMode','auto'); clf;
     
    for n=3:length(files) % skip 1+2, they are are "." and ".."
        filename = fullfile(localDir,currentFileList{m},files(n).name);
        dot = strfind(files(n).name,'.');
        ext = files(n).name(dot(end)+1:end);  % dot(end) because spikes files have two dots
        disp(['processing ' currentFileList{m} '/' files(n).name])
        switch ext
            case 'continuous'
                [data, timestamps, info] = load_open_ephys_data_faster(filename);
                chan = str2double(info.header.channel(3:end));
                und = strfind(files(n).name,'_');
                if length(und)==1
                    run = 1;
                else
                    run = str2double(files(n).name(und(end)+1:dot-1));
                end
                figure((m-1)*100+run); set(gcf,'Color',[1 1 1],'Position',[300+20*m 300+20*run 800 1000],'PaperPositionMode','auto');
                try
                    subplot(24,1,chan-8); plot(data(1:30000));
                catch
                    subplot(24,1,chan); plot(data(1:30000));
                end
                set(gca,'yticklabel',[],'box','off');
                ylabel(num2str(chan));
                if chan==9
                    fileForTitle = currentFileList{m};
                    fileForTitle(strfind(fileForTitle,'_')) = '.';
                    title([fileForTitle ' - run' num2str(run)]);
                end
             
            case 'events'
                [data, timestamps, info] = load_open_ephys_data_faster(filename);

            
%             case 'spikes'
%                 [data, timestamps, info] = load_open_ephys_data_faster(filename);
%                 here = strfind(info.header.electrode,'0 n');
%                 chan = str2double(info.header.electrode(here+3:end)) + 1; % +1 because electrodes start with 0
                

        end
        
    end
    
    
end

    

% 
%     
%         
%     
% try
%     load([localDir remoteFiles{n}],'-mat');
% catch me
%     fprintf(' Could not load! Check file. Skipping...\n');
%     return
% end
% 
% try PDS = rmfield(PDS,'initialParameters'); catch; end
% try PDS = rmfield(PDS,'initialParameterNames'); catch; end
% try PDS = rmfield(PDS,'initialParametersMerged'); catch; end
% try PDS = rmfield(PDS,'functionHandles'); catch; end
% try PDS = rmfield(PDS,'conditionNames'); catch; end
% 
% for t = 1:length(PDS.data) % loop over trials for this file
% 
%     if removeAnalogData % also removes analog-derived vars
%         try PDS.data{t}.datapixx = rmfield(PDS.data{t}.datapixx,'adc'); catch; end
%         try PDS.data{t}.behavior = rmfield(PDS.data{t}.behavior,'fixFP'); catch; end
%         try PDS.data{t}.behavior = rmfield(PDS.data{t}.behavior,'eyeXYs'); catch; end
%         try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'eyeXYs'); catch; end
%     else
%         try PDS.data{t}.behavior.eyeXYs = PDS.data{t}.stimulus.eyeXYs; catch; end
%         try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'eyeXYs'); catch; end
%     end
%     
%     if removeTimingData
%         try PDS.data{t} = rmfield(PDS.data{t},'timing'); catch; end
%     end
% 
%     if removeMotionTrackingData
%         try PDS.data{t} = rmfield(PDS.data{t},'imu'); catch; end
%         try PDS.data{t} = rmfield(PDS.data{t},'nexonar'); catch; end
%         try PDS.data{t} = rmfield(PDS.data{t},'mp'); catch; end
%     end
%     
%     if removeDotPositionData
%         try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'dotX'); catch; end
%         try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'dotY'); catch; end
%         try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'dotZ'); catch; end
%         try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'dotPos'); catch; end        
%         try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'dotSize'); catch; end
%     end
%     
%     try PDS.data{t}.behavior = rmfield(PDS.data{t}.behavior,'reward'); catch; end
%     try PDS.data{t} = rmfield(PDS.data{t},'mouse'); catch; end
%     try PDS.data{t}.behavior.goodtrial = PDS.data{t}.pldaps.goodtrial; catch;  end
%     try PDS.data{t} = rmfield(PDS.data{t},'pldaps'); catch; end                               
%     try PDS.data{t} = rmfield(PDS.data{t},'plexon'); catch; end
%     
%     try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'timeFpOn'); catch; end
%     try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'timeFpOff'); catch; end
%     try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'timeStimOn'); catch; end
%     try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'timeStimOff'); catch; end
%     try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'timeChoice'); catch; end
%     try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'frameFpOn'); catch; end
%     try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'frameFpOff'); catch; end
%     try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'frameTargetOn'); catch; end
%     try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'frameTargetOff'); catch; end
%     try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'frameFpOn'); catch; end
%     try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'frameFpOff'); catch; end
%     try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'frameStimOn'); catch; end
%     try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'frameStimOff'); catch; end
%     
%     try PDS.data{t}.behavior.timeTargDisappears = PDS.data{t}.task.timeTargDisappears; catch; end
%     try PDS.data{t}.behavior.probOfMemorySaccade = PDS.data{t}.task.probOfMemorySaccade; catch; end
%     try PDS.data{t} = rmfield(PDS.data{t},'task'); catch; end
%     try PDS.data{t} = rmfield(PDS.data{t},'postTarget'); catch; end                               
%     try PDS.data{t} = rmfield(PDS.data{t},'flagNextTrial'); catch; end                               
%     try PDS.data{t} = rmfield(PDS.data{t},'iFrame'); catch; end                               
%     try PDS.data{t} = rmfield(PDS.data{t},'state'); catch; end                               
%     try PDS.data{t} = rmfield(PDS.data{t},'keyboard'); catch; end                               
%     try PDS.data{t} = rmfield(PDS.data{t},'currentFrameState'); catch; end                               
%     try PDS.data{t} = rmfield(PDS.data{t},'eyeX'); catch; end                               
%     try PDS.data{t} = rmfield(PDS.data{t},'eyeY'); catch; end                               
%     try PDS.data{t} = rmfield(PDS.data{t},'pupil'); catch; end                               
%     try PDS.data{t} = rmfield(PDS.data{t},'remainingFrameTime'); catch; end                               
%     try PDS.data{t} = rmfield(PDS.data{t},'framePreLastDrawIdleCount'); catch; end                               
%     try PDS.data{t} = rmfield(PDS.data{t},'framePostLastDrawIdleCount'); catch; end                               
% 
%     % rename fields (dots3DMP)
%     try PDS.conditions{t}.stimulus = renameStructField(PDS.conditions{t}.stimulus,'conflictAngle','delta'); catch; end
%     try PDS.conditions{t}.stimulus = renameStructField(PDS.conditions{t}.stimulus,'stimCondition','modality'); catch; end
%     try
%         if PDS.conditions{t}.stimulus.modality==1 % don't need a bunch of nans for vestib trials, make it just one nan
%             PDS.data{t}.stimulus.dotX = NaN;
%             PDS.data{t}.stimulus.dotY = NaN;
%             PDS.data{t}.stimulus.dotZ = NaN;
%             PDS.data{t}.stimulus.dotSize = NaN;
%         end
%     catch    
%     end
% 
% end
% 
% save([localDir remoteFiles{n}],'PDS','-v7.3');
% fprintf(' done\n');

