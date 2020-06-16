function [output_op_info]=find_noise_free_electrodes(datafile,subject_op_info)
    if ~isempty(subject_op_info)
    %ecog_ch=subject_op_info.op_info.ecog_ch;
    ref_ch=subject_op_info.op_info.Ref;
    gnd_ch=subject_op_info.op_info.GND;
    bad_ch=subject_op_info.op_info.bad_channels;
    else 
    ecog_ch=[];
    ref_ch=[];
    gnd_ch=[];
    bad_ch=[];
    end 

ecog.param.filter_car  = 1;     % 0 == off % 1 == on
ecog.param.filter_type = 'IIR'; % IIR | FIR                             


% --- highpass filter --- 
ecog.param.highpass.Wp = 0.50; % Hz
ecog.param.highpass.Ws = 0.05; % Hz
ecog.param.highpass.Rp = 3;    % dB
ecog.param.highpass.Rs = 30;   % dB

% --- lowpass filter --- 
ecog.param.lowpass.Wp = 100;     % Hz
ecog.param.lowpass.Ws = 110;     % Hz
ecog.param.lowpass.Rp = 3;     % dB
ecog.param.lowpass.Rs = 30;    % dB

%---- downsampling rate ----

param=ecog.param;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CHECK DATA FILE DEPENDENCIES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 fprintf(1, '> Checking data file dependencies \n');
% 
 % go through all data files
for idx=1:length(datafile)
    
    % load just the header of this data file
    [ ~, ~, parameter ] = load_bcidat(datafile{idx},[0 0]); %#ok<*SAGROW>.
    parameters{idx}=parameter;

end


% % go through all combinations of data files
% for idxA=1:length(parameters),
%     for idxB=1:length(parameters),
%         
%         % compare these parameters between all files
%         dependency.SamplingRate        =  dependency.SamplingRate       &      (parameters{idxA}.SamplingRate.NumericValue   == parameters{idxB}.SamplingRate.NumericValue);
%         dependency.SourceChList        =  dependency.SourceChList       & (mean(parameters{idxA}.SourceChList.NumericValue   == parameters{idxB}.SourceChList.NumericValue)   == 1);
%        
%     end
% end
sampling_rates=cellfun(@(x) x.SamplingRate.NumericValue,parameters);
SourceChList=cellfun(@(x) x.SourceChList.NumericValue,parameters,'UniformOutput',false);
SourceCh=cell2mat(cellfun(@(x) x.SourceCh.NumericValue,parameters,'UniformOutput',false));
transmit_channels=cellfun(@(x) x.TransmitChList.NumericValue,parameters,'UniformOutput',false);


if length(unique(sampling_rates))~=1
    warning('error: Data files dont have the same SamplingRate!'); %#ok<WNTAG>
else 
    sampling_rate=unique(sampling_rates);
end

if length(unique(cellfun(@length,SourceChList)))~=1 
    warning('error: Data files dont have the same SourceChList!'); %#ok<WNTAG>
end

if length(unique(SourceCh))~=1 
    error('error: Data files dont have the same SourceNumber!'); %#ok<WNTAG>
end

if length(unique(cellfun(@length,SourceChList)))~=1    || length(unique(sampling_rates))~=1,
   error('error: Data files are too different to be analyzed together!'); 
end
if length(unique(cellfun(@length,transmit_channels)))==1 
    
    param.channels = transmit_channels{1};
end
if unique(SourceCh)~=length(param.channels)
    warning('error: channel ids and number dont match up! changing the param.channels to reflect SourceCh'); 
    param.channels=[1:(unique(SourceCh))]';
%    keyboard
end 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE HIGH-PASS FILTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define passband, stopband and attenuation
highpass{1}.Wp = param.highpass.Wp/(sampling_rate/2); 
highpass{1}.Ws = param.highpass.Ws/(sampling_rate/2);
highpass{1}.Rp = param.highpass.Rp; 
highpass{1}.Rs = param.highpass.Rs;
% calculate the minimum filter order for butterworth filter
[highpass{1}.n,highpass{1}.Wn] = buttord(highpass{1}.Wp,highpass{1}.Ws,highpass{1}.Rp,highpass{1}.Rs);
highpass{1}.n = highpass{1}.n + rem(highpass{1}.n,2);
% caclulate the filter coefficients in Zero-Pole-Gain design
[highpass{1}.z,highpass{1}.p,highpass{1}.k] = butter(highpass{1}.n,highpass{1}.Wn,'high');
[highpass{1}.sos,highpass{1}.g]=zp2sos(highpass{1}.z,highpass{1}.p,highpass{1}.k);
highpass{1}.h=dfilt.df2sos(highpass{1}.sos,highpass{1}.g);





%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE LINE-NOISE FILTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
peak.fcenter = 60;
peak.bw      = 0.001;

% calculate the IIR-peak filter coefficients in a,b format 
peak.wo = peak.fcenter/(sampling_rate/2);  
peak.bw = peak.bw;
[peak.b,peak.a] = iirpeak(peak.wo,peak.bw);  

% define the harmonics of line noise frequency
param.filter.notch.fcenter = [60,120,180,240];
param.filter.notch.bw      = ones(1,length(param.filter.notch.fcenter)).*0.001;

% calculate the IIR-peak filter coefficients in a,b format 
for idx = 1:length(param.filter.notch.fcenter),
    notch{idx}.wo = param.filter.notch.fcenter(idx)/(sampling_rate/2);  
    notch{idx}.bw = param.filter.notch.bw(idx);
    [notch{idx}.b,notch{idx}.a] = iirnotch(notch{idx}.wo,notch{idx}.bw);  
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD AND CONCATENATE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal = single([]);
session_length=[];
% go through all data files
for idx_file = 1:length(datafile)

    fprintf(1, '> Processing data file %d/%d \n', idx_file,length(datafile));
    
    %% load data file
    fprintf(1, '>> Loading data file %s \n',datafile{idx_file});

    [ signal_loop, states_loop, parameters_loop ] = load_bcidat(datafile{idx_file});
    signal_loop=signal_loop(:,param.channels);
    signal_broadband=signal_loop;
    session_length=[session_length;size(signal_loop,1)];
    %% highpass filter
    fprintf(1, '>> Highpass filtering signal \n');

    fprintf(1,'[');
    parfor idx_channel=1:size(signal_loop,2) 
        warning('off', 'signal:filtfilt:ParseSOS');
        signal_loop(:,idx_channel) = single(filtfilt(highpass{1}.sos,highpass{1}.g,double(signal_loop(:,idx_channel)))); %#ok<PFBNS>
        fprintf(1,'.');    
    end
    fprintf(1,'] done\n');
    %
    
    buffer_before = sampling_rate * 1; % 1 sec 
    buffer_after  = sampling_rate * 2; % 2 sec
    
    
    
    % concatenate states and signals from data files
    signal              = [signal;signal_loop]; %#ok<AGROW>
    clear signal_loop
end
%% plot the dataset throughout 
session_length=[1;session_length];
session_start=cumsum(session_length);
% get sections 
t_length=25e4;
kk_total=ceil(size(signal,1)/t_length);
x_cell=mat2cell(signal',ones(1,size(signal,2)));
min_max=mean(cell2mat(cellfun(@(y) [min(y),max(y)],x_cell,'UniformOutput',false)));
x_norm_cell=cellfun(@(x) (x-min_max(1))./(min_max(2)-min_max(1)),x_cell,'UniformOutput',false);
% 
figure(1);
%colors=inferno(size(x_norm_cell,1));
col_inf=inferno(floor(.8*size(x_norm_cell,1)));
col_vir=viridis(floor(.8*size(x_norm_cell,1)));
colors=[col_vir(1:floor(size(x_norm_cell,1)/2),:);col_inf(1:(floor(size(x_norm_cell,1)/2)+1),:)];
for kk=1:kk_total
    clf;
    set(gcf,'color',[.7,.7,.7],'position',[17,1,2454,1337])
    ax=axes('position',[.02,.02,.95,.95]);
    hold on
    session_end=min([(kk*t_length),size(signal,1)]);
    t_window=((kk-1)*t_length+1):2:(session_end);
    
    arrayfun(@(x) plot(t_window,x_norm_cell{x}(t_window)+x-.5,'color',colors(x,:)),[1:size(x_norm_cell,1)])
    
    set(ax,'ytick',[1:size(x_norm_cell,1)]);
    set(ax,'yticklabel','')
    arrayfun(@(x) text(t_window(1),x,num2str(x),'Color',colors(x,:),'HorizontalAlignment','right','VerticalAlignment','middle'),[1:size(x_norm_cell,1)]);
    set(ax,'ylim',[0,size(x_norm_cell,1)+1])
    if any(ismember(t_window,session_start))
        arrayfun(@(x) plot([t_window(x),t_window(x)],[0,size(x_norm_cell,1)+1],'k-','LineWidth',2),find(ismember(t_window,session_start)))
        session_name=num2str(find(ismember(session_start,t_window)));
        arrayfun(@(x) text(t_window(x),size(x_norm_cell,1)+1,['sess: ',session_name],'VerticalAlignment','bottom','fontsize',20),find(ismember(t_window,session_start)))
    end 
    shg;
    title({'pre common source averaging',[num2str(kk),' of ', num2str(kk_total)]})
    pause
end 
close all; 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEASSURE LINE-NOISE POWER BEFORE SIGNAL PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, '> Meassuring 60 Hz noise power before signal processing \n');
fprintf(1,'[');
parfor idx_channel=1:size(signal,2),
    % calculate average root-mean-square of the line-noise
    signal_noise_before(idx_channel) = mean(sqrt(filtfilt(peak.b,peak.a,double(signal(:,idx_channel))).^2),1); %#ok<PFBNS>
    fprintf(1,'.');    
end
fprintf(1,'] done\n');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND DATA CHANNELS WITH SIGNIFICANT LINE-NOISE POWER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1, '> Searching for channels with significant 60 Hz noise \n');

fprintf(1,'[');
if isfield(parameters{1},'DeviceIDMaster')
    signal_noise=zeros(size(signal,2),1);
    % determine the number of g.USBamps where each amp has 16 channels
    num_amps = round(size(signal,2) / 16);

    % for each of these amps
    for idx_amp = 1:num_amps,
        
        % the the set of channels
        idx_low  = (idx_amp-1)*16+1;
        idx_high = (idx_amp-0)*16+0;

        % get the channels that actually had signals 
        list_channels = intersect(idx_low:idx_high,param.channels);
        % check if any channels are left and 
        if ~isempty(list_channels),
            % calculate the common average reference signal 
            signal_mean = mean(signal(:,list_channels),2);
        else
            % if no channels is left then the common average reference signal is zero
            signal_mean = zeros(size(signal,1),1);
        end
            
        % for each channel on this amp
        for idx_ch=1:length(list_channels),
            % subtract the common average signal from each channel of this amp          
            signal_preliminary = double(signal(:,list_channels(idx_ch)));% - double(signal_mean);
            
            % calculate the residual line-noise root-mean-square power
            v=mean(sqrt(filtfilt(peak.b,peak.a,signal_preliminary).^2),1);
            signal_noise(list_channels(idx_ch)) = v;
            fprintf(1,'.');
            
        end
    end  
else
    warning('error: common average reference filtering is only supported for g.USBamp data');  %#ok<WNTAG>
end

clear signal_mean
fprintf(1,'] done\n');
%%
% find those channels for which the line-noise power is 1.5 standard deviations higher than the average 

param.channels_noise = param.channels( signal_noise > (mean(signal_noise)+1.5*std(signal_noise)));
%param.channels_deselect = sort(unique([param.channels_noise param.channel_ground param.channel_ref param.channel_seizure]));

    
param.channels_deselect = sort(unique([reshape(param.channels_noise,1,[]),unique([ref_ch,gnd_ch,bad_ch]) ]));
param.channels_selected = setdiff(param.channels,param.channels_deselect);

fprintf(1, '> Not using %d channels for common source averaging: ',length(param.channels_deselect));
fprintf(1, '%d ',param.channels_deselect);
fprintf(1, '\n');

%% 
fprintf(1, '> Common average filtering signal \n');

% now calculate the common-average reference signal again, without those channels that have significant line-noise
fprintf(1,'[');
if isfield(parameters{1},'DeviceIDMaster') && param.filter_car > 0,
    
    % determine the number of g.USBamps where each amp has 16 channels
    num_amps = round(size(signal,2) / 16);

    % for each of these amps
    for idx_amp = 1:num_amps,
        idx_low  = (idx_amp-1)*16+1;
        idx_high = (idx_amp-0)*16+0;

        % exclude the channels that had signifiant line-noise
        list_channels = intersect(idx_low:idx_high,param.channels_selected);
        % check if any channels are left and 
        if ~isempty(list_channels) && length(list_channels)>1 ,   
          
            % calculate the common average reference signal 
            signal_mean = mean(signal(:,intersect(idx_low:idx_high,param.channels_selected)),2);
        
            % subtract the common average signal from each channel of this amp
            for idx_ch=idx_low:idx_high,
                signal(:,idx_ch) = signal(:,idx_ch) - signal_mean;
                fprintf(1,'.');
            end
        elseif ~isempty(list_channels) && length(list_channels)==1 ,
            param.channels_selected(param.channels_selected==list_channels)=[];
            param.channels_deselect=[param.channels_deselect,list_channels];
            fprintf('remove channel with no additional pair for common reference averaging')
        else
            % if no channel is left then the signal is not filtered
            for idx_ch=idx_low:idx_high,
                signal(:,idx_ch) = signal(:,idx_ch);
                fprintf(1,'.');
            end
        end
    end  
end 
fprintf(1,'] done\n');
%% 
fprintf(1, '> Meassuring 60 Hz noise power after signal processing \n');

fprintf(1,'[');
% for each channels calculate the root-mean-square line-noise power
parfor idx_channel=1:size(signal,2),
    signal_noise_after(idx_channel) = mean(sqrt(filtfilt(peak.b,peak.a,double(signal(:,idx_channel))).^2),1); %#ok<PFBNS>
    fprintf(1,'.');    
end
fprintf(1,'] done\n');
%% 
fprintf(1, '> Reduced 60 Hz noise from %.2f to %.2f uV',mean(signal_noise_before(param.channels_selected)),mean(signal_noise_after(param.channels_selected)));
fprintf(1, '\n');


%% 
fprintf(1, '> Notch filtering signal \n');
fprintf(1,'[');
% for each channel
parfor idx_channel=1:size(signal,2),
    
    % get the signal for this channel
    signal_preliminary = double(signal(:,idx_channel));
    
    % remove all harmonics of line-noise
    for idx = 1:length(param.filter.notch.fcenter), %#ok<PFBNS>
        signal_preliminary = filtfilt(notch{idx}.b,notch{idx}.a,signal_preliminary); %#ok<PFBNS>
    end 
    
    % return the signal
    signal(:,idx_channel) = single(signal_preliminary);
    
    fprintf(1,'.');
end
fprintf(1,'] done\n');
%% plot post filtering 
t_length=25e4;
kk_total=ceil(size(signal,1)/t_length);
x_cell=mat2cell(signal',ones(1,size(signal,2)));
valid_channels=ones(1,size(signal,2));
valid_channels(param.channels_deselect)=nan;
x_cell=arrayfun(@(x) x_cell{x}*valid_channels(x),1:size(signal,2),'UniformOutput',false);
x_cell=x_cell';
min_max=nanmean(cell2mat(cellfun(@(y) [min(y),max(y)],x_cell,'UniformOutput',false)));
x_norm_cell=cellfun(@(x) (x-min_max(1))./(.5*(min_max(2)-min_max(1))),x_cell,'UniformOutput',false);
% 
figure(1);

col_inf=inferno(floor(.8*size(x_norm_cell,1)));
col_vir=viridis(floor(.8*size(x_norm_cell,1)));
colors=[col_vir(1:floor(size(x_norm_cell,1)/2),:);col_inf(1:(floor(size(x_norm_cell,1)/2)+1),:)];
for kk=1:kk_total
    clf;
    set(gcf,'color',[.7,.7,.7],'position',[17,1,2454,1337])
    ax=axes('position',[.02,.02,.95,.95]);
    hold on
    session_end=min([(kk*t_length),size(signal,1)]);
    t_window=((kk-1)*t_length+1):2:(session_end);
    arrayfun(@(x) plot(t_window,x_norm_cell{x}(t_window)+x-1,'color',colors(x,:)),param.channels)
    
    set(ax,'ytick',[1:size(x_norm_cell,1)]);
    set(ax,'yticklabel','')
    arrayfun(@(x) text(t_window(1),x,num2str(x),'Color',colors(x,:),'HorizontalAlignment','right','VerticalAlignment','middle'),param.channels);
    set(ax,'ylim',[0,size(x_norm_cell,1)+1]);
    if any(ismember(t_window,session_start))
        arrayfun(@(x) plot([t_window(x),t_window(x)],[0,size(x_norm_cell,1)+1],'k-','LineWidth',2),find(ismember(t_window,session_start)))
        session_name=num2str(find(ismember(session_start,t_window)));
        arrayfun(@(x) text(t_window(x),size(x_norm_cell,1)+1,['sess: ',session_name],'VerticalAlignment','bottom','fontsize',20),find(ismember(t_window,session_start)))
    end 
    shg;
    title({'post common source averaging',[num2str(kk),' of ', num2str(kk_total)]})
    pause
end 
close all; 
%% 
prompt='additional channels to remove? format :[1,2] ';
x=input(prompt);
if ~isempty(x)
    param.channels_deselect=sort([param.channels_deselect,x]);
    param.channels_selected = setdiff(param.channels,param.channels_deselect);
    t_length=25e4;
    kk_total=ceil(size(signal,1)/t_length);
    x_cell=mat2cell(signal',ones(1,size(signal,2)));
    valid_channels=ones(1,size(signal,2));
    valid_channels(param.channels_deselect)=nan;
    x_cell=arrayfun(@(x) x_cell{x}*valid_channels(x),1:size(signal,2),'UniformOutput',false);
    x_cell=x_cell';
    min_max=nanmean(cell2mat(cellfun(@(y) [min(y),max(y)],x_cell,'UniformOutput',false)));
    x_norm_cell=cellfun(@(x) (x-min_max(1))./(.5*(min_max(2)-min_max(1))),x_cell,'UniformOutput',false);
    %
    figure(1);
    
    col_inf=inferno(floor(.8*size(x_norm_cell,1)));
    col_vir=viridis(floor(.8*size(x_norm_cell,1)));
    colors=[col_vir(1:floor(size(x_norm_cell,1)/2),:);col_inf(1:(floor(size(x_norm_cell,1)/2)+1),:)];
    for kk=1:kk_total
        clf;
        set(gcf,'color',[.7,.7,.7],'position',[17,1,2454,1337])
        ax=axes('position',[.02,.02,.95,.95]);
        hold on
        session_end=min([(kk*t_length),size(signal,1)]);
        t_window=((kk-1)*t_length+1):2:(session_end);
        arrayfun(@(x) plot(t_window,x_norm_cell{x}(t_window)+x-1,'color',colors(x,:)),param.channels)
        
        set(ax,'ytick',[1:size(x_norm_cell,1)]);
        set(ax,'yticklabel','')
        arrayfun(@(x) text(t_window(1),x,num2str(x),'Color',colors(x,:),'HorizontalAlignment','right','VerticalAlignment','middle'),param.channels);
        set(ax,'ylim',[0,size(x_norm_cell,1)+1]);
        if any(ismember(t_window,session_start))
        arrayfun(@(x) plot([t_window(x),t_window(x)],[0,size(x_norm_cell,1)+1],'k-','LineWidth',2),find(ismember(t_window,session_start)))
        session_name=num2str(find(ismember(session_start,t_window)));
        arrayfun(@(x) text(t_window(x),size(x_norm_cell,1)+1,['sess: ',session_name],'VerticalAlignment','bottom','fontsize',20),find(ismember(t_window,session_start)))
        end 
        shg;
        title({'post common source averaging',[num2str(kk),' of ', num2str(kk_total)]})
        pause
    end
    close all;
end

%% add the results to subject_op_info
subject_op_info.op_info.channel_removed_by_user=x;
subject_op_info.op_info.transmit_chan=param.channels;
subject_op_info.op_info.unselected_channels=param.channels_deselect;
subject_op_info.op_info.clean_channels=param.channels_selected;
subject_op_info.op_info.sampling_rate=sampling_rate;
subject_op_info.op_info.channel_noise_across_all_sess=signal_noise_before;
subject_op_info.op_info.channel_denoise_across_all_sess=signal_noise_after;
output_op_info=subject_op_info;


end 
