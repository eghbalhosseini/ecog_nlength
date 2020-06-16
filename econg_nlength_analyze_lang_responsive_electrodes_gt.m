%% clean up workspace
clear all
close all
home 
%% specify where the data is
% main_path='~/MyData/ecog_nlength/';
% data_path='~/MyData/ecog_nlength/crunched/';

addpath(genpath('\GitHub\evlab_ecog_tools\'));
addpath(genpath('\GitHub\evlab_ecog_tools\ecog-filters\'));
addpath(genpath('\GitHub\evlab_ecog_tools\albany_mex_files'));
addpath(genpath('\GitHub\evlab_matlab_tools\beeswarm\'));

data_path='C:\Users\greta\Dropbox (MIT)\ECoG_data\crunched\';
main_path='C:\Users\greta\Dropbox (MIT)\ECoG_data\crunched\';
analysis_path=strcat(main_path,'analysis/distributed_oscilatory_power/');
subject_id='AMC082';
d_data= dir(strcat(data_path,'/',subject_id,'*_crunched_v2.mat'));
fprintf(' %d .mat files were found \n', length(d_data));
d_data=arrayfun(@(x) strcat(d_data(x).folder,'/',d_data(x).name),[1:length(d_data)]','uni',false);

%% combine responses from all sessions for a given condition. (S= sentence,....)
s3w8i_word_hilb_ave_tensor_all=[];
s3w8s_word_hilb_ave_tensor_all=[];

for k=1:length(d_data)
    
    subj=load(d_data{k});
    
%     subj_sess_id=fieldnames(subj);
%     subj=subj.(subj_sess_id{1});
%     data=subj.data;
%     info=subj.info;
%     
%     trial_type=info.word_type;
%     trial_type_str = cellfun(@string,trial_type,'UniformOutput',false);
%     trial_type_u = cellfun(@unique,trial_type_str,'UniformOutput',false);
%     
%     if k == 1
%         fprintf('Conditions available: %s\n', string(trial_type_u))
%     end
    stim='word'

    cond='3sents_8words_intact';
    stim_intact_zs_ave_tensor = extract_cond(subj, cond, stim);
    
    cond='3sents_8words_scrambled';
    stim_scrambled_zs_ave_tensor = extract_cond(subj, cond, stim);
    
    % combine trials for the condition
    s3w8i_word_hilb_ave_tensor_all=cat(3,s3w8i_word_hilb_ave_tensor_all,stim_intact_zs_ave_tensor);
    s3w8s_word_hilb_ave_tensor_all=cat(3,s3w8s_word_hilb_ave_tensor_all,stim_scrambled_zs_ave_tensor);
    
    fprintf('added %s\n', d_data{k})
end 


%% code for plotting : (not fully tested here)

close all 
num_rows=5;
num_columns=2;
total_plots=num_rows*num_columns;
p=0;
for i=1:size(s3w8i_word_hilb_ave_tensor_all,1)
    if i > 6
        break
    end
    electrode_resp=squeeze(s3w8i_word_hilb_ave_tensor_all(i,:,:));
    word_pos=repmat(1:size(electrode_resp,1),size(electrode_resp,2),1)';
    x=word_pos(:);
    y=electrode_resp(:);
    f=figure(fix((i-1)/total_plots)+1);
    set(f,'position',[2937 -115 1266 1275])
%     sub_title=sprintf('subj: %s,\n electrode: %s ', info.subject ,strrep(info.ecog_channels_labels{ valid_chan_id(i)},' ','')); 
    sub_title=sprintf('subj: %s,\n electrode: %s ', subject_id, string(i));
    ax=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
    x_swarm=beeswarm(x,y,'sort_style','rand','corral_style','omit','dot_size',.1);
    hold off
    bl=scatter(x_swarm,y,10,[1,.2,.2],'filled');
    %bl.MarkerEdgeColor=[1,1,1];
    bl.MarkerEdgeAlpha=.5;bl.MarkerFaceAlpha=1;
    hAnnotation = get(bl,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
    hold on
    
    plot(1:size(electrode_resp,1),nanmean(electrode_resp,2),'k-','LineWidth',3)
    h1=plot(get(gca,'xlim'),0*get(gca,'ylim'),'k--');hAnnotation = get(h1,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
    ax.XAxis.Visible = 'on'; 
    ax.XTick=1:size(electrode_resp,1);
    y_quantile=quantile(y,90);
    set(ax,'ydir', 'normal','box','off','ylim',[y_quantile(1),y_quantile(end)]);
    title(sub_title);
    if ~mod(i,total_plots) | i==size(s3w8i_word_hilb_ave_tensor_all,1)
        %legend('show','Location','northeastoutside')
        p=p+1;
        if ~exist(strcat(analysis_path,info.subject))
            mkdir(strcat(analysis_path,info.subject))
        end 
        xlabel('word position');
        set(gcf,'PaperPosition',[.25 .25 8 6])
        %print(gcf,'-fillpage', '-dpdf', strcat(analysis_path,info.subject,'/',info.subject,'_electrode_res_',num2str(p),'_new.pdf')); 
        close(gcf)
    end 
    
end


%% 
function [stim_zs_ave_tensor] = extract_cond(subj, cond, stim)
    subj_sess_id=fieldnames(subj);
    subj=subj.(subj_sess_id{1});
    data=subj.data;
    info=subj.info;
    
    trial_type=info.word_type;
    trial_type_str = cellfun(@string,trial_type,'UniformOutput',false);
    trial_type_u = cellfun(@unique,trial_type_str,'UniformOutput',false);

    cond_id=find(cellfun(@(x) x==cond, trial_type_u));
    data_cond=data(cond_id);

    stim_loc=cellfun(@(x) find(strcmp(x.stimuli_type,stim)), data_cond,'uni',false);
    hilb_zs_ave_cell=cellfun(@(x) x.signal_ave_hilbert_zs_downsample_parsed, data_cond,'uni',false);
    stim_zs_ave_cell=arrayfun(@(x) hilb_zs_ave_cell{x}(stim_loc{x}),[1:size(hilb_zs_ave_cell,1)],'uni',false );
    % reshape stim to correct form if it is not 
    stim_zs_ave_cell=cellfun(@(x) cell2mat(reshape(x,1,[])),stim_zs_ave_cell,'uni',false );
    % tensor format (elec*stim(words)*trial)
    stim_zs_ave_tensor=double(cell2mat(permute(stim_zs_ave_cell,[1,3,2])));
end