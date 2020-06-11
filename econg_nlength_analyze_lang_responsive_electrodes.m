%% clean up workspace
clear all
close all
home 
%% specify where the data is
main_path='~/MyData/ecog_nlength/';
data_path='~/MyData/ecog_nlength/crunched/';
analysis_path=strcat(main_path,'analysis/distributed_oscilatory_power/');
subject_id='AMC026';
d_data= dir(strcat(data_path,'/',subject_id,'*_crunched_v2.mat'));
fprintf(' %d .mat files were found \n', length(d_sn));
d_data=arrayfun(@(x) strcat(d_data(x).folder,'/',d_data(x).name),[1:length(d_data)]','uni',false);
%% combine responses from all sessions for a given condition. (S= sentence,....)
sent_word_hilb_ave_tensor_all=[];
for k=1:length(d_data)
    subj=load(d_data{k});
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    % sentences 
    cond='S';
    trial_type=info.word_type;
    cond_id=find(cellfun(@(x) x==cond,trial_type));
    data_cond=data(cond_id );
    % stimuli to find in the condition trial 
    stim='word';
    stim_loc=cellfun(@(x) find(strcmp(x.stimuli_type,stim)), data_cond,'uni',false);
    hilb_zs_ave_cell=cellfun(@(x) x.signal_ave_hilbert_zs_downsample_parsed, data_cond,'uni',false);
    stim_zs_ave_cell=arrayfun(@(x) hilb_zs_ave_cell{x}(stim_loc{x}),[1:size(hilb_zs_ave_cell,1)],'uni',false );
    % reshape stim to correct form if it is not 
    stim_zs_ave_cell=cellfun(@(x) cell2mat(reshape(x,1,[])),stim_zs_ave_cell,'uni',false );
    % tensor format (elec*stim(words)*trial)
    stim_zs_ave_tensor=double(cell2mat(permute(stim_zs_ave_cell,[1,3,2])));
    % combine trials for the condition
    sent_word_hilb_ave_tensor_all=cat(3,sent_word_hilb_ave_tensor_all,stim_zs_ave_tensor);
    fprintf('added %s\n', d_data{k})
end 