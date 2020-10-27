% selection of language specific electrodes by comparing electrode response
% during pretrial to the electrode response 4 words into any sentence
% stimulus

%requirements: a pretrial duration, stimuli of more than 4 words per trial

%for each electrode, 
%1. get the response amplitudes of all of the pretrials
%2. get the response amplitudes to the fourth word of the non-scrambled
%trials
%3. conduct a t-test comparing the pretrial amplitudes to the fourth word
%amplitudes -- if significant, add this electrode to language responsive
%electrodes

%% clean up workspace
clear all
close all

%% parameters
subject_id = 'AMC096';
experiment_name = 'MITNLengthSentences';
condition1 = [{'3sents_8words_intact'},{'6sents_4words_intact'}, {'1sent_24words_intact'}];
condition2 = [{'3sents_8words_intact'},{'6sents_4words_intact'}, {'1sent_24words_intact'}];
%condition2 = [{'3sents_8words_scrambled'},{'6sents_4words_scrambled'}, {'1sent_24words_scrambled'}];
cond1_use_pretrial = false;
cond2_use_pretrial = true;
cond1_target_word_idx = 8;
cond2_target_word_idx = 8;
data_to_use = "hilbert_zs"; %options: hilbert, hilbert_zs;

global print_statements
print_statements = true;

%% specify where the data is
addpath(genpath('~/GitHub/evlab_ecog_tools'));
addpath(genpath('~/GitHub/evlab_ecog_tools/ecog-filters/'));
addpath(genpath('~/GitHub/evlab_ecog_tools/albany_mex_files'));
addpath(genpath('~/GitHub/evlab_matlab_tools/Colormaps'));
code_path='~/GitHub/evlab_ecog_tools';    
ecog_path = '~/Desktop/ECOG';

crunched_data_path = [ecog_path filesep 'crunched' filesep experiment_name filesep]; %save it into an experiment specific folder

%path specifically for expt sub_op_info
expt_sub_op_info_path = [crunched_data_path 'sub_op_info_' experiment_name filesep];

%% load the subject's crunched data for the given experiment
d_data= dir(strcat(crunched_data_path,filesep,subject_id,'*_crunched.mat'));
fprintf(' %d .mat files were found \n', length(d_data));
d_data=arrayfun(@(x) strcat(d_data(x).folder,filesep,d_data(x).name),[1:length(d_data)]','uni',false);

if(print_statements)fprintf("Using subject %s's %s data from %s for language electrode selection \n", subject_id,data_to_use, experiment_name);end

[cond1,cond2] = extract_subj_data(d_data, ...
                                  data_to_use,... 
                                  condition1,condition2,...
                                  cond1_use_pretrial, cond2_use_pretrial,...
                                  cond1_target_word_idx, cond2_target_word_idx);
%% get language electrodes
[lang_electrodes_list, electrodes] = determine_lang_electrodes(cond1,cond2);




%%
function [list, lang_electrodes] = determine_lang_electrodes(cond1,cond2)
global print_statements
%each row of cond1 and cond2 contain one electrode's avg response
%amplitudes over all the relevant trials

%t-test comparing the rows of cond1 and cond2
num_electrodes = length(cond1);
list = [];
if(print_statements)fprintf("Conducting t-tests comparing condition1 and condition2 responses for %d electrodes\n", num_electrodes);end
t_test_results = [1:num_electrodes];
for i=1:num_electrodes
    cond1_vector = cond1(i,:);
    cond2_vector = cond2(i,:);
    result = ttest2(cond1_vector,cond2_vector);
    t_test_results(i) = result;
    if (result==1)
        list = [list i]; %save the electrode number in the list
    end
end

lang_electrodes = t_test_results;

end

%%
function [cond1_amplitudes, cond2_amplitudes] = extract_subj_data(d_data,data_to_use, condition1, condition2, cond1_use_pretrial, cond2_use_pretrial, cond1_target_word_idx, cond2_target_word_idx)
global print_statements
cond1_amplitudes = [];
cond2_amplitudes = [];
for k=1:length(d_data)
    if(print_statements)fprintf('Session %d \n', k);end
    subj=load(d_data{k});
    
    amplitudes1 = extract_condition_amplitudes(subj,data_to_use, condition1, cond1_use_pretrial,cond1_target_word_idx);
    cond1_amplitudes = cat(1,cond1_amplitudes, amplitudes1);
    amplitudes2 = extract_condition_amplitudes(subj,data_to_use, condition2, cond2_use_pretrial, cond2_target_word_idx);
    cond2_amplitudes = cat(1,cond2_amplitudes, amplitudes2);
end

%put data into easy to use format
cond1_amplitudes = cell2mat(reshape(cond1_amplitudes,1,[]));
cond2_amplitudes = cell2mat(reshape(cond2_amplitudes,1,[]));


end

function [avg_amplitudes] = extract_condition_amplitudes(subj,data_to_use, conditions,use_pretrial, target_word_idx)
global print_statements
avg_amplitudes = [];
if(print_statements)fprintf('extracting average response amplitudes to... \n');end
for i=1:length(conditions)
    extracted_amplitudes = get_avg_amplitudes(subj,data_to_use, conditions{i}, use_pretrial, target_word_idx);
    avg_amplitudes = cat(1,avg_amplitudes, extracted_amplitudes);
end
end

%% extract relevant conditions
% subj: struct, a subj's crunched materials, containing info and data structs
% data_to_use: string, which type of data from the data struct to use
% cond: string, the condition to extract avg amplitudes from
% use_pretrial: boolean, specifies whether you should extract the amplitude from the pretrials in the relevant trials
% target_word_idx: int, only used if use_pretrial is false, specifies the word to collect avg amplitude from
function [avg_amplitudes] = get_avg_amplitudes(subj, data_to_use, cond, use_pretrial,target_word_idx)
global print_statements
subj_sess_id=fieldnames(subj);
subj=subj.(subj_sess_id{1}); %get rid of top layer of struct
info = subj.info;
data = subj.data;

%get strings of the condition names on each trial
trial_type_str = cellfun(@string, info.condition_name,'UniformOutput',false); 
%convert the logical arrays (one per word) into 1 dimension (one per trial)
trial_type_unique = cellfun(@unique,trial_type_str,'UniformOutput',false);
%get the indices of the relevant trials
relevant_trials = cellfun(@(x) x==cond, trial_type_unique);

relevant_trial_data=data(relevant_trials);
if(use_pretrial)
    %get the average pretrial response amplitude for all relevant trials
    data_name = strcat("signal_ave_pre_trial_",data_to_use, "_downsample");
    if(print_statements)fprintf('pretrial before %s from %s \n',cond,data_name);end
    ave_electrodes=cellfun(@(x) x.(data_name), relevant_trial_data,'UniformOutput',false);
else
    %get the average response amplitude to the target word for all relevant
    %trials
    data_name = strcat("signal_ave_", data_to_use, "_downsample_parsed");
    if(print_statements)fprintf('word %d in %s from  %s \n', target_word_idx, cond, data_name);end
    ave_electrodes=cellfun(@(x) x.(data_name){target_word_idx,1}, relevant_trial_data,'UniformOutput',false);
end

avg_amplitudes = ave_electrodes;
end




      