% selection of language specific electrodes by comparing electrode response
% of the pretrial of a condition and a target word index in that condition
% or the target word index of two different conditions

%for each electrode the script, 
%1. gets the response amplitudes to condition1 (or pretrial)
%2. get the response amplitudes to condition2 (or pretrial)
%3. conducts a two sample t-test comparing the two conditions (or pretrials) 
%if significant, adds the electrode to language responsive electrodes

%% clean up workspace
clear all
close all

%% specify parameters
% subject_id = 'AMC096';
% experiment_name = 'MITNLengthSentences';
% condition1 = [{'3sents_8words_intact'},{'6sents_4words_intact'}, {'1sent_24words_intact'}];
% %condition2 = [{'3sents_8words_intact'},{'6sents_4words_intact'}, {'1sent_24words_intact'}];
% %condition1 = [{'3sents_8words_scrambled'},{'6sents_4words_scrambled'}, {'1sent_24words_scrambled'}];
% condition2 = [{'3sents_8words_scrambled'},{'6sents_4words_scrambled'}, {'1sent_24words_scrambled'}];
% cond1_use_pretrial = true;
% cond2_use_pretrial = true;
% cond1_target_word_idx = 0;
% cond2_target_word_idx = 8;
% data_to_use = "hilbert_zs"; %options: hilbert, hilbert_zs;

subject_id = 'AMC096'; %must be string
experiment_name = 'MITLangloc'; %must be string
condition1 = [{'Sentences'}]; %must be a list
condition2 = [{'Jabberwocky'}]; %must be a list
cond1_use_pretrial = false; %must be boolean
cond2_use_pretrial = false; %must be boolean
cond1_target_word_idx = [1:12]; %word(s) to get the average response amplitude to, must be a list
cond2_target_word_idx = [1:12]; %word(s) to get the average response amplitude to, must be a list
data_to_use = "hilbert_zs"; %must be string
verbose = true;


num_permutations = 1000;
p_threshold = 0.01;
%TODO- convert parameters to a struct? pass the struct through
%would get rid of global verbose variable


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

%% extract subject's data for the given experiment
d_data= dir(strcat(crunched_data_path,filesep,subject_id,'*_crunched.mat'));
fprintf(' %d .mat files were found \n', length(d_data));
d_data=arrayfun(@(x) strcat(d_data(x).folder,filesep,d_data(x).name),[1:length(d_data)]','uni',false);

if(verbose)fprintf("Using subject %s's %s data from %s for language electrode selection \n", subject_id,data_to_use, experiment_name);end


[cond1] = extract_subj_data(d_data, ...
                            data_to_use,... 
                            condition1,cond1_use_pretrial,cond1_target_word_idx,...
                            verbose);
[cond2] = extract_subj_data(d_data, ...
                            data_to_use,... 
                            condition2,cond2_use_pretrial,cond2_target_word_idx,...
                            verbose);
%% get language electrodes
[lang_electrodes_list, electrodes] = determine_lang_electrodes(cond1,cond2,verbose);


%%
function [list, lang_electrodes] = determine_lang_electrodes(cond1,cond2,verbose)
%each row of cond1 and cond2 contain one electrode's avg response
%amplitudes over all the relevant trials

%t-test comparing the rows of cond1 and cond2
num_electrodes = length(cond1);
list = [];
if(verbose)fprintf("Conducting t-tests comparing condition1 and condition2 responses for %d electrodes\n", num_electrodes);end
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

function [cond_trial_amplitudes] = extract_subj_data(d_data,data_to_use, condition, use_pretrial, target_word_idx,verbose)

cond_amplitudes = [];
for k=1:length(d_data)
    if(verbose)fprintf('Session %d \n', k);end
    subj=load(d_data{k});
    amplitudes = extract_condition_amplitudes(subj,data_to_use, condition, use_pretrial,target_word_idx, verbose);
    cond_amplitudes = cat(1,cond_amplitudes, amplitudes);
end

%put data into easy to use format
cond_word_amplitudes = cell2mat(reshape(cond_amplitudes,1,[]));

num_words = length(target_word_idx);
num_trials = size(cond_word_amplitudes,2)/num_words;

%compute the mean across the word positions in each trial. 
cond_trial_amplitudes = zeros(length(cond_word_amplitudes),num_trials);

for i=1:length(cond_word_amplitudes)
    for j=1:num_trials
        start_idx=1;
        if(j>1)
            start_idx = start_idx + (j-1)*num_words;
        end
        end_idx = start_idx+num_words-1;
        cond_trial_amplitudes(i,j) = mean(cond_word_amplitudes(i, start_idx:end_idx));
    end
end

end
function [avg_amplitudes] = extract_condition_amplitudes(subj,data_to_use, conditions,use_pretrial, target_word_idx,verbose)
avg_amplitudes = [];
avg_word_amplitudes = [];
if(verbose)fprintf('extracting average response amplitudes to... \n');end
for i=1:length(conditions)
    for j=1:length(target_word_idx)
        extracted_amplitudes = get_avg_amplitudes(subj,data_to_use, conditions{i}, use_pretrial, target_word_idx(j), verbose);
        avg_word_amplitudes = cat(1,avg_word_amplitudes, extracted_amplitudes);
    end
    avg_amplitudes = cat(1,avg_amplitudes, avg_word_amplitudes);
end
end
% subj: struct, a subj's crunched materials, containing info and data structs
% data_to_use: string, which type of data from the data struct to use
% cond: string, the condition to extract avg amplitudes from
% use_pretrial: boolean, specifies whether you should extract the amplitude from the pretrials in the relevant trials
% target_word_idx: int, only used if use_pretrial is false, specifies the word to collect avg amplitude from
function [avg_amplitudes] = get_avg_amplitudes(subj, data_to_use, cond, use_pretrial,target_word_idx, verbose)
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
    if(verbose)fprintf('pretrial before %s\n',cond);end
    ave_electrodes=cellfun(@(x) x.signal_ave_pre_trial_hilbert_zs_downsample, relevant_trial_data,'UniformOutput',false);
else
    %get the average response amplitude to the target word for all relevant
    %trials
    data_name = strcat("signal_ave_", data_to_use, "_downsample_parsed");
    if(verbose)fprintf('word %d in %s\n', target_word_idx, cond);end
    ave_electrodes=cellfun(@(x) x.signal_ave_hilbert_zs_downsample_parsed{target_word_idx,1}, relevant_trial_data,'UniformOutput',false);
end

avg_amplitudes = ave_electrodes;
end