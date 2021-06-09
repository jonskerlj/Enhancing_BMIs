%% Code used to preprocess the data for classification of 8 targets - CO
% Task
clear all; close all;
%% Load the data
fpath = '../raw_data/';
fname = '._MM_S1_processed.mat';

load([fpath fname])
td = trial_data;
%% Remove unsuccesful trials
N = size(trial_data,2);
idx = [];
for i = 1:N
    if trial_data(i).result == 'I' || trial_data(i).result == 'F'
        idx = [idx, i];         
    end
end
td(idx) = [];

%% Split train and test data
M =740;
td_train = td(1:M);
%td_test = td(M+1:end);

%%
y = [];
for row = 1:M
    angle = determine_angle(td_train(row));
    M1_trial = td_train(row).M1_spikes';
    PMd_trial = td_train(row).PMd_spikes';
        
    M1_spikes(:,row) = sum(M1_trial,2);
    PMd_spikes(:,row) = sum(PMd_trial,2);
    
    % only cue data
    M1_spikes_test(:,row) = sum(M1_trial(:,td_train(row).idx_target_on:td_train(row).idx_go_cue),2);
    PMd_spikes_test(:,row) = sum(PMd_trial(:,td_train(row).idx_target_on:td_train(row).idx_go_cue),2);
    y = [y; angle];
end

%% Build features vector
F_M1 = []; F_test_M1 = [];
F_PMd = []; F_test_PMd = [];

for i = 1: size(M1_spikes,2)
    % total M1 spikes
    total_M1_spikes = sum(M1_spikes(:,i));
    total_M1_spikes_test = sum(M1_spikes_test(:,i));
    % total PMd spikes
    total_PMd_spikes = sum(PMd_spikes(:,i));
    total_PMd_spikes_test = sum(PMd_spikes_test(:,i));
    
    f_M1 = M1_spikes(:,i)';
    f_test_M1 = M1_spikes_test(:,i)';
    
    f_PMd = PMd_spikes(:,i)';
    f_test_PMd = PMd_spikes_test(:,i)';
    % Average firing rates for full trial
    f_M1 = f_M1./total_M1_spikes;
    F_M1 = [F_M1;f_M1]; 
    f_PMd = f_PMd./total_PMd_spikes;
    F_PMd = [F_PMd;f_PMd]; 

    % Average firing rates only for t: idx_target_on - idx_go_cue
    f_test_M1 = f_test_M1./total_M1_spikes_test;
    F_test_M1 = [F_test_M1;f_test_M1]; 
    
    f_test_PMd = f_test_PMd./total_PMd_spikes_test;
    F_test_PMd = [F_test_PMd;f_test_PMd];
   
end
%%
F_test_M1_PMd = [F_test_M1,F_test_PMd];

%% Features vectors
F_final = [[F_test_M1,F_test_PMd,y];[F_M1, F_PMd, y]];
F_final_PMd = [[F_test_PMd,y];[F_PMd, y]];

