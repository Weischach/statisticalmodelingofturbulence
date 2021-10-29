clear all; close all; clc;


%% loading the 500 data point model to make comparisons between estimate and models%%



load logmldata500.mat

%% loading ensemble model parameters

cd minicurve

load ens_mods.mat

a1 = kf;
b1 = sz;
% a1 = 0;
% b1 = 25;
verfcn = 0;


%% loading file by file to create an ensemble model and model deviation for every individual models %%
t2 = 1;
for shuff = 2:2:8
    for k2 = 1:length(a1)
        kf = a1(k2);
        t1 = 1;
        sz = b1(k2);
        m1 = strcat(num2str(shuff),'-kf-pt', num2str(sz),'logpww');
        for j = 0:4
            kf1 = kf + j;
            m = strcat(m1,num2str(kf1),'.mat');
            load(m)
            for l = 1:size(testset,2)
                samp = testset(l)+1;
                nnorder(t1) = samp;
                for l1 = 1:(500/sz)
                    pwwtest((l1-1)*sz+1:l1*sz,t1) = test_y((l-1)*(500/sz)+l1,:);
                    pwwmodel((l1-1)*sz+1:l1*sz,t1) = yhat((l-1)*(500/sz)+l1,:);
                end
                gtnum(t2,t1) = samp;
                gt(t1) = samplist(samp);
                t1 = t1+1;
            end
        end
        [samtemp, samorder] = sort(nnorder); %% sort the NN shuffle
        pwwtest = pwwtest(:,[samorder]);
        pwwmodel = pwwmodel(:,[samorder]);
        pm1(t2,:,:) = pwwmodel;
        pwwver = reshape(pwwseq(:,23,:),size(pwwtest));
        pwwpred = reshape(pwwseq(:,22,:),size(pwwtest));
        verfn = mean(mean((pwwtest-pwwver).^2));
        verfcn = verfcn + verfn;    % verfcn = 0 implies the computation was model data was calculated correctly
        ind_model_devn(t2,:) = mean((pwwtest-pwwmodel).^2);
        act_devn = mean((pwwtest-pwwpred).^2);
        ind_mod_perf(t2,:) = (act_devn-ind_model_devn(t2,:));
        t2 = t2+1;
    end
end


ens_mod = reshape(mean(pm1,1),size(pwwmodel)); %for the ensemble average



model_devn = mean((pwwtest-ens_mod).^2);
mod_perf = (act_devn-model_devn);
[compar, corder] = sort(mod_perf); % performance of the ensemble
[ind_compar, ind_corder] = sort(ind_mod_perf,2); % performance of individual 48 models


%% list the models performing worse than the estimate



for k3 = 1:size(ind_compar,1)
    temp = ind_compar(k3,:);
    mod_lim(k3) = length(temp(temp<0));
end

smax = max(mod_lim);
bad_mods = zeros(size(ind_compar,1),smax);




for k2 = 1:size(ind_compar,1)
    for k3 = 1:mod_lim(k2)
        bad_mods(k2,k3) = ind_corder(k2,k3);
    end
    bad_mods(k2,:) = sort(bad_mods(k2,:));
end


%% count the sample recurring as performing worse than estimate in all the models 

nm = 1;
for k2 = min(min(bad_mods)):max(max(bad_mods))
    cnt1=0;
    for k3 = 1:size(ind_compar,1)
        cnt = 1;
        temp = bad_mods(k3,cnt);
        while temp< k2 + 1 && cnt<smax
            if temp == k2 && temp>0
                cnt1 = cnt1+1;
            end
            cnt = cnt+1;
            if cnt<smax+1
                temp = bad_mods(k3,cnt);
            end
        end
    end
    if cnt1>0
        sim(nm, :) = [k2, cnt1];
        nm = nm+1;
    end
end





%% plot the individual model performance with respect to the estimate for every sample in a sorted order %%



m = length(compar(compar<0));
plot(compar)
hold on
plot(zeros(249))
xlabel('samples')
ylabel('Difference between model and estimate')
title('Comparison between model error to best estimate error. +ve implies model performed better')




cd ../