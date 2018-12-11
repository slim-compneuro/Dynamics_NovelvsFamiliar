close all;clc;clear all;
%% Experimental data
load('Data_Sheinberg.mat')
% From the data published in Woloszyn, L., & Sheinberg, D. L. (2012). 
% Effects of long-term visual experience on responses of distinct classes of 
% single units in inferior temporal cortex. Neuron, 74(1), 193-205. doi:10.1016/j.neuron.2012.01.032
% 
% Data_Sheinberg.mat has firing rates for novel and familiar stimuli - 
% Firing rates in each neuron were normalized using its mean and std of time-averaged firing rates
% and stimuli were rank-ordered according to time-averaged rates among the
% sets of novel and familiar stimuli, respectively.
% At each rank of stimuli, normalized rates were averaged over neurons
%
% dt = 5; 
% T_exp, RE_fam_exp and RE_nov_exp are time and firing rates for novel and familiar stimuli.
% Mean_NormalizedExcR_Nov and Input_current for input-output transfer
% function - Mean_NormalizedExcR_Nov is normalized activities at each rank
% of novel stimuli averaged over neurons and Input_current is drawn from
% Gaussian statistics with mean 0 and std 1 (see S Lim, NN 2015)

%% Inference of post-synaptic dependence of recurrent and feedforward connections
% Parameters
w = 0;              % uniform connectivity strength before learning
k = 1.8;            % Strength of adaptation
tau_A = 200;        % time constant for adaptation current (ms)
alpha_A = dt/tau_A;

index_stim = 1:125; % 125 novel and familiar stimuli each
Nstim = 125;
%%
% Linear or nonlinear dynamics (flag_nonlinear = 0 if linear or
% flag_nonlinear = 1 if nonlinear)
flag_nonlinear = 0;
if (flag_nonlinear ==0)
    Mean_NormalizedExcR_Nov = Input_Current;
end

tic
% larger weights lambda in the late phase to capture the rebound better (lambda = 5 from 230 ms after the stimulus onset, and otherwise 1)
lambda = 5;              

% Fitting between 80 and 300 ms after the stimulus onset
Tinit = 80;Tend = 300; 
index_time = round((Tinit-T_exp(1))/dt)+1:round((Tend-T_exp(1))/dt);  

T = T_exp(index_time);   % time bwt 80 and 200
NT = length(T);

RE_fam = RE_fam_exp(:,index_time);
RE_nov = RE_nov_exp(:,index_time);

for i = 1:Nstim
    IE_nov(i,:) = interp1(Mean_NormalizedExcR_Nov,Input_Current,RE_nov(i,:),'linear','extrap');
    IE_fam(i,:) = interp1(Mean_NormalizedExcR_Nov,Input_Current,RE_fam(i,:),'linear','extrap');
end

sE_nov_A = RE_nov;
sE_fam_A = RE_fam;
for j = 2:NT
    sE_nov_A(:,j)  = alpha_A*RE_nov(:,j-1)+(1-alpha_A)*sE_nov_A(:,j-1);
    sE_fam_A(:,j)  = alpha_A*RE_fam(:,j-1)+(1-alpha_A)*sE_fam_A(:,j-1);
end

% Inferring post-synaptic dependence from maximum rate (stimulus rank
% corresponding to i_max = 125)
fit_start = 2;

% External input for mean (mean_I) and max (max_I) inferred from firing
% rate for novel stimuli
mean_I = mean(IE_nov(:,fit_start:end))-w*(mean(RE_nov(:,fit_start-1:end-1),1)-mean(RE_nov(:,fit_start-1),1))+k*(mean(sE_nov_A(:,fit_start-1:end-1),1)-mean(sE_nov_A(:,fit_start-1),1));
i_max = 125;
max_I = mean(IE_nov(i_max,fit_start:end),1)-w*(mean(RE_nov(:,fit_start-1:end-1),1)-mean(RE_nov(:,fit_start-1),1))+k*(mean(sE_nov_A(i_max,fit_start-1:end-1),1)-mean(sE_nov_A(i_max,fit_start-1,1)));

max_fpost_Rec_r = -0.1:0.01:0.5;
max_fpost_FF_r = -1:0.02:1;
error = zeros(length(max_fpost_Rec_r),length(max_fpost_Rec_r));
for p = 1:length(max_fpost_Rec_r)
    for q = 1:length(max_fpost_FF_r)
        max_fpost_Rec = max_fpost_Rec_r(p);
        max_fpost_FF = max_fpost_FF_r(q);

        % validation
        max_RE_fam_Sim = mean(RE_fam(i_max,:),1);
        max_sE_fam_A_Sim = mean(sE_fam_A(i_max,:),1);

        for j = 1:NT-1
            max_Input = - k*(max_sE_fam_A_Sim(j)-max_sE_fam_A_Sim(1))...
                + max_fpost_Rec*(max_RE_fam_Sim(j)-max_RE_fam_Sim(1)) + max_I(j) + max_fpost_FF*(max_I(j)-max_I(1));

            max_sE_fam_A_Sim(j+1)  = alpha_A*max_RE_fam_Sim(j)+(1-alpha_A)*max_sE_fam_A_Sim(j);
            max_RE_fam_Sim(j+1) = interp1(Input_Current,Mean_NormalizedExcR_Nov,max_Input,'linear','extrap');
        end

        error(p,q) = sum((max_RE_fam_Sim(1:30)-RE_fam(i_max,1:30)).^2)+lambda*sum((max_RE_fam_Sim(31:end)-RE_fam(i_max,31:end)).^2);

    end
end

[M1,I1] = min(error);
[M2,I2] = min(M1);

min_max_fpost_Rec = max_fpost_Rec_r(I1(I2));
min_max_fpost_FF = max_fpost_FF_r(I2);

imagesc(max_fpost_FF_r,max_fpost_Rec_r,error)
set(gca,'YDir','normal')

max_fpost_Rec = min_max_fpost_Rec;
max_fpost_FF = min_max_fpost_FF;

% validation
max_RE_fam_Sim = mean(RE_fam(i_max,:),1);
max_sE_fam_A_Sim = mean(sE_fam_A(i_max,:),1);

for j = 1:NT-1
    max_Input = - k*(max_sE_fam_A_Sim(j)-max_sE_fam_A_Sim(1))...
        + max_fpost_Rec*(max_RE_fam_Sim(j)-max_RE_fam_Sim(1)) + max_I(j) + max_fpost_FF*(max_I(j)-max_I(1));

    max_sE_fam_A_Sim(j+1)  = alpha_A*max_RE_fam_Sim(j)+(1-alpha_A)*max_sE_fam_A_Sim(j);
    max_RE_fam_Sim(j+1) = interp1(Input_Current,Mean_NormalizedExcR_Nov,max_Input,'linear','extrap');
end

figure;plot(T,mean(RE_nov(i_max,:),1),'r');hold on
plot(T,mean(RE_fam(i_max,:),1),'b')
plot(T,max_RE_fam_Sim,'k','LineWidth',2);
legend('Emp. Nov.','Emp. Fam.','Sim. Fam.')
xlabel('Time (ms)');ylabel('Normalized rate')
title('Maximum firing rate')
xlim([T(1) T(end)])
%% Fitting the response to the remaining rank of stimuli assuming that the maximum is the same as the experimental data
max_RE_fam_Exp = RE_fam(i_max,:);

fpost_Rec =zeros(Nstim,1);
fpost_FF = zeros(Nstim,1);
Iext = zeros(Nstim,43);

fpost_Rec_r = -0.1:0.01:0.5;
fpost_FF_r = -1:0.02:1;

RE_fam_Sim_o = zeros(size(RE_fam));
for i = 1:124
    index = i;
    error = zeros(length(fpost_Rec_r),length(fpost_FF_r));
    for p = 1:length(fpost_Rec_r)
        for q = 1:length(fpost_FF_r)
            i_fpost_Rec = fpost_Rec_r(p);
            i_fpost_FF = fpost_FF_r(q);

            I_index = mean(IE_nov(index,2:end),1)+k*(mean(sE_nov_A(index,1:end-1),1)-mean(sE_nov_A(index,1),1))...
                    +mean(IE_fam(index,1))-mean(IE_nov(index,1));

    % validation
    RE_fam_Sim = RE_fam(i,:);
    sE_fam_adap_Sim = sE_fam_A(i,:);

    for j = 1:NT-1
        Input = - k*(sE_fam_adap_Sim(j)-sE_fam_adap_Sim(1))...
            + i_fpost_Rec*(max_RE_fam_Exp(j)-max_RE_fam_Exp(1)) + I_index(j) + i_fpost_FF*(max_I(j)-max_I(1));

        sE_fam_adap_Sim(j+1)  = alpha_A*RE_fam_Sim(j)+(1-alpha_A)*sE_fam_adap_Sim(j);
        RE_fam_Sim(j+1) = interp1(Input_Current,Mean_NormalizedExcR_Nov,Input,'linear','extrap');
    end

    error(p,q) = sum((RE_fam_Sim(1:30)-RE_fam(i,1:30)).^2)+lambda*sum((RE_fam_Sim(31:end)-RE_fam(i,31:end)).^2);

        end
    end

    [M1,I1] = min(error);
    [M2,I2] = min(M1);

    min_i_fpost_Rec = fpost_Rec_r(I1(I2));
    min_i_fpost_FF = fpost_FF_r(I2);

    i_fpost_Rec = min_i_fpost_Rec;
    i_fpost_FF = min_i_fpost_FF;

    % validation
    RE_fam_Sim = RE_fam(i,:);
    sE_fam_adap_Sim = sE_fam_A(i,:);

    for j = 1:NT-1
        Input = - k*(sE_fam_adap_Sim(j)-sE_fam_adap_Sim(1))...
            + i_fpost_Rec*(max_RE_fam_Exp(j)-max_RE_fam_Exp(1)) + I_index(j) + i_fpost_FF*(max_I(j)-max_I(1));

        sE_fam_adap_Sim(j+1)  = alpha_A*RE_fam_Sim(j)+(1-alpha_A)*sE_fam_adap_Sim(j);
        RE_fam_Sim(j+1) = interp1(Input_Current,Mean_NormalizedExcR_Nov,Input,'linear','extrap');
    end

    fpost_Rec(i) = i_fpost_Rec;
    fpost_FF(i) = i_fpost_FF;
    Iext(i,:) = I_index;
    
    RE_fam_Sim_o(i,:) = RE_fam_Sim;
end

Iext(Nstim,:) = max_I;
fpost_Rec(Nstim) = max_fpost_Rec;
fpost_FF(Nstim) = max_fpost_FF;
RE_fam_Sim_o(Nstim,:) = max_RE_fam_Sim;

figure;plot(T,mean(RE_nov,1),'r');hold on
plot(T,mean(RE_fam,1),'b')
plot(T,mean(RE_fam_Sim_o),'k','LineWidth',2);
legend('Emp. Nov.','Emp. Fam.','Sim. Fam.')
xlabel('Time (ms)');ylabel('Normalized rate')
title('Mean firing rate')
xlim([T(1) T(end)])

figure;
subplot(1,2,1);plot(fpost_Rec,'ko');xlim([0 125])
ylabel('fpost_{Rec}')
subplot(1,2,2);plot(fpost_FF,'ko');xlim([0 125])
ylabel('fpost_{FF}')