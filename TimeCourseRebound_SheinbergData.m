clear all;clc;close all;
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

%% Time Course
% T_exp = 25:5:320;
NT = length(T_exp);
figure; hold on
plot(T_exp,mean(RE_nov_exp),'r','LineWidth',2);
plot(T_exp,mean(RE_fam_exp),'b','LineWidth',2);
hold off;
xlim([50 320])
legend('Nov','Fam')
xlabel('Time (ms)');ylabel('Normalized firing rate')
title('Mean excitatory response')

figure; hold on
plot(T_exp,RE_nov_exp(end,:),'r','LineWidth',2);
plot(T_exp,RE_fam_exp(end,:),'b','LineWidth',2);
hold off;
xlim([50 320])
legend('Nov','Fam')
xlabel('Time (ms)');ylabel('Normalized firing rate')
title('Max excitatory response')


%% Rebound strength measured by the slope of activity changes btw 230 ms and 320 ms
T_Slope = 230:5:320;
Tinit = 230; Tend = 320;
index_time = round((Tinit-T_exp(1))/dt)+1:round((Tend-T_exp(1))/dt)+1;  

NT_Slope = length(T_Slope);

RE_fam_Slope = RE_fam_exp(:,index_time);
RE_nov_Slope = RE_nov_exp(:,index_time);

for i = 1:125
    X = RE_nov_Slope(i,:);
    f = fit(T_Slope',X','poly1');
    Slope_Exc_nov(i)= f.p1;

    X = RE_fam_Slope(i,:);
    f = fit(T_Slope',X','poly1');
    Slope_Exc_fam(i)= f.p1;
    
end

figure;
x = 1:125;
plot(x,Slope_Exc_nov,'r',x,Slope_Exc_fam,'b','LineWidth',1)
xlim([0 125]);
legend('Nov','Fam')
xlabel('Neuronal Index');ylabel('Slope btw 200 and 320ms')
title('Exc Neuron')
legend('Nov','Fam')