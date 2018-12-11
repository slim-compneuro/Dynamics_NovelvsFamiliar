clc;clear all;close all;
%% Experimental data
load('Data_Freedman.mat')
% Data_Freedman.mat has mean firing rates for nov and fam stimuli presented
% in "McKee, J. L., Thomas, S. L., & Freedman, D. J. (2013). Neuronal
% representations of novel and familiar visual stimuli in macaque inferior temporal, perirhinal and prefrontal cortices. Paper presented at the Society for Neuroscience, San Diego"

% t: time
% R_nov: mean rate for nov stim
% R_fam: mean rate for fam stim
dt  = 5;             % time step in the data (ms)

%% Simulation for one stimulus presentation
% Parameters
tau = 5;             % time constant for rate dynamics (ms)
k = 1.8;            % Strength of adaptation
w = 0;              % uniform connectivity strength before learning
ws = .9;            % Potentiation in high rate neurons after learning = strength of positive fdbk in m dynamics

% Adapation variables which are linearly filtered firing rates
tau_A = 200;        % time constant for adaptation current (ms)
alpha_A = dt/tau_A;

A_nov = R_nov;
A_fam = R_fam;
for i = 2:length(t)
    A_nov(i) = (1-alpha_A)*A_nov(i-1) + alpha_A*R_nov(i-1);
    A_fam(i) = (1-alpha_A)*A_fam(i-1) + alpha_A*R_fam(i-1);
end

t_fit = 575:dt:1100;
index_fit = round(t_fit(1)/dt):round(t_fit(end)/dt);
Inov_fit = R_nov(index_fit);
Ifam_fit = R_fam(index_fit);
Rnov_fit = R_nov(index_fit-1);
Rfam_fit = R_fam(index_fit-1);
Anov_fit = A_nov(index_fit-1);
Afam_fit = A_fam(index_fit-1);

% % External input for mean (IextMean for r dynamics) and max(IextMax for m dynamics) rate
Iext = Inov_fit - w*Rnov_fit + k*Anov_fit;
% DiffExpEqn = 'a*exp(-x/t1)+b*exp(-x/t2)-a-b';
% xdata = t_fit-t_fit(1);
% ydata = (Iext-Iext(1))';
% IextFit = fit(xdata',ydata',DiffExpEqn,'Start',[-6 5 25 200]);
IextFit.a = -5;
IextFit.b = 6;
IextFit.t1 = 40;
IextFit.t2 = 700;

dt_sim = .1;
t_sim = t(index_fit(1)):dt_sim:t(index_fit(end));
IextMean = IextFit.a*(exp(-(t_sim-t_sim(1))/IextFit.t1)-1)+IextFit.b*(exp(-(t_sim-t_sim(1))/IextFit.t2)-1)+Iext(1);
IextMax  = IextFit.a*(exp(-(t_sim-t_sim(1))/20)-1)-IextFit.a*(exp(-(t_sim-t_sim(1))/400)-1)+Iext(1);    % shorter rise and decay time constant for faster dynamics after learning
IextMean = IextMean';
IextMax  = IextMax';

% Simulation of firing rate before learning, Rnov_sim
Rnov_sim(1) = Rnov_fit(1);
Anov_sim(1) = Rnov_fit(1);

for i = 2:length(t_sim)
    Rnov_sim(i) = Rnov_sim(i-1) + dt_sim/dt*(-Rnov_sim(i-1)+w*Rnov_sim(i-1)-k*(Anov_sim(i-1)-Anov_sim(1))+ (IextMean(i-1)-IextMean(1)) +Rnov_sim(1));
    Anov_sim(i) = Anov_sim(i-1) + dt_sim/tau_A*(-Anov_sim(i-1) + Rnov_sim(i-1));
end
figure;
plot(t(index_fit),Rnov_fit,'Color',[0.5 0.5 0.5]);hold on
plot(t_sim,Rnov_sim,'r');hold off
xlim([t(index_fit(1)) t(index_fit(end))]);
xlabel('Time (ms)');ylabel('Firing Rate (Hz)')
legend('Emp. Nov.','Sim. Nov.')

%  Simulation of firing rate after learning, Rfam_sim with m_sim
m_sim = 0;
n_sim = 0;

for i = 2:length(t_sim)+1
    m_sim(i) = m_sim(i-1) + dt_sim/dt*(-m_sim(i-1)+ws*(m_sim(i-1)) - k*(n_sim(i-1)) + (IextMax(i-1)-Iext(1)));
    n_sim(i) = n_sim(i-1) + dt_sim/tau_A*(-n_sim(i-1) + m_sim(i-1));
end
% figure;
% plot(t_sim,m_sim(2:end)/max(m_sim),'k')
% xlim([t(index_fit(1)) t(index_fit(end))]);
% xlabel('Time (ms)');ylabel('Normalized m')

fpost_FF = -0.7;
fpost_Rec = 0.3;

Rfam_sim(1) = Rfam_fit(1);
Afam_sim(1) = Rfam_fit(1);

for i = 2:length(t_sim)
    Rfam_sim(i) = Rfam_sim(i-1) + dt_sim/dt*(-Rfam_sim(i-1)+w*Rfam_sim(i-1)-k*(Afam_sim(i-1)-Afam_sim(1))...
        + fpost_Rec*m_sim(i) + (IextMean(i-1)-IextMean(1)) + fpost_FF*(IextMax(i-1)-Iext(1))+Rfam_sim(1));
    Afam_sim(i) = Afam_sim(i-1) + dt_sim/tau_A*(-Afam_sim(i-1) + Rfam_sim(i-1));
end

figure;
plot(t(index_fit),Rnov_fit,'Color',[0.5 0.5 0.5]);hold on
plot(t(index_fit),Rfam_fit,'k');
plot(t_sim,Rnov_sim,'r','LineWidth',1)
plot(t_sim,Rfam_sim,'b','LineWidth',1)
xlim([t(index_fit(1)) t(index_fit(end))]);
xlabel('Time (ms)');ylabel('Firing Rate (Hz)')
legend('Emp. Nov','Emp. Fam','Sim. Nov','Sim. Fam')


%% Simulation for multiple stimulus presentation 
StimPeriod = 150;
NT_Period = StimPeriod/dt_sim;
T_multiple = t(index_fit(1)):dt_sim:t(index_fit(1))+StimPeriod*3;
NT_multiple = length(T_multiple);

Iext_multiple = zeros(NT_multiple,1);
I_const = IextMean(1);
Iext_multiple(1:NT_Period) = IextMean(1:NT_Period)-I_const;
DecayRate = exp(-dt_sim/50);
for k = 1:2;
    for j = k*NT_Period+1:(k+1)*NT_Period
        Iext_multiple(j) = Iext_multiple(j-1)*DecayRate;
    end
    Iext_multiple(k*NT_Period+1:(k+1)*NT_Period) = Iext_multiple(k*NT_Period+1:(k+1)*NT_Period)+IextMean(1:NT_Period)-I_const;
end
Iext_multiple = Iext_multiple+I_const;

IextMax_multiple = zeros(NT_multiple,1);
I_const = IextMean(1);
IextMax_multiple(1:NT_Period) = IextMax(1:NT_Period)-I_const;
for k = 1:2;
    for j = k*NT_Period+1:(k+1)*NT_Period
        IextMax_multiple(j) = IextMax_multiple(j-1)*DecayRate;
    end
    IextMax_multiple(k*NT_Period+1:(k+1)*NT_Period) = IextMax_multiple(k*NT_Period+1:(k+1)*NT_Period)+IextMax(1:NT_Period)-I_const;
end
IextMax_multiple = IextMax_multiple+I_const;

Rnov_multiple(1) = Rnov_fit(1);
Anov_multiple(1) = Rnov_fit(1);
Rfam_multiple(1) = Rfam_fit(1);
Afam_multiple(1) = Rfam_fit(1);
m_multiple = 0;
n_multiple = 0;

for i = 2:NT_multiple+1
    m_multiple(i) = m_multiple(i-1) + dt_sim/dt*(-m_multiple(i-1)+ws*(m_multiple(i-1)) - k*(n_multiple(i-1)) + (IextMax_multiple(i-1)-Iext(1)));
    n_multiple(i) = n_multiple(i-1) + dt_sim/tau_A*(-n_multiple(i-1) + m_multiple(i-1));
end

for i = 2:NT_multiple
    Rnov_multiple(i) = Rnov_multiple(i-1) + dt_sim/dt*(-Rnov_multiple(i-1)+w*Rnov_multiple(i-1)-k*(Anov_multiple(i-1)-Anov_multiple(1)) + (Iext_multiple(i-1)-IextMean(1)) +Rnov_multiple(1));
    Anov_multiple(i) = Anov_multiple(i-1) + dt_sim/tau_A*(-Anov_multiple(i-1) + Rnov_multiple(i-1));

    Rfam_multiple(i) = Rfam_multiple(i-1) + dt_sim/dt*(-Rfam_multiple(i-1)+w*Rfam_multiple(i-1)-k*(Afam_multiple(i-1)-Afam_multiple(1))...
        + fpost_Rec*m_multiple(i) + (Iext_multiple(i-1)-IextMean(1)) + fpost_FF*(IextMax_multiple(i-1)-Iext(1))+Rfam_multiple(1));
    Afam_multiple(i) = Afam_multiple(i-1) + dt_sim/tau_A*(-Afam_multiple(i-1) + Rfam_multiple(i-1));
end

figure;
plot(T_multiple,Rnov_multiple,'r','LineWidth',1)
hold on;
plot(T_multiple,Rfam_multiple,'b','LineWidth',1)
hold off
xlim([T_multiple(1) T_multiple(end)]);
xlabel('Time (ms)');ylabel('Firing Rate (Hz)')
legend('Sim. Nov','Sim. Fam')
title(['Mean Rate for serial presentation with T=',num2str(StimPeriod)])