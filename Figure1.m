clear all;clc;close all;
%% Experimental data 
load('DataFigure1.mat')
% From the data published in Woloszyn, L., & Sheinberg, D. L. (2012). 
% Effects of long-term visual experience on responses of distinct classes of 
% single units in inferior temporal cortex. Neuron, 74(1), 193-205. doi:10.1016/j.neuron.2012.01.032
% 
% Data_Sheinberg.mat has firing rates for novel and familiar stimuli - 
% Firing rates in each neuron were normalized using its mean and std of time-averaged firing rates
% and stimuli were rank-ordered according to time-averaged rates among the
% sets of novel and familiar stimuli, respectively.
% 
% Norm_RE_fam, Norm_RE_nov, Norm_RI_fam, Norm_RI_nov are normalized firing
% rates for fam and nov stimuli in selected neurons (19 exc and 9 inh)
% between 80 ms and 320 ms with time step dt = 5 ms
% size: # of stimuli by # of time steps by # of neurons

%%
NExc = size(Norm_RE_fam,3);
NInh = size(Norm_RI_fam,3);

% T_exp = 25:5:320;
NT = length(T_exp);
%% Time course
% smoothing the activity 
nt = 5;
for k = 1:NExc
    for i = 1:125;
        Norm_RE_fam_smooth(i,:,k) = smooth(squeeze(Norm_RE_fam(i,:,k)),nt);
        Norm_RE_nov_smooth(i,:,k) = smooth(squeeze(Norm_RE_nov(i,:,k)),nt);
    end
end

for k = 1:NInh
    for i = 1:125;
        Norm_RI_fam_smooth(i,:,k) = smooth(squeeze(Norm_RI_fam(i,:,k)),nt);
        Norm_RI_nov_smooth(i,:,k) = smooth(squeeze(Norm_RI_nov(i,:,k)),nt);
    end
end

mE_nov = mean(mean(Norm_RE_nov_smooth,1),3);
mE_fam = mean(mean(Norm_RE_fam_smooth,1),3);
mI_nov = mean(mean(Norm_RI_nov_smooth,1),3);
mI_fam = mean(mean(Norm_RI_fam_smooth,1),3);

sE_nov = std(mean(Norm_RE_nov_smooth,1),0,3)/sqrt(NExc);
sE_fam = std(mean(Norm_RE_fam_smooth,1),0,3)/sqrt(NExc);
sI_nov = std(mean(Norm_RI_nov_smooth,1),0,3)/sqrt(NInh);
sI_fam = std(mean(Norm_RI_fam_smooth,1),0,3)/sqrt(NInh);

figure; 
X = [T_exp fliplr(T_exp)];
Y = [sE_nov + mE_nov fliplr(-sE_nov+mE_nov)]; 
h = fill(X,Y,[0.5 0 0],'FaceAlpha',0.3);
    set(h,'EdgeColor','None');hold on;
plot(T_exp,mE_nov,'r','LineWidth',2);

X = [T_exp fliplr(T_exp)];
Y = [sE_fam + mE_fam fliplr(-sE_fam+mE_fam)]; 
h = fill(X,Y,[0 0 0.5],'FaceAlpha',0.3);
    set(h,'EdgeColor','None');
plot(T_exp,mE_fam,'b','LineWidth',2);
hold off;
xlim([50 320])

figure; 
X = [T_exp fliplr(T_exp)];
Y = [sI_nov + mI_nov fliplr(-sI_nov+mI_nov)]; 
h = fill(X,Y,[0.5 0 0],'FaceAlpha',0.3);
    set(h,'EdgeColor','None');hold on;
plot(T_exp,mI_nov,'r','LineWidth',2);

X = [T_exp fliplr(T_exp)];
Y = [sI_fam + mI_fam fliplr(-sI_fam+mI_fam)]; 
h = fill(X,Y,[0 0 0.5],'FaceAlpha',0.3);
    set(h,'EdgeColor','None');
plot(T_exp,mI_fam,'b','LineWidth',2);
hold off;
xlim([50 320])

mE_nov = mean(Norm_RE_nov_smooth(125,:,:),3);
mE_fam = mean(Norm_RE_fam_smooth(125,:,:),3);
mI_nov = mean(Norm_RI_nov_smooth(125,:,:),3);
mI_fam = mean(Norm_RI_fam_smooth(125,:,:),3);

sE_nov = std(Norm_RE_nov_smooth(125,:,:),0,3)/sqrt(NExc);
sE_fam = std(Norm_RE_fam_smooth(125,:,:),0,3)/sqrt(NExc);
sI_nov = std(Norm_RI_nov_smooth(125,:,:),0,3)/sqrt(NInh);
sI_fam = std(Norm_RI_fam_smooth(125,:,:),0,3)/sqrt(NInh);

figure; 
X = [T_exp fliplr(T_exp)];
Y = [sE_nov + mE_nov fliplr(-sE_nov+mE_nov)]; 
h = fill(X,Y,[0.5 0 0],'FaceAlpha',0.3);
    set(h,'EdgeColor','None');hold on;
plot(T_exp,mE_nov,'r','LineWidth',2);

X = [T_exp fliplr(T_exp)];
Y = [sE_fam + mE_fam fliplr(-sE_fam+mE_fam)]; 
h = fill(X,Y,[0 0 0.5],'FaceAlpha',0.3);
    set(h,'EdgeColor','None');
plot(T_exp,mE_fam,'b','LineWidth',2);
hold off;
xlim([50 320])

figure; 
X = [T_exp fliplr(T_exp)];
Y = [sI_nov + mI_nov fliplr(-sI_nov+mI_nov)]; 
h = fill(X,Y,[0.5 0 0],'FaceAlpha',0.3);
    set(h,'EdgeColor','None');hold on;
plot(T_exp,mI_nov,'r','LineWidth',2);

X = [T_exp fliplr(T_exp)];
Y = [sI_fam + mI_fam fliplr(-sI_fam+mI_fam)]; 
h = fill(X,Y,[0 0 0.5],'FaceAlpha',0.3);
    set(h,'EdgeColor','None');
plot(T_exp,mI_fam,'b','LineWidth',2);
hold off;
xlim([50 320])

%% Rebound strength measured by the slope of activity changes btw 230 ms and 320 ms
T_Slope = 230:5:320;
Tinit = 230; Tend = 320;
index_time = round((Tinit-T_exp(1))/dt)+1:round((Tend-T_exp(1))/dt)+1;  

NT_Slope = length(T_Slope);

RE_fam_Slope = Norm_RE_fam(:,index_time,:);
RE_nov_Slope = Norm_RE_nov(:,index_time,:);

RI_fam_Slope = Norm_RI_fam(:,index_time,:);
RI_nov_Slope = Norm_RI_nov(:,index_time,:);

Slope_Exc_nov = zeros(125,NExc);
Slope_Exc_fam = zeros(125,NExc);

for j = 1:NExc
    for i = 1:125
        X = RE_nov_Slope(i,:,j);
        f = fit(T_Slope',X','poly1');
        Slope_Exc_nov(i,j)= f.p1;
        
        X = RE_fam_Slope(i,:,j);
        f = fit(T_Slope',X','poly1');
        Slope_Exc_fam(i,j)= f.p1;
    end
end
m_Slope_Exc_nov = mean(Slope_Exc_nov,2);
s_Slope_Exc_nov = std(Slope_Exc_nov,0,2)/sqrt(NExc);
m_Slope_Exc_fam = mean(Slope_Exc_fam,2);
s_Slope_Exc_fam = std(Slope_Exc_fam,0,2)/sqrt(NExc);

figure;hold on
x = 1:125;
plot(x,m_Slope_Exc_nov,'r',x,m_Slope_Exc_fam,'b','LineWidth',1)

y = m_Slope_Exc_nov';
z = s_Slope_Exc_nov';
f  = fill([x flip(x)],[y+z flip(-z+y)],'r');
set(f,'EdgeColor','none','FaceAlpha',0.2)

y = m_Slope_Exc_fam';
z = s_Slope_Exc_fam';
f  = fill([x flip(x)],[y+z flip(-z+y)],'b');
set(f,'EdgeColor','none','FaceAlpha',0.2) 
hold off
xlim([0 125]);
xlabel('Neuronal Index');ylabel('Slope btw 200 and 320ms')
title('Exc Neuron')
legend('Nov','Fam')

Slope_Inh_nov = zeros(125,NInh);
Slope_Inh_fam = zeros(125,NInh);

for j = 1:NInh
    for i = 1:125
        X = RI_nov_Slope(i,:,j);
        f = fit(T_Slope',X','poly1');
        Slope_Inh_nov(i,j)= f.p1;

        X = RI_fam_Slope(i,:,j);
        f = fit(T_Slope',X','poly1');
        Slope_Inh_fam(i,j)= f.p1;
    end
end

m_Slope_Inh_nov = mean(Slope_Inh_nov,2);
s_Slope_Inh_nov = std(Slope_Inh_nov,0,2)/sqrt(NInh);
m_Slope_Inh_fam = mean(Slope_Inh_fam,2);
s_Slope_Inh_fam = std(Slope_Inh_fam,0,2)/sqrt(NInh);

figure;hold on
x = 1:125;
plot(x,m_Slope_Inh_nov,'r',x,m_Slope_Inh_fam,'b','LineWidth',1)
legend('Nov','Fam')

y = m_Slope_Inh_nov';
z = s_Slope_Inh_nov';
f  = fill([x flip(x)],[y+z flip(-z+y)],'r');
set(f,'EdgeColor','none','FaceAlpha',0.2)

y = m_Slope_Inh_fam';
z = s_Slope_Inh_fam';
f  = fill([x flip(x)],[y+z flip(-z+y)],'b');
set(f,'EdgeColor','none','FaceAlpha',0.2) 
hold off
xlim([0 125]); ylim([-0.1 0.1])
xlabel('Neuronal Index');ylabel('Slope btw 200 and 320ms')
title('Inh Neuron')