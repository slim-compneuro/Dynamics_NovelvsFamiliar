close all;clear all;
%%  Parameters
w = 0;          % uniform connectivity before learning
tauR = 5;       % time constant for firing rate (ms)
tauA = 200;     % time constant for adaptation (ms)
kA = 2;         % strength of adaptation
    
dt = 0.05;      % numerical time step
t = 0:dt:600;   % time vector
NT = length(t);

Nneuron = 2000; % number of neurons

alpha = 0.9;    % strength of Hebbian learning rule for recurrent connections
FFdep = 0.4;    % uniform scaling factor for feedforward connection
 
%% Input-output transfer function f and its inverse
f = @(x) x;
g = @(x) x;
%% External input current = Constant Input  + InputFactor*Temporal Profile 
r_base = 5;     % baseline firing rate (Hz)
InputConst = feval(g,r_base)-(w-kA)*r_base; % Constant Input generating baseline firing rate

% Temporal profile which is the sum of two exponential
tRise = 50; tDecay = 150;
InputTemporal = -3*exp(-t/tRise)+3*exp(-t/tDecay);
m = max(InputTemporal);
InputTemporal = InputTemporal/m;

% Different input strength following gamma distribution with parameter lambda
lambda = 3;
InputFactor = sort(.8*randg(lambda,Nneuron,1));
FFTemporal = InputFactor*InputTemporal;

%% Simulation for one stimulus presentation
% Hebbian learning rule of the form of r_i*(r_j-mean(r))
DelW = alpha/Nneuron/var(InputFactor)*InputFactor*(InputFactor-mean(InputFactor))';

rE_nov = r_base*ones(Nneuron,NT);
aE_nov = rE_nov;

rE_fam = r_base*ones(Nneuron,NT);
aE_fam = rE_fam;

for j = 1:NT-1
    Input = w*mean(rE_nov(:,j)) - kA*aE_nov(:,j)+ InputConst+ FFTemporal(:,j) ;
    rE_nov(:,j+1) = rE_nov(:,j) + dt/tauR*(-rE_nov(:,j)+feval(f,Input));
    aE_nov(:,j+1) = aE_nov(:,j) + dt/tauA*(-aE_nov(:,j) + rE_nov(:,j));
    
    Input = w*mean(rE_fam(:,j)) + DelW*rE_fam(:,j) - kA*aE_fam(:,j)+ InputConst + FFdep*FFTemporal(:,j);
    rE_fam(:,j+1) = rE_fam(:,j) + dt/tauR*(-rE_fam(:,j)+feval(f,Input));
    aE_fam(:,j+1) = aE_fam(:,j) + dt/tauA*(-aE_fam(:,j) + rE_fam(:,j));
end

figure;
plot(t,mean(rE_nov),t,mean(rE_fam));
xlim([0 300])
Index_rep = Nneuron/8:Nneuron/8:Nneuron;
figure;
plot(t,rE_nov(Index_rep,:))
xlim([0 300]);ylim([2 15])
figure;
plot(t,rE_fam(Index_rep,:))
xlim([0 300]);ylim([2 15])

%% Simulation for successive presentation of two stimuli
rng('default')

% Permutation of the stimuli for another stimuli
IndexNew = randperm(Nneuron);
DelW = alpha/Nneuron/var(InputFactor)*InputFactor*(InputFactor-mean(InputFactor))'...
      +alpha/Nneuron/var(InputFactor)*InputFactor(IndexNew)*(InputFactor(IndexNew)-mean(InputFactor))';

% period of one stimulus presentation with linear decay after the offset of one stimulus
period = 150;
decayWindow = 100;
TempWeight1 = (t<=period)-(t>period).*(t<period+decayWindow).*(t-(period+decayWindow))/decayWindow;
TempWeight2 = (t>period).*(t<=2*period) -(t>2*period).*(t<2*period+decayWindow).*(t-(2*period+decayWindow))/decayWindow;
TempWeight3 = (t>2*period).*(t<=3*period) -(t>3*period).*(t<3*period+decayWindow).*(t-(3*period+decayWindow))/decayWindow;
TempWeight4 = (t>3*period).*(t<=4*period) -(t>4*period).*(t<4*period+decayWindow).*(t-(4*period+decayWindow))/decayWindow;

for i = 1:Nneuron
    FFTemporal(i,:) = InputFactor(i)*TempWeight1.*(-3*exp(-t/tRise)+3*exp(-t/tDecay))/m...
        + InputFactor(IndexNew(i))*TempWeight2.*(-3*exp(-(t-period)/tRise)+3*exp(-(t-period)/tDecay))/m...
        + InputFactor(i)*TempWeight3.*(-3*exp(-(t-2*period)/tRise)+3*exp(-(t-2*period)/tDecay))/m...
        + InputFactor(IndexNew(i))*TempWeight4.*(-3*exp(-(t-3*period)/tRise)+3*exp(-(t-3*period)/tDecay))/m;
end

rE_nov_rep = r_base*ones(Nneuron,NT);
aE_nov_rep = rE_nov;
rE_fam_rep = r_base*ones(Nneuron,NT);
aE_fam_rep = rE_fam_rep;

for j = 1:NT-1
    Input = w*mean(rE_nov_rep(:,j)) + InputConst + FFTemporal(:,j) - kA*aE_nov_rep(:,j);
    rE_nov_rep(:,j+1) = rE_nov_rep(:,j) + dt/tauR*(-rE_nov_rep(:,j)+feval(f,Input));
    aE_nov_rep(:,j+1) = aE_nov_rep(:,j) + dt/tauA*(-aE_nov_rep(:,j) + rE_nov_rep(:,j));
    
    Input = w*mean(rE_fam_rep(:,j)) + DelW*rE_fam_rep(:,j)+ InputConst + FFdep*FFTemporal(:,j) - kA*aE_fam_rep(:,j);
    rE_fam_rep(:,j+1) = rE_fam_rep(:,j) + dt/tauR*(-rE_fam_rep(:,j)+feval(f,Input));
    aE_fam_rep(:,j+1) = aE_fam_rep(:,j) + dt/tauA*(-aE_fam_rep(:,j) + rE_fam_rep(:,j));
end

figure;
plot(t,mean(rE_nov_rep),t,mean(rE_fam_rep));

IndexNewSample1 = find(IndexNew == Nneuron);
IndexNewSample2 = find(IndexNew == Nneuron*3/4);
IndexNewSample3 = find(IndexNew == Nneuron*2/4);
IndexNewSample4 = find(IndexNew == Nneuron*1/4);
Index_rep = [Nneuron/4:Nneuron/4:Nneuron IndexNewSample1 IndexNewSample2 IndexNewSample3 IndexNewSample4];

figure;
plot(t,rE_nov_rep(Index_rep,:))

figure;
plot(t,rE_fam_rep(Index_rep,:))