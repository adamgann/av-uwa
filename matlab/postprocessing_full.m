% postprocessing_full.m
% 
% Generate Fig 5, showing the SINR and chosen rate over time. 
% 
% Adam Gannnon, SUNY Buffalo, 2018. 

clear variables;
close all;
clc;

%% Parameters 

fixed24Name = 'workspace_fixedL24.mat';
adaptive_name = 'workspace_adaptiveL.mat';


%% Load the 

load(fixed24Name)
ber24 = berAv;
av24 = sinrAv;
dur24 = packetDurStore;
rate24  = rateStore;
errFreeData24 = totalDataTransfer/8;

load(adaptive_name)
bera = berAv;
ava = sinrAv;
prefilta = sinrPrefilt;
tveca = sinrTvec;
dura = packetDurStore;
ratea = rateStore;
errFreeDataa = totalDataTransfer/8;


%% Plot Data


rateFig = figure;

% First subplot is postfilt SINR
subplot(2,1,1);
plot(sinrTvec,ava);

% Add the refline target of 10db
targetSinr = 10;
refline = targetSinr.*ones(size(ava));
hold on
plot(sinrTvec,refline,'r--','LineWidth',2);

xlim([0 400])
ylim([-1 20])

xH =xlabel('Time (s)');
yH =ylabel({'Pre-Detection','SINR (dB)'});
lH = legend({'Instantaneous','Target'},'Orientation','horizontal','location','SouthEast');

set(lH,'Interpreter','latex','Fontsize',12);
set(xH,'Interpreter','latex','Fontsize',14);
set(yH,'Interpreter','latex','Fontsize',14);

% Second subplot is the data rate
subplot(2,1,2);
plot(sinrTvec,ratea);
hold on
plot(sinrTvec,rate24,'r:','LineWidth',3);

ylim([0 7])
xlim([0 400])

xH = xlabel('Time (s)');
yH = ylabel('Data Rate (kbps)');
lH = legend({'Proposed Adaptive  ','Fixed (L=24)  '},'Orientation','horizontal');

set(lH,'Interpreter','latex','Fontsize',12);
set(xH,'Interpreter','latex','Fontsize',14);
set(yH,'Interpreter','latex','Fontsize',14);


