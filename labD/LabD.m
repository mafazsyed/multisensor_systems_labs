
clc
clear all
close all

%% CUSUM test algorithm on strain gauge sensor
% Declaration of algorithm parameters
load('straingauge.mat')
theta_0 = 3000;
sigma_0 = 1;
threshold_pos = 20;
threshold_neg = -20;
n_alarm_pos=0;
n_alarm_neg=0;
k_alarm_pos=[];
k_alarm_neg=[];

% Two-Sided CUSUM Test with threshold = 20
g_pos = 0*straingauge;
g_neg = 0*straingauge;
s = (straingauge - theta_0)/sigma_0;
for k = 1:size(straingauge,2)-1
 g_pos(k+1) = g_pos(k) + s(k); %neglect leakage term
 g_neg(k+1) = g_neg(k) + s(k); %neglect leakage term

 %positive test
 if g_pos(k+1) <0
     g_pos(k+1)=0;
 end
 if g_pos(k+1) > threshold_pos
    n_alarm_pos=n_alarm_pos+1;
    k_alarm_pos=[k_alarm_pos;k+1];
    % g_pos(k+1)=0; %reset
 end

 %negative test
 if g_neg(k+1) >0
     g_neg(k+1)=0;
 end
 if g_neg(k+1) < threshold_neg
    n_alarm_neg=n_alarm_neg+1;
    k_alarm_neg=[k_alarm_neg;k+1];
    % g_neg(k+1)=0;%reset
 end
end

figure
p1=plot(g_pos);
hold on
p2=plot(g_neg);
for i=1:length(k_alarm_pos)
    p3=plot([k_alarm_pos(i) k_alarm_pos(i)],[-20 20],'r--');
end
for i=1:length(k_alarm_neg)
    p4=plot([k_alarm_neg(i) k_alarm_neg(i)],[-20 20],'r-.');
end
legend([p1,p2,p3,p4],'Positive Test', 'Negative Test','Alarms for positive test','Alarms for negative test')
xlabel('Step')
ylabel('g_t')