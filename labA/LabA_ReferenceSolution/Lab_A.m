clear all;
close all;
clc;

% Define noise statistics
noise_mean = 0;
noise_var = 9;

%% Task 1.1 Load the data
encoder_x = load('encoder.mat').encoder;
n_measurements_setup1=length(encoder_x); %full measurement
n_measurements_setup2=3; %first 3 measurement

%% Task 1.2 Compute the mean
mean_data_full=mean(encoder_x);
mean_data_setup2=mean(encoder_x(1:n_measurements_setup2));

%% Task 1.3 Compute the truncated mean
% Define prior distribution
prior_min = 21;
prior_max = 30;

% Data truncation
encoder_x_truncated=encoder_x;
encoder_x_truncated(encoder_x_truncated>prior_max)=prior_max;
encoder_x_truncated(encoder_x_truncated<prior_min)=prior_min;

% Mean value for the truncated data
mean_data_truncated_full=mean(encoder_x_truncated);
mean_data_truncated_setup2=mean(encoder_x_truncated(1:5));

%% Task 1.4 MMSE estimation with uniform prior
fun1 = @(x) (x/sqrt(2*pi*(noise_var)/n_measurements_setup1)).*exp(-n_measurements_setup1/(2*noise_var)*(x-sum(encoder_x)/n_measurements_setup1).^2);
fun2 = @(x) (1/sqrt(2*pi*(noise_var)/n_measurements_setup1)).*exp(-n_measurements_setup1/(2*noise_var)*(x-sum(encoder_x)/n_measurements_setup1).^2);
mmse_estimation_uniprior_full = integral(fun1,prior_min,prior_max)/integral(fun2,prior_min,prior_max);

fun1 = @(x) (x/sqrt(2*pi*(noise_var)/n_measurements_setup2)).*exp(-n_measurements_setup2/(2*noise_var)*(x-sum(encoder_x(1:n_measurements_setup2))/n_measurements_setup2).^2);
fun2 = @(x) (1/sqrt(2*pi*(noise_var)/n_measurements_setup2)).*exp(-n_measurements_setup2/(2*noise_var)*(x-sum(encoder_x(1:n_measurements_setup2))/n_measurements_setup2).^2);
mmse_estimation_uniprior_setup2 = integral(fun1,prior_min,prior_max)/integral(fun2,prior_min,prior_max);

%% Task 1.5 MMSE estimation with Gaussian prior
prior_mean = 25;
prior_var = 2;
data = (prior_var/(prior_var + noise_var/n_measurements_setup1))*sum(encoder_x)/n_measurements_setup1;
prior = (noise_var/n_measurements_setup1/(prior_var + noise_var/n_measurements_setup1))*prior_mean;
mmse_estimation_Gauprior_full = data + prior;

data = (prior_var/(prior_var + noise_var/n_measurements_setup2))*sum(encoder_x(1:n_measurements_setup2))/n_measurements_setup2;
prior = (noise_var/n_measurements_setup2/(prior_var + noise_var/n_measurements_setup2))*prior_mean;
mmse_estimation_Gauprior_setup2 = data + prior;


%% Task 2 Sequential MMSE estimation with Gaussian prior
prior_mean = 25;
prior_var_1 = 0.1;
prior_var_2 = 4;

mmse_estimation_Gauprior_seq=zeros(length(encoder_x),1);
mmse_estimation_Gauprior_seq(1)=prior_mean;
mmse_estimation_Gauprior_seq2=zeros(length(encoder_x),1);
mmse_estimation_Gauprior_seq2(1)=prior_mean;
moving_average=zeros(length(encoder_x),1);
moving_average(1)=encoder_x(1);

for i=1:n_measurements_setup1-1
    mmse_estimation_Gauprior_seq(i+1)=mmse_estimation_Gauprior_seq(i)+prior_var_1/((i+1)*prior_var_1+noise_var)*(encoder_x(i)-mmse_estimation_Gauprior_seq(i));
    mmse_estimation_Gauprior_seq2(i+1)=mmse_estimation_Gauprior_seq2(i)+prior_var_2/((i+1)*prior_var_2+noise_var)*(encoder_x(i)-mmse_estimation_Gauprior_seq2(i));
    moving_average(i+1)=mean(encoder_x(1:i+1));
end

% Plot measured angles and MMSE estimators
figure;
hold on;
scatter(1:n_measurements_setup1, encoder_x)
plot(1:n_measurements_setup1, mmse_estimation_Gauprior_seq, '-', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(1:n_measurements_setup1, mmse_estimation_Gauprior_seq2, '.-', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(1:n_measurements_setup1, moving_average, '--', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on
xlabel('Measurement index');
ylabel('Angle (degrees)');
legend('Raw measurement data','Sequential MMSE estimator with prior knowledge (25,0.1)', 'Sequential MMSE estimator with prior knowledge (25,4)','Moving average');

