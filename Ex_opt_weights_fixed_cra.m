clc; clear; close all;
addpath(genpath('functions'));
rng(2023);

%% Initialization
params.c = 340; % speed of sound
params.fband = [0, 8e3]; % frequency range
params.K = 75; % num of freq bins
params.Fs = 16000; % sample frequency
params.L = 32; % filter order (L+1 coefficeints), type 1 - symmetric
params.theta_d = 0;
params.theta_BW = 30; % desired beamwidth
params.A_BW = db2mag(-3); % amplitude at beamwidth
params.alpha = 1; % trade off between DF and WNG of cost obj
params.M = 5;  % num of rings

% uniform spacing
% params.Rm = linspace(0, 25, params.M).'*1e-2;
% params.Nm = ceil(4*pi*params.Rm*max(params.fband)/params.c);
% params.Nm(params.Nm == 0) = 1;

% nonuniform optimized spacing
params.Rm = [2, 4.8, 8.1, 13.9, 25].'*1e-2;
params.Nm = [6, 15, 17, 19, 19].';

params = update_params(params);

% calc FIR coefficients
[~, coeff] = calc_proposed_FIR_beamformer(params);

% calc filters frequency response
H = calc_freq_rep(coeff, params);

% calc beampattern
bp = B(params.T_normalized*H, params.d, params.f_grid, params.theta_grid);

% calc directivity factor as function of frequency
df = DF(params.T_normalized*H, params);

% calc white noise gain as function of frequency
wng = WNG(params.T_normalized*H, params);

% plot array geometry
plot_array_geometry("elem_pos", params.r, "partition", params.T);

% plot beampattern as function of frequency and theta-elevation
plot_beampattern_vs_f_vs_theta("bp", bp, "f", params.f_grid, "theta",...
    params.theta_grid, "normalize", false, "plt_contour", true);

% plot weights in freq domain
plot_weights(H, params)

% plot DF and WNG
figure;
subplot(211)
semilogx(params.f_grid/1e3, 10*log10(df));
grid on;
ylabel('Directivity');
xlabel('f[Hz]');
ylim([0 15]);
xlim([0, 8]);

subplot(212)
semilogx(params.f_grid/1e3, 10*log10(wng));
grid on;
ylabel('WNG');
xlabel('f[Hz]')
ylim([0 25]);
xlim([0, 8]);


