clc; clear; close all;
addpath(genpath('functions'));
rng(282);

%% CRA configuration
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
params.Rmin = 0;
params.Rmax = 25e-2;

% uniform spacing
params.Rm = linspace(0, 25, params.M).'*1e-2;
params.Nm = ceil(4*pi*params.Rm*max(params.fband)/params.c);
params.Nm(params.Nm == 0) = 1;

% params.Rm = linspace(0, 25, params.M).'*1e-2;
% params.Nm = [6, 15, 17, 19, 19].';

% % nonuniform optimized spacing
% params.Rm = [2, 4.8, 8.1, 13.9, 25].'*1e-2;
% params.Nm = [6, 15, 17, 19, 19].';

params = update_params(params);

% plot array geometry
% plot_array_geometry("elem_pos", params.r, "partition", params.T);

% GA parameters

ga_params.iter_max = 100; % max number of generations
ga_params.p_mu = 0.2; % probability of a chromosome to be mutated
ga_params.n_pop = 100;
ga_params.n_genes = params.M; 
ga_params.n_elite = 0.1*ga_params.n_pop;
ga_params.n_crossover = (ga_params.n_pop-ga_params.n_elite);

ga_params.n_mutation = ga_params.n_pop - ga_params.n_elite - ga_params.n_crossover;
ga_params.gamma = 2; % mutation shape parameter
ga_params.p_mu_genes = 0.2; % mutation probability of a gene

ga_params.mutation = @mutation_dynanic_nonuniform; % mutation_gaussian, mutation_nonuniform, mutation_dynanic_nonuniform
ga_params.crossover = @crossover_intermediate; % crossover_single_point, crossover_uniform, crossover_intermediate


params = genetic_algorithm_CRA(params, ga_params);
save('M_5_alpha_1_Seed_282.mat')

%%

load('M_5_alpha_1_Seed_282.mat')
% calc FIR coefficients
[~, coeff] = calc_proposed_FIR_beamformer(params);

% calc filters frequency response
H = calc_freq_rep(coeff, params);

% find minimum num of mic in each ring
mask = double(H > 0.02);
f_max = [];
for i = 1:params.M
    idx = find(mask(i, :), 1, 'last');
    f_max = [f_max; params.f_grid(idx)];
end
params.Nm = ceil(4*pi*params.Rm.*f_max/params.c);
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

