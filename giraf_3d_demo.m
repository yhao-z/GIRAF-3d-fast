clear;
clc;
%% load data
addpath('./data','./etc','./fft_selfdefined');
load acc=10_192_radialray.mat
load data25.mat
x_origin=data;
% x0=fftn(x_origin)/numel(x_origin);
x0=fft3(x_origin);
sampmask=mask0;
res=size(x0);
%% define index sets
ind_samples = find(sampmask);
[S,St] = defSSt(ind_samples,res);
b = S(ifft_t(x0));
xinit = fft_t(St(b));
%% global problem settings
settings.filter_siz = [15 15 7];
settings.res = size(x0);
settings.exit_tol = 5e-4;  %exit algorithm if relative change in NRMSE between iterates less than this
settings.lambda = 0.01; %regularization parameter (0=equality constrained)0.01
settings.p = 0; %Schatten p penalty (0 <= p <= 1)
%% GIRAF parameters
param.iter = 30; %number of IRLS iterations
param.eps0 = 0; %inital epsilon (0=auto-init) 
param.eta = 15; %epsilon decrease factor (typically between 1.1-1.5);
param.epsmin = 1e-7;
param.ADMM_iter = 30;
param.ADMM_tol = 1e-9;
param.delta1 = 25; %ADMM conditioning parameter (typically between 10-100);
param.delta2 = 25; %ADMM conditioning parameter (typically between 10-100);
%% run GIRAF
[x,cost] = giraf_3d(x_origin,xinit,b,S,St,sampmask,param,settings);
x_data = ifft3(x);
SNR = -20*log10(norm(x_data(:)-x_origin(:))/norm(x_origin(:)));