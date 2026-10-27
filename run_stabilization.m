%% RUN_STABILIZATION  Figures and data for the Stabilization chapter.
%  Feedback stabilization instances (CT, DT, impulsive), the SOS/grid
%  minimum-dwell-time synthesis, and the guaranteed-rate (alpha,rho)
%  dwell-time synthesis trade-off.
%  Outputs: vector PDFs to figures/, data macros to figures/data/.
%  Runtime: order of a couple of minutes (the SOS dwell-time synthesis and the
%  sample-path simulations dominate; SOS needs an SDP solver on the path).
close all
% requires mvst_setup (library and solvers on the path)
mvst_runscripts({ ...
    'ex_stabz_ct', ...          % continuous-time: Brownian+Poissonian
    'ex_stabz_levy', ...        % continuous-time: Levy noise
    'ex_stabz_delay', ...       % continuous-time: diffusions with delays
    'ex_stabz_mjls', ...        % continuous-time: Markov jump linear systems
    'ex_stabz_poissample', ...  % continuous-time: Poissonian sampling
    'ex_stabz_dt', ...          % discrete-time: random parameters
    'ex_stabz_dtdelay', ...     % discrete-time: delays
    'ex_stabz_dtmjls', ...      % discrete-time: Markov jump
    'ex_stabz_imp', ...         % impulsive: deterministic impulses
    'ex_stabz_switched', ...    % impulsive: switched systems
    'ex_stabz_sd', ...          % impulsive: sampled-data (periodic)
    'ex_stabz_sdrange', ...     % impulsive: range sampled-data, scheduled vs constant refresh
    'ex_stabz_msi', ...         % impulsive: stochastic impulses and switching
    'ex_stabz_ncs', ...         % impulsive: networked control, packet dropouts
    'ex_stabz_sos', ...         % dwell-time synthesis via SOSTOOLS (kept)
    'ex_stabz_sd_sos', ...      % sampled-data synthesis via SOSTOOLS
    'ex_stabz_switched_sos', ... % switched minimum dwell-time synthesis via SOSTOOLS
    'ex_stabz_rate'});          % guaranteed hybrid rate
