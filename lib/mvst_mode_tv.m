function md = mvst_mode_tv(n, m, varargin)
%MVST_MODE_TV  Build a TIME-VARYING controlled-mode data structure.
%   MD = MVST_MODE_TV(N,M,'A',@(t)...,'B',...,...) describes one
%   continuous-time controlled diffusion mode whose coefficients may be
%   constant matrices OR function handles of time t:
%
%       dx = (A(t)x + B(t)u) dt + sum_l (Ac_l(t)x + Bc_l(t)u) dW_l,
%
%   with cost [x;u]'*[Q(t) S(t); S(t)' R(t)]*[x;u].  Name/value pairs (all
%   optional, defaulting to the right-size constant zeros / identity):
%       'A','B','Q','R','S'  : constants or handles @(t) of the usual sizes
%       'Ac','Bc'            : constants or handles @(t) returning n x n x Nc
%                              and n x m x Nc arrays (multiplicative noise)
%       'W'                  : n x n additive process-noise covariance
%       'X0'                 : n x n initial state covariance (or handle)
%   The returned struct is consumed by the *_tv solvers (MVST_RIC_CT_TV,
%   MVST_DRE_CT_TV, MVST_SOLVE_FINITE_TV, MVST_SOLVE_INFINITE_TV), which
%   evaluate every field at the current time through MVST_EVALC.
%
%   See also MVST_MODE, MVST_EVALC, MVST_RIC_CT_TV, MVST_SOLVE_FINITE_TV.
p = inputParser;
p.addParameter('A',  zeros(n));
p.addParameter('B',  zeros(n,m));
p.addParameter('Q',  zeros(n));
p.addParameter('R',  eye(m));
p.addParameter('Ac', zeros(n,n,0));
p.addParameter('Bc', zeros(n,m,0));
p.addParameter('S',  zeros(n,m));
p.addParameter('W',  zeros(n));
p.addParameter('X0', eye(n));
p.parse(varargin{:});
o = p.Results;

md.A = o.A; md.B = o.B; md.Q = o.Q; md.R = o.R;
md.Ac = o.Ac; md.Bc = o.Bc; md.S = o.S; md.W = o.W; md.X0 = o.X0;
md.n = n; md.m = m;
md.Nc = size(mvst_evalc(o.Ac,0), 3);          % number of noise channels (at t=0)
md.tv = true;
end
