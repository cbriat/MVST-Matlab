function md = mvst_mode(A, B, Q, R, varargin)
%MVST_MODE  Build a (normalized) controlled-mode data structure.
%   MD = MVST_MODE(A,B,Q,R) describes one continuous-time controlled
%   diffusion mode
%
%       dx = (A x + B u) dt + sum_l (Ac(:,:,l) x + Bc(:,:,l) u) dW_l,
%
%   with quadratic cost  [x;u]'*[Q S; S' R]*[x;u].  Optional name/value
%   pairs (all default to the right-size zeros):
%       'Ac' : n x n x Nc array of state noise coefficients
%       'Bc' : n x m x Nc array of input noise coefficients
%       'S'  : n x m cross-cost block
%       'W'  : n x n additive process-noise covariance (average-cost case)
%       'X0' : n x n initial state covariance E[x x']
%   The returned struct is consumed by MVST_RIC_CT and all solvers.
%
%   UNIFIED-MODEL EXTENSIONS (optional; see lib/MODEL_CONVENTION.md).  The
%   following jump/reset fields let the SAME mode struct describe the
%   impulsive / Poisson / discrete classes handled by the SOS layer
%   (MVST_STAB_SOS, MVST_STABZ_SOS, MVST_DISS_SOS, ...).  They map onto the
%   rich SOS naming through MVST_TORICH (A -> A0, B -> B0):
%       'J'    : n x n [x p]  SCHEDULED (dwell-time) state-jump / reset maps
%       'BJ'   : n x m [x p]  scheduled-jump input selectors (paired with J)
%       'Jf'   : n x n x pf   SPONTANEOUS (Poisson) jump maps
%       'lamf' : 1 x pf       Poisson rates (paired with Jf)
%       'BJf'  : n x m x pf   spontaneous-jump input selectors (paired with Jf)
%   Every extension DEFAULTS TO ABSENT: when it is not supplied it is not
%   added to the struct, so a mode built without these arguments is
%   byte-identical to the pre-extension output and every existing caller
%   (MVST_RIC_CT and the Riccati solvers, which read fields by name) is
%   unaffected.  The Riccati layer ignores any extension it does not read.
%
%   See also MVST_SYSTEM, MVST_MODE_TV, MVST_TORICH, MVST_TOTHIN, MVST_RIC_CT.

if ~(ismatrix(A) && size(A,1) == size(A,2))
    error('mvst:mode:notSquare', 'A must be a square matrix.');
end
if size(B,1) ~= size(A,1)
    error('mvst:mode:Bdim', 'B must have %d rows to match A.', size(A,1));
end
if ~(ismatrix(Q) && size(Q,1) == size(Q,2))
    error('mvst:mode:Qdim', 'Q must be a square matrix.');
end
if ~(ismatrix(R) && size(R,1) == size(R,2))
    error('mvst:mode:Rdim', 'R must be a square matrix.');
end

p = inputParser;
p.addParameter('Ac', []); p.addParameter('Bc', []);
p.addParameter('S',  []); p.addParameter('W',  []);
p.addParameter('X0', []);
p.addParameter('J',   []); p.addParameter('BJ',  []);
p.addParameter('Jf',  []); p.addParameter('lamf',[]); p.addParameter('BJf', []);
p.parse(varargin{:});
opt = p.Results;

n = size(A,1); m = size(B,2);
md.A = A;  md.B = B;  md.Q = (Q+Q.')/2;  md.R = (R+R.')/2;
md.n = n;  md.m = m;

if isempty(opt.Ac); md.Ac = zeros(n,n,0); else; md.Ac = opt.Ac; end
if isempty(opt.Bc); md.Bc = zeros(n,m,0); else; md.Bc = opt.Bc; end
% One channel array may be omitted: pad the ABSENT one with zeros to the
% channel count of the one that was supplied (e.g. supplying only 'Ac' for a
% state-noise-only diffusion; the input-noise Bc is then zero).  Only both
% present with DIFFERING counts is an error.
nAc = size(md.Ac,3); nBc = size(md.Bc,3);
if nAc > 0 && nBc == 0
    md.Bc = zeros(n,m,nAc);
elseif nBc > 0 && nAc == 0
    md.Ac = zeros(n,n,nBc);
end
md.Nc = size(md.Ac,3);
if size(md.Bc,3) ~= md.Nc
    error('mvst:mode:noise','Ac and Bc must have the same number of channels.');
end
if isempty(opt.S);  md.S  = zeros(n,m); else; md.S  = opt.S;  end
if isempty(opt.W);  md.W  = zeros(n);   else; md.W  = (opt.W+opt.W.')/2; end
if isempty(opt.X0); md.X0 = eye(n);     else; md.X0 = (opt.X0+opt.X0.')/2; end

% Optional jump/reset extensions: attach ONLY when supplied, so the struct
% stays byte-identical for callers that do not use them.
if ~isempty(opt.J);    md.J    = opt.J;    end
if ~isempty(opt.BJ);   md.BJ   = opt.BJ;   end
if ~isempty(opt.Jf);   md.Jf   = opt.Jf;   end
if ~isempty(opt.lamf); md.lamf = opt.lamf; end
if ~isempty(opt.BJf);  md.BJf  = opt.BJf;  end
end
