function sys = mvst_system(modes, varargin)
%MVST_SYSTEM  Assemble a (possibly multi-mode) system for the solvers.
%   SYS = MVST_SYSTEM(MD) wraps a single mode MD (from MVST_MODE).
%   SYS = MVST_SYSTEM({MD1,...,MDM}, 'Pi',PI, 'pi0',P0, 'ZT',ZT) assembles a
%   multi-mode (Markov-jump / switched) system, where
%       Pi  : M x M Markov generator coupling the modes (default 0),
%       pi0 : 1 x M initial mode distribution (default e_1),
%       ZT  : terminal weight for the finite-horizon solver (default 0).
%
%   UNIFIED-MODEL EXTENSIONS (optional; see lib/MODEL_CONVENTION.md).  For
%   the SOS layer (MVST_STAB_SOS_MODES, MVST_DISS_SOS_MODES) two more
%   coupling data may be supplied and are carried through MVST_TORICH:
%       'P'  : M x M discrete-time Markov TRANSITION matrix (rows sum to 1),
%              for the DT-MJLS class (SYS.Pi is the CT generator);
%       'J'  : n x n common switch / impulse RESET applied at a mode change.
%   Both DEFAULT TO ABSENT (not added to the struct when omitted), so a
%   system assembled without them is byte-identical to the pre-extension
%   output and the Riccati solvers, which read only .modes/.Pi/.pi0/.ZT,
%   are unaffected.
%
%   See also MVST_MODE, MVST_TORICH, MVST_TOTHIN, MVST_SOLVE_INFINITE.
if ~iscell(modes); modes = {modes}; end
M = numel(modes);
if M < 1
    error('mvst:system:noModes', 'At least one mode is required.');
end
for i = 1:M
    if ~isstruct(modes{i})
        error('mvst:system:mode', ...
            'Mode %d must be a struct (as built by MVST_MODE).', i);
    end
end
if ~all(isfield(modes{1}, {'n','m'}))
    error('mvst:system:mode', ...
        'Mode structs must carry .n and .m (build them with MVST_MODE).');
end
p = inputParser;
p.addParameter('Pi',  zeros(M));
p.addParameter('pi0', [1 zeros(1,M-1)]);
p.addParameter('ZT',  []);
p.addParameter('P',   []);
p.addParameter('J',   []);
p.parse(varargin{:});
if ~isequal(size(p.Results.Pi), [M M])
    error('mvst:system:generator', 'Pi must be %d-by-%d for %d modes.', M, M, M);
end
if numel(p.Results.pi0) ~= M
    error('mvst:system:pi0', 'pi0 must have %d entries for %d modes.', M, M);
end
sys.modes = modes;
sys.Pi    = p.Results.Pi;
sys.pi0   = p.Results.pi0;
if isempty(p.Results.ZT); sys.ZT = zeros(modes{1}.n); else; sys.ZT = p.Results.ZT; end
sys.M = M; sys.n = modes{1}.n; sys.m = modes{1}.m;
% Optional SOS-layer coupling extensions: attach ONLY when supplied.
if ~isempty(p.Results.P); sys.P = p.Results.P; end
if ~isempty(p.Results.J); sys.J = p.Results.J; end
end
