function rich = mvst_torich(sys)
%MVST_TORICH  Convert a thin mode/system struct to the rich SOS-layer struct.
%   RICH = MVST_TORICH(SYS) maps the thin OC/Riccati data convention
%   (MVST_MODE / MVST_SYSTEM, drift .A, input .B, noise .Ac/.Bc, plus the
%   optional jump/reset extensions .J/.BJ/.Jf/.lamf/.BJf) onto the rich
%   struct consumed by the SOS layer (MVST_STAB_SOS, MVST_STABZ_SOS,
%   MVST_DISS_SOS and their *_MODES counterparts), whose naming uses .A0 for
%   the drift and .B0 for the input.  See lib/MODEL_CONVENTION.md.
%
%   The routine is IDEMPOTENT and permissive:
%     * a struct that already carries .A0 is returned unchanged (already
%       rich), so the SOS entry points may call MVST_TORICH unconditionally;
%     * a thin single mode (has .A, from MVST_MODE) becomes a single-mode
%       rich struct: A->A0, B->B0, and .Ac/.Bc/.J/.BJ/.Jf/.lamf/.BJf plus
%       any dissipativity I/O fields (.Bw/.Dwc/.C/.Dw/.Bdf/.Bd/.CJ/.DJ)
%       already present are copied through verbatim;
%     * a thin multi-mode system (has .modes, from MVST_SYSTEM) becomes a
%       mode-indexed rich struct: A0 is n x n x M, Ac is n x n x q x M
%       (channels zero-padded to the common maximum q), B0 / Bc are stacked
%       likewise, the CT generator .Pi is exposed under BOTH names .Pi
%       (MVST_STAB_SOS_MODES) and .Pig (MVST_DISS_SOS_MODES), and the
%       optional .P (DT transition) / .J (reset) are carried through.
%
%   Cost data (.Q/.R/.S/.W/.X0/.ZT) has no SOS-layer counterpart and is not
%   copied.  Dissipativity I/O data (disturbance / performance-output blocks)
%   has no home in a thin mode built by MVST_MODE, but is passed through when
%   the caller has attached it to the input struct.
%
%   See also MVST_TOTHIN, MVST_MODE, MVST_SYSTEM, MVST_STAB_SOS.

if ~isstruct(sys)
    error('mvst:torich:input','SYS must be a struct.');
end

% Already rich: nothing to do (idempotent).
if isfield(sys,'A0')
    rich = sys; return;
end

if isfield(sys,'modes')                       % thin multi-mode system
    rich = system_to_rich(sys);
elseif isfield(sys,'A')                        % thin single mode
    rich = mode_to_rich(sys);
else
    error('mvst:torich:input', ...
        ['SYS is neither rich (no .A0), a thin mode (no .A), nor a thin ', ...
         'system (no .modes).']);
end
end

% =========================================================================
function rich = mode_to_rich(md)
%MODE_TO_RICH  Single thin mode -> single-mode rich struct.
rich.A0 = md.A;
if isfield(md,'B') && ~isempty(md.B); rich.B0 = md.B; end
% Noise + jump/reset + dissipativity passthrough (same field names on both
% sides except A/B); copy each only when present and non-empty.
passthru = {'Ac','Bc','J','BJ','Jf','lamf','BJf', ...
            'Bw','Dwc','C','Dw','Bdf','Bd','CJ','DJ'};
for f = passthru
    if isfield(md,f{1}) && ~isempty(md.(f{1})); rich.(f{1}) = md.(f{1}); end
end
end

% =========================================================================
function rich = system_to_rich(sys)
%SYSTEM_TO_RICH  Thin multi-mode system -> mode-indexed rich struct.
modes = sys.modes; M = numel(modes);
n = size(modes{1}.A,1);
m = size(modes{1}.B,2);

% --- drift A0 : n x n x M
A0 = zeros(n,n,M);
for i = 1:M; A0(:,:,i) = modes{i}.A; end
rich.A0 = A0;

% --- state noise Ac : n x n x q x M (zero-pad channels to common max q)
q = 0;
for i = 1:M
    if isfield(modes{i},'Ac'); q = max(q, size(modes{i}.Ac,3)); end
end
if q > 0
    Ac = zeros(n,n,q,M);
    for i = 1:M
        if isfield(modes{i},'Ac')
            qi = size(modes{i}.Ac,3);
            Ac(:,:,1:qi,i) = modes{i}.Ac;
        end
    end
    rich.Ac = Ac;
end

% --- input B0 : n x m x M (only if any mode carries a non-empty B)
if m > 0
    B0 = zeros(n,m,M);
    for i = 1:M; B0(:,:,i) = modes{i}.B; end
    rich.B0 = B0;
    if q > 0                                   % input noise Bc : n x m x q x M
        Bc = zeros(n,m,q,M);
        for i = 1:M
            if isfield(modes{i},'Bc')
                qi = size(modes{i}.Bc,3);
                Bc(:,:,1:qi,i) = modes{i}.Bc;
            end
        end
        rich.Bc = Bc;
    end
end

% --- mode coupling: expose the CT generator under both consumer names.
if isfield(sys,'Pi') && ~isempty(sys.Pi)
    rich.Pi  = sys.Pi;      % MVST_STAB_SOS_MODES (generator / MSI coupling)
    rich.Pig = sys.Pi;      % MVST_DISS_SOS_MODES (generator)
end
if isfield(sys,'P') && ~isempty(sys.P); rich.P = sys.P; end   % DT transition
if isfield(sys,'J') && ~isempty(sys.J); rich.J = sys.J; end   % switch reset
end
