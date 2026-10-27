function res = mvst_game_gamma(sys, opts)
%MVST_GAME_GAMMA  Class-dispatched H-infinity game level (thin router).
%   RES = MVST_GAME_GAMMA(SYS) computes the critical disturbance level
%   gamma* of the state-feedback zero-sum H-infinity game described by the
%   struct SYS, routing to the matching per-class exact backend and
%   returning a UNIFORM result struct.  The backends are used unchanged;
%   this function only normalises the fields of SYS into the positional
%   arguments each backend expects.
%
%   RES = MVST_GAME_GAMMA(SYS, OPTS) forces the class and passes tuning
%   parameters through OPTS (see below).
%
%   OUTPUT (uniform):
%     RES.gammastar  the critical level gamma* returned by the backend;
%     RES.P          the settled dual(s): a matrix (ct, sd, poisson) or a
%                    1 x M cell of per-mode duals (mjls); absent (dwell);
%     RES.gains      the certified gains cell (dwell only);
%     RES.method     the resolved class tag.
%
%   SYS FIELDS (drift accepts both the thin and the rich naming):
%     .A  or .A0    flow drift  (n x n, or n x n x M; may be a function
%                                handle @(t) for the ct/sd time-varying case)
%     .Bu           control input
%     .Bw           disturbance input
%     .Q            state weight
%     .Ru           control weight
%   plus, per class:
%     mjls    : .Pig (or .Pi) generator, multi-mode .A0
%     dwell   : .Td dwell-time, multi-mode .A0
%     poisson : .Ac diffusion channels, .Jf spontaneous jump maps,
%               .lamf jump rates
%
%   OPTS FIELDS (all optional):
%     .class   explicit class tag, one of
%              'ct'|'sd' | 'mjls' | 'dwell'|'switched' | 'poisson'
%     .Tlong   long horizon for the boundedness shooting (ct, mjls, poisson)
%     .Td      dwell-time (dwell; overrides SYS.Td)
%     .solver  SDP backend name, carried coherently for callers that pass a
%              shared options struct (or an MVST.OPTIONS); the exact backends
%              below are purely numerical, so the field is normalized (via
%              MVST_OPTS, without a path probe) and does not affect the result
%
%   CLASS INFERENCE (when OPTS.class is absent; first match wins):
%     1. .Jf or .lamf present                    -> 'poisson'
%     2. multi-mode (.A0 3-D / cell) and .Pig|.Pi -> 'mjls'
%     3. multi-mode and .Td (SYS or OPTS)        -> 'dwell'
%     4. otherwise                               -> 'ct'
%   (The sampled-data game shares the 'ct' backend; use OPTS.class='ct'.)
%
%   Backends: MVST_GAME_GAMMA_CT, MVST_GAME_GAMMA_CT_MJLS,
%   MVST_GAME_GAMMA_CT_DWELL, MVST_GAME_GAMMA_CT_POISSON.
%
%   See also MVST_L2GAIN.

if nargin < 2 || isempty(opts); opts = struct(); end
opts = mvst_opts(opts, false);   % accept mvst.Options / .solver; no path probe

A  = getfld(sys, {'A0','A'}, []);
Bu = getfld(sys, {'Bu'}, []);
Bw = getfld(sys, {'Bw'}, []);
Q  = getfld(sys, {'Q'},  []);
Ru = getfld(sys, {'Ru'}, []);
if isempty(A);  error('mvst:game_gamma:noA',  'SYS must carry a drift .A0 (or .A).'); end
if isempty(Bu); error('mvst:game_gamma:missing', 'SYS must carry .Bu.'); end
if isempty(Bw); error('mvst:game_gamma:missing', 'SYS must carry .Bw.'); end
if isempty(Q);  error('mvst:game_gamma:missing', 'SYS must carry .Q.');  end
if isempty(Ru); error('mvst:game_gamma:missing', 'SYS must carry .Ru.'); end

cls = resolve_class(sys, opts);
Tlong = getopt(opts, 'Tlong', []);

res = struct('gammastar', [], 'method', cls);

Cn  = tocell3(getfld(sys, {'Ac'},  {}));   % state-mult diffusion (single mode)
Dwc = tocell3(getfld(sys, {'Dwc'}, {}));   % disturbance-mult diffusion (shared)

switch cls
    case 'ct'
        Afun = A; if ~isa(A,'function_handle'); Afun = A(:,:,1); end
        [gstar, Pinf] = mvst_game_gamma_ct(Afun, Bu, Bw, Q, Ru, Tlong, Cn, Dwc);
        res.gammastar = gstar; res.P = Pinf;

    case 'mjls'
        Pig = getfld(sys, {'Pig','Pi'}, []);
        if isempty(Pig)
            error('mvst:game_gamma:mjls', 'MJLS class needs a generator .Pig (or .Pi).');
        end
        Ac = tocell3(A);
        Chan = permodechan(getfld(sys, {'Ac'}, []), numel(Ac));   % per-mode diffusion
        if isempty(Tlong); Tlong = []; end
        [gstar, P0] = mvst_game_gamma_ct_mjls(Ac, Bu, Bw, Q, Ru, Pig, Tlong, Chan, Dwc);
        res.gammastar = gstar; res.P = P0;

    case 'dwell'
        Td = getopt(opts, 'Td', getfld(sys, {'Td'}, []));
        if isempty(Td)
            error('mvst:game_gamma:dwell', 'Dwell class needs a dwell-time (.Td).');
        end
        Ac = tocell3(A);
        Chan = permodechan(getfld(sys, {'Ac'}, []), numel(Ac));   % per-mode diffusion
        [gam, Kc] = mvst_game_gamma_ct_dwell(Ac, Bu, Bw, Q, Ru, Td, Chan, Dwc);
        res.gammastar = gam; res.gains = Kc;

    case 'poisson'
        Cc  = tocell3(getfld(sys, {'Ac'}, {}));    % diffusion channels
        Jc  = tocell3(getfld(sys, {'Jf'}, {}));    % spontaneous jump maps
        lam = getfld(sys, {'lamf'}, []);           % jump rates
        if isempty(Tlong)
            [gstar, Pinf] = mvst_game_gamma_ct_poisson(A(:,:,1), Bu, Bw, Cc, Jc, lam, Q, Ru);
        else
            [gstar, Pinf] = mvst_game_gamma_ct_poisson(A(:,:,1), Bu, Bw, Cc, Jc, lam, Q, Ru, Tlong);
        end
        res.gammastar = gstar; res.P = Pinf;

    otherwise
        error('mvst:game_gamma:class', 'Unknown class tag ''%s''.', cls);
end
end

% ------------------------------------------------------------------------
function cls = resolve_class(sys, opts)
if isfield(opts, 'class') && ~isempty(opts.class)
    cls = lower(opts.class);
    if strcmp(cls,'switched'); cls = 'dwell'; end
    if strcmp(cls,'sd');       cls = 'ct';    end
    return;
end
multi = ismode_multi(sys);
if isfield(sys,'Jf') || isfield(sys,'lamf')
    cls = 'poisson';
elseif multi && (isfield(sys,'Pig') || isfield(sys,'Pi'))
    cls = 'mjls';
elseif multi && (isfield(sys,'Td') || (isfield(opts,'Td') && ~isempty(opts.Td)))
    cls = 'dwell';
else
    cls = 'ct';
end
end

% ------------------------------------------------------------------------
function tf = ismode_multi(sys)
A = getfld(sys, {'A0','A'}, []);
if isa(A,'function_handle'); tf = false; return; end
tf = (iscell(A) && numel(A) > 1) || (~isempty(A) && size(A,3) > 1) || isfield(sys,'modes');
end

% ------------------------------------------------------------------------
function v = getfld(s, names, default)
v = default;
for k = 1:numel(names)
    if isfield(s, names{k}) && ~isempty(s.(names{k}))
        v = s.(names{k}); return;
    end
end
end

% ------------------------------------------------------------------------
function v = getopt(opts, name, default)
if isfield(opts, name) && ~isempty(opts.(name)); v = opts.(name); else; v = default; end
end

% ------------------------------------------------------------------------
function c = tocell3(X)
%TOCELL3  Slice a 3-D array along its 3rd dimension into a row cell.
if iscell(X)
    c = X;
elseif isempty(X)
    c = {};
else
    p = size(X, 3);
    c = cell(1, p);
    for k = 1:p; c{k} = X(:,:,k); end
end
end

% ------------------------------------------------------------------------
function C = permodechan(Ac4, M)
%PERMODECHAN  Multi-mode diffusion array -> 1xM cell of 1xq channel cells.
if isempty(Ac4); C = repmat({{}}, 1, M); return; end
q = size(Ac4, 3); Mi = size(Ac4, 4);
C = cell(1, M);
for i = 1:M
    idx = i; if Mi == 1; idx = 1; end
    ci = cell(1, q);
    for l = 1:q; ci{l} = Ac4(:, :, l, idx); end
    C{i} = ci;
end
end
