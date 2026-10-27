function res = mvst_l2gain(sys, opts)
%MVST_L2GAIN  Class-dispatched L2-gain (thin router over the exact backends).
%   RES = MVST_L2GAIN(SYS) computes the (mean-square) L2-gain of the system
%   described by the struct SYS, routing to the matching per-class exact
%   backend and returning a UNIFORM result struct.  The backends are used
%   unchanged; this function only normalises the fields of SYS into the
%   positional arguments each backend expects.
%
%   RES = MVST_L2GAIN(SYS, OPTS) forces the class and passes tuning
%   parameters through OPTS (see below).
%
%   OUTPUT (uniform):
%     RES.gamma   the L2-gain returned by the backend;
%     RES.method  the resolved class tag (one of the strings below);
%     RES.info    backend diagnostics (struct; sd: Nper used, increment,
%                 converged; dtmjls: certificate cell S; otherwise empty).
%
%   SYS FIELDS (drift/input accept both the thin and the rich naming):
%     .A  or .A0    flow drift            (n x n, or n x n x M for multi-mode)
%     .Bw           disturbance input     (n x m)
%     .C            performance output    (nz x n, or nz x n x M)
%     .Dw           feedthrough           (nz x m; default 0)
%   plus, per class:
%     ct      : .Ac state-multiplicative diffusion, .Dwc disturbance-diffusion
%               (when present the MEAN-SQUARE gain is returned; absent = the
%               deterministic H-infinity norm, unchanged)
%     mjls    : .Pig (or .Pi) generator, multi-mode .A0/.C, per-mode .Ac
%               (n x n x q x M) diffusion, shared .Dwc disturbance-diffusion
%     dwell   : .Td dwell-time, multi-mode .A0/.C, per-mode .Ac diffusion,
%               shared .Dwc disturbance-diffusion
%     poisson : .Ac diffusion channels, .Dwc disturbance-diffusion channels,
%               .Jf spontaneous jump maps, .lamf jump rates, .Bdf jump
%               disturbance channels
%     dt      : pure discrete time x+ = A x + Bw w (A Schur): exact bounded-real
%               test, no extra fields (see MVST_L2GAIN_DT)
%     dtmjls  : discrete-time MJLS: multi-mode .A0/.C and a row-stochastic
%               transition matrix .Ptr (or .P) (see MVST_L2GAIN_DT_MJLS)
%     sd      : .Tbar sampling period (hybrid sampled-output gain; no .Dw).
%               BEHAVIOR NOTE: a sampled-data struct carrying .Ac (flow state-
%               multiplicative) and/or .Dwc (flow disturbance-multiplicative)
%               now has that flow diffusion FOLDED into the hybrid Riccati,
%               returning the mean-square hybrid gain; absent = the
%               deterministic hybrid gain, unchanged.
%
%   OPTS FIELDS (all optional):
%     .class   explicit class tag, one of
%              'ct' | 'mjls' | 'dwell'|'switched' | 'poisson' | 'sd'|'hybrid'
%              | 'dt' | 'dtmjls'
%     .Td      dwell-time  (dwell; overrides SYS.Td)
%     .Tbar    sampling period (sd; overrides SYS.Tbar)
%     .Thoriz  shooting horizon (mjls, poisson)
%     .Nper    number of periods (sd); ABSENT/EMPTY = adaptive with a convergence
%              check (see MVST_L2GAIN_SD_HYBRID); a scalar fixes the horizon
%     .Ntol    relative convergence tolerance of the adaptive Nper (sd)
%     .solver  SDP backend name, carried coherently for callers that pass a
%              shared options struct (or an MVST.OPTIONS); the exact backends
%              below are purely numerical, so the field is normalized (via
%              MVST_OPTS, without a path probe) and does not affect the result
%
%   CLASS INFERENCE (when OPTS.class is absent; first match wins):
%     0. SYS.domain = 'dt' (or SYS.dt = true)    -> 'dt' (single mode) or
%        'dtmjls' (multi-mode)
%     1. .Jf or .lamf present                    -> 'poisson'
%     2. .Tbar present and single-mode           -> 'sd'
%     3. multi-mode (.A0 3-D / cell) and .Pig|.Pi -> 'mjls'
%     4. multi-mode and .Td (SYS or OPTS)        -> 'dwell'
%     5. otherwise                               -> 'ct'
%
%   Backends: MVST_L2GAIN_DT, MVST_L2GAIN_DT_MJLS, MVST_L2GAIN_CT, MVST_L2GAIN_CT_MJLS, MVST_L2GAIN_CT_DWELL,
%   MVST_L2GAIN_CT_POISSON, MVST_L2GAIN_SD_HYBRID.
%
%   See also MVST_GAME_GAMMA.

if nargin < 2 || isempty(opts); opts = struct(); end
opts = mvst_opts(opts, false);   % accept mvst.Options / .solver; no path probe

A  = getfld(sys, {'A0','A'}, []);
Bw = getfld(sys, {'Bw'}, []);
C  = getfld(sys, {'C'},  []);
if isempty(A);  error('mvst:l2gain:noA',  'SYS must carry a drift .A0 (or .A).'); end
if isempty(Bw); error('mvst:l2gain:noBw', 'SYS must carry a disturbance input .Bw.'); end

cls = resolve_class(sys, opts);

m  = size(Bw, 2);
nz = maxdim(C, 1);
Dw = getfld(sys, {'Dw'}, zeros(nz, m));

info = struct();
switch cls
    case 'dt'
        gam = mvst_l2gain_dt(A(:,:,1), Bw, C(:,:,1), Dw);

    case 'dtmjls'
        Ptr = getfld(sys, {'Ptr','P'}, []);
        if isempty(Ptr)
            error('mvst:l2gain:dtmjls', 'DT-MJLS class needs a transition matrix .Ptr (or .P).');
        end
        Nit = []; if isfield(opts,'Thoriz') && ~isempty(opts.Thoriz); Nit = opts.Thoriz; end
        Ad = tocell3(A); Cd = tocell3(C);
        if numel(Cd) == 1; Cd = repmat(Cd, 1, numel(Ad)); end    % common output
        [gam, Scert] = mvst_l2gain_dt_mjls(Ad, Bw, Cd, Dw, Ptr, Nit);
        info.S = Scert;

    case 'ct'
        Amat = A(:,:,1);
        Cn  = tocell3(getfld(sys, {'Ac'},  {}));   % state-mult diffusion
        Dwc = tocell3(getfld(sys, {'Dwc'}, {}));   % disturbance-mult diffusion
        gam = mvst_l2gain_ct(Amat, Bw, C(:,:,1), Dw, Cn, Dwc);

    case 'mjls'
        Pig = getfld(sys, {'Pig','Pi'}, []);
        if isempty(Pig)
            error('mvst:l2gain:mjls', 'MJLS class needs a generator .Pig (or .Pi).');
        end
        Ac = tocell3(A); Cc = tocell3(C);
        if numel(Cc) == 1; Cc = repmat(Cc, 1, numel(Ac)); end   % common output
        Chan = permodechan(getfld(sys, {'Ac'}, []), numel(Ac));   % per-mode diffusion
        Dwc  = tocell3(getfld(sys, {'Dwc'}, {}));                 % shared dist-diffusion
        if isfield(opts,'Thoriz') && ~isempty(opts.Thoriz)
            gam = mvst_l2gain_ct_mjls(Ac, Bw, Cc, Dw, Pig, opts.Thoriz, Chan, Dwc);
        else
            gam = mvst_l2gain_ct_mjls(Ac, Bw, Cc, Dw, Pig, [], Chan, Dwc);
        end

    case 'dwell'
        Td = getopt(opts, 'Td', getfld(sys, {'Td'}, []));
        if isempty(Td)
            error('mvst:l2gain:dwell', 'Dwell class needs a dwell-time (.Td).');
        end
        Ac = tocell3(A); Cc = tocell3(C);
        if numel(Cc) == 1; Cc = repmat(Cc, 1, numel(Ac)); end   % common output
        Chan = permodechan(getfld(sys, {'Ac'}, []), numel(Ac));   % per-mode diffusion
        Dwc  = tocell3(getfld(sys, {'Dwc'}, {}));                 % shared dist-diffusion
        gam = mvst_l2gain_ct_dwell(Ac, Bw, Cc, Dw, Td, Chan, Dwc);

    case 'poisson'
        Cc  = tocell3(getfld(sys, {'Ac'},  {}));            % diffusion channels
        Dwc = tocell3(getfld(sys, {'Dwc'}, {}));            % disturbance-diffusion
        Jc  = tocell3(getfld(sys, {'Jf'},  {}));            % spontaneous jump maps
        lam = getfld(sys, {'lamf'}, []);                    % jump rates
        Dr  = tocell3(getfld(sys, {'Bdf'}, {}));            % jump-disturbance
        n = size(A,1);
        % Match the disturbance-diffusion cell to the diffusion channels.
        if numel(Dwc) < numel(Cc)
            Dwc(end+1:numel(Cc)) = {zeros(n, m)};
        end
        % Match the jump-disturbance cell to the jump maps.
        if numel(Dr) < numel(Jc)
            Dr(end+1:numel(Jc)) = {zeros(n, m)};
        end
        if isfield(opts,'Thoriz') && ~isempty(opts.Thoriz)
            gam = mvst_l2gain_ct_poisson(A(:,:,1), Bw, C(:,:,1), Dw, Cc, Dwc, Jc, lam, Dr, opts.Thoriz);
        else
            gam = mvst_l2gain_ct_poisson(A(:,:,1), Bw, C(:,:,1), Dw, Cc, Dwc, Jc, lam, Dr);
        end

    case 'sd'
        Tbar = getopt(opts, 'Tbar', getfld(sys, {'Tbar'}, []));
        if isempty(Tbar)
            error('mvst:l2gain:sd', 'Sampled-data class needs a period (.Tbar).');
        end
        Cn  = tocell3(getfld(sys, {'Ac'},  {}));   % flow state-mult diffusion
        Dwc = tocell3(getfld(sys, {'Dwc'}, {}));   % flow disturbance-mult diffusion
        Nper = getopt(opts, 'Nper', []);
        Ntol = getopt(opts, 'Ntol', []);
        [gam, info] = mvst_l2gain_sd_hybrid(A(:,:,1), Bw, C(:,:,1), Tbar, Nper, Cn, Dwc, [], Ntol);

    otherwise
        error('mvst:l2gain:class', 'Unknown class tag ''%s''.', cls);
end

res = struct('gamma', gam, 'method', cls, 'info', info);
end

% ------------------------------------------------------------------------
function cls = resolve_class(sys, opts)
if isfield(opts, 'class') && ~isempty(opts.class)
    cls = lower(opts.class);
    if strcmp(cls,'switched'); cls = 'dwell'; end
    if strcmp(cls,'hybrid');   cls = 'sd';    end
    return;
end
multi = ismode_multi(sys);
isdt = (isfield(sys,'domain') && strcmpi(sys.domain,'dt')) || ...
       (isfield(sys,'dt') && isequal(sys.dt, true));
if isdt
    if multi; cls = 'dtmjls'; else; cls = 'dt'; end
elseif isfield(sys,'Jf') || isfield(sys,'lamf')
    cls = 'poisson';
elseif isfield(sys,'Tbar') && ~multi
    cls = 'sd';
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
function d = maxdim(X, k)
if isempty(X); d = 0; elseif iscell(X); d = size(X{1}, k); else; d = size(X, k); end
end

% ------------------------------------------------------------------------
function c = tocell3(X)
%TOCELL3  Slice a 3-D array along its 3rd dimension into a row cell.
%   A cell is returned unchanged; a 2-D matrix becomes a 1-element cell; an
%   empty input becomes an empty cell.
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
%   AC4 is n x n x q x M (rich multi-mode state noise, as MVST_TORICH stacks
%   it); slice out mode i's q channels.  An n x n x q array (shared across
%   modes) is broadcast to every mode; empty gives all-empty per-mode cells.
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
