function mvst_stabtable(figdir, key)
%MVST_STABTABLE  Emit one row of the Stabilization-chapter summary table.
%   MVST_STABTABLE(FIGDIR, KEY) writes (or refreshes) a single row macro
%   \StabzRow<Key> into FIGDIR/data/stabz_summary.tex, the machine-generated
%   counterpart of the hand-built Table~\ref{tab:stabz:summary}.  KEY is the
%   same short identifier the example script passes to MVST_STABFIG (one of
%   ct, levy, delay, mjls, poissample, dt, dtdelay, dtmjls, imp, switched,
%   sd, msi, ncs).
%
%   The row body references the per-example data macros (\Stabz*K,
%   \Stabz*Open, \Stabz*Closed, ...) already emitted by the same scripts, so
%   the generated file, once \input alongside those data files, reproduces
%   the table rows from live run results.  Column layout matches the book
%   table: class & synthesized gain & measure & open$\to$closed & dwell/margin.
%
%   The function is idempotent and order-independent: each call rewrites the
%   file keeping every row previously registered (canonical order below), so
%   the thirteen scripts may be run in any order or repeatedly.
%
%   The book does NOT yet switch over to this file; a source comment in the
%   hand-built table points here for the eventual switch.
%
%   See also MVST_STABFIG, MVST_TEXDATA.

outfile = fullfile(figdir, 'data', 'stabz_summary.tex');

% ------- canonical row registry (key -> {macroSuffix, rowBody}) ------------
% rowBody carries the trailing \\ so each macro is a drop-in tabular row.
reg = { ...
 'ct',         'Ct',    'Brownian and Poisson (CT) & $\StabzCtK$ & $\alpha$ & $\StabzCtOpen\to\StabzCtClosed$ & --\\'; ...
 'levy',       'Levy',  'L\''evy noise (CT) & $\StabzLevyK$ & $\alpha$ & $\StabzLevyOpen\to\StabzLevyClosed$ & --\\'; ...
 'delay',      'Delay', 'Delayed diffusion (CT) & $\StabzDelayK$ & $\alpha$ & $\StabzDelayOpen\to\StabzDelayClosed$ & delay-indep.\\'; ...
 'mjls',       'Mjls',  'Markov jump (CT) & \shortstack[l]{$K_1=\StabzMjlsKone$\\$K_2=\StabzMjlsKtwo$} & $\alpha$ & $\StabzMjlsOpen\to\StabzMjlsClosed$ & --\\'; ...
 'poissample', 'Pois',  'Poissonian sampling (CT) & $\StabzPoisK$ & $\alpha$ & $\StabzPoisOpen\to\StabzPoisClosed$ & $\lambda=\StabzPoisLam$\\'; ...
 'dt',         'Dt',    'Random parameters (DT) & $\StabzDtK$ & $\rho$ & $\StabzDtOpen\to\StabzDtClosed$ & --\\'; ...
 'dtdelay',    'Dtd',   'Delayed (DT) & $\StabzDtdK$ & $\rho$ & $\StabzDtdOpen\to\StabzDtdClosed$ & delay-indep.\\'; ...
 'dtmjls',     'Dtm',   'Markov jump (DT) & \shortstack[l]{$K_1=\StabzDtmKone$\\[2pt]$K_2=\StabzDtmKtwo$} & $\rho$ & $\StabzDtmOpen\to\StabzDtmClosed$ & --\\'; ...
 'imp',        'Imp',   'Deterministic impulses & \shortstack[l]{$K(0)=\StabzImpKzero$\\$K(\bar T)=\StabzImpKTbar$} & $\rho$ & $\StabzImpOpen\to\StabzImpClosed$ & $\bar T=\StabzImpTbar$\\'; ...
 'switched',   'Sw',    'Switched & \shortstack[l]{$K_1=\StabzSwKone$\\$K_2=\StabzSwKtwo$} & $\rho$ & $\StabzSwOpen\to\StabzSwClosed$ & $T_{\mathrm{sw}}=\StabzSwTsw$\\'; ...
 'sd',         'Sd',    'Sampled-data & $\StabzSdK$ & $\rho$ & $\StabzSdOpen\to\StabzSdClosed$ & $\bar T=\StabzSdTbar$\\'; ...
 'msi',        'Msi',   'Impulses and switching & \shortstack[l]{$K_1=\StabzMsiKone$\\$K_2=\StabzMsiKtwo$} & $\alpha$ & $\StabzMsiOpen\to\StabzMsiClosed$ & --\\'; ...
 'ncs',        'Ncs',   'Networked, packet dropouts & $\StabzNcsK$ & $\rho$ & $\StabzNcsOpen\to\StabzNcsClosed$ & \shortstack[l]{$\bar T=\StabzNcsTbar$\\$p=\StabzNcsP$}\\'; ...
};
keys  = reg(:,1);
idx   = find(strcmp(key, keys), 1);
if isempty(idx)
    error('mvst:stabtable:key', 'Unknown table key "%s".', key);
end
thisMacro = ['\StabzRow' reg{idx,2}];

% ------- load rows already registered (line-based, one \def per line) ------
present = false(size(keys));
bodies  = reg(:,3);            % default bodies (overwritten by file if newer)
d = fileparts(outfile);
if ~isempty(d) && ~exist(d, 'dir'); mkdir(d); end
if exist(outfile, 'file')
    fid = fopen(outfile, 'r');
    while true
        L = fgetl(fid);
        if ~ischar(L); break; end
        tok = regexp(L, '^\\def\\StabzRow([A-Za-z]+)\{(.*)\}$', 'tokens', 'once');
        if isempty(tok); continue; end
        j = find(strcmp(tok{1}, reg(:,2)), 1);
        if ~isempty(j); present(j) = true; bodies{j} = tok{2}; end
    end
    fclose(fid);
end
present(idx) = true;                 % this call registers/refreshes its row
bodies{idx} = reg{idx,3};

% ------- rewrite the file in canonical order ------------------------------
fid = fopen(outfile, 'w');
if fid < 0; error('mvst:stabtable:open', 'Cannot open %s', outfile); end
fprintf(fid, '%% Auto-generated Stabilization summary-table rows. Do not edit.\n');
fprintf(fid, '%% Generated counterpart of Table~\\ref{tab:stabz:summary}; one row macro per example.\n');
for j = 1:numel(keys)
    if present(j)
        fwrite(fid, ['\def\StabzRow' reg{j,2} '{' bodies{j} '}' newline]);
    end
end
fclose(fid);
fprintf('  stabz_summary.tex: registered row %s\n', thisMacro);
end
