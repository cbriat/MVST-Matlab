function o = mvst_sosopts(x)
%MVST_SOSOPTS  Options struct handed to MVST_SOSSOLVE, built from normalized options.
%   O = MVST_SOSOPTS(X) accepts a solver name or the struct returned by MVST_OPTS
%   and returns the struct for MVST_SOSSOLVE: .solver, the forwarded .reltol,
%   .verbosity and .nosdp when present, and .params.fid = 0 (silent backend).
%
%   See also MVST_OPTS, MVST_SOSSOLVE.
if ischar(x) || isstring(x)
    o = struct('solver', char(x));
else
    o = struct('solver', x.solver);
    for f = {'reltol','verbosity','nosdp'}
        if isfield(x, f{1}); o.(f{1}) = x.(f{1}); end
    end
end
o.params.fid = 0;
end
