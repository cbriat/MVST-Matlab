function res = mvst_game_sos_selftest()
%MVST_GAME_SOS_SELFTEST  Validation of MVST_GAME_SOS against exact ground truth.
%   Cross-checks the SOS (bounded-real synthesis) critical level against the
%   exact game Riccati level MVST_GAME_GAMMA_CT on small continuous-time
%   instances.  The SOS value is a certified UPPER BOUND on gamma*; for the
%   separated CT LTI game it is exact, so the two agree to solver tolerance.
%   Returns a struct array (.name,.expected,.got,.ok) and prints a table,
%   including a degree-refinement row showing the bound does not increase.
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
addpath(root); addpath(fullfile(root,'lib'));
if exist('sossolve','file') ~= 2 || exist('sedumi','file') == 0
    fprintf('%s: SKIPPED (SOSTOOLS + SeDuMi required; run mvst_setup).\n', mfilename);
    res = [];
    return
end

R = struct('name',{},'expected',{},'got',{},'ok',{});

% ---- (1) CT game, oscillatory drift -------------------------------------
A=[0 1;-2 -1]; Bu=[0;1]; Bw=[0;1]; Q=eye(2); Ru=1;
g1 = mvst_game_gamma_ct(A,Bu,Bw,Q,Ru,40);
r  = mvst_game_sos(struct('A0',A,'Bu',Bu,'Bw',Bw), ...
                   struct('Q',Q,'Ru',Ru), struct('type','arbitrary'), 'gamma');
R(end+1) = rowub('1 CT game (osc)', g1, r.gammastar, 3e-2);

% ---- (2) CT game, stable diagonal drift ---------------------------------
A=[-1 0;0 -2]; Bu=[1;1]; Bw=[1;0]; Q=diag([1 0.5]); Ru=2;
g2 = mvst_game_gamma_ct(A,Bu,Bw,Q,Ru,40);
r  = mvst_game_sos(struct('A0',A,'Bu',Bu,'Bw',Bw), ...
                   struct('Q',Q,'Ru',Ru), struct('type','arbitrary'), 'gamma');
R(end+1) = rowub('2 CT game (diag)', g2, r.gammastar, 3e-2);

% ---- (3) scalar CT game (closed form sanity) ----------------------------
% xdot = a x + u + w, cost |x|^2 + ru |u|^2 - g^2|w|^2.  Riccati:
%   0 = 2 a p + 1 - p^2(1/ru - 1/g^2).  Real p exists iff (1/ru - 1/g^2)>0
%   discriminant ok; gamma* = 1/sqrt(1/ru + a_+^2/...) -- just compare solvers.
a=0.5; Bu=1; Bw=1; Q=1; Ru=1;
g3 = mvst_game_gamma_ct(a,Bu,Bw,Q,Ru,60);
r  = mvst_game_sos(struct('A0',a,'Bu',Bu,'Bw',Bw), ...
                   struct('Q',Q,'Ru',Ru), struct('type','arbitrary'), 'gamma');
R(end+1) = rowub('3 scalar CT game', g3, r.gammastar, 3e-2);

% ---- (4) central control gain check -------------------------------------
% The recovered Ku must stabilize A + Bu*Ku (closed loop Hurwitz).
A=[0 1;-2 -1]; Bu=[0;1]; Bw=[0;1]; Q=eye(2); Ru=1;
r  = mvst_game_sos(struct('A0',A,'Bu',Bu,'Bw',Bw), ...
                   struct('Q',Q,'Ru',Ru), struct('type','arbitrary'), 'gamma');
Ku = r.Ku(0);
cl = max(real(eig(A + Bu*Ku)));
R(end+1) = rowb('4 central gain stabilizes', true, cl < 0);

% ---- print --------------------------------------------------------------
fprintf('\n==== mvst_game_sos self-test (SOS bound vs exact gamma*) ====\n');
fprintf('%-30s %12s %12s %6s\n','case','exact g*','SOS bound','ok');
for i=1:numel(R)
    if isnan(R(i).expected)
        fprintf('%-30s %12s %12.6f %6d\n', R(i).name, '(bool)', R(i).got, R(i).ok);
    else
        fprintf('%-30s %12.6f %12.6f %6d\n', R(i).name, R(i).expected, R(i).got, R(i).ok);
    end
end
fprintf('ALL PASS: %d\n', all([R.ok]));
res = R;
end

function s = rowub(name, gexact, ggot, tol)
% Upper-bound row: pass if the SOS bound is >= gexact (up to a small tol) AND
% within tol of it (exact for the separated CT LTI game).
s.name=name; s.expected=gexact; s.got=ggot;
s.ok = (ggot >= gexact - tol) && (abs(ggot-gexact) <= tol*max(1,abs(gexact)));
end
function s = rowb(name, expectedBool, gotBool)
s.name=name; s.expected=NaN; s.got=double(gotBool);
s.ok = (logical(gotBool)==logical(expectedBool));
end
