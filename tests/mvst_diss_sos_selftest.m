function res = mvst_diss_sos_selftest()
%MVST_DISS_SOS_SELFTEST  Validation of MVST_DISS_SOS against exact ground truth.
%   Runs the cases documented in the MVST_DISS_SOS header and returns a struct
%   array with fields .name, .expected, .got, .ok.  Prints a table.
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
addpath(root); addpath(fullfile(root,'lib'));
if exist('sossolve','file') ~= 2 || exist('sedumi','file') == 0
    fprintf('%s: SKIPPED (SOSTOOLS + SeDuMi required; run mvst_setup).\n', mfilename);
    res = [];
    return
end

R = struct('name',{},'expected',{},'got',{},'ok',{});

% ---- (1a) CT L2-gain, stable diagonal drift -----------------------------
A=[-3 0;0 -3]; Bw=[0;1]; C=[0.5 1]; Dw=0.5;
exp1 = mvst_l2gain_ct(A,Bw,C,Dw);
r = mvst_diss_sos(struct('A0',A,'Bw',Bw,'C',C,'Dw',Dw), ...
                  struct('type','l2gain'), struct('type','arbitrary'), 'gamma');
R(end+1) = row('1a CT L2-gain (diag)', exp1, r.gamma, 1.5e-2);

% ---- (1b) CT L2-gain, oscillatory stable drift --------------------------
A=[0 1;-2 -1]; Bw=[0;1]; C=[1 0]; Dw=0.5;
exp2 = mvst_l2gain_ct(A,Bw,C,Dw);
r = mvst_diss_sos(struct('A0',A,'Bw',Bw,'C',C,'Dw',Dw), ...
                  struct('type','l2gain'), struct('type','arbitrary'), 'gamma');
R(end+1) = row('1b CT L2-gain (osc)', exp2, r.gamma, 1.5e-2);

% ---- (2a) passivity of a positive-real system (feasible) ----------------
A=[0 1;-1 -1]; Bw=[0;1]; C=[1 0]; Dw=0.5;
r = mvst_diss_sos(struct('A0',A,'Bw',Bw,'C',C,'Dw',Dw), ...
                  struct('type','passivity'), struct('type','arbitrary'), 'feas');
R(end+1) = rowb('2a passivity (PR system)', true, r.feas);

% ---- (2b) non-passive system (infeasible) -------------------------------
% C = [0 1] with Dw = 0: relative degree, strictly non-passive at high freq.
A=[0 1;-1 -1]; Bw=[0;1]; C=[0 1]; Dw=0;
r = mvst_diss_sos(struct('A0',A,'Bw',Bw,'C',C,'Dw',Dw), ...
                  struct('type','passivity'), struct('type','arbitrary'), 'feas');
R(end+1) = rowb('2b non-passive system', false, r.feas);

% ---- (3) periodic sampled-output L2-gain --------------------------------
A=[0 1;-1 -1]; Bw=[0;1]; C=[1 0]; Tbar=0.3;
exp3 = mvst_l2gain_sd_hybrid(A,Bw,C,Tbar);
% flow: pure disturbance credit (no flow output); output sampled at the jump
sys3 = struct('A0',A,'Bw',Bw,'J',eye(2),'CJ',C);
r = mvst_diss_sos(sys3, struct('type','l2gain'), ...
                  struct('type','periodic','Tbar',Tbar,'deg',4), 'gamma');
R(end+1) = row('3 periodic sampled gain', exp3, r.gamma, 3e-2);

% ---- (4) Poisson mean-square L2-gain (ex_diss_poisson data) -------------
A=[0 1;-2 -1]; Bw=[0;1]; C=[1 0]; Dw=0.5;
C1=0.15*eye(2); Dw1=[0.1;0]; J1=1.15*eye(2); lam=0.5;
exp4 = mvst_l2gain_ct_poisson(A,Bw,C,Dw,{C1},{Dw1},{J1},lam);
sys4 = struct('A0',A,'Bw',Bw,'C',C,'Dw',Dw, ...
              'Ac',C1,'Dwc',Dw1,'Jf',J1,'lamf',lam);
r = mvst_diss_sos(sys4, struct('type','l2gain'), struct('type','arbitrary'), 'gamma');
R(end+1) = row('4 Poisson ms gain', exp4, r.gamma, 3e-2);

% ---- print --------------------------------------------------------------
fprintf('\n==== mvst_diss_sos self-test ====\n');
fprintf('%-28s %12s %12s %6s\n','case','expected','got','ok');
for i=1:numel(R)
    if isnan(R(i).expected)
        fprintf('%-28s %12s %12.6f %6d\n', R(i).name, '(bool)', R(i).got, R(i).ok);
    else
        fprintf('%-28s %12.6f %12.6f %6d\n', R(i).name, R(i).expected, R(i).got, R(i).ok);
    end
end
fprintf('ALL PASS: %d\n', all([R.ok]));
res = R;
end

function s = row(name, expected, got, tol)
s.name=name; s.expected=expected; s.got=got;
s.ok = abs(got-expected) <= tol*max(1,abs(expected));
end
function s = rowb(name, expectedBool, gotBool)
s.name=name; s.expected=NaN; s.got=double(gotBool);
s.ok = (logical(gotBool)==logical(expectedBool));
end
