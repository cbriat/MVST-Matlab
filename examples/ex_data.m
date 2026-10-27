function D = ex_data()
%EX_DATA  Example systems used by the optimal-control example scripts.
%   D = EX_DATA() returns a struct with the systems of the chapter:
%     D.ct   single-mode continuous-time LQ with multiplicative noise
%            (used by the finite/infinite/average figures),
%     D.mjls Markov-jump-with-diffusion system (2 modes),
%     D.sd   sampled-data data (A,B,Ac,Bc,Q,R,Tbar),
%     D.sw   switched system (2 modes) with a dwell-time.

% ----- single-mode continuous-time LQ with one multiplicative-noise channel
A  = [ 0 1 ; -1 0.2 ];                 % unstable (Re lambda = +0.1): needs control
B  = [ 0 ; 1 ];
Ac = 0.30*[ 0 0 ; 1 0 ];               % state-multiplicative noise
Bc = [ 0 ; 0.20 ];                     % input-multiplicative noise
Q  = diag([2.0 , 1.0]);
R  = 1.0;
W  = diag([0.30 , 0.30]);              % additive process noise (average cost)
md = mvst_mode(A,B,Q,R,'Ac',Ac,'Bc',Bc,'W',W,'X0',eye(2));
D.ct = mvst_system(md);

% ----- Markov-jump-with-diffusion system: two modes
A1 = [ 0 1 ; -1  0.2 ];   A2 = [ 0 1 ; -2 -0.1 ];
B1 = [ 0 ; 1 ];           B2 = [ 0 ; 1.2 ];
Ac1 = 0.25*[0 0;1 0];     Ac2 = 0.35*[0 0;1 0];
Bc1 = [0;0.15];           Bc2 = [0;0.20];
md1 = mvst_mode(A1,B1,diag([2 1]),1.0,'Ac',Ac1,'Bc',Bc1,'W',diag([.3 .3]),'X0',eye(2));
md2 = mvst_mode(A2,B2,diag([1 1]),1.0,'Ac',Ac2,'Bc',Bc2,'W',diag([.2 .2]),'X0',eye(2));
Pi = [ -1.0  1.0 ;        % generator (rows sum to 0)
        2.0 -2.0 ];
D.mjls = mvst_system({md1,md2}, 'Pi', Pi, 'pi0', [0.5 0.5]);

% ----- sampled-data data (double integrator with small noise + forcing)
D.sd = struct('A',[0 1;0 0], 'B',[0;1], ...
              'Ac',0.10*[0 0;1 0], 'Bc',[0;0.10], ...
              'Q',diag([1 1]), 'R',1.0, 'Tbar',0.30, ...
              'W',diag([0.20 0.20]), 'X0',eye(2));

% ----- switched system: two modes, minimum/range dwell-time
sw1 = mvst_mode([0 1;-1 0.1], [0;1], diag([2 1]), 1.0, 'Ac',0.2*[0 0;1 0],'Bc',[0;0.1],'W',diag([.2 .2]),'X0',eye(2));
sw2 = mvst_mode([0 1;-3 0.2], [0;1], diag([1 2]), 1.0, 'Ac',0.2*[0 0;1 0],'Bc',[0;0.1],'W',diag([.2 .2]),'X0',eye(2));
D.sw = struct('modes',{{sw1,sw2}}, 'Tdwell',0.5);

% ----- finite-time covariance steering instance (minimum-energy, CT) -----
%   steer EX(0)E* = X0 to EX(T)E* = X1f; Q = 0, R = I (minimum control energy)
stA = [0 1;-1 0.1]; stB = [0;1]; stAc = 0.15*[0 0;1 0]; stBc = [0;0.10];
D.steer = struct( ...
   'md',  mvst_mode(stA,stB,zeros(2),1.0,'Ac',stAc,'Bc',stBc,'X0',[3 0;0 0.5]), ...
   'X1f', [1.0 0;0 1.5], 'T', 2.0);
end
