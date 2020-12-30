function [res, success, et, niter] = runpf_simp(mpc, mpopt)
%% Simplified Newton-Raphson Power Flow, as an alternative for
%% standard matpower runpf(), to save time and improve readability.
%% input: matpower case mpc, matpower option mpopt (compatible)
%% output: solved case res, convergence flag success, elasped time
%% et, and number of iterations

%% Dealing with argins
if nargin < 2 || isempty(mpopt)
    mpopt.pf.nr.max_it = 30;
    mpopt.pf.nr.tol = 1e-8;
    mpopt.pf.nr.lin_solver = '';
    mpopt.verbose = 1;
end

max_itFlag = ...
    ~isfield(mpopt, 'pf') || isempty(mpopt.pf) || ~isfield(mpopt.pf, 'nr') ||...
    isempty(mpopt.pf.nr) || ~isfield(mpopt.pf.nr, 'max_it') || ...
    isempty(mpopt.pf.nr.max_it);
tolFlag = ...
    ~isfield(mpopt, 'pf') ||isempty(mpopt.pf) || ~isfield(mpopt.pf, 'nr') ||...
    isempty(mpopt.pf.nr) || ~isfield(mpopt.pf.nr, 'tol') || ...
    isempty(mpopt.pf.nr.tol);
lin_solverFlag = ...
    ~isfield(mpopt, 'pf') || isempty(mpopt.pf) || ~isfield(mpopt.pf, 'nr') ||...
    isempty(mpopt.pf.nr) || ~isfield(mpopt.pf.nr, 'lin_solver') || ...
    isempty(mpopt.pf.nr.lin_solver);
verboseFlag = isempty(mpopt.verbose) || ~isfield(mpopt, 'verbose');

if max_itFlag
    mpopt.pf.nr.max_it = 30;
end

if tolFlag
    mpopt.pf.nr.tol = 1e-8;
end

if lin_solverFlag
    mpopt.pf.nr.lin_solver = '';
end

if verboseFlag
    mpopt.verbose = 1;
end


%% Specify parameters
tol = mpopt.pf.nr.tol; % Convergence tolerance
max_iter = mpopt.pf.nr.max_it;
lin_solver  = mpopt.pf.nr.lin_solver;

%% Load the test case
[isampc, errInfo] = ismpc(mpc);
if ~isampc
    error(['this may not be a standard matpower case, ', errInfo]);
end

%% Store placeholders for various columns in the Matpower data format.
% See Appendix B of the Matpower manual.
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, VM, VA] = idx_bus;
[GEN_BUS, ~, ~, ~, ~, VG, ~,GEN_STATUS] = idx_gen;
[~, ~, ~, ~, ~, ~, ~] = idx_cost;

[Sbase, bus, gen, branch] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);

nb = size(bus,1);

%% Compute the admittance matrix
[Ybus, Yf, Yt] = makeYbus(mpc);

[ref, pv, pq] = bustypes(bus, gen);
on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
gbus = gen(on, GEN_BUS);                %% what buses are they at?

% Store voltage initialization
V0 = bus(:,VM) .* exp(1j * bus(:,VA) * pi / 180); % Complex voltage phasor
Vm0 = bus(:,VM);
Va0 = bus(:,VA);
iter = 0;

vcb = ones(size(V0));           %% create mask of voltage-controlled buses
vcb(pq) = 0;                    %% exclude PQ buses
k = find(vcb(gbus));            %% in-service gens at v-c buses
V0(gbus(k)) = gen(on(k), VG) ./ abs(V0(gbus(k))).* V0(gbus(k));

npv = length(pv);
npq = length(pq);
j1 = 1;         j2 = npv;           %% j1:j2 - V angle of pv buses
j3 = j2 + 1;    j4 = j2 + npq;      %% j3:j4 - V angle of pq buses
j5 = j4 + 1;    j6 = j4 + npq;      %% j5:j6 - V mag of pq buses

mis = V0 .* conj(Ybus * V0) - makeSbus(Sbase, bus, gen, [], Vm0);
F = [   real(mis([pv; pq]));
    imag(mis(pq))   ];

normF = norm(F, inf);
converged = 0;
if mpopt.verbose > 1
    fprintf('\n it    max P & Q mismatch (p.u.)');
    fprintf('\n----  ---------------------------');
    fprintf('\n%3d        %10.3e', iter, normF);
end
if normF < tol
    converged = 1;
    if mpopt.verbose > 1
        fprintf('\nConverged!\n');
    end
end

if isempty(lin_solver)
    nx = length(F);
    if nx <= 10
        lin_solver = '\';       %% default \ operator
    else    %% MATLAB and nx > 10 or Octave and nx > 2000
        lin_solver = 'LU3';     %% LU decomp with 3 output args, AMD ordering
    end
end

%% do Newton iterations
tic;
V = V0;
Vm = abs(V);
Va = angle(V);
while (~converged && iter < max_iter)
    %% update iteration counter
    iter = iter + 1;
    
    %% evaluate Jacobian
    [dSbus_dVa, dSbus_dVm] = dSbus_dV(Ybus, V);
    [~, neg_dSd_dVm] =  makeSbus(Sbase, bus, gen, [], Vm);
    dSbus_dVm = dSbus_dVm - neg_dSd_dVm;
    
    j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
    j12 = real(dSbus_dVm([pv; pq], pq));
    j21 = imag(dSbus_dVa(pq, [pv; pq]));
    j22 = imag(dSbus_dVm(pq, pq));
    
    J = [   j11 j12;
        j21 j22;    ];
    
    %% compute update step
    dx = mplinsolve(J, -F, lin_solver);
    
    %% update voltage
    if npv
        Va(pv) = Va(pv) + dx(j1:j2);
    end
    if npq
        Va(pq) = Va(pq) + dx(j3:j4);
        Vm(pq) = Vm(pq) + dx(j5:j6);
    end
    V = Vm .* exp(1j * Va);
    Vm = abs(V);            %% update Vm and Va again in case
    Va = angle(V);          %% we wrapped around with a negative Vm
    
    %% evalute F(x)
    mis = V .* conj(Ybus * V) - makeSbus(Sbase, bus, gen, [], Vm);
    F = [   real(mis([pv; pq]));
        imag(mis(pq))   ];
    
    %% check for convergence
    normF = norm(F, inf);
    if mpopt.verbose > 1
        fprintf('\n%3d        %10.3e', iter, normF);
    end
    converged = normF < tol;
end
et = toc;
if mpopt.verbose
    fprintf('\nNewton''s method power flow (power balance, polar) converged in %d iterations.\n', iter);
    fprintf('\nElasped time is %d seconds.\n', et);
end

if mpopt.verbose
    if ~converged
        fprintf('\nNewton''s method power flow (power balance, polar) did not converge in %d iterations.\n', iter);
        fprintf('\n>>>>> Did NOT converge (%d seconds).\n', etc);
    end
end

[bus, gen, branch] = pfsoln(Sbase, bus, gen, branch, Ybus, Yf, Yt, V, ref, pv, pq);
[res.bus, res.gen, res.branch] = deal(bus, gen, branch);

%% Update the success flag
success = logical(converged);
niter = iter;
end
