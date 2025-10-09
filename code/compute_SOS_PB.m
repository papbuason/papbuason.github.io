function Avoltage = compute_SOS_PB(Jinv,Ybus,cfg,mpc)
%Paprapee Buason
%Georgia Tech and Los Alamos National Laboratory
%Main File for Second Order Sensitivity Matrix for voltage
%Spring 2024
%Jinv = inverse of Jacobian, cfg = PQ bus number

define_constants;

pq = find(mpc.bus(:,BUS_TYPE) == 1); 
pv = find(mpc.bus(:,BUS_TYPE) == 2); 

idx_Jinv = find(abs(Jinv(cfg + length(pv) + length(pq),:)) >= 1e-4); %find non-zero indices
a = -Jinv(cfg + length(pv) + length(pq),:);
a = a(idx_Jinv);

H = cell(length(pv) + 2*length(pq),1);
for ii = 1:length(pv) + 2*length(pq)
    if ii <= length(pv) + length(pq)
        H{ii} = makeHes_PB('t', ii, Ybus, mpc);
    else 
        H{ii} = makeHes_PB('vmag', ii - length(pq) - length(pv), Ybus, mpc);
    end
end

sz_H = 0; reps = zeros(length(pv) + 2*length(pq),1);
for ii = 1:length(pv) + 2*length(pq)
   H{ii} = H{ii}(idx_Jinv, :);
   sz_H = sz_H + nnz(H{ii});
   reps(ii) = nnz(H{ii});
end

big_H = zeros(sz_H,4); 
for ii = 1:length(pv) + 2*length(pq)
    [kk,ll,v] = find(H{ii});
    big_H(sum(reps(1:ii-1)) + 1:sum(reps(1:ii)),:) = [repmat(ii,nnz(H{ii}),1), kk, ll, v];
end


Avoltage = zeros(length(pv) + 2*length(pq));
for n = 1:length(pv)+2*length(pq) %n is index of power injections
    temp_Ahes = big_H(:,4).*repelem(Jinv(:,n),reps);
    Ahes = sparse(big_H(:,2),big_H(:,3),temp_Ahes,length(a),2*length(pq) + length(pv));    
    Avoltage(n,:) = a*Ahes;
end 

Avoltage = Avoltage*full(Jinv);

%% remove some sensitivity if there is no load.
P_demand_0 = find(mpc.bus(:,PD) == 0);
Q_demand_0 = find(mpc.bus(:,QD) == 0);
pvpq = find(mpc.bus(:,BUS_TYPE) ~= 3);
ind_remove_p = find(ismember(pvpq,P_demand_0) == 1);
ind_remove_q = length(pv) + length(pq) + find(ismember(pq,Q_demand_0) == 1);    
ind_remove = union(ind_remove_p,ind_remove_q);

Avoltage(ind_remove,:) = [];
Avoltage(:,ind_remove) = [];

function J = makeHes_PB(inp, m, Ybus, baseMVA, bus, branch, gen, fullJac)
%MAKEJAC  Forms the power flow Jacobian.
%   J = MAKEJAC(MPC)
%   J = MAKEJAC(MPC, FULLJAC)
%   J = MAKEJAC(BASEMVA, BUS, BRANCH, GEN)
%   J = MAKEJAC(BASEMVA, BUS, BRANCH, GEN, FULLJAC)
%   [J, YBUS, YF, YT] = MAKEJAC(MPC)


if nargin < 8
    mpc     = baseMVA;
    if nargin > 5
        fullJac = bus;
    else
        fullJac = 0;
    end
    baseMVA = mpc.baseMVA;
    bus     = mpc.bus;
    branch  = mpc.branch;
    gen     = mpc.gen;
elseif nargin < 9
    fullJac = 0;
end

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% build Ybus
% [Ybus, ~, ~] = makeYbus(baseMVA, bus, branch);

%% extract voltage
V = bus(:, VM) .* exp(1j * pi/180 * bus(:, VA));

%% make sure we use generator setpoint voltage for PV and slack buses
on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
gbus = gen(on, GEN_BUS);                %% what buses are they at?
k = find(bus(gbus, BUS_TYPE) == PV | bus(gbus, BUS_TYPE) == REF);
V(gbus(k)) = gen(on(k), VG) ./ abs(V(gbus(k))).* V(gbus(k));

%% build Jacobian
[dSbus_dVa, dSbus_dVm] = d2Sbus_dV_DV_PB(Ybus, V, m, inp, bus);

%% get bus index lists of each type of bus
[~, pv, pq] = bustypes(bus, gen);

j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
j12 = real(dSbus_dVm([pv; pq], pq));
j21 = imag(dSbus_dVa(pq, [pv; pq])); 
j22 = imag(dSbus_dVm(pq, pq));


J = [   j11 j12;
        j21 j22;    ];

function [dSbus_dVa, dSbus_dVm] = d2Sbus_dV_DV_PB(Ybus, V, m ,inp, bus)
%% default input args
n = length(V);
Ibus = Ybus * V;
pq = find(bus(:,2) == 1);
pvpq = [find(bus(:,2) == 2); find(bus(:,2) == 1)];
if inp == 'vmag'
    m = pq(m);
elseif inp == 't'
    m = pvpq(m);
end

diagV       = sparse(1:n, 1:n, V, n, n);
diagIbus    = sparse(1:n, 1:n, Ibus, n, n);
diagVnorm   = sparse(1:n, 1:n, V./abs(V), n, n);
Vm = sparse(m,1,V(m),n,1);
diagVm = diag(Vm); 

tm = sparse(m,1,exp(1j*angle(V(m))),n,1); 
diagtm = diag(tm);

if inp == 't'    
    dSbus_dVa = diagVm*conj(-diagIbus + Ybus*diagV) + diagV*conj(diag(Ybus*Vm) - Ybus*diagVm);                     %% d2Sbus/dVa_dVa
    dSbus_dVm = 1j*(-diagVnorm*conj(diag(Ybus*Vm)) + diagtm*conj(diagIbus) - diagV*conj(Ybus*diagtm) + diagVm*conj(Ybus*diagVnorm)); %% d2Sbus/dVm_dVa
elseif inp == 'vmag'
    dSbus_dVa =  1j * (diagtm * conj(diagIbus - Ybus * diagV) + diagV*(conj(diag(Ybus*tm)) - conj(Ybus*diagtm)));    %% d2Sbus/dVa_dVm
    dSbus_dVm = diagVnorm*conj(diag(Ybus*tm)) + diagtm*conj(Ybus*diagVnorm);    %% d2Sbus/dVm_dVm
end

