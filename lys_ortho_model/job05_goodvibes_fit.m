%%
% load reference data and lattice model

load proc/goodvibes_model.mat ENM ResidueGroups
load proc/reference_halos.mat I sigma
load proc/reference_calc.mat Gk ind LDT

% initialize variables to store refinement info (4 stages)
refineLog = struct('seconds',{},'history',{},'fitinfo',{},'kfit',{},'springType',{},'parameterization',{},'pfit',{});

useParallel = true;

%% Stage 1
% fit overall Gaussian spring constant

parameterization = 'global';
springType = 'Gaussian';

[param2k,p0] = ENM.parameterize(parameterization,0.5);

Vfun = @(p) ENM.Hessian(springType,param2k(p));

tic
[pfit,fitinfo,history] = LDT.fitHessianToHalos(I,sigma,Gk,ind,Vfun,p0,0,Inf,...
    'MaxFunctionEvaluations',1000,'UseParallel',useParallel,'Display','iter');

refineLog(1) = struct(...
    'seconds',toc,...
    'history',history,...
    'fitinfo',fitinfo,...
    'kfit',param2k(pfit),...
    'springType',springType,...
    'parameterization',parameterization,...
    'pfit',pfit);

%% Stage 2
% fit Gaussian spring constant for each interface

parameterization = 'interface';
springType = 'Gaussian';

[param2k,p0] = ENM.parameterize(parameterization,refineLog(1).kfit);

Vfun = @(p) ENM.Hessian(springType,param2k(p));

tic
[pfit,fitinfo,history] = LDT.fitHessianToHalos(I,sigma,Gk,ind,Vfun,p0,0*p0,Inf*p0,...
    'MaxFunctionEvaluations',1000,'UseParallel',useParallel,'Display','iter');

refineLog(2) = struct(...
    'seconds',toc,...
    'history',history,...
    'fitinfo',fitinfo,...
    'kfit',param2k(pfit),...
    'springType',springType,...
    'parameterization',parameterization,...
    'pfit',pfit);

%% Stage 3
% fit hybrid spring constants for each interface

parameterization = 'interface';
springType = 'hybrid';

k0 = repmat(refineLog(2).kfit,[1,2]);
[param2k,p0] = ENM.parameterize(parameterization,k0);

Vfun = @(p) ENM.Hessian(springType,param2k(p));
tic
[pfit,fitinfo,history] = LDT.fitHessianToHalos(I,sigma,Gk,ind,Vfun,p0,0*p0,Inf*p0,...
    'MaxFunctionEvaluations',1000,'UseParallel',useParallel,'Display','iter');

refineLog(3) = struct(...
    'seconds',toc,...
    'history',history,...
    'fitinfo',fitinfo,...
    'kfit',param2k(pfit),...
    'springType',springType,...
    'parameterization',parameterization,...
    'pfit',pfit);
%% Stage 4
% fit individual hybrid springs

parameterization = 'uniquegrouped';
groupAssignment = ResidueGroups.mdxGroupByResidue;
springType = 'hybrid';

[param2k,p0] = ENM.parameterize(parameterization,refineLog(3).kfit,groupAssignment);

Vfun = @(p) ENM.Hessian(springType,param2k(p));

tic
[pfit,fitinfo,history] = LDT.fitHessianToHalos(I,sigma,Gk,ind,Vfun,p0,0*p0,Inf*p0,...
    'MaxFunctionEvaluations',10000,'UseParallel',useParallel,'Display','iter');

refineLog(4) = struct(...
    'seconds',toc,...
    'history',history,...
    'fitinfo',fitinfo,...
    'kfit',param2k(pfit),...
    'springType',springType,...
    'parameterization',parameterization,...
    'pfit',pfit);

%% save result

LD = nm.LatticeDynamics('V',Vfun(pfit),'supercell',LDT.supercell,'M',LDT.M);

save proc/goodvibes_fit.mat LD refineLog

%% export joint adps

Tsymm = proc.script.LatticeDynamicsTools.calc_symmetrized_covariances(ENM.Cell,LD);
Tasu = proc.script.LatticeDynamicsTools.calc_asu_covariances(ENM.Cell,LD);

writetable(Tsymm,'export/lys_ortho_goodvibes_joint_adps.csv');
writetable(Tasu,'export/lys_ortho_goodvibes_com_covariances.csv');