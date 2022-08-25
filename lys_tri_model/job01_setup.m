%%

mkdir proc
mkdir export

addpath('~/Documents/GitHub/ando-lab/mdx-lib/')

%% get atomic model and structure factor from REFMAC output
%
% 1) load the REFMAC reflections file
% 2) load the AIMLESS intensities file (Imerge)
% 3) figure out what ADP scaling was applied by REFMAC's Fobs
% 4) undo the ADP scaling to F columns
% 5) output structure factors in a standard table
% 6) read in the PDB file
% 7) undo the ADP scaling of anisou parameters
% 8) apply extinction correction
% 9) save everything: Atoms, AtomFF, Basis, SpaceGroup, hklTable


clear opts
opts.mtzREFMAC =   'import/lys_nitr_ubatch4_refmac.mtz';
opts.mtzMerged =   'import/june2017_nitrate_ubatch4_aimless_truncate1.mtz';
opts.pdbFileName = 'import/lys_nitr_ubatch4_refmac.pdb';

opts.refmac_filter = struct('H','h','K','k','L','l','FREER','isFree','FP','Fobs','SIGFP','sigmaFobs',...
    'FC_ALL_LS','Fcalc','PHIC_ALL_LS','phiFcalc');

opts.merged_filter = struct('H','h','K','k','L','l','IMEAN','Imerge','SIGIMEAN','sigmaImerge');

%% load refmac reflections file
[hklREFMAC,Basis,SpaceGroup] = proc.script.ImportMTZ('mtz2mdx',opts.refmac_filter).run(opts.mtzREFMAC);

% reformat the table (--> T1)
T1 = hklREFMAC(:,{'h','k','l','isFree','Fobs','sigmaFobs','Fcalc'});
T1.Fcalc = T1.Fcalc.*exp(1i*(pi/180).*hklREFMAC.phiFcalc); % add phases

% load the aimless reflections file
hklMerged = proc.script.ImportMTZ('mtz2mdx',opts.merged_filter).run(opts.mtzMerged); % Imerge

% reformat the table (--> T2)
T2 = hklMerged(:,{'h','k','l','Imerge','sigmaImerge'});
T2 = T2(~isnan(T2.Imerge),:);

% combine tables (T <-- T1, T2)
T = outerjoin(T1,T2,'Keys',{'h','k','l'},'MergeKeys',true);

%%
% figure out what ADP scaling was applied by REFMAC to Fobs

[params,anisoscales] = anisoscaling(T,Basis,SpaceGroup);
Uadd = params.U;
disp(Uadd)

%     0.0065   -0.0021    0.0006
%    -0.0021    0.0166   -0.0017
%     0.0006   -0.0017   -0.0160

%%
% undo the ADP scaling to F columns

hklTable = T1;
hklTable.Fobs = hklTable.Fobs.*anisoscales;
hklTable.sigmaFobs = hklTable.sigmaFobs.*anisoscales;
hklTable.Fcalc = hklTable.Fcalc.*anisoscales;

%%
% apply extinction correction

[sx,sy,sz] = Basis.invert.frac2lab(hklTable.h,hklTable.k,hklTable.l);
s = sqrt(sx.^2 + sy.^2 + sz.^2);

lambda = 0.97;
sin2theta = lambda*s/2;
x = 0.002; % determined previously by fitting Fobs vs. Fcalc

p = (1 + 0.001*x*abs(hklTable.Fcalc).^2*lambda^3./sin2theta).^(-1/4);

hklTable.Fobs = hklTable.Fobs./p;
%%
% read in the PDB file
Atoms = proc.script.ImportPDB().run(opts.pdbFileName);

AtomFF = proc.script.CoordinateTools.to_xray_structure(Atoms);

% apply anisotropic ADP to atomic form factors
AtomFF.mdxFormFactor = AtomFF.mdxFormFactor.addU(Uadd);

%%
% save everything: Atoms, AtomFF, Basis, SpaceGroup, hklTable
save proc/atomic_model.mat Atoms AtomFF Basis SpaceGroup hklTable opts Uadd

%% helper functions

function [params,anisoscales] = anisoscaling(T,Basis,SpaceGroup)

SM = proc.script.ScaleModelToFobs('Basis',Basis.orient(),'SpaceGroup',SpaceGroup);

[T.sx,T.sy,T.sz] = SM.Basis.invert.frac2lab(T.h,T.k,T.l);

tt = T;
tt = tt(~isnan(tt.Imerge) & ~isnan(tt.Fobs) & tt.sigmaImerge>0,:); % get rid of NaNs

[numParams,param2struct] = SM.parameterize_model();
p2s = @(v) param2struct([v(1),0,0,v(2:(numParams-2))]); % HACK to remove ksol, Bsol from fitting model

Taniso = @(s) latt.Blob(1,0).addU(s.U).rescale(s.ktot);
residfun = @(s,t) (t.Imerge - abs(Taniso(s).scatteringAmplitude(t.sx,t.sy,t.sz).*t.Fobs).^2)./t.sigmaImerge;

solvFit = lsqnonlin(@(v) residfun(p2s(v/1000),tt),1000*ones(1,numParams-2),[],[]);

params = p2s(solvFit/1000);

anisoscales = latt.Blob(1,0).addU(params.U).scatteringAmplitude(T.sx,T.sy,T.sz);

end