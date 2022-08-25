load proc/atomic_model.mat AtomFF Basis SpaceGroup Atoms hklTable

UnitCellGrid = proc.script.GridDesigner(...
    'Basis',Basis,...
    'SpaceGroup',SpaceGroup,...
    'dgrid',1/3).run();

MT = proc.script.MapTools('Grid',UnitCellGrid.PeriodicGrid,...
    'Basis',UnitCellGrid.Basis.orient,...
    'SpaceGroup',SpaceGroup,...
    'type','density',...
    'isPeriodic',true);

% calculate electron density from Fobs using calculated phases
hklTable.F = hklTable.Fobs.*exp(1i*angle(hklTable.Fcalc));
hklTable.F(isnan(hklTable.F)) = hklTable.Fcalc(isnan(hklTable.F));

% add Friedel symmetry mates
tmp = hklTable;
tmp.F = conj(tmp.F);
tmp.h = -tmp.h;
tmp.k = -tmp.k;
tmp.l = -tmp.l;
hklTable = [hklTable;tmp];

F = MT.invert.table2array(hklTable(:,{'h','k','l','F'}),'symexpand');
rho0 = MT.inverse_fourier_transform(F);

A = Atoms;
A.mdxAltGroup = true(size(Atoms,1),1);

BSM = proc.script.BulkSolventMask(...
    'LatticeGrid',latt.LatticeGrid(MT.Grid,MT.Basis),...
    'Atoms',A,...
    'SpaceGroup',SpaceGroup,...
    'rMax',5);

[indMapASU,probDist] = BSM.distance_map();
indMaps = BSM.run();

%% calculate asu mask, dilate, and weight by multiplicity

is_inside = indMaps{1}>0;
is_inside = imdilate(is_inside,strel('sphere',1));
index_extended = indMapASU.*is_inside;

% compute the multiplicity of the ASU mask
Tasu = MT.array2table(double(is_inside),is_inside);
mult = MT.table2array(Tasu,'symexpand','sum');
ind = MT.frac2ind(Tasu.x,Tasu.y,Tasu.z);
Tasu.mult = 1./mult(ind);
msk_asu = MT.table2array(Tasu(:,{'x','y','z','mult'}));

% compute the solvent mask (unit cell)
tmp = mult;
tmp(isnan(tmp)) = 0;
is_solv = not(logical(tmp));
rho_avg_solv = mean(rho0(is_solv));

msk_asu_blur = blur_density(MT,msk_asu,50);

% compute the weighted electron density of the ASU
rho_asu = msk_asu_blur.*(rho0-rho_avg_solv);

%% convert to table, add atom coordinates, and unwrap

T = MT.array2table(rho_asu,logical(indMapASU));
tmp = MT.array2table(indMapASU,logical(indMapASU));
T.atom_index = tmp.(4);

[xa,ya,za] = MT.Basis.lab2frac(A.x,A.y,A.z);
T.x = T.x - round(T.x-xa(T.atom_index));
T.y = T.y - round(T.y-ya(T.atom_index));
T.z = T.z - round(T.z-za(T.atom_index));

% make supercell grid (temporarily)
MTc = MT.resize('roi',[min(T.x),min(T.y),min(T.z)],[max(T.x),max(T.y),max(T.z)]);

rho_asu_unwrap = MTc.table2array(T);

%% export

% resample on a coarser grid by Fourier cropping
[rho_asu_unwrap2,MTc2] = MTc.fourier_interpolate(rho_asu_unwrap,2/3);

%save the masked density
h5Out = 'export/lys_tet_edens.h5';

% save
rho = rho0;
MT.export(h5Out,'unitcell',rho,'datatype','single');

% save the unwrapped density
rho = rho_asu_unwrap;
MTc.export(h5Out,'unwrapped',rho,'newmap',true);

% save the unwrapped cropped density
rho = rho_asu_unwrap2;
MTc2.export(h5Out,'unwrapped_cropped',rho,'newmap',true);
%% FUNCTIONS

function rhoBlur = blur_density(MT,rho,Badd)
[F,MTf] = MT.fourier_transform(rho);
[sx,sy,sz] = MTf.Grid.grid();
[sx,sy,sz] = MTf.Basis.frac2lab(sx,sy,sz);
Fblur = latt.Blob(1,0).addB(Badd).scatteringAmplitude(sx,sy,sz);
rhoBlur = MT.inverse_fourier_transform(F.*Fblur);
end
