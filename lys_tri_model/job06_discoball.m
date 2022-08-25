%% read in peaks from the 3D-deltaPDF

mapFileName = '../lys_tri_map/export/lys_tri_delta_pdf.h5';
rmax = 4;

[MTp,M] = proc.script.MapTools.fromfile(mapFileName,'delta_pdf');
MTp.Basis = MTp.Basis.orient(); % fix the orientation in real space
MT_peak = MTp.resize('radius',rmax,[0,0,0]);
DT = proc.script.DeltaPDFTools('MT_peak',MT_peak,'rmax',rmax,'supercell',MTp.Grid.P);

delta_pdf_peaks = DT.read_peak_data(mapFileName,'delta_pdf','P');

%% Patterson origin peak

[MT,rho] = proc.script.MapTools.import('export/lys_tri_edens.h5','unitcell','rho');
MT.Basis = MT.Basis.orient(); % fix the orientation in real space

[F,MTf] = MT.fourier_transform(rho);

% compute the patterson map
MTi = MTf;
MTi.type = 'intensity';
MTi.isPeriodic = false;
[P0,MTp] = MTi.fourier_transform(F.*conj(F));

% crop out the central peak
[MTp,resizefun] = MTp.resize('radius',DT.rmax + 1);
P0 = resizefun(P0);

% interpolate on the grid
[x0,y0,z0] = MTp.Grid.grid();
GI = griddedInterpolant(x0,y0,z0,P0,'makima');

[x,y,z] = DT.MT_peak.Grid.grid();
patterson_peak = GI(x,y,z);

patterson_peak(~DT.peak_mask) = NaN;

%% deconvolution

[fit_info,delta_pdf_fit] = DT.deconvolve_peaks(delta_pdf_peaks,patterson_peak);

%% save the results
save proc/discoball_fit.mat DT delta_pdf_peaks patterson_peak delta_pdf_fit fit_info

%% export table of joint adps

T = DT.export_joint_adps(fit_info);
writetable(T,'export/lys_tri_discoball_joint_adps.csv');

