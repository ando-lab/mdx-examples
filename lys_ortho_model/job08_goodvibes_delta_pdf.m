%% delta pdf map from goodvibes

mapFileIn = '../lys_ortho_map/export/lys_ortho_delta_pdf.h5';
MTp = proc.script.MapTools.import(mapFileIn,'delta_pdf');

% load goodvibes simulation
[MT,I] = proc.script.MapTools.import('export/lys_ortho_goodvibes.h5','supercell','I');
tmp = MT.array2table(I);
MT = MTp.invert;
I = MT.table2array(tmp,'symexpand','replace',NaN);

% punch out Bragg peaks and fill using nearest neighbors
[h,k,l] = MT.Grid.grid();
isBragg = abs(mod(h,1))<1E-6 & abs(mod(k,1))<1E-6 & abs(mod(l,1))<1E-6;
clear h k l

I(isBragg) = NaN;
I = fillin_nans(I);

% calculate the dpdf

[P,MTp] = MT.fourier_transform(I);

%% export

chunksize = MTp.Grid.P; 

MTp.export('export/lys_ortho_goodvibes_delta_pdf.h5','delta_pdf',P,'Datatype','single','ChunkSize',chunksize);

%% FUNCTIONS

function I = fillin_nans(I)

isOK = ~isnan(I);

ker = zeros(3,3,3);
ker([1,3],2,2) = 1;
ker(2,[1,3],2) = 1;
ker(2,2,[1,3]) = 1;

A0 = I;
A0(~isOK) = 0;
n0 = convn(isOK,ker,'same');
A0 = convn(A0,ker,'same')./n0;
isnew = ~isOK & n0>0;
I(isnew) = A0(isnew);

end
