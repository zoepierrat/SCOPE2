function bin_to_csv(Output_dir,V,vmax)

f.pars_file               = fopen([Output_dir,'pars_and_input_short.bin'],'r');
f.veg_file                = fopen([Output_dir,'vegetation.bin'],'r');
f.fluor_file              = fopen([Output_dir,'fluorescence_scalars.bin'],'r');
f.fluor_spectrum_file     = fopen([Output_dir,'fluorescence.bin'],'r');
f.sigmaF_file             = fopen([Output_dir,'sigmaF.bin'],'r');
f.fhemis_file             = fopen([Output_dir,'fluorescence_hemis.bin'],'r');
f.r_file                  = fopen([Output_dir,'reflectance.bin'],'r');
f.Eout_file               = fopen([Output_dir,'Eout_spectrum.bin'],'r');
f.Lo_file                 = fopen([Output_dir,'Lo_spectrum.bin'],'r');
f.Lo2_file                = fopen([Output_dir,'Lo_spectrum_inclF.bin'],'r');
f.Esun_file               = fopen([Output_dir,'Esun.bin'], 'r');
f.Esky_file               = fopen([Output_dir,'Esky.bin'],'r');


g.pars_file               = fopen([Output_dir,'pars_and_input_short.csv'],'w');
g.veg_file                = fopen([Output_dir,'vegetation.csv'],'w');
g.r_file                  = fopen([Output_dir,'reflectance.csv'],'w');
g.Lo_file                 = fopen([Output_dir,'Lo_spectrum.csv'],'w');
g.Esun_file               = fopen([Output_dir,'Esun.csv'], 'w');
g.Esky_file               = fopen([Output_dir,'Esky.csv'],'w');
g.Eout_file               = fopen([Output_dir,'Eout_spectrum.csv'],'w');
if f.fluor_file>0
    g.fluor_file              = fopen([Output_dir,'fluorescence_scalars.csv'],'w');
    g.fluor_spectrum_file     = fopen([Output_dir,'fluorescence.csv'],'w');
    g.sigmaF_file             = fopen([Output_dir,'sigmaF.csv'],'w');
    g.fhemis_file             = fopen([Output_dir,'fluorescence_hemis.csv'],'w');
    g.Lo2_file                = fopen([Output_dir,'Lo_spectrum_inclF.csv'],'w');
end

pars            = fread(f.pars_file, 'double');
npars           = pars(1);
pars2           = reshape(pars,[npars+1, (length(pars))/(npars+1)]);
pars2           = pars2';
ns              = size(pars2,1);

for j = find(vmax>1)
    fprintf(g.pars_file,'%s,',V(vmax>1).Name);
end
fprintf(g.pars_file,' \n');
for j = 1:size(pars2,1)
    fprintf(g.pars_file,'%d,',pars2(j,2:end-1));
    fprintf(g.pars_file,'%d',pars2(j,end));
    fprintf(g.pars_file,'\n');
end
clear('pars2')

fprintf(g.veg_file,'%s','simulation_number, year, DoY, aPAR, aPARbyCab, aPARbyCab(energyunits), Photosynthesis, Electron_transport, NPQ_energy, LST');
fprintf(g.veg_file,'\n');
fprintf(g.veg_file,'%s',', , , umol m-2 s-1, umol m-2 s-1, Wm-2, umol m-2 s-1, umol m-2 s-1, Wm-2, K');
fprintf(g.veg_file,'\n');
veg_out = fread(f.veg_file,'double');
for k = 1:ns
    fprintf(g.veg_file,'%d,',veg_out(1+10*(k-1):10*k));
    fprintf(g.veg_file,'\n');
end
clear('veg_out')
if f.fluor_file>0
    fprintf(g.fluor_file,'%s','F_1stpeak, wl_1stpeak, F_2ndpeak, wl_2ndpeak, F687, F760, LFtot, EFtot, EFtot_RC');
    fprintf(g.fluor_file,'\n');
    fprintf(g.fluor_file,'%s','Wm-2um-1sr-1,nm,Wm-2um-1sr-1,nm,Wm-2um-1sr-1,Wm-2um-1sr-1,Wm-2sr-1,Wm-2,Wm-2');
    fprintf(g.fluor_file,'\n');
    fluor_out = fread(f.fluor_file,'double');
    for k = 1:ns
        fprintf(g.fluor_file,'%d,',fluor_out(1+9*(k-1):9*k));
        fprintf(g.fluor_file,'\n');
    end
    fprintf(g.fluor_spectrum_file,'%s','fluorescence_spectrum 640:1:850 nm');
    fprintf(g.fluor_spectrum_file,'\n');
    fprintf(g.fluor_spectrum_file,'%s','Wm-2um-1sr-1');
    fprintf(g.fluor_spectrum_file,'\n');
    fluor_spectrum_out = fread(f.fluor_spectrum_file,'double');
    for k = 1:ns
        fprintf(g.fluor_spectrum_file,'%d,',fluor_spectrum_out(1+211*(k-1):211*k-1));
        fprintf(g.fluor_spectrum_file,'%d',fluor_spectrum_out(211*k));
        fprintf(g.fluor_spectrum_file,'\n');
    end
    clear('fluor_spectrum_out')
    
    
    fprintf(g.sigmaF_file,'%s','escape probability 640:1:850 nm');
    fprintf(g.sigmaF_file,'\n');
    fprintf(g.sigmaF_file,'%s','sr-1');
    fprintf(g.sigmaF_file,'\n');
    sigmaF_out = fread(f.sigmaF_file,'double');
    for k = 1:ns
        fprintf(g.sigmaF_file,'%d,',sigmaF_out(1+211*(k-1):211*k-1));
        fprintf(g.sigmaF_file,'%d',sigmaF_out(211*k));
        fprintf(g.sigmaF_file,'\n');
    end
    clear('sigmaF_out');
    
    fprintf(g.fhemis_file,'%s','fluorescence_spectrum 640:1:850 nm hemispherically integrated');
    fprintf(g.fhemis_file,'\n');
    fprintf(g.fhemis_file,'%s','Wm-2um-1');
    fprintf(g.fhemis_file,'\n');
    fhemis_out = fread(f.fhemis_file,'double');
    for k = 1:ns
        fprintf(g.fhemis_file,'%d,',fhemis_out(1+211*(k-1):211*k-1));
        fprintf(g.fhemis_file,'%d',fhemis_out(211*k));
        fprintf(g.fhemis_file,'\n');
    end
    clear('fhemis_out');
    
    fprintf(g.Lo2_file,'%s','upwelling radiance including fluorescence');
    fprintf(g.Lo2_file,'\n');
    fprintf(g.Lo2_file,'%s','W m-2 um-1 sr-1');
    fprintf(g.Lo2_file,'\n');
    Lo2_out = fread(f.Lo2_file,'double');
    for k = 1:ns
        fprintf(g.Lo2_file,'%d,',Lo2_out(1+2162*(k-1):2162*k-1));
        fprintf(g.Lo2_file,'%d',Lo2_out(2162*k));
        fprintf(g.Lo2_file,'\n');
    end
    clear('Lo2_file');
end


fprintf(g.r_file,'%s','reflectance');
fprintf(g.r_file,'\n');
fprintf(g.r_file,'%s','pi*upwelling radiance/irradiance');
fprintf(g.r_file,'\n');
r_out = fread(f.r_file,'double');
for k = 1:ns
    fprintf(g.r_file,'%d,',r_out(1+2162*(k-1):2162*k-1));
    fprintf(g.r_file,'%d',r_out(2162*k));
    fprintf(g.r_file,'\n');
end
clear('r_out');

fprintf(g.Eout_file,'%s','hemispherically integrated upwelling radiance');
fprintf(g.Eout_file,'\n');
fprintf(g.Eout_file,'%s','W m-2 um-1');
fprintf(g.Eout_file,'\n');
Eout_out = fread(f.Eout_file,'double');
for k = 1:ns
    fprintf(g.Eout_file,'%d,',Eout_out(1+2162*(k-1):2162*k-1));
    fprintf(g.Eout_file,'%d',Eout_out(2162*k));
    fprintf(g.Eout_file,'\n');
end
clear('Eout_out');

fprintf(g.Lo_file,'%s','upwelling radiance excluding fluorescence');
fprintf(g.Lo_file,'\n');
fprintf(g.Lo_file,'%s','W m-2 um-1 sr-1');
fprintf(g.Lo_file,'\n');
Lo_out = fread(f.Lo_file,'double');
for k = 1:ns
    fprintf(g.Lo_file,'%d,',Lo_out(1+2162*(k-1):2162*k-1));
    fprintf(g.Lo_file,'%d',Lo_out(2162*k));
    fprintf(g.Lo_file,'\n');
end
clear('Lo_out');


fprintf(g.Esun_file,'%s','direct solar irradiance');
fprintf(g.Esun_file,'\n');
fprintf(g.Esun_file,'%s','W m-2 um-1 sr-1');
fprintf(g.Esun_file,'\n');
Esun_out = fread(f.Esun_file,'double');
for k = 1:ns
    fprintf(g.Esun_file,'%d,',Esun_out(1+2162*(k-1):2162*k-1));
    fprintf(g.Esun_file,'%d',Esun_out(2162*k));
    fprintf(g.Esun_file,'\n');
end
clear('Esun_out');

fprintf(g.Esky_file,'%s','direct solar irradiance');
fprintf(g.Esky_file,'\n');
fprintf(g.Esky_file,'%s','W m-2 um-1 sr-1');
fprintf(g.Esky_file,'\n');
Esky_out = fread(f.Esky_file,'double');
for k = 1:ns
    fprintf(g.Esky_file,'%d,',Esky_out(1+2162*(k-1):2162*k-1));
    fprintf(g.Esky_file,'%d',Esky_out(2162*k));
    fprintf(g.Esky_file,'\n');
end
clear('Esky_out');
fclose('all');

%% deleting .bin
delete(fullfile(Output_dir, '*.bin'))
end
