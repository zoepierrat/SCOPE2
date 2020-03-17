function [Output_dir,f] = create_output_files_binary(parameter_file, F, path_of_code,input_path,spectral,options)
%% Create Output dir
string          = clock;
Outdir_Name     = char(F(1).FileName);
Output_dir      = sprintf(['output/',Outdir_Name,'_%4.0f-%02.0f-%02.0f-%02.0f%02.0f/'],[string(1) string(2) string(3) string(4) string(5)]);
warning('off','MATLAB:DELETE:FileNotFound')
if any(~exist(Output_dir,'dir'))
    mkdir(Output_dir)
    mkdir([Output_dir,'Parameters' filesep])
end


%% Log File
for i = 1:length(parameter_file{1})
    copyfile([input_path parameter_file{1}{i}],[Output_dir,'Parameters/', parameter_file{1}{i}],'f')
end
fidpath          = fopen([Output_dir,'Parameters/SCOPEversion.txt'],'w');      % complete path of the SCOPE code
fprintf(fidpath,'%s', path_of_code);

%% Open files for writing
f.pars_file               = fopen([Output_dir,'pars_and_input_short.bin'],'w');
f.veg_file                = fopen([Output_dir,'vegetation.bin'],'w');
if options.calc_fluor
    f.fluor_file              = fopen([Output_dir,'fluorescence_scalars.bin'],'w');
    f.fluor_spectrum_file     = fopen([Output_dir,'fluorescence.bin'],'w');
    f.sigmaF_file             = fopen([Output_dir,'sigmaF.bin'],'w');
    f.fhemis_file             = fopen([Output_dir,'fluorescence_hemis.bin'],'w');
    f.Lo2_file                = fopen([Output_dir,'Lo_spectrum_inclF.bin'],'w');
end
f.r_file                  = fopen([Output_dir,'reflectance.bin'],'w');
f.Eout_file               = fopen([Output_dir,'Eout_spectrum.bin'],'w');
f.Lo_file                 = fopen([Output_dir,'Lo_spectrum.bin'],'w');
f.Esun_file               = fopen([Output_dir,'Esun.bin'], 'w');
f.Esky_file               = fopen([Output_dir,'Esky.bin'],'w');


wlS = spectral.wlS; %#ok<*NASGU>
wlF = spectral.wlF;

save([Output_dir 'wlS.txt'], 'wlS', '-ascii');
save([Output_dir 'wlF.txt'], 'wlF', '-ascii');
