function output_data_binary(f,k, xyt, rad,  canopy,V, vi, vmax,options)
%% OUTPUT DATA
% author C. Van der Tol
% date:      30 Nov 2019 

%%
if isdatetime(xyt.t)
    get_doy = @(x) juliandate(x) - juliandate(datetime(year(x), 1, 0));
    V(46).Val = get_doy(io.timestamp2datetime(xyt.startDOY));
    V(47).Val = get_doy(io.timestamp2datetime(xyt.endDOY));
    xyt.t = get_doy(xyt.t);
end

%% Vegetation products
veg_out = [k xyt.year(k) xyt.t(k) canopy.Pntot canopy.Pntot_Cab canopy.Rntot_Cab canopy.A canopy.Ja canopy.ENPQ  canopy.LST];
fwrite(f.veg_file,veg_out,'double');

%% Fluorescence scalar outputs
if options.calc_fluor
    fluor_out = [rad.F685  rad.wl685 rad.F740 rad.wl740 rad.F687 rad.F760 ...
        rad.LoutF rad.EoutF rad.EoutFrc];
    fwrite(f.fluor_file,fluor_out,'double');
    
    %% Fluorescence spectral outputs
    % fluorescence radiance (L) in observation direction [mW m-2 nm-1 sr-1]
    fwrite(f.fluor_spectrum_file, rad.LoF_, 'double');
    fwrite(f.sigmaF_file, rad.sigmaF, 'double');
    fwrite(f.fhemis_file,rad.EoutF_, 'double');
end

%% Other radiance and reflectance outputs
fwrite(f.r_file,canopy.reflectance,'double');
fwrite(f.Eout_file,rad.Eout_,'double');
fwrite(f.Lo_file, rad.Lotot_, 'double');
if options.calc_fluor
    fwrite(f.Lo2_file, rad.Lototf_,'double');
end
fwrite(f.Esun_file, rad.Esun_, 'double');
fwrite(f.Esky_file, rad.Esky_, 'double');

k2 = find(vmax>1);  % really needed for the first one, later vi > 1
V_short = nan(1,length(k2)+1);
V_short(1) = length(k2);
for i = 1:length(k2)
    V_short(i+1) = V(k2(i)).Val(vi(k2(i)));
end
fwrite(f.pars_file, V_short,'double');