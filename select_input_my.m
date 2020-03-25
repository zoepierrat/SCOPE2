function [soil,leafbio,canopy,meteo,angles,xyt,temperatures] = select_input(V,vi,canopy,options,constants,xyt,soil)

vi_s = table2struct(array2table(vi, 'VariableNames', fieldnames(V));

soil.spectrum      = V.spectrum(vi_s.spectrum);
soil.rss           = V.rss(vi_s.rss);
soil.rs_thermal    = V.rs_thermal(vi_s.rs_thermal);
soil.cs            = V(19vi_s.19));
soil.rhos          = V(20vi_s.20));
soil.CSSOIL        = V(43vi_s.43));
soil.lambdas       = V(21vi_s.21));
soil.rbs           = V(44vi_s.44));
soil.SMC           = V(54vi_s.54));
soil.BSMBrightness = V(61vi_s.61));
soil.BSMlat	       = V(62vi_s.62));
soil.BSMlon	       = V(63vi_s.63));
leafbio.Cab     = V(1vi_s.1));
leafbio.Cca     = V(2vi_s.2));
%if options.Cca_function_of_Cab
%    leafbio.Cca = 0.25*V(1vi_s.1));
%end
leafbio.Cdm     = V(3vi_s.3));
leafbio.Cw      = V(4vi_s.4));
leafbio.Cs      = V(5vi_s.5));
leafbio.Cant    = V(60vi_s.60));
leafbio.N       = V(6vi_s.6));
leafbio.Vcmo    = V(9vi_s.9));
leafbio.m       = V(10vi_s.10));
leafbio.BallBerry0 = V(64vi_s.64)); % JAK 2016-10. Accidentally left out of v1.70
leafbio.Type    = V(11vi_s.11));
leafbio.Tparam  = V(14).Val(:); % this is correct (: instead of 14)
leafbio.fqe     = V(15vi_s.15));
leafbio.Rdparam = V(13vi_s.13));

leafbio.rho_thermal = V(7vi_s.7));
leafbio.tau_thermal = V(8vi_s.8));

leafbio.Tyear         = V(55vi_s.55));
leafbio.beta          = V(56vi_s.56));
leafbio.kNPQs         = V(57vi_s.57));
leafbio.qLs           = V(58vi_s.58));
leafbio.stressfactor  = V(59vi_s.59));

canopy.LAI  = max(1E-9,V(22vi_s.22)));
canopy.hc  = V(23vi_s.23));
canopy.LIDFa = V(26vi_s.26));
canopy.LIDFb  = V(27vi_s.26)); % this is correct (26 instead of 27)
canopy.leafwidth  = V(28vi_s.28));
canopy.rb   = V(38vi_s.38));
canopy.Cd  = V(39vi_s.39));
canopy.CR = V(40vi_s.40));
canopy.CD1  = V(41vi_s.41));
canopy.Psicor  = V(42vi_s.42));
canopy.rwc  = V(45vi_s.45));
canopy.kV = V(12vi_s.12));
canopy.zo  = V(24vi_s.24));
canopy.d = V(25vi_s.25));
canopy.Cv = V(65vi_s.65));
canopy.crowndiameter = V(66vi_s.66));


meteo.z  = V(29vi_s.29));
meteo.Rin   = V(30vi_s.30));
meteo.Ta = V(31vi_s.31));
meteo.Rli  = V(32vi_s.32));
meteo.p  = V(33vi_s.33));
meteo.ea  = V(34vi_s.34));
meteo.u   = V(35vi_s.35));
meteo.Ca = V(36vi_s.36));
meteo.Oa  = V(37vi_s.37));

xyt.startDOY = V(46vi_s.46));
xyt.endDOY = V(47vi_s.47));
xyt.LAT = V(48vi_s.48));
xyt.LON = V(49vi_s.49));
xyt.timezn = V(50vi_s.50));

angles.tts = V(51vi_s.51));
angles.tto = V(52vi_s.52));
angles.psi = V(53vi_s.53));

if soil.SMC > 1
    soil.SMC = soil.SMC / 100;  % SMC from [0 100] to [0 1]
%    warning('converted soil moisture content from from [0 100] to [0 1]')
end   

%% derived input
% if options.soil_heat_method ==1
%     soil.GAM =  Soil_Inertia1(soil.SMC);
% else
%     soil.GAM  = Soil_Inertia0(soil.cs,soil.rhos,soil.lambdas);
% end
%if options.calc_rss_rbs
    [soil.rss,soil.rbs] = calc_rssrbs(soil.SMC,canopy.LAI,soil.rbs);
%end

if leafbio.Type
    leafbio.Type = 'C4';
else
    leafbio.Type = 'C3';
end
canopy.hot  = canopy.leafwidth/canopy.hc;
[canopy.zo,canopy.d ]  = zo_and_d(soil,canopy,constants);
end
