function [soil,leafbio,canopy,meteo,angles,xyt] = select_input_my2(V,vi,canopy,options,constants,xyt,soil)

vi = table2struct(array2table(vi', 'VariableNames', fieldnames(V)));

soil.spectrum      = V.spectrum(vi.spectrum);
soil.rss           = V.rss(vi.rss);
soil.rs_thermal    = V.rs_thermal(vi.rs_thermal);
soil.cs            = V.cs(vi.cs);
soil.rhos          = V.rhos(vi.rhos);
soil.CSSOIL        = V.CSSOIL(vi.CSSOIL);
soil.lambdas       = V.lambdas(vi.lambdas);
soil.rbs           = V.rbs(vi.rbs);
soil.SMC           = V.SMC(vi.SMC);
soil.BSMBrightness = V.BSMBrightness(vi.BSMBrightness);
soil.BSMlat	       = V.BSMlat(vi.BSMlat);
soil.BSMlon	       = V.BSMlon(vi.BSMlon);
leafbio.Cab     = V.Cab(vi.Cab);
leafbio.Cca     = V.Cca(vi.Cca);
%if options.Cca_function_of_Cab
%    leafbio.Cca = 0.25*V(1vi_s.1));
%end
leafbio.Cdm     = V.Cdm(vi.Cdm);
leafbio.Cw      = V.Cw(vi.Cw);
leafbio.Cs      = V.Cs(vi.Cs);
leafbio.Cant    = V.Cant(vi.Cant);
leafbio.N       = V.N(vi.N);
leafbio.Vcmo    = V.Vcmo(vi.Vcmo);
leafbio.m       = V.m(vi.m);
leafbio.BallBerry0 = V.BallBerry0(vi.BallBerry0);
leafbio.Type    = V.Type(vi.Type);
leafbio.Tparam  = V.Tparam;  % everything
leafbio.fqe     = V.fqe(vi.fqe);
leafbio.Rdparam = V.Rdparam(vi.Rdparam);

leafbio.rho_thermal = V.rho_thermal(vi.rho_thermal);
leafbio.tau_thermal = V.tau_thermal(vi.tau_thermal);

leafbio.Tyear         = V.Tyear(vi.Tyear);
leafbio.beta          = V.beta(vi.beta);
leafbio.kNPQs         = V.kNPQs(vi.kNPQs);
leafbio.qLs           = V.qLs(vi.qLs);
leafbio.stressfactor  = V.stressfactor(vi.stressfactor);

canopy.LAI  = max(1E-9, V.LAI(vi.LAI));
canopy.hc  = V.hc(vi.hc);
canopy.LIDFa = V.LIDFa(vi.LIDFa);
canopy.LIDFb  = V.LIDFb(vi.LIDFb);
canopy.leafwidth  = V.leafwidth(vi.leafwidth);
canopy.rb   = V.rb(vi.rb);
canopy.Cd  = V.Cd(vi.Cd);
canopy.CR = V.CR(vi.CR);
canopy.CD1  = V.CD1(vi.CD1);
canopy.Psicor  = V.Psicor(vi.Psicor);
canopy.rwc  = V.rwc(vi.rwc);
canopy.kV = V.kV(vi.kV);
canopy.zo  = V.zo(vi.zo);
canopy.d = V.d(vi.d);
canopy.Cv = V.Cv(vi.Cv);
canopy.crowndiameter = V.crowndiameter(vi.crowndiameter);


meteo.z  = V.z(vi.z);
meteo.Rin   = V.Rin(vi.Rin);
meteo.Ta = V.Ta(vi.Ta);
meteo.Rli  = V.Rli(vi.Rli);
meteo.p  = V.p(vi.p);
meteo.ea  = V.ea(vi.ea);
meteo.u   = V.u(vi.u);
meteo.Ca = V.Ca(vi.Ca);
meteo.Oa  = V.Oa(vi.Oa);

xyt.startDOY = V.startDOY(vi.startDOY);
xyt.endDOY = V.endDOY(vi.endDOY);
xyt.LAT = V.LAT(vi.LAT);
xyt.LON = V.LON(vi.LON);
xyt.timezn = V.timezn(vi.timezn);

angles.tts = V.tts(vi.tts);
angles.tto = V.tto(vi.tto);
angles.psi = V.psi(vi.psi);

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
