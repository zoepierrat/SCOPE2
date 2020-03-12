function [rad,gap] = RTMo(spectral,atmo,soil,leafopt,canopy,angles,constants,meteo,options)

%% 0. Preparations
deg2rad = constants.deg2rad;    % degree to rad
wl      = spectral.wlS;         % SCOPE wavelengths as a column-vector
nwl     = length(wl);           %
wlPAR   = spectral.wlPAR;       % PAR wavelength range
minPAR  = min(wlPAR);           % min PAR
maxPAR  = max(wlPAR);           % max PAR
Ipar    = find(wl>=minPAR & wl<=maxPAR); % Indices for PAR wavelenghts within wl
tts     = angles.tts;           % solar zenith angle
tto     = angles.tto;           % observer zenith angle
psi     = angles.psi;           % relative azimuth anglee

nl      = canopy.nlayers;       % number of canopy layers (nl)
litab   = canopy.litab;         % SAIL leaf inclibation angles % leaf inclination angles PY
lazitab = canopy.lazitab;       % leaf azimuth angles relative to the sun
%nli     = canopy.nlincl;        % numler of leaf inclinations (13)
nlazi   = canopy.nlazi;         % number of azimuth angles (36)
LAI     = canopy.LAI;           % leaf area index
lidf    = canopy.lidf;          % leaf inclination distribution function
x       = canopy.x;             % all levels except for the top
dx      = 1/nl;

kChlrel =   leafopt.kChlrel;
rho     =   leafopt.refl;
tau     =   leafopt.tran;
rs      =   soil.refl;          % [nwl,nsoils] soil reflectance spectra
epsc    =   1-rho-tau;          % [nl,nwl]        emissivity of leaves
epss    =   1-rs;              % [nwl]        emissivity of soil
iLAI    =   LAI/nl;             % [1]          LAI of elementary layer
xl      =   [0; x];             % [nl+1]       all levels + soil

% 0.2 initialisations (allocation of memory)
[Es_,Emin_,Eplu_]           = deal(zeros(nl+1,nwl));       % [nl+1,nwl]     direct, up and down diff. rad.

%% 1. Geometric quantities
% 1.1 general geometric quantities. these variables are scalars
cos_tts     = cos(tts*deg2rad);             %           cos solar       angle
tan_tto     = tan(tto*deg2rad);             %           tan observation angle

cos_tto     = cos(tto*deg2rad);             %           cos observation angle
sin_tts     = sin(tts*deg2rad);             %           sin solar       angle
tan_tts     = tan(tts*deg2rad);             %           tan observation angle

psi         = abs(psi-360*round(psi/360));  %           (to ensure that volscatt is symmetric for psi=90 and psi=270)
dso         = sqrt(tan_tts.^2 + tan_tto.^2 - 2*tan_tts.*tan_tto.*cos(psi*deg2rad));

% 1.2 geometric factors associated with extinction and scattering
[chi_s,chi_o,frho,ftau]=volscat(tts,tto,psi,litab);   % volume scattering

cos_ttlo    = cos(lazitab*deg2rad);         % [1,36]    cos leaf azimuth angles

cos_ttli    = cos(litab*deg2rad);           % [13]      cos leaf angles
sin_ttli    = sin(litab*deg2rad);           % [13]      sinus leaf angles

ksli        = chi_s./cos_tts;               % [13]      p306{1} extinction coefficient in direction of sun        per leaf angle
koli        = chi_o./cos_tto;               % [13]      p307{1} extinction coefficient in direction of observer   per leaf angle

sobli       = frho*pi/(cos_tts*cos_tto);    % [13]      pag 309{1} area scattering coefficient fractions
sofli       = ftau*pi/(cos_tts*cos_tto);    % [13]      pag 309{1}
bfli        = cos_ttli.^2;                  % [13]

%integration over angles (using a vector inproduct) -> scalars
k           = ksli'*lidf;                   %           pag 306{1}    extinction coefficient in direction of sun.
K           = koli'*lidf;                   %           pag 307{1}    extinction coefficient in direction of observer
bf          = bfli'*lidf;                   %
sob         = sobli'*lidf;                  %           weight of specular2directional back    scatter coefficient
sof         = sofli'*lidf;                  %           weight of specular2directional forward scatter coefficient
% 1.3 geometric factors to be used later with rho and tau, f1 f2 of pag 304:
% these variables are scalars
sdb         = 0.5*(k+bf);                   % fs*f1
sdf         = 0.5*(k-bf);                   % fs*f2     weight of specular2diffuse     foward  scatter coefficient
ddb         = 0.5*(1+bf);                   % f1^2+f2^2 weight of diffuse2diffuse      back    scatter coefficient
ddf         = 0.5*(1-bf);                   % 2*f1*f2   weight of diffuse2diffuse      forward scatter coefficient
dob         = 0.5*(K+bf);                   % fo*f1     weight of diffuse2directional  back    scatter coefficient
dof         = 0.5*(K-bf);                   % fo*f2     weight of diffuse2directional  forward scatter coefficient

% 1.4 solar irradiance factor for all leaf orientations
Cs          = cos_ttli*cos_tts;             % [nli]     pag 305 modified by Joris
Ss          = sin_ttli*sin_tts;             % [nli]     pag 305 modified by Joris

cos_deltas  = Cs*ones(1,nlazi) + Ss*cos_ttlo;  % [nli,nlazi]
fs          = abs(cos_deltas/cos_tts);         % [nli,nlazi] pag 305

% 1.5 probabilities Ps, Po, Pso
Ps          =   exp(k*xl*LAI);                                              % [nl+1]  p154{1} probability of viewing a leaf in solar dir
Po          =   exp(K*xl*LAI);                                              % [nl+1]  p154{1} probability of viewing a leaf in observation dir

Ps(1:nl)    =   Ps(1:nl) *(1-exp(-k*LAI*dx))/(k*LAI*dx);                                      % Correct Ps/Po for finite dx
Po(1:nl)    =   Po(1:nl) *(1-exp(-K*LAI*dx))/(K*LAI*dx);  % Correct Ps/Po for finite dx
q           =   canopy.hot;
Pso         =   zeros(size(Po));
for j=1:length(xl)
    Pso(j,:)=   quad(@(y)Psofunction(K,k,LAI,q,dso,y),xl(j)-dx,xl(j))/dx; %#ok<FREMO>
end

Pso(Pso>Po)= min([Po(Pso>Po),Ps(Pso>Po)],[],2);    %takes care of rounding error
Pso(Pso>Ps)= min([Po(Pso>Ps),Ps(Pso>Ps)],[],2);    %takes care of rounding error
gap.Pso      = Pso;

%% 2. Calculation of reflectance
% 2.1  reflectance, transmittance factors in a thin layer
% the following are vectors with lenght [nl,nwl]
sigb        = ddb*rho + ddf*tau;            % [nl,nwl]     sigmab, p305{1} diffuse     backscatter scattering coefficient for diffuse  incidence
sigf        = ddf*rho + ddb*tau;            % [nl,nwl]     sigmaf, p305{1} diffuse     forward     scattering coefficient for forward  incidence
sb          = sdb*rho + sdf*tau;            % [nl,nwl]     sb,     p305{1} diffuse     backscatter scattering coefficient for specular incidence
sf          = sdf*rho + sdb*tau;            % [nl,nwl]     sf,     p305{1} diffuse     forward     scattering coefficient for specular incidence
vb          = dob*rho + dof*tau;            % [nl,nwl]     vb,     p305{1} directional backscatter scattering coefficient for diffuse  incidence
vf          = dof*rho + dob*tau;            % [nl,nwl]     vf,     p305{1} directional forward     scattering coefficient for diffuse  incidence
w           = sob*rho + sof*tau;            % [nl,nwl]     w,      p309{1} bidirectional scattering coefficent (directional-directional)
a           = 1-sigf;                       % [nl,nwl]     attenuation

% the following six lines may be removed in a next version when RTMt is revised
m           = sqrt(a.*a-sigb.*sigb);% 
rinf        = (a-m)./sigb;         %               Reflection coefficient for infinite thick canopy    
rad.fHs     = (1-rinf.*rinf).*(1-rs(end))./(1-rinf.*rs);
rad.fHc     = iLAI*m.*(1-rinf);
rad.fbottom = (rs-rinf)./(1-rs.*rinf);
rad.rinf    = rinf;
rad.m       = m;

%2.2. convert scattering and extinction coefficients into thin layer reflectances and transmittances
% Eq. 17 in mSCOPE paper
tau_ss = 1-k.*iLAI;
tau_dd = 1-a.*iLAI;
tau_sd = sf.*iLAI;
rho_sd = sb.*iLAI;
rho_dd = sigb.*iLAI;

%2.3. reflectance calcualtion
% Eq. 18 in mSCOPE paper
[R_sd,R_dd] =   deal(zeros(nl+1,nwl));      % surface reflectance
[Xsd,Xdd]   =   deal(zeros(nl,nwl));        % Effective transmittance
Xss         =   zeros(1,nl);                % Effective transmittance
R_sd(nl+1,:)   =   rs;
R_dd(nl+1,:)   =   rs;
for j=nl:-1:1 % from bottom to top. note nl+1 the background. 1 is the top of canopy.
    Xss      = tau_ss;
    dnorm       = 1-rho_dd'.*R_dd(j+1,:);
    Xsd(j,:)    = (tau_sd'+tau_ss.*R_sd(j+1,:).*rho_dd')./dnorm;
    Xdd(j,:)    = tau_dd'./dnorm;
    R_sd(j,:)   = rho_sd'+tau_dd'.*(R_sd(j+1,:).*Xss+R_dd(j+1,:).*Xsd(j,:));
    R_dd(j,:)   = rho_dd'+tau_dd'.*R_dd(j+1,:).*Xdd(j,:);
end

%% 3. Flux calculation
% 3.1 calculation of incoming light Esun_ and Esky_
% Extract MODTRAN atmosphere parameters at the SCOPE wavelengths

Fd      = zeros(nwl,1);
Ls      = Planck(wl,meteo.Ta+273.15);
% Solar and sky irradiance using 6 atmosperic functions
rdd     = R_dd(1,:)';
rsd     = R_sd(1,:)';

if ~isfield(atmo,'Esun_')
    
    t1  = atmo.M(:,1);
    t3  = atmo.M(:,2);
    t4  = atmo.M(:,3);
    t5  = atmo.M(:,4);
    t12 = atmo.M(:,5);
    t16 = atmo.M(:,6);
    
    % radiation fluxes, downward and upward (these all have dimenstion [nwl]
    % first calculate hemispherical reflectances rsd and rdd according to SAIL
    % these are assumed for the reflectance of the surroundings
    % rdo is computed with SAIL as well
    % assume Fd of surroundings = 0 for the momemnt
    % initial guess of temperature of surroundings from Ta;

    Esun_   = pi*t1.*t4;
    Esky_   = pi./(1-t3.*rdd).*(t1.*(t5+t12.*rsd)+Fd+(1-rdd).*Ls.*t3+t16);
    
    % fractional contributions of Esun and Esky to total incident radiation in
    % optical and thermal parts of the spectrum
    if meteo.Rin ~= -999
        % fractional contributions of Esun and Esky to total incident radiation in
        % optical and thermal parts of the spectrum
        
        [fEsuno,fEskyo,fEsunt,fEskyt]          = deal(0*Esun_);   %initialization
        
        J_o             = wl<3000;                          %find optical spectrum
        Esunto          = 0.001 * Sint(Esun_(J_o),wl(J_o)); %Calculate optical sun fluxes (by Integration), including conversion mW >> W
        Eskyto          = 0.001 * Sint(Esky_(J_o),wl(J_o)); %Calculate optical sun fluxes (by Integration)
        Etoto           = Esunto + Eskyto;                  %Calculate total fluxes
        fEsuno(J_o)     = Esun_(J_o)/Etoto;                 %fraction of contribution of Sun fluxes to total light
        fEskyo(J_o)     = Esky_(J_o)/Etoto;                 %fraction of contribution of Sky fluxes to total light
        
        J_t             = wl>=3000;                         %find thermal spectrum
        Esuntt          = 0.001 * Sint(Esun_(J_t),wl(J_t)); %Themal solar fluxes
        Eskytt          = 0.001 * Sint(Esky_(J_t),wl(J_t)); %Thermal Sky fluxes
        Etott           = Eskytt + Esuntt;                  %Total
        fEsunt(J_t)     = Esun_(J_t)/Etott;                 %fraction from Esun
        fEskyt(J_t)     = Esky_(J_t)/Etott;                 %fraction from Esky
        
        Esun_(J_o) = fEsuno(J_o)*meteo.Rin;
        Esky_(J_o) = fEskyo(J_o)*meteo.Rin;
        Esun_(J_t) = fEsunt(J_t)*meteo.Rli;
        Esky_(J_t) = fEskyt(J_t)*meteo.Rli;
    end
    
else
    Esun_ = atmo.Esun_;
    Esky_ = atmo.Esky_;
end

% 3.2 flux profile calculation
% Eq. 19 in mSCOPE paper
Es_(1,:)       = Esun_;
Emin_(1,:)     = Esky_;

for j=1:nl % from top to bottom
    Es_(j+1,:)    =   Xss.*Es_(j,:);
    Emin_(j+1,:)  =   Xsd(j,:).*Es_(j,:)+Xdd(j,:).*Emin_(j,:);
    Eplu_(j,:)    =   R_sd(j,:).*Es_(j,:)+R_dd(j,:).*Emin_(j,:);
end
Eplu_(j+1,:)    =   rs'.*(Es_(j+1,:)+Emin_(j+1,:)); 
%CvdT added calculation of Eplu_(j+1,:)

% 3.3 outgoing fluxes, hemispherical and in viewing direction, spectrum
% hemispherical, spectral
Eout_   = Eplu_(1,:)';                  %           [nwl]
% in viewing direction, spectral
piLoc_      = (vb.*(Emin_(1:nl,:)'*Po(1:nl)) +...
    vf.*(Eplu_(1:nl,:)'*Po(1:nl)) +...
    w.*Esun_*sum(Pso(1:nl)))*iLAI;
piLos_      = rs.*Emin_(nl+1,:)'*Po(nl+1) + rs.*Esun_*Pso(nl+1);
piLo_       = piLoc_ + piLos_;              %           [nwl]
Lo_         = piLo_/pi;

Refl        = piLo_./(Esky_+Esun_); % rso and rdo are not computed separately
           
% up and down and hemispherical out, cumulative over wavelenght
Eouto       = 0.001 * Sint(Eout_(spectral.IwlP),spectral.wlP);   %     [1] hemispherical out, in optical range (W m-2)
Eoutt       = 0.001 * Sint(Eout_(spectral.IwlT),spectral.wlT);   %     [1] hemispherical out, in thermal range (W m-2)
Lot         = 0.001 * Sint(Lo_(spectral.IwlT),spectral.wlT);   %     [1] hemispherical out, in thermal range (W m-2)

%% 4. net fluxes, spectral and total, and incoming fluxes
%4.1 incident PAR at the top of canopy, spectral and spectrally integrated
P_          = e2phot(wl(Ipar)*1E-9,(Esun_(Ipar)+Esky_(Ipar)),constants);      %
P           = .001 * Sint(P_,wlPAR);                                % mol m-2s-1
%Psun        = 0.001 * Sint(e2phot(wlPAR*1E-9,Esun_(Ipar),constants),wlPAR);   % Incident solar PAR in PAR units
% Incident and absorbed solar radiation

%4.2 Absorbed radiation
%    absorbed radiation in Wm-2         (Asun)
%    absorbed PAR in mol m-2s-1         (Pnsun)
%    absorbed PAR in Wm-2               (Rnsun_PAR)
%    absorbed PAR by Chl in mol m-2s-1  (Pnsun_Cab)

%Psky        = 0.001 * Sint(e2phot(wlPAR*1E-9,Esky_(Ipar)),wlPAR);
Asun        = 0.001 * Sint(Esun_.*epsc,wl);                         % Total absorbed solar radiation
Pnsun       = 0.001 * Sint(e2phot(wlPAR*1E-9,Esun_(Ipar).*epsc(Ipar),constants),wlPAR);  % Absorbed solar radiation  in PAR range in moles m-2 s-1
Rnsun_Cab   = 0.001 * Sint(Esun_(Ipar).*epsc(Ipar).*kChlrel(Ipar),wlPAR);
Pnsun_Cab   = 0.001 * Sint(e2phot(wlPAR*1E-9,kChlrel(Ipar).*Esun_(Ipar).*epsc(Ipar),constants),wlPAR);
Rnsun_PAR   = 0.001 * Sint(Esun_(Ipar).*epsc(Ipar),wlPAR);

%4.3 total direct radiation (incident and net) per leaf area (W m-2 leaf)

% total direct radiation (incident and net) per leaf area (W m-2 leaf)
%Pdir        = fs * Psun;                        % [13 x 36]   incident
if options.lite
    fs      = lidf'*mean(fs,2);%
end
Rndir       = fs * Asun;                        %    net
Pndir       = fs * Pnsun;                       %    net PAR
Pndir_Cab   = fs * Pnsun_Cab;                   %    net PAR Cab
Rndir_Cab   = fs * Rnsun_Cab;                   %   net PAR energy units
Rndir_PAR   = fs * Rnsun_PAR;                   %  net PAR energy units

%4.4 total diffuse radiation (net) per leaf area (W m-2 leaf)
% diffuse incident radiation for the present layer 'j' (mW m-2 um-1)
E_          = .5*(Emin_(1:end-1,:) + Emin_(2:end,:)+ Eplu_(1:end-1,:)+ Eplu_(2:end,:));
%Pdif        = .001 * Sint(e2phot(wlPAR*1E-9,E_(:,Ipar)',constants),wlPAR);
% canopy layers, diffuse radiation

% net radiation (mW m-2 um-1) and net PAR (moles m-2 s-1 um-1), per wavelength
Rndif_        = E_.*epsc';                                                   % [nl,nwl]  Net (absorbed) radiation by leaves
Pndif_         = .001 *(e2phot(wlPAR*1E-9, Rndif_(:,Ipar)',constants))';                % [nl,nwl]  Net (absorbed) as PAR photons
Pndif_Cab_     = .001 *(e2phot(wlPAR*1E-9, kChlrel(Ipar).*Rndif_(:,Ipar)',constants))';  % [nl,nwl]  Net (absorbed) as PAR photons by Cab
Rndif_Cab_     = kChlrel(Ipar)'.*Rndif_(:,Ipar);  % [nl,nwlPAR]  Net (absorbed) as PAR energy
Rndif_PAR_     = Rndif_(:,Ipar);  % [nl,nwlPAR]  Net (absorbed) as PAR energy

% net radiation (W m-2) and net PAR (moles m-2 s-1), integrated over all wavelengths
Rndif            = .001 * Sint(Rndif_',wl);              % [nl]  Full spectrum net diffuse flux
Pndif           =        Sint(Pndif_(:,Ipar)',wlPAR);        % [nl]  Absorbed PAR
Pndif_Cab        =        Sint(Pndif_Cab_(:,Ipar)',wlPAR);    % [nl]  Absorbed PAR by Cab integrated
Rndif_Cab        = .001 * Sint(Rndif_Cab_(:,Ipar)',wlPAR);      % [nl]  Absorbed PAR by Cab integrated
Rndif_PAR        = .001 * Sint(Rndif_PAR_(:,Ipar)',wlPAR);    % [nl]  Absorbed PAR by Cab integrated

% soil layer, direct and diffuse radiation
Rndirsoil   = .001 * Sint(Esun_.*epss,wl);          % [1] Absorbed solar flux by the soil
Rndifsoil   = .001 * Sint(Emin_(nl+1,:).*epss',wl); % [1] Absorbed diffuse downward flux by the soil (W m-2)

% net (n) radiation R and net PAR P per component: sunlit (u), shaded (h) soil(s) and canopy (c),
% [W m-2 leaf or soil surface um-1]
Rnhc        = Rndif;            % [nl] shaded leaves or needles
Pnhc        = Pndif;            % [nl] shaded leaves or needles
Pnhc_Cab    = Pndif_Cab;        % [nl] shaded leaves or needles
Rnhc_Cab    = Rndif_Cab;        % [nl] shaded leaves or needles
Rnhc_PAR    = Rndif_PAR;        % [nl] shaded leaves or needles

if ~options.lite
    [Rnuc,Pnuc,Pnuc_Cab,Rnuc_PAR,Rnuc_Cab] = deal(0*Rndir);
    for j = 1:nl
        %Puc(:,:,j)  = Pdir      + Pdif(j);      % [13,36,nl] Total fluxes on sunlit leaves or needles
        Rnuc(:,:,j) = Rndir     + Rndif(j);     % [13,36,nl] Total fluxes on sunlit leaves or needles
        Pnuc(:,:,j)  = Pndir     + Pndif(j);     % [13,36,nl] Total fluxes on sunlit leaves or needles
        Pnuc_Cab(:,:,j)  = Pndir_Cab  + Pndif_Cab(j);% [13,36,nl] Total fluxes on sunlit leaves or needles
        Rnuc_PAR(:,:,j)  = Rndir_PAR  + Rndif_PAR(j);% [13,36,nl] Total fluxes on sunlit leaves or needles
        Rnuc_Cab(:,:,j)  = Rndir_Cab  + Rndif_Cab(j);% [13,36,nl] Total fluxes on sunlit leaves or needles   
    end
else   
    %    Puc(:,:,j)  = Pdir      + Pdif(j);      % [13,36,nl] Total fluxes on sunlit leaves or needles
    Rnuc = Rndir     + Rndif;     % [13,36,nl] Total fluxes on sunlit leaves or needles
    Pnuc  = Pndir     + Pndif;     % [13,36,nl] Total fluxes on sunlit leaves or needles
    Pnuc_Cab  = Pndir_Cab  + Pndif_Cab;% [13,36,nl] Total fluxes on sunlit leaves or needles
    Rnuc_Cab  = Rndir_Cab  + Rndif_Cab;% [13,36,nl] Total fluxes on sunlit leaves or needles
    Rnuc_PAR  = Rndir_PAR  + Rndif_PAR;% [13,36,nl] Total fluxes on sunlit leaves or needles       
end
Rnus        = Rndifsoil + Rndirsoil; % [1] sunlit soil 
Rnhs        = Rndifsoil;  % [1] shaded soil
%% 5 place output in structure rad
gap.k       = k;        % extinction cofficient in the solar direction
gap.K       = K;        % extinction cofficient in the viewing direction
gap.Ps      = Ps;       % gap fraction in the solar direction
gap.Po      = Po;       % gap fraction in the viewing direction
rad.rsd     = rsd;      % TOC directional-hemispherical reflectance
rad.rdd     = rdd;      % TOC hemispherical-hemispherical reflectance
%rad.rdo     = Refl;     % TOC hemispherical-directional reflectance %not computed
rad.rso     = Refl;     % TOC directional-directional reflectance   %not computed
%rad.refl    = Refl;     % TOC reflectance
rad.rho_dd  = rho_dd;   % diffuse-diffuse reflectance for the thin layers
rad.tau_dd  = tau_dd;   % diffuse-diffuse transmittance for the thin layers
rad.rho_sd  = rho_sd;   % direct-diffuse reflectance for the thin layers
rad.tau_ss  = tau_ss;   % direct-direct transmittance for the thin layers
rad.tau_sd  = tau_sd;   % direct-diffuse transmittance for the thin layers
rad.R_sd    = R_sd;
rad.R_dd    = R_dd;

rad.vb      = vb;       % directional backscatter coefficient for diffuse incidence
rad.vf      = vf;       % directional forward scatter coefficient for diffuse incidence
rad.sigf    = sigf;     % forward scatter coefficient for specular flux
rad.sigb    = sigb;     % backscatter coefficient for specular flux

rad.Esun_   = Esun_;    % [2162x1 double]   incident solar spectrum (mW m-2 um-1)
rad.Esky_   = Esky_;    % [2162x1 double]   incident sky spectrum (mW m-2 um-1)
rad.PAR     = P;        % [1 double]        incident spectrally integrated PAR (moles m-2 s-1)

rad.Eplu_   = Eplu_;   % [61x2162 double]  upward diffuse radiation in the canopy (mW m-2 um-1)
rad.Emin_   = Emin_;   % [61x2162 double]  downward diffuse radiation in the canopy (mW m-2 um-1)

rad.Lo_     = Lo_;      % [2162x1 double]   TOC radiance in observation direction (mW m-2 um-1 sr-1)
rad.Eout_   = Eout_;    % [2162x1 double]   TOC upward radiation (mW m-2 um-1)
rad.Eouto   = Eouto;        % [1 double]        TOC spectrally integrated upward optical ratiation (W m-2)
rad.Eoutt   = Eoutt;        % [1 double]        TOC spectrally integrated upward thermal ratiation (W m-2)
rad.Lot     = Lot;

rad.Rnhs    = Rnhs;     % [1 double]        net radiation (W m-2) of shaded soil
rad.Rnus    = Rnus;     % [1 double]        net radiation (W m-2) of sunlit soil
rad.Rnhc    = Rnhc;     % [60x1 double]     net radiation (W m-2) of shaded leaves
rad.Rnuc    = Rnuc;     % [13x36x60 double] net radiation (W m-2) of sunlit leaves

rad.Pnh     = 1E6*Pnhc;     % [60x1 double]     net PAR (moles m-2 s-1) of shaded leaves
rad.Pnu     = 1E6*Pnuc;     % [13x36x60 double] net PAR (moles m-2 s-1) of sunlit leaves
rad.Pnh_Cab = 1E6*Pnhc_Cab;% [60x1 double]      net PAR absorbed by Cab (moles m-2 s-1) of shaded leaves
rad.Pnu_Cab = 1E6*Pnuc_Cab; % [13x36x60 double] net PAR absorbed by Cab (moles m-2 s-1) of sunlit leaves
rad.Rnh_Cab = Rnhc_Cab; % [60x1 double]    net PAR absorbed by Cab (W m-2) of shaded leaves
rad.Rnu_Cab = Rnuc_Cab; % [13x36x60 double] net PAR absorbed by Cab (W m-2) of sunlit leaves
rad.Rnh_PAR = Rnhc_PAR;     % [60x1 double]     net PAR absorbed by Cab (W m-2) of shaded leaves
rad.Rnu_PAR = Rnuc_PAR;     % [13x36x60 double] net PAR absorbed (W m-2) of sunlit
rad.Xdd     =   Xdd;

%% APPENDIX I function volscat

function [chi_s,chi_o,frho,ftau]    =   volscat(tts,tto,psi,ttli)

%Volscatt version 2.
%created by W. Verhoef
%edited by Joris Timmermans to matlab nomenclature.
% date: 11 February 2008
%tts    [1]         Sun            zenith angle in degrees
%tto    [1]         Observation    zenith angle in degrees
%psi    [1]         Difference of  azimuth angle between solar and viewing position
%ttli   [ttli]      leaf inclination array

deg2rad = pi/180;
nli     = length(ttli);

psi_rad         = psi*deg2rad*ones(nli,1);

cos_psi         = cos(psi*deg2rad);                 %   cosine of relative azimuth angle

cos_ttli        = cos(ttli*deg2rad);                %   cosine of normal of upperside of leaf
sin_ttli        = sin(ttli*deg2rad);                %   sine   of normal of upperside of leaf

cos_tts         = cos(tts*deg2rad);                 %   cosine of sun zenith angle
sin_tts         = sin(tts*deg2rad);                 %   sine   of sun zenith angle

cos_tto         = cos(tto*deg2rad);                 %   cosine of observer zenith angle
sin_tto         = sin(tto*deg2rad);                 %   sine   of observer zenith angle

Cs              = cos_ttli*cos_tts;                 %   p305{1}
Ss              = sin_ttli*sin_tts;                 %   p305{1}

Co              = cos_ttli*cos_tto;                 %   p305{1}
So              = sin_ttli*sin_tto;                 %   p305{1}

As              = max([Ss,Cs],[],2);
Ao              = max([So,Co],[],2);

bts             = acos(-Cs./As);                    %   p305{1}
bto             = acos(-Co./Ao);                    %   p305{2}

chi_o           = 2/pi*((bto-pi/2).*Co + sin(bto).*So);
chi_s           = 2/pi*((bts-pi/2).*Cs + sin(bts).*Ss);

delta1          = abs(bts-bto);                     %   p308{1}
delta2          = pi-abs(bts + bto - pi);           %   p308{1}

Tot             = psi_rad + delta1 + delta2;        %   pag 130{1}

bt1             = min([psi_rad,delta1],[],2);
bt3             = max([psi_rad,delta2],[],2);
bt2             = Tot - bt1 - bt3;

T1              = 2.*Cs.*Co + Ss.*So.*cos_psi;
T2              = sin(bt2).*(2*As.*Ao + Ss.*So.*cos(bt1).*cos(bt3));

Jmin            = (   bt2).*T1 - T2;
Jplus           = (pi-bt2).*T1 + T2;

frho            =  Jplus/(2*pi^2);
ftau            = -Jmin /(2*pi^2);

% pag.309 wl-> pag 135{1}
frho            = max([zeros(nli,1),frho],[],2);
ftau            = max([zeros(nli,1),ftau],[],2);
return

%% APPENDIX II function e2phot

function molphotons = e2phot(lambda,E,constants)
%molphotons = e2phot(lambda,E) calculates the number of moles of photons
%corresponding to E Joules of energy of wavelength lambda (m)

e           = ephoton(lambda,constants);
photons     = E./e;
molphotons  = photons./constants.A;
return;

function E = ephoton(lambda,constants)
%E = phot2e(lambda) calculates the energy content (J) of 1 photon of
%wavelength lambda (m)

h       = constants.h;           % [J s]         Planck's constant
c       = constants.c;           % [m s-1]       speed of light
E       = h*c./lambda;           % [J]           energy of 1 photon
return;

%% APPENDIX III function Pso

function pso    =   Psofunction(K,k,LAI,q,dso,xl)
if dso~=0
    alf         =   (dso/q) *2/(k+K);
    pso         =   exp((K+k)*LAI*xl + sqrt(K*k)*LAI/(alf  )*(1-exp(xl*(alf  ))));% [nl+1]  factor for correlation of Ps and Po
else
    pso         =   exp((K+k)*LAI*xl - sqrt(K*k)*LAI*xl);% [nl+1]  factor for correlation of Ps and Po
    
end