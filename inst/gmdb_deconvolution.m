function [OutputRecord] = gmdb_deconvolution(Record,OutcropSite,ToplayerSite)
%% ========================================================================
% Copyright SRK/FIUBA (C) 2019
% Coded By: p. Barbieri (pbarbieri@fi.uba.ar)
% Version: V101
% Date: 15-05-2019
% -------------------------------------------------------------------------
% USAGE: 
%   gmdb_deconvolution:  
%       For every record stored as a FourierSeries in a
%       Record struct, the function deconvolves the record from the outcrop
%       to the bedrock throgh the OutcropSite and then convolves the
%       bedrock record to the toplayer through the ToplayerSite
% -------------------------------------------------------------------------
% INPUT:
%       Record       struct with the records stored as FourierSeries
%       OutcropSite  struct with the following fields
%                  - h:     Vector of strata height (top to bottom) [m]
%                  - rho:   Vector of density of each strata (top to
%                           bottom) [kg/m3]
%                  - GSm    Initial shear stiffness of each strata (top to
%                           bottom) [N/m2]
%                  - XCH    Initial damage of each strata (top to bottom)
%                           [ ]
%                  - SH     Initial damping of each strata (top to bottom)
%                           [ ]
%       ToplayerSite  struct with the following fields
%                  - h:     Vector of strata height (top to bottom) [m]
%                  - rho:   Vector of density of each strata (top to
%                           bottom) [kg/m3]
%                  - GSm    Initial shear stiffness of each strata (top to
%                           bottom) [N/m2]
%                  - XCH    Initial damage of each strata (top to bottom)
%                           [ ]
%                  - SH     Initial damping of each strata (top to bottom)
%                           [ ]
%                  - DamageModel Damage model to apply in each strata for a
%                                bedrock to toplayer deconvolution. Models
%                                supported are: 'fixed', 'seed', 'rollings'
%                               ,'vrock', 'vrock', 'sedd & idris sand mean'
% -------------------------------------------------------------------------
% OUTPUT:
%       OutputRecord  struct with the timeseries and FourierSeries of the
%                     processed records
% -------------------------------------------------------------------------
% EXAMPLE:
% % Outcrop site
%   OutcropSite.h = [30;1];
%   OutcropSite.rho = [2400;2400];
%   OutcropSite.GSm = [9.4;9.4];
%   OutcropSite.XCH = [1;1];
%   OutcropSite.SH = [0.03;0.03];
% % Toplayer site of two strata with underlaying rcok
%   ToplayerSite.h = [10,20,1];
%   ToplayerSite.rho = [2100;2300;2400];
%   ToplayerSite.GSm = [9;9.5;12];
%   ToplayerSite.XCH = [1;1;1];
%   ToplayerSite.SH = [0.03;0.03;0.03];
%   ToplayerSite.DamageModel = {'seed';'seed';'fixed'};
% [OutputRecord] = gmdb_deconvolution(Record,OutcropSite,ToplayerSite);
% -------------------------------------------------------------------------
% BIBLIO:
%       
% -------------------------------------------------------------------------
% VALIDATE:
% Version:
% Date:
% Validated by:
% -------------------------------------------------------------------------
% LOG
%   V101    27/05/2019      First version
% =========================================================================


global Options
SetOptions;

BSR = Options.BSR;
DH = Options.DH;

h_o = OutcropSite.h; if size(h_o,2)>size(h_o,1), h_o=h_o.'; end
rho_o = OutcropSite.rho; if size(rho_o,2)>size(rho_o,1), rho_o=rho_o.'; end
GSm_o = OutcropSite.GSm; if size(GSm_o,2)>size(GSm_o,1), GSm_o=GSm_o.'; end
XCH_o = OutcropSite.XCH; if size(XCH_o,2)>size(XCH_o,1), XCH_o=XCH_o.'; end
SH_o = OutcropSite.SH; if size(SH_o,2)>size(SH_o,1), SH_o=SH_o.'; end

h = ToplayerSite.h;     if size(h,2)>size(h,1), h=h.'; end
rho = ToplayerSite.rho; if size(rho,2)>size(rho,1), rho=rho.'; end
GSm = ToplayerSite.GSm; if size(GSm,2)>size(GSm,1), GSm=GSm.'; end
SH = ToplayerSite.SH;   if size(SH,2)>size(SH,1), SH=SH.'; end
DamageModel = ToplayerSite.DamageModel; if size(DamageModel,2)>size(DamageModel,1), DamageModel=DamageModel.'; end

NR = height(Record.Index);
for r = 1:NR
    % Read record
    FourierSeries = Record.FourierSeries{1,r};
    f = seconds(FourierSeries.Time);
    % Outcrop to bedcrok deconvolution
    if ~isempty(find(strcmpi(FourierSeries.Properties.VariableNames,'DWo'),1)) && Options.TryDisp 
        DWo = FourierSeries.DWo;
        % Trim simmetric Fourier transform
        if Options.TrimSymmetricFourierTransform
            DWo = DWo(1:ceil(numel(DWo)/2));
        end
        % set units to meters
        Units = FourierSeries.Properties.VariableUnits{find(strcmpi(FourierSeries.Properties.VariableNames,'DWo'),1)};
        [DWo] = Set_Units(DWo,Units);
        % Deconvolve
        [DWb] = S2B(DWo,f,h_o,rho_o,GSm_o,XCH_o,SH_o);
    else
        AWb = FourierSeries.AWo;
        % Trim simmetric Fourier transform
        if Options.TrimSymmetricFourierTransform
            f = f(1:ceil(numel(AWb)/2));
            AWb = AWb(1:ceil(numel(AWb)/2));
        end
        % set units to meters
        Units = FourierSeries.Properties.VariableUnits{find(strcmpi(FourierSeries.Properties.VariableNames,'AWo'),1)};
        [AWb] = Set_Units(AWb,Units);
        [AWb] = S2B(AWb,f,h_o,rho_o,GSm_o,XCH_o,SH_o);
        DWb = -AWb./(2*pi*f).^2;
    end
    
    % Bedcrock to toplayer deconvolution
    [RUF,~] = SRA(DWb,f,h,rho,GSm,SH,DH,BSR,DamageModel);
    
    DW = RUF(:,1);
    VW = 1i.*DW*2*pi.*f;
    AW = -DW.*(2*pi*f).^2;
    
    [AT,t] = Get_TS(AW,f);
    [VT,~] = Get_TS(VW,f);
    [DT,~] = Get_TS(DW,f);
    
    if Options.PEER_Procesing
        [AT,VT,DT] = PEER_Procesing(AT,t,hpf);
        [AW,~] = Get_FS(AT,t);
        [VW,~] = Get_FS(VT,t);
        [DW,~] = Get_FS(DT,t);
    end
    
    OutputRecord.TimeSeries{1,r} = timetable(AT,VT,DT,'TimeStep',seconds(t(2)-t(1)));
    OutputRecord.TimeSeries{1,r}.Properties.VariableUnits = {'m/s/s','m/s','m'};
    OutputRecord.FourierSeries{1,r} = timetable(AW,VW,DW,'TimeStep',seconds(f(2)-f(1)));
    OutputRecord.FourierSeries{1,r}.Properties.VariableUnits = {'m/s/s','m/s','m'};
    OutputRecord.Intesity = table(); % TO COMPLETE!
end

end

function [] = SetOptions()
global Options 

% Trims forueir transform in half if true
Options.TrimSymmetricFourierTransform = true;

% Reads displacement Fourier Series (if available) for input.
Options.TryDisp = false;

% Max height of layer for shear computation in SRA with damage
Options.DH = 1; % [m]

% Silva's factor
Options.BSR = 0.75; % [ ]

% Highpass Buitterwood half power frequency
Options.hpf = .01; % [1/s]
Options.PEER_Procesing = false;

end

function [BF] = S2B(SF,f,h,rho,GSm,XCH,SH)
% SRA SURFACE TO BEDROCK FOR DISPLACEMENT RECORD(S)
%   SF:  Frequency content of surface aceleration, velocity or displacement record(s) by column.
%   f:   Frequency vector by column.
%   h:   Height of each strata.
%   rho: Density of each strata
%   GSm: Shear modulus of each strata.
%   SH: Hysteretic damping of each strata.
%   XCH: Damage of each strata
%   BF:  Frequency content of bedrock aceleration, velocity or displacement record(s) by column.

[a,b,~] = Get_SRA_ab(f,h,rho,GSm,XCH,SH);
NUP = size(SF,1);
NR = size(SF,2);
NL = size(h,1);
BF = zeros(NUP,NR);
for r  = 1:NR
    BF(:,r) = (a(:,NL)+b(:,NL))/2.*SF(:,r);
end % r
end

function [RUF,h] = SRA(BUF,f,h,rho,GSm,SH,DH,BSR,DamageModel)
% SRA BEDROCK TO TOPLAYER WITH DAMAGE
%   BUF: Frequency content of bedrock displacement record [m]
%   f:   Frequency vector by column [1/s]
%   h:   Height of each strata [m]
%   rho: Density of each strata [kg/m3]
%   GSm: Shear modulus of each strata [N/m2]
%   SH:  Hysteretic damping of each strata [ ]
%   DH:  Max height of layer for shear computation [m]
%   BSR: Silva's factor. usually 75% Reduction A,B [ ]
%   DamageModel: Damage model of each layer [ ]. Options: 'fixed', 'seed',
%                'rollings', 'vrock', 'vrock', 'sedd & idris sand mean'
%   RUF: Frequency content of displacement records for each depth by column.

% Division of layer in sublayers
h(end) = DH;
NL = size(h,1);
NH = [0;cumsum(ceil(h/DH))];
NUP = size(BUF,1);
DamageModel = categorical(DamageModel);
for j = 1:NL
    k = NL-j+1;
    h(NH(k)+1:NH(k+1),1) = h(k)/ceil(h(k)/DH);
    rho(NH(k)+1:NH(k+1),1) = rho(k);
    GSm(NH(k)+1:NH(k+1),1) = GSm(k);
    SH(NH(k)+1:NH(k+1),1) = SH(k);
    DamageModel(NH(k)+1:NH(k+1)) = DamageModel(k);
end
NL = size(h,1);
XCH = ones(NL,1);
% Shear deformation deconvolution of each layer with damage
RF = repmat(BUF,1,NL);
for j = 1:6
    [a,b,D] = Get_SRA_ab(f,h,rho,GSm,XCH,SH);
    K = repmat(2*pi*f,1,NL).*repmat(sqrt(rho./D).',NUP,1);
    GF = sqrt(BSR)*1i*K.*(-a+b)./repmat(a(:,NL)+b(:,NL),1,NL).*RF;
    [GT,~] = Get_TS(GF,f);
    Geff = max(GT,[],1).';
    [SH,XCH] = arrayfun(@(G,xi,model) get_damage(G,xi,model),Geff,SH,DamageModel);
end
RUF = sqrt(BSR)*RF.*(a+b)./repmat(a(:,NL)+b(:,NL),1,NL);
end

function [SH,XCH] = get_damage(Geff,SH,DamageModel)
switch lower(DamageModel)
    case 'fixed'
        XCH = 1;
    case 'seed'
        SH = 0.23*(Geff/(Geff+1.56*10^-3))+0.05;
        XCH = 1.56*10^-3/(Geff+1.56*10^-3);
    case 'rollings'
        SH = 0.008+0.18*(1+0.15*(100*Geff)^-0.9)^-0.75;
        XCH = 1/(1.2+16*100*Geff*(1+10^(-20*100*Geff)));
    case 'vrock'
        SH = (3+(Geff/0.007)/(1/3+(Geff/0.007)/(23-3)))/100;
        XCH = 1/(1+(Geff/0.005)^1.8)^0.35;  
    case 'sedd & idris sand mean'
        Geff = 100*Geff;
        strain = [1.00e-4, 3.16e-4, 1.00e-3, 3.16e-3, 1.00e-2, 3.16e-2, 1.00e-1, 3.16e-1, 1.00e+0];
        damage = [1.000, 0.990, 0.960, 0.880, 0.740, 0.520, 0.290, 0.150, 0.060];
        damping = [0.570, 0.860, 1.700, 3.100, 5.500, 9.500, 15.500, 21.100, 24.600]/100;
        if Geff < strain(1)
            XCH = damage(1);
            SH = damping(1);
        elseif Geff > strain(end)
            XCH = damage(end);
            SH = damping(end);
        else
            XCH = interp1(strain, damage, Geff);
            SH = interp1(strain, damping, Geff);
        end
end
end

function [a,b,D] = Get_SRA_ab(f,h,rho,GSm,XCH,SH)
% SRA TRANSFER FUNCTION FOR EACH STRATA 
%   f:  Frequency vector by column.
%   h:   Height of each strata.
%   rho: Density of each strata
%   GSm: Shear modulus of each strata.
%   SH: Hysteretic damping of each strata.
%   XCH: Damage of each strata
%   a,b: Site filter amplitudes for each frequency and strata.
%   D:   Complex rigidity of each layer

Wo = 2*pi*f;
NUP = size(f,1);
NL = size(h,1);

%D = XCH.*GSm.*(1+2*1i*SH);
D =  XCH.*GSm.*(1-2*SH.^2+2*1i*SH.*sqrt(1-SH.^2));

c = h.*sqrt(rho./D);
a = zeros(NUP,NL); a(:,1) = 1;
b = zeros(NUP,NL); b(:,1) = 1;
for j=1:(NL-1) %#ok<*FXUP>
    alpha = D(j)/D(j+1)*sqrt(rho(j)*D(j+1)/(rho(j+1)*D(j)));
    ap = 1 + alpha;
    an = 1 - alpha;
    a(:,j+1) = 0.5*ap.*exp(1i*c(j)*Wo).*a(:,j) + 0.5*an.*exp(-1i*c(j)*Wo).*b(:,j);% AMPLITUD NORMALIZADA ONDAS ASCENTENDES Y DESCENDENTES
    b(:,j+1) = 0.5*an.*exp(1i*c(j)*Wo).*a(:,j) + 0.5*ap.*exp(-1i*c(j)*Wo).*b(:,j);
end

end

function [AT,VT,UT] = PEER_Procesing(AT,t,hpf)
% PEER PROCESSING OF ACELERATION TIMESERIE(s) TO OBTAIN VELOCITY AND
% DISPLACEMENTS
%   AT: Timeseries of record(s) by column.
%   t:  Time vector by column.
%   VT: Velocity timeseries of record(s) by column.
%   UT: Displacement timeseries of record(s) by column.

AT = AT - mean(AT);
dt = t(2)-t(1);
[HPFilterObj,LPFilterObj] = Build_Filters(dt,hpf);
% Filtered acceleration
[AT] = Filter_AT(AT,t,HPFilterObj,LPFilterObj);
% Baseline correction
[AT,VT,UT] = BaselineCorrection(AT,t);
end

function [HPFilterObj,LPFilterObj] = Build_Filters(dt,hpf)
% BUILDS FILTERS FOR RECORD PROCESSING
%   dt: Time step of the record
%   HPFilterObj: High-pass filter object
%   LPFilter: : Low-pass filter object

HPFilterObj = designfilt('highpassiir', 'FilterOrder', 6, ...
                 'HalfPowerFrequency',hpf, 'SampleRate',1/dt, ...
                 'DesignMethod', 'butter');
if 1/(2*dt)<50
    lpf = 1/(2*dt);
else
    lpf = 50;
end
LPFilterObj =  designfilt('lowpassiir', 'PassbandFrequency', lpf-5, ...
            'StopbandFrequency', lpf, 'PassbandRipple', 1, ...
            'StopbandAttenuation', lpf-10, 'SampleRate', 1/dt, 'MatchExactly', 'passband');
end

function [AT] = Filter_AT(AT,t,HPFilterObj,LPFilterObj)
% FILTER OF ACELERATION TIMESERIE(s)
%   AT: Timeseries of record(s) by column.
%   t:  Time vector by column.
%   HPFilterObj: High-pass filter object
%   LPFilter: : Low-pass filter object

NP =  size(AT,1);
dt = t(2) - t(1);
fhp = HPFilterObj.HalfPowerFrequency;
nfilt = HPFilterObj.FilterOrder;
% Padding 
NPad = ceil(1.5/2*nfilt/fhp*1/dt);
AT(NPad+1:NPad+NP,:) = AT(:,:);
AT(1:NPad,:) = 0;
AT(NPad+NP+1:2*NPad+NP,:) = 0;
% Filtering
for c = 1:size(AT,2)
    AT(:,c) = filtfilt(HPFilterObj,AT(:,c));
    AT(:,c) = filtfilt(LPFilterObj,AT(:,c));
end
% Unpaddig
AT = AT(NPad+1:NPad+NP);
end

function [AT,VT,UT] = BaselineCorrection(AT,t)
% Newmark beta integration
[~,UT] = Get_VUT(AT,t);
% Base-line correction
for r = 1:size(AT,2)
    % 6th order polynom fit to fit displacement with ao=a1=0
    p = PolyFit(t,UT(:,r),[1 1 1 1 0 0]);
    % Acceleration correction
    AT(:,r) = AT(:,r) - polyval(polyder(polyder(p)),t);
end
% Newmark beta integration
[VT,UT] = Get_VUT(AT,t);
end

function [p] = PolyFit(x,y,coefs)
% returns the coefficients for a polynomial p(x) that is a ...
% best fit (in a least-squares sense) for the data in y. coefs is a vector
% of 1 or 0 such that the polynom is defined by:
%
% n = numel(coef)
% p(x) = coef(1)*p(1)*x^(n-1)+coef(2)*p(2)*x^(n-2)+....

NP = numel(x);
N = numel(coefs); % order of the polynom
n = sum(coefs); % number of nonzero coef
M = zeros(n,n);
b = zeros(n,1);
p = zeros(N,1); % coeficients

fil = 0;
for r = 1:N
if coefs(r)==1
    fil = fil+1;
    for j = 1:NP
        b(fil) = b(fil) + y(j) * x(j)^(N-fil);
    end
    col = 0;
    for k = 1:N
    if coefs(k)==1
       col = col+1;
       pwr_fil = x.^(N-k);
       pwr_col = x.^(N-r);
       M(fil,col) = sum(pwr_fil.*pwr_col);
    end
    end
end
end

a = M\b;
k = 0;
for j = 1:N
    if coefs(j)==1
        k = k+1;
        p(j) = a(k);
    end
end
end

function [AF,f] = Get_FS(AT,t)
% FOURIER TRANSFORM OF RECORD(S)
%   AT: Timeseries of record(s) by column.
%   t:  Time vector by column.
%   AF: Frequency content of record(s) by column.
%   f:  Frequency vector by column.

NP = size(AT,1);
dt = t(2)-t(1);
NFFT = pow2(nextpow2(NP)+1);
NUP = NFFT/2+1;
f = 1/(2*dt)*linspace(0,1,NUP).';
AF = fft(AT,NFFT,1);
AF = AF(1:NUP,:);
end

function [AT,t] = Get_TS(AF,f)
% INVERSE FOURIER TRANSFORM OF RECORD(S)
%   AF: Frequency content of record(s) by column.
%   f:  Frequency vector by column.
%   AT: Timeseries of record(s) by column.
%   t:  Time vector by column.

df = f(2)-f(1);
NUP = size(AF,1);
dt = 1/(2*df*(NUP-1));
t = linspace(0,(NUP-1)*dt,NUP).';
AF = cat(1,AF,fliplr(AF));
AT = ifft(AF,'symmetric');
AT = AT(1:NUP,:);
end

function [AT]=Set_Units(AT,units)
% SE UNIT OF RECORD(S) TO METERS
%   AT: Timeseries of record(s) by column.
%   units:  units of the record
switch lower(units)
    case {'mm','mm/s','mm/sec','mm/s2','mm/s/s','mm/s^2','mm/sec/sec'}
        SF=1/1000;
    case {'cm','cm/s','cm/sec','cm/s2','cm/s/s','cm/s^2','cm/sec/sec'}
        SF=1/100;
    case {'m','m/s','m/sec','m/s2','m/s/s','m/s^2','m/sec/sec'}
        SF=1;
    case 'g'
        SF=9.81;
    case 'g/10'
        SF=0.981;
    otherwise
        error('** UNKNOWN RECORD UNITS %s',units);
end
AT = AT*SF;
end
