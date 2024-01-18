function a = DropEvol
% Sta³e fundamentalne
global q R e  n c eps0 M Cw Delta_c Delta_t Dv CpAir s   Tinf K Ro   Ro_air alfa_c alfa_t Ma;
% 
R = 8.31447;
e = 1.602e-19;
c = 299792458;
eps0 = 1e7/(4*pi*c^2);
% 
%Sta³e fizyczne
% 
M  = 18.015268e-3;
q  = 226e4;
Cw = 4200;
% 
%Wielkoœci zalezne od temperatury
% 
Tinf = 273.15+15;
P  = 100600;
s  = 0;
Dv =  2.11e-5*( (Tinf/273.15)^1.94 )*( 101325/P );
CpAzot  = 1040;
Mazot   = 28.013e-3;
Ro_azot = @(x)  Mazot *  ( P - s*Psat( x ) ) / ( R * x );
Ro_air  =  Ro_azot( Tinf ) + ( M * s * Psat( Tinf ) ) /( R * Tinf );
% 
% Ciep³o w³asciwe mieszaniny gazów
% 
ffa  = @(z) z / 647.096;
alfa = @(x) ( -1135.905627715 - 5.65134998e-8 * ffa(x)^-19 + 2690.66631 * ffa(x)  + 127.287297*ffa(x)^4.5 -...
             135.003439*ffa(x)^5 + 0.981825814*ffa(x)^54.5)*1000;
    
CpparaWodna = ( alfa( Tinf ) + ( ( R * Tinf^2 )/( M * Psat( Tinf ) ) ) * dPsat( Tinf ) - 2502000 )/(Tinf -273.15);

CpAir = ( CpparaWodna * M * s*Psat( Tinf )+CpAzot*Mazot*(P-s*Psat( Tinf )))/...
        ( Mazot *( P - s * Psat( Tinf ) ) + Mazot * s * Psat( Tinf ) );
%     
% MASA MOLOWA MIESZANINY AZOTU I PARY WODNEJ. GESTOSC WLASCIWA CIEKLEJ WODY
% 
Ma = (M - Mazot)*(s*Psat(Tinf)/P) + Mazot;
Ro = 1e3 * ( 1 - ( ( Tinf + 15.7914 )/( 508929.2*( Tinf - 205.0204 ) ) )*(Tinf-277.1363)^2 ); 
% 
% PRZEWODNOSC CIEPLNA MIESZANINY GAZOW
% 
Kazot = -1.5748e-3 + ( 1.14388e-4*Tinf ) - ( 7.65412e-8 * Tinf^2 ) + 1.18325e-3*( Ro_azot(Tinf)/( Mazot*1e3 ) )+...
         2.90297e-5*( Ro_azot(Tinf)/(Mazot*1e3))^2 + 8.53486e-7*(Ro_azot(Tinf)/(Mazot*1e3))^3 +...
         8.50205e-8*(Ro_azot(Tinf)/(Mazot*1e3))^4;
KparaWodna = -0.35376e-2 + 0.654755e-4*Tinf + 0.17446e-7*Tinf^2;
AzotWoda   = ( sqrt(2)/4 )*sqrt( M/(M + Mazot))*(1 + sqrt( Kazot/KparaWodna )*(M/Mazot)^(1/4) )^2;
WodaWoda   = ( sqrt(2)/4 )*sqrt( Mazot/(M + Mazot))*(1 + sqrt( Kazot\KparaWodna )*(M\Mazot)^(1/4) )^2;

K = Kazot/(1 + AzotWoda*(s*Psat(Tinf)/(P - s*Psat(Tinf)))) + (KparaWodna*s*Psat(Tinf)/...
    (s*Psat(Tinf) + WodaWoda*(P - s*Psat(Tinf))));

Delta_c = 4*( Dv/sqrt( ( 8 * R * Tinf /( pi * M ) ) ) );
Delta_t = Delta_c + 4*( K/(CpAir*Ro_air*sqrt( ( 8 * R * Tinf /( pi * M ) ) ) ) );
% 
% Generacja ewolucji promienia kropli
% 
a = 15950e-9;
dzielnik = 10;
krok = 0.04/dzielnik;
s = 0.94;
alfa_c = 0.13;
alfa_t = 1;
n = 0.45e6;
stop = ( (e*n)^2/(64*pi^2*eps0*Gamma(Tinf) ) )^(1/3);
% 
wb = waitbar(0,'I am thinking...' );
ii = 1;
key = 1;
while (a>=stop)&(ii<7500)
    waitbar(stop/a,wb );

    out.a(ii) = a;
    out.t(ii) = ttP(a);
    out.time(ii) = (ii-1)*krok;
    out.da(ii) = daP(a);
    a = a + out.da(ii)*krok;
    
    ii = ii + 1;
end
close(wb);

figure;
plot(out.time,out.a);
a = out;

% =========================================================================
function temper = ttP ( a )
global q R M s  Tinf ;
ft = @(t) (((q*M*dyf(a,t))/(R*hcond(a,t)) )*( Psat(Tinf)/Tinf )*(s-(Tinf/t)*Kelvin(a,t)*exp( ((q*M)/(R*t*Tinf))*(t - Tinf) ) ))-(t-Tinf);
temper =  fzero(ft,Tinf);
%--------------------------------------------------------------------------  
function Ps = Psat(x)
ff = 1 - (x/647.096);
Ps =  22.064e6 * exp( ( 647.096/x )*( -7.85951783 * ff + 1.84408259 * ff^1.5 - 11.7866497*ff^3 +...
      22.6807411*ff^3.5 -15.9618719*ff^4 + 1.80122502*ff^7.5 ) );
%--------------------------------------------------------------------------                       
function dPs = dPsat( x )
ffe =  1 - (0.0015453657571674064*x);
expon =  exp( (647.096*(-7.85951783*ffe+1.84408259*ffe^1.5 - 11.7866497*ffe^3 +...
            22.6807411*ffe^3.5 - 15.9618719*ffe^4+1.80122502*ffe^7.5))/x );
dPs =  2.2064e7 * expon*( (-1 / x^2 )*(647.096*(-7.85951783*ffe +1.84408259*ffe^1.5 -11.7866497*ffe^3 + 22.6807411*ffe^3.5 - ...
       15.9618719*ffe^4 +1.80122502*ffe^7.5) )+( 1 / x )*(647.096*(0.01214582972232868 - 0.0042746731319618725*ffe^0.5 +0.05464405451432244*ffe^2 -...
       0.12267514225091794*ffe^2.5 +0.09866772101821059*ffe^3 -0.020876636001458827*ffe^6.5 )) );
%--------------------------------------------------------------------------
function Gm = Gamma(x)
ff = 1 - (x/647.096);
Gm =  235.8e-3*ff^1.256 * (1 - 0.625*ff);
%--------------------------------------------------------------------------
function df = dyf(x,y)
global  Dv Delta_c alfa_c M R;
df = Dv/( (x/( x+Delta_c ) )+( Dv/(x*alfa_c) )*sqrt( (2*pi*M)/(R*y) ) );
%--------------------------------------------------------------------------
function hc = hcond(x,y)
global  K  Delta_t alfa_t CpAir R Ro_air  Ma ;
hc =  K/( (x/(x+Delta_t) )+( K/(CpAir*Ro_air*alfa_t*x))*sqrt(2*pi*Ma/(R*y)) );
%--------------------------------------------------------------------------
function Kl = Kelvin (x,y)
global M R Ro e n eps0;
Kl =  exp(M/(R*Ro*y)*(2*Gamma(y)/x - ((e*n)^2)/(32*pi^2*eps0*x^4) ) );
%--------------------------------------------------------------------------
function da = daP(x)
global M R Ro s q  Tinf;
da =  ( M*dyf(x,ttP(x))/(R*x*Ro) )*(Psat(Tinf)/Tinf)*(s-(Tinf/ttP(x))*Kelvin(x,ttP(x))*exp( (q*M)/(Tinf*R*ttP(x))*(ttP(x)-Tinf) ) );
