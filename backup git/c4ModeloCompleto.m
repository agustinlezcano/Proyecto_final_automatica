%Parámetros del sistema
%Carga mecánica
Jl = 0.2520; %+/- 0.1260 [kg.m^2]
bl = 0; %+/-0.0630 [N.m/rad/s]
Tl = 0; %+/- 6.28 [N.m] Asumir funcion escalon

%Tren de transmision
r = 314.3008;   %reduccion total 314.3008:1
wl_nom = 2.2;   %[rad/s] Velocidad nominal (salida)->21 rpm
Tq_nom = 7.26;  %[N.m]  Torque nominal
Tq_max = 29.42; %[N.m]  Torque maximo

%Maquina electrica
Jm = 3.1e-6;    %+/- 1% [kg.m^2]
bm = 1.5e-5;    %+/- 1% [N.m/rad/s]
Pp = 3;         %pares de polos
lambda_m = 0.01546; %+/- 1% [V/rad/s]
Lq = 5.8;       %+/- 1% [mH]
Ld = 6.6;       %+/- 1% [mH]
Lls = 0.8;      %+/- 1% [mH]
Rs_40 = 1.02;   %+/- 1% [ohm]
Rs_115 = 1.32;  %+/- 1% [ohm]
alfa_cu = 3.9e-3;   %+/- 1% [1/°C]
Cts = .818;     %+/- 1% [W.s/°C]
Rts_amb = 146.7;    %+/- 1% [°C/W]
Tau_tsamb = 120;    %+/- 1% [s] 
wm_nom = 691.15;    %[rad/s] velocidad nominal del motor
Vsl_nom = 24;       %[V]
Is_nom = 0.4;       %[A]
Is_max = 2;         %[A]
Ts_max = 115;       %[°C]. Rango de operacion -15<=Tam<=40

