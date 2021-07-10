%clc; clear; close all;
%Parametros del sistema
%Carga mecanica
Jl = 0.2520; %+/- 0.1260 [kg.m^2]
bl = 0; %+/-0.0630 [N.m/rad/s]
Tl = 0; %+/- 6.28 [N.m] Asumir funcion escalon
%Descomentar para variacion de parametros de carga
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
Lq = 5.8e-3;       %+/- 1% [H]
Ld = 6.6e-3;       %+/- 1% [H]
Lls = 0.8e-3;      %+/- 1% [H]
R_Ref = 1.02;   %+/- 1% [ohm]
Rs_115 = 1.32;  %+/- 1% [ohm]
alfa_cu = 3.9e-3;   %+/- 1% [1/C]
Cts = 0.818;     %+/- 1% [W.s/C]
Rts_amb = 146.7;    %+/- 1% [C/W]
Tau_tsamb = 120;    %+/- 1% [s] 
wm_nom = 691.15;    %[rad/s] velocidad nominal del motor
Vsl_nom = 24;       %[V]
Is_nom = 0.4;       %[A]
Is_max = 2;         %[A]
Ts_max = 115;       %[C]. Rango de operacion -15<=Tam<=40

Beq = bm+(bl)/r^2;
Jeq = Jm+(Jl/r^2);

TREF = 40;
id0 = +0.5; 

%Resistencias para control proporcional del modulador de torque
Rq=29;
Rd=33;
R0=4;

w_pos = 800;
n = 2.5;
Kd = Jeq*n*w_pos;
Kp = Kd*w_pos;
Ki = Kp*w_pos/n;

% Valores de Kw y de Ktheta sin tener en cuenta el integrador en el
% observador
% K_omega_o = 10223018.9;
% K_tita_o =  6397.3456;

% Valores de Ki, Ktita, Komega agregando integrador en observador
K_w_o = 3.2768e10;
K_tita_o = 9600;
K_i_o = 3.072e7;

%Curvas parametricas
% SCOPE_T = ScopeData1toeque_velocidad.signals(2).values;
% figure()
% plot(SCOPE_T,SCOPE_omega)
% grid on
