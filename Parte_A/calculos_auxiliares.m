%% Determinacion de polos y ceros del sistema (sistema LTI equivalente
% aumentado)
%Funcion de transferencia 1
Rs = R_Ref;
NUM1 = [0 0 0 (3/2)*Pp*lambda_m];
DEN1 = [(Jeq*Lq) (Jeq*Rs+Beq*Lq) (Beq*Rs+(3/2)*Pp^2*lambda_m^2) 0];

figure(1)
[Z1,P1,K1] = tf2zp(NUM1,DEN1)
subplot(1,2,1)
H1=tf(NUM1, DEN1);
pzmap(H1)
grid on
title('Polos y ceros para Gvqs')
%funcion de transferencia 2

NUM2 = [0 0 Lq Rs];
DEN2 = r*DEN1;

[Z2,P2,K2] = tf2zp(NUM2,DEN2)

H2=tf(NUM2, DEN2);
subplot(1,2,2)
pzmap(H2)
grid on
title('Polos y ceros para GTl')

%La figura de abajo es exclusivamente para el informe
figure(2)
pzmap(H2)
grid on
title('Polos y ceros del sistema')

%La unica diferencia entre las dos funciones de transferencia es que GTl
%aporta un cero en -Rs/Lq. 
%La funcion de transferencia aporta un polo en el origen que equivale a
%un integrador puro. Los integradores aportan inestabilidad 
%El grafico representativo es el de GTl

[A,B,C,D] = tf2ss(NUM2,DEN2)

%MODOS DE VIBRAR: Autovalores de la matriz  A

e = eig(A)  %Se verifica que los autovalores de A son los polos de la funci�n de transferencia. Los autovalores son los modos de oscilacion
    
%Frecuencia Natural
Wn = sqrt((1/(Lq*Jeq)*(Beq*Rs+(3/2)*Pp^2*lambda_m^2)))   %Frecuencia natural del sistema

%Amortiguamiento
Zitta = (1/(Lq*Jeq)*(Jeq*Rs+Beq*Lq)/(2*Wn)) 

%De lo anterior se puede concluir que el sistema esta subamortiguado y
%entonces es estable.

%VER DE AGREGAR UN ANALISIS DE TEMPERATURA (CONSULTAR)

%An�lisis de settling time y otros par�metros
% figure(3)
% step(H2)
% S = stepinfo(H2) %Error: ver step

%A medida que aumenta la parte imaginaria->la frecuencia natural de la respuesta aumenta

%% Matrices de Obsevabilidad y Controlabilidad
Matriz_A=[0 1 0; 
          0 (-Beq/Jeq) ((3*Pp*lambda_m)/(2*Jeq));
          0 ((Pp*lambda_m)/(Jeq)) -Rs/Lq];
Matriz_B=[0 0;
          (-1)/(r*Jeq) 0;
          0 1/Lq];
Matriz_C=[1 0 0];

Controlabilidad = [Matriz_B Matriz_A*Matriz_B (Matriz_A*Matriz_A)*(Matriz_B)];
fprintf('\nrango de Matriz de controlabilidad\n')
rank(Controlabilidad)
Observabilidad = [Matriz_C;
                  Matriz_C*Matriz_A;
                  Matriz_C*(Matriz_A*Matriz_A)];
fprintf('\nrango de Matriz de observabilidad\n')
rank(Observabilidad)

syms b_eq j_eq P_p lam_m R_s R_r L_q L_d

M_A=[0 1 0; 
          0 (-b_eq/j_eq) ((3*P_p*lam_m)/(2*j_eq));
          0 ((P_p*lam_m)/(j_eq)) -R_s/L_q];
%Duda: usar solo variable controlada para analizar controlabilidad
% M_B=[0 0;
%           (-1)/(R_r*j_eq) 0;
%           0 1/L_q];
M_B=[0;
     0;
     1/L_q];

      
%M_C=[1 0 0];
M_C=[0 1 0];    %Para salida omega

M_cont = simplify([M_B M_A*M_B (M_A*M_A)*(M_B)])
M_obs = simplify([M_C;
         M_C*M_A;
         M_C*(M_A*M_A)])
M_A_2=[0 1 0 0; 
          0 (-b_eq/j_eq) ((3*P_p*lam_m)/(2*j_eq)) 0;
          0 0 0 (-R_s/L_d);
          0 ((P_p*lam_m)/(j_eq)) -R_s/L_q 0];
M_B_2=[0;
         0;
     1/L_q;
     0];
 M_C_2=[1 0 0 1];
fprintf('Matriz de observabilidad para LTI aumentado')
M_obs_aumentado = simplify([M_C_2; M_C_2*M_A_2; M_C_2*M_A_2^2; M_C_2*M_A_2^3])
if (rank(M_obs_aumentado)<4)
    fprintf('Sistema LTI aumentado NO completamente observable')
    rank(M_obs_aumentado)
else
    fprintf('Sistema LTI aumentado completamente observable')
    rank(M_obs_aumentado)
end


fprintf('\nMatriz de controlabilidad para LTI aumentado')

M_cont_2 = simplify([M_B_2 M_A_2*M_B_2 (M_A_2*M_A_2)*(M_B_2)])
 
if rank(M_cont_2)~=4
    fprintf('Matriz LTI Aumentada no controlable completamente\n')
else
    fprintf('Sistema completamente controlable para Matriz LTI Aumentada \n')
end

M_B_3 = [0 0;
         0 0;
         1/L_q 0;
         0 1/L_d];
     
M_cont_3 = simplify([M_B_3 M_A_2*M_B_3 (M_A_2*M_A_2)*(M_B_3) (M_A_2*M_A_2*M_A_2)*(M_B_3)])

if rank(M_cont_3)~=4
    fprintf('Matriz LTI Aumentada no controlable completamente\n')
else
    fprintf('Sistema completamente controlable para Matriz LTI Aumentada \n')
end

%% Parametros: settling time, rise time, etc
t1 = 0.1:.2:.9;
%Valores de Iq
iq1 = [0 0.0907 0.376 -0.1950 -0.2858];   %Para perturbacion de torque
iq2 = [0.0907 0.376 -0.1950 -0.2858 0];	%Para perturbacion de torque
iq10 = zeros(1,5);
iq90 = zeros(1,5);
iq_p = zeros(1,5);
iq_n = zeros(1,5);
DT_iq = [101.725e-3 45.80e-3  45.77e-3]    %Tiempo de establecimiento iq(+-1%)
sobrepico_iq = [7.404];
%Valores de w
w1 = [0 420 414.2 426.8 6.28];
w2 = [420 414.2 426.8 6.28 0];
w10 = zeros(1,5);
w90 = zeros(1,5);

w_p = zeros(1,5);
w_n = zeros(1,5);

for i=1:5
    iq_p(i) = iq2(i)+ abs(iq2(i) - iq1(i)) * 0.01;
    iq_n(i) = iq2(i)- abs(iq2(i) - iq1(i)) * 0.01;
    iq10(i) = (iq2(i) - iq1(i)) * 0.1 + iq1(i);
    iq90(i) = (iq2(i) - iq1(i)) * 0.9 + iq1(i);

    w_p(i) = w2(i)+ abs(w2(i) - w1(i)) * 0.01;
    w_n(i) = w2(i)- abs(w2(i) - w1(i)) * 0.01;
    w10(i) = (w2(i) - w1(i)) * 0.1 + w1(i);
    w90(i) = (w2(i) - w1(i)) * 0.9 + w1(i);
end
iq_p
iq_n
iq10
iq90
w_p
w_n
w10
w90

%% Cálculo de Kio para observador

syms Ki_o K_o Kw_o s

M_A_o = [0 1 0;
        0 0 1;
        0 0 0];
    
M_C_o = [1 0 0];

K = [K_o;  Kw_o; Ki_o;];

EC_obs = collect((det(s*eye(3) - (M_A_o-K*M_C_o))),s)
expand((s+3200)^3)


% vars = [Ki_o K_o Kw_o];
% eq1=solve(eqns, vars)
% eq1.double
%s^3+3*s^2*3200+3*s*3200^2+3200^3