%% Parte 1- 3A Determinacion de polos y ceros del sistema (sistema LTI equivalente
% aumentado)
Rs = R_Ref;
%Funcion de transferencia 1
NUM1 = [0 0 0 (3/2)*Pp*lambda_m];
DEN1 = [(Jeq*Lq) (Jeq*Rs+Beq*Lq) (Beq*Rs+(3/2)*Pp^2*lambda_m^2) 0];

figure()
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

syms b_eq j_eq P_p lam_m R_s R_r L_q

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
%% Parte 2- Ejercicio 2- Polos del controlador
format longg
Jl = [0.2520 (0.2520+0.1260) (0.2520-0.1260)];
bl = [0 0.0630 -0.0630];
num_2tita = zeros(3,4);
den_2=zeros(3,4);
num_2TL=zeros(3,2);

Beq = bm+(bl)/r^2;
Jeq = Jm+(Jl/r^2);

w_pos = 800;
n = 2.5;
Kd = Jeq*n*w_pos;
Kp = Kd*w_pos;
Ki = Kp*w_pos/n;


for i=1:3
    num_2tita(i,:) = [0 Kd(i) Kp(i) Ki(i)];
    den_2(i,:) = [Jeq(i) Kd(i) Kp(i) Ki(i)];
    num_2TL(i,:) = [1 0];
end

[Z1,P1,K1] = tf2zp(num_2tita(1,:),den_2(1,:))
H2a=tf(num_2tita(1,:), den_2(1,:));

[Z2,P2,K2] = tf2zp(num_2tita(2,:),den_2(2,:))
H2b=tf(num_2tita(2,:), den_2(2,:));

[Z3,P3,K3] = tf2zp(num_2tita(3,:),den_2(3,:))
H2c=tf(num_2tita(3,:), den_2(3,:));

figure()
pzmap(H2a,H2c)
grid on
% num_2tita = [0 Kd Kp Ki];
% den_2 = [Jeq Kd Kp Ki];
% num_2TL = [1 0];
% [Z1,P1,K1] = tf2zp(num_2tita,den_2)
% H2a=tf(num_2tita, den_2);
% 
% num_2re_corriente = [0 1];
% den_2re_corriente = [1 5000];
% [Z2,P2,K2] = tf2zp(num_2re_corriente,den_2re_corriente)
% H2b=tf(num_2re_corriente, den_2re_corriente);
% 
% figure()
% pzmap(H2a,H2b,H2)
% grid on
% % hold on
% % pzmap(H2b,'LineWidth',2)
% % grid on

