%Codigo para comparar NL completo desacoplado vs LTI eq aumentado
%Parte 1- ejercicio 6
%% Inciso A
close all; 
LTI_tita = LTI_tita_w.signals(1).values;
LTI_omega = LTI_tita_w.signals(2).values;
LTI_iq = LTI_i.signals(1).values; 
LTI_id = LTI_i.signals(2).values;  
time_LTI = LTI_tita_w.time;

NL_tita = NL_tita_w.signals(1).values;
NL_omega = NL_tita_w.signals(2).values;
NL_iq = NL_i.signals(1).values;
NL_id = NL_i.signals(2).values;
NL_i0 = NL_i.signals(3).values;
time_NL = NL_tita_w.time;

%tita y omega
figure();
%subplot(2,1,1);
mi_plot(time_LTI, LTI_tita, time_NL, NL_tita, 'tita\_LTI','tita\_NL');
%subplot(2,1,2);
mi_plot(time_LTI, LTI_omega, time_NL, NL_omega, 'w\_LTI','w\_NL');

%corrientes qd0
figure()
subplot(2,1,1);
mi_plot(time_LTI, LTI_iq, time_NL, NL_iq, 'iq\_LTI', 'iq\_NL');
subplot(2,1,2);
mi_plot(time_LTI, LTI_id, time_NL, NL_id, 'id\_LTI', 'id\_NL');

function mi_plot(x1,y1,x2,y2,text1,text2)
    plot(x1, y1, 'r');
    hold on
    plot(x2, y2, 'b');
    grid on
    legend(text1,text2);
end
