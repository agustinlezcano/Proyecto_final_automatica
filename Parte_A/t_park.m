function park_dir = t_park(fa,fb,fc,tita_r)

park_MATRIZ = [cos(tita_r) cos(tita_r-2*pi/3) cos(tita_r+2*pi/3); 
               sin(tita_r) sin(tita_r-2*pi/3) sin(tita_r+2*pi/3);
                  0.5           0.5            0.5];
                 
park_dir = (2/3) * park_MATRIZ * [fa;fb;fc];

end