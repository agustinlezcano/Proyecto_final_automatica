function park_inv = t_park_inv(fq,fd,f0,tita_r)


park_MATRIZ = [cos(tita_r)        sin(tita_r)        1; 
               cos(tita_r-2*pi/3) sin(tita_r-2*pi/3) 1;
               cos(tita_r+2*pi/3) sin(tita_r+2*pi/3) 1];
                 
park_inv = park_MATRIZ * [fq;fd;f0];

end