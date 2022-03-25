function [lv_p] = f_propagation(lv_nj,lv_f,lv_c,lv_d)

lv_p = exp(-1i*(lv_f*2*pi/lv_c)*lv_d.*lv_nj);