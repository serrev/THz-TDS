function [lv_FP] = lv_FP(lv_nj,lv_nk,lv_nl,lv_f,lv_c,lv_dk)

lv_FP = 1/(1+lv_reflect(lv_nj,lv_nk)*lv_reflect(lv_nk,lv_nl)*lv_propagation(lv_nk,lv_f,lv_c,lv_dk)^2);