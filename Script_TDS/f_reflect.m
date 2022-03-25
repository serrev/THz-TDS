function [lv_r] = f_reflect(lv_nj,lv_nk)

lv_r = (lv_nj-lv_nk)./(lv_nj+lv_nk);