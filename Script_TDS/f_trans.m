function [lv_t] = f_trans(lv_nj,lv_nk)

lv_t = (2*lv_nj)./(lv_nj+lv_nk);