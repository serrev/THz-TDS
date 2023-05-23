function [m_field] = f_getcomplexfield(lv_n0,lv_n1,lv_n2,lv_n3,lv_n4,c_c,lv_d0,lv_d1,lv_d2,lv_d3,lv_d4,lv_freq)
% 
%With FP
% m_field = f_trans(lv_n0,lv_n1).*f_propagation(lv_n1,lv_freq,c_c,lv_d1).*f_trans(lv_n1,lv_n2).*f_propagation(lv_n2,lv_freq,c_c,lv_d2).*f_FP(lv_n0,lv_n1,lv_n2,lv_freq,c_c,lv_d1).*...
%           f_trans(lv_n2,lv_n3).*f_propagation(lv_n3,lv_freq,c_c,lv_d3).*f_FP(lv_n1,lv_n2,lv_n3,lv_freq,c_c,lv_d2).*f_trans(lv_n3,lv_n4)*f_propagation(lv_n4,lv_freq,c_c,lv_d4).*...
%           f_FP(lv_n2,lv_n3,lv_n4,lv_freq,c_c,lv_d3);

% m_field = f_trans(lv_n1,lv_n2).*f_propagation(lv_n2,lv_freq,c_c,lv_d2).*f_trans(lv_n2,lv_n3).*f_FP(lv_n1,lv_n2,lv_n3,lv_freq,c_c,lv_d2);

% % Without FP
% m_field = f_trans(lv_n0,lv_n1).*f_propagation(lv_n1,lv_freq,c_c,lv_d1).*f_trans(lv_n1,lv_n2).*f_propagation(lv_n2,lv_freq,c_c,lv_d2).*...
%           f_trans(lv_n2,lv_n3).*f_propagation(lv_n3,lv_freq,c_c,lv_d3).*f_trans(lv_n3,lv_n4).*f_propagation(lv_n4,lv_freq,c_c,lv_d4);

% % Without FP AMA
m_field = f_trans(lv_n1,lv_n2).*f_propagation(lv_n2,lv_freq,c_c,lv_d2).*...
          f_trans(lv_n2,lv_n3);

% With FP AMA
% m_field = f_trans(lv_n1,lv_n2).*f_propagation(lv_n2,lv_freq,c_c,lv_d2).*...
%           f_trans(lv_n2,lv_n3).*f_FP(lv_n1,lv_n2,lv_n3,lv_freq,c_c,lv_d2);
