function la_hampulse = f_hamming(lv_timestep,la_pulse)

[lv_maxpeakref,lv_indpeakref] = max(abs(la_pulse));
lv_NHam = 4e-12/lv_timestep;
lv_Ham = hamming(lv_NHam);
la_Window = [(zeros(lv_indpeakref - ceil(lv_NHam/2),1))',lv_Ham', (zeros(size(la_pulse,1)-ceil(lv_NHam/2)-lv_indpeakref,1))'];
la_hampulse = la_pulse.*la_Window';