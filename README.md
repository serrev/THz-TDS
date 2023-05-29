# THz-TDS
THz Time-Domain Software for retrieve the refractive index, permittivity and conductivity

%% Important things  %%

Developed with love by SerRev 

In the script's name ->

A -> Air

S -> Substrate

M -> Sample

S_TDS_AAMAA -> Script for resolve numerically the complex refractive index in a scheeme of having a free-standing sample (3 layer system)

S_TDS_ASASA -> Script for resolve numerically the complex refractive index in a scheeme of having Air+Substrate+Air+Substrate+Air (5 layer system)

Inside the script (good practices all of us should follow to write codes that can be understand by anyone at any time) ->

a_* -> array

la_* -> local array

m_* -> matrix

v_* -> variable

lv_* -> local variable

f_File* -> file

lv_File* -> local file

f_* -> MATLAB function

lv_f -> Figure number

---------------

To measure the TDS of a free-standing sample, use code S_TDS_AAMAA

To measure the TDS of a sample inside two substrates, use firstly S_TDS_ASASA, then compute the script S_TDS_ASMSA loading the solutions of the substrate obained before

--------------

For any information regarding the code, free email can be send to ser.revuelta@gmail.com
