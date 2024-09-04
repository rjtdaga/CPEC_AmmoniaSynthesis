As compared to version Version_Dybeck_Orig_Mine, this code will have the following additional functions:
1) Changed the Dybeck approach by decoupling isQE calculation from the change of superbasin due to execution of nonQE reaction step.
2) Added feature of varying parameters with charge, charge function currently defined by square function waves.
3) Added an analysis input to analyse a particular quantity and vary it over different runs, running on different nodes using parallel implementation. 
The period can be provided with inputs being x, x+d,.... x is the input value and d is the difference.
4) The diffusion probability is changed from exp(dE/RT) to exp(dE/RT)/(exp(dE/RT)+1)
5) H adsorption desorption step is assumed to be equilibrated
6) The effective bond order calculation method is changed to the one using UBI-QEP theory instead of simple averaging.
