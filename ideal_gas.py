from __future__ import print_function
import espressomd as es
from espressomd import analyze
from espressomd import grand_canonical
import numpy as np
import func as func
import sys, os, math
#-------------------------------------------------------------------------#
int_steps=10;
cycles=10000;
box_l = 10
volume = np.power(box_l, 3.)
n_H= 100


# Simulation Parameters..
print ("features:",es.features())
types = {'H':0}
charges = {'H': 0}
kT=1.0; gamma=1.0;
time_step=0.01; skin=0.4

#----------------------------------------#
system=es.System()
system.box_l = [box_l, box_l, box_l]
system.periodicity = [1, 1, 1]
system.thermostat.set_langevin(kT=kT, gamma=gamma)
system.time_step=time_step
system.cell_system.skin = skin

func.add_H(system, n_H=n_H, box_l= box_l, H_type=types['H'],q_H=charges['H'], verbose=0)
func.write_xyz(system, filename="initial.xyz", mode="w", pids="all")
r,rdf_sA =system.analysis.rdf(rdf_type='rdf',
        type_list_a=[types["H"]],
        type_list_b=[types["H"]],
        r_min=0, r_max=box_l/2,  r_bins=5)
r_bins = 50
avg_gR = 0
rdf_avg = 0
avg_num_density = 0
for i in range(0,cycles):
    system.integrator.run(int_steps);
    r,rdf=system.analysis.rdf(rdf_type='rdf',
            type_list_a=[types["H"]],
            type_list_b=[types["H"]],
            r_min=0, r_max=box_l/2,  r_bins=50)
    gR, slice, rho =func.calculate_rdf(system,r_bins=50, r_min=0, r_max=box_l/2,  box_l=box_l, type_a=types["H"], type_b=types["H"])
    gr_value = np.array(gR)
    rdf_value = np.array(rdf)
    num_density = gr_value*rho
    avg_gR = avg_gR + gr_value
    rdf_avg =  rdf_avg + rdf_value
    avg_num_density = avg_num_density + num_density
    #print (avg_num_density, rdf_avg, avg_gR)
avg_my_rdf = np.array(avg_gR)/cycles
avg_rdf_es = np.array(rdf_avg)/cycles
avg_nD = np.array(avg_num_density)/cycles
save_data = open("rdf.dat", "w")
save_data.write("r      es_rdf      my_rdf         density\n")
length = len(slice)
for i in range (0, length):
     save_data.write("{}        {}      {}      {}\n".format(slice[i], avg_rdf_es[i], avg_my_rdf[i], avg_nD[i]))
