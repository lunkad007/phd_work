from __future__ import print_function
import numpy as np
import random, math
from random import randrange, uniform
pi = math.pi
sin = math.sin
cos = math.cos
import sys
from math import sqrt, ceil, floor, pow, pi
def calc_sphere_coord(radius, theta, phi):
    theta = theta*pi/180
    phi= phi*pi/180
    x= radius * sin(phi) * cos(theta)
    y= radius * sin(phi) * sin(theta)
    z= radius * cos(phi)
    return(x,y,z)

def create_protein(system,n_pA, n_pB, pA_type, pB_type, q_pA, q_pB, box_l, n_prot, verbose):
    Direction = (1,0,0)
    for i in range(0,n_prot):
        pos=np.random.randint(box_l*0.8, size=3)
        for j in range(0, n_pA):
            if(verbose):
                print ("pid: %d"%pid);
            pid = len(system.part)
            if j==0:
                system.part.add(id=pid, pos=pos, type=pA_type, q=q_pA)
            else:
                bond_partner= pid-1
                system.part.add(id=pid, pos=pos, type=pA_type, q=q_pA)
                system.part[pid].add_bond((system.bonded_inter[0],bond_partner));
            pid+=1
            pos+= Direction
        pid= len(system.part)
        for k in range (0,n_pB):
            system.part.add(id=pid, pos=pos, type=pB_type,q=q_pB)
            bond_partner= pid-1
            system.part[pid].add_bond((system.bonded_inter[0], bond_partner))
            pid+=1
            pos+=Direction


def create_star_pe(system, n_stars, n_arms, mpc, q_sC, q_A, sC_type, A_type, box_l, center_pos, verbose):
    for i in range(0, n_stars):
        if n_stars ==1:
            centeral_pos=center_pos*system.box_l
        else:
            centeral_pos=random.random()*system.box_l
        pid = len(system.part)
        centeral_id = pid
        system.part.add(id=pid, pos=centeral_pos, type=sC_type, q=q_sC)
        for arm in range(0, n_arms):
            pos=centeral_pos.copy()
            radius=1
            theta= 360*arm/n_arms
            phi= 180*arm/n_arms
            for i in range(0, mpc):
                if(verbose):
                    print ("pid: %d"%pid);

                pid+=1
                pos+=calc_sphere_coord(radius, theta, phi)
                if (i==0):
                    bond_partner=centeral_id
                else:
                    bond_partner = pid-1
                system.part.add(id=pid,pos=pos, type= A_type, q= q_A)
                system.part[pid].add_bond((system.bonded_inter[0], bond_partner))

def add_H(system, n_H, box_l, H_type, q_H, verbose):
    pid=len(system.part)
    for i in range(0, n_H):
        pos=box_l**np.random.random(3)
        if(verbose):
            print ("pid: %d"%pid);
        system.part.add(id=pid, pos=pos, type = H_type, q= q_H)
        pid+=1

def write_xyz(system, filename, mode, pids):
    part=system.part
    if(pids=="all"):
        pids=[];
        for p in part:
            pids.append(p.id);
    npart = len(pids)
    f=open(filename, mode)
    f.write("%d\n#box"%(npart)+str(system.box_l)+"\n");
    for p in part:
        pos=p.pos
        f.write("C{0:.0f} {1:.6g} {2:.6g} {3:.6g}\n".format(p.type,  pos[0],pos[1], pos[2]))


def add_wca_all(system,types,elj,slj,cutlj,shift="auto"):
    for t1 in types:
        for t2 in types:
            i=types[t1]; j=types[t2];
            if (i<=j):
                system.non_bonded_inter[i,j].lennard_jones.                set_params(epsilon=elj, sigma=slj, cutoff=cutlj, shift="auto");
                print ("Add WCA for pair %d %d "%(i,j),t1,t2);

def run_warm(system,warm_step, warm_n_time, lj_cap, min_dist):
    print("Start warm-up integration (capped LJ-interactions)")
    system.time = 0
    print("""
              Start warmup integration:
               At maximum {} times {} steps
                Stop if minimal distance is larger than {}
                 """.strip().format(warm_n_time, warm_step, min_dist))
    i=0
    tmp_cap = lj_cap
    system.non_bonded_inter.set_force_cap(tmp_cap)
    act_min_dist = system.analysis.mindist()
    while i < warm_n_time and act_min_dist < min_dist :
        system.integrator.run(warm_step)
        act_min_dist = system.analysis.mindist()
        energies = system.analysis.energy()
        print("run {} at time = {} (LJ cap= {} ) min dist = {}  E={}".strip().format(i, system.time, lj_cap, act_min_dist, energies['total']))
        i+=1
        system.non_bonded_inter.set_force_cap(lj_cap)
        lj_cap += 10
    system.non_bonded_inter.set_force_cap(0)
    print("\nWarm up finished\n")




def distance(system, n_prot,  cs_id, prot, p_start_id):
    mid_prot =(p_start_id +prot)/2
    t_prot_particle= n_prot*(prot)
    part=system.part
    Cs_pos=part[cs_id].pos_folded
    prot_pos=[]
    for i in range(mid_prot,t_prot_particle,prot):
        prot_pos.append(part[i].pos_folded)
    each_pdist=[]
    for p in prot_pos:
        x_axis = Cs_pos[0]-p[0]
        y_axis = Cs_pos[1]-p[1]
        z_axis = Cs_pos[2]-p[2]
        dist = (x_axis**2+ y_axis**2+ z_axis**2)**0.5
        each_pdist.append(dist)
    return (each_pdist)


def calculate_rdf(system,r_bins, r_min, r_max, box_l, type_a, type_b):
    """
    r_bins -- Number of bins
    r_max -- Maxium length
    r_min -- Minimum length
    box_l -- system box length
    """
    #print (r_bins, r_min, r_max)
    part= system.part
    pos_a=[]
    pos_b=[]
    xx=[]; yy=[]; zz=[]
    hist=np.zeros(r_bins+1)
    gR = np.zeros(r_bins)
    bin_width = (r_max - r_min)/float(r_bins)
    cnt = 0
    xx=[];yy=[];zz=[]
    for i in part:
        if i.type==type_a:
            pos_a.append(i.pos_folded)
        if i.type==type_b:
            pos_b.append(i.pos_folded)
    length_a = len(pos_a)
    length_b = len(pos_b)
    if type_a == type_b:
        npart=length_a
    else:
        npart = length_b
    for j in range (0, length_a):
        xj=pos_a[j][0]
        yj=pos_a[j][1]
        zj=pos_a[j][2]
        if type_a == type_b:
            for k in range (j+1, length_a):
                xx.append(pos_a[k][0]-xj)
                yy.append(pos_a[k][1]-yj)
                zz.append(pos_a[k][2]-zj)
                cnt+=1
        else:
            for k in range (0, length_b):
                xx.append(pos_a[k][0]-xj)
                yy.append(pos_a[k][1]-yj)
                zz.append(pos_a[k][2]-zj)
                cnt+=1
    length_box=len(xx)
    for i in range(0, length_box):
        #minimum image
            if (xx[i] < -r_max):   xx[i] = xx[i] + box_l
            if (xx[i] > r_max):   xx[i] = xx[i] - box_l
            if (yy[i] < -r_max):   yy[i] = yy[i] + box_l
            if (yy[i] > r_max):   yy[i] = yy[i] - box_l
            if (zz[i] < -r_max):   zz[i] = zz[i] + box_l
            if (zz[i] > r_max):   zz[i] = zz[i] - box_l
            #distance between i and j
            rd = xx[i] * xx[i] + yy[i] * yy[i] + zz[i] * zz[i]
            rij = sqrt(rd)
            #determine in which bin the distance falls
            bin = int(ceil(rij/bin_width))
            if (bin <= r_bins):
                hist[bin] +=1
    volume = pow(box_l, 3.)
    rho = npart/volume
    norm = 2.0*math.pi*rho*npart*bin_width
    gR=[]; slice=[]
    density = []
    for i in range(1, r_bins+1):
        r = (i - 0.5) * bin_width
        i#bin_vol = (4.0/3.0) * math.pi * (pow(r+(bin_width/2.0),3.0) - pow(r-(bin_width/2.0), 3.0))
        val = hist[i]/norm/((r*r)+(bin_width*bin_width)/12)
        gR.append(val)
        slice.append(r)
    return gR, slice, rho
