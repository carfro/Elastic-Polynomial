#!/usr/bin/python
# -*- coding: utf-8 -*-
""" This script calculates the bulk modulus of a solid from Quantum Espresso outputs.
It is called with 8 arguments "python modifiedApprEos_V1.py 'output_dir' 'output_name' 'atom_ref' 'sym_grp' 'show_plot' 'verbosity' 'tex_table_out' 'plot_out' "

    Parameters
    ----------
    :arg output_dir:          TODO 
    :arg output_name:       TODO
    :arg atom_ref:          TODO
    :arg sym_grp:           TODO
    :arg show_plot:         TODO
    :arg verbosity:         TODO
    :arg tex_table_out:     TODO
    :arg plot_out:          TODO

    Returns
    -------
    :returns: 
"""
#______________________________IMPORTS_______________________________________________

import os, getopt, sys, math, numpy as np, glob, argparse, re

import matplotlib.pyplot as plt
from data_parsing import read_input
from fit_and_plot import pdf_plot_oneD,data_plot_oneD,data_plot_twoD,bmfit_1D,one_var_fit,two_var_fit,E_der_2D,pdf_violin_plot_oneD,pdf_bar_plot_oneD,pdf_bar_plot_oneD_OLD,pdf_plot_twoD
from table_generator import generate_tex_table,generate_tex_table_light,print_ascii_output_table,generate_tex_table_light_transposed,exp_val
from numpy.linalg import inv 

#______________________________OLD IMPORTS_______________________________________________
#from   scipy import *
#from   scipy.optimize import leastsq,curve_fit
#from   lmfit import Model

#from matplotlib.widgets import LassoSelector
#from matplotlib.path import Path


#from   prettytable import PrettyTable

#import matplotlib
#from   matplotlib import cm
#from   mpl_toolkits.mplot3d import Axes3D

#matplotlib.use('TkAgg')
#matplotlib.rcParams['text.usetex'] = True

#from   pylatex import LongTable, MultiColumn


#from elastic_properties import *
#from plot_generator import data_plot_twoD
#from fitting import bmfit_1D,one_var_fit,two_var_fit

#from   matplotlib import rc

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

# ------------------  ALL THINGS BULK MODULUS  --------------------


def inv_nbr_atms(sym):
    if (sym == "bcc") or (sym == "hcp") or (sym=="halite"):
        vrat =1./2.
    elif (sym == "fcc"):
        vrat =1./4.
    elif sym == "dia":
        vrat =1./4.
    else:
        vrat=1
        print("Vrat =",vrat)
    return vrat

def nbr_atms_comp(sym):
    if (sym == "bcc") or (sym == "hcp") or (sym=="fcc"):
        vrat =1./1.
    elif (sym == "dia") or (sym=="halite"):
        vrat =1./2.
    else:
        vrat=1
        print("Vrat =",vrat)
    return vrat


def one_dim_bulkmod(name,functional,param_list,v_list,e_list,eRy_list,geo,vrat,nat,nat_type):
    """TODO: Docstring for one_dim_bulkmod.

    :param_list: TODO
    :v_list: TODO
    :e_list: TODO
    :returns: TODO

    """
    B0=0
    BP=0
    print("NAT",nat)
    print("vrat",vrat)
    #if nat_type==2:
    #    nat=vrat#*nat#*nat_type
    
    if geo!="3L-c": 
        oneD_poly = one_var_fit(param_list,e_list)
        BM_poly=bmfit_1D(v_list,e_list)
        mInv = np.array([[1/oneD_poly.params['m1'].value]]) 
        gam1 = oneD_poly.params['gam1'].value

        [V0,g,h] = ana_expr(oneD_poly.params['a0'].value,geo,vrat)

        B0 = V0 / np.matmul( g,np.matmul(mInv,np.transpose(g)) )#*1.60219*10**(-19 + 30 - 9)

        Bp = V0 * (3 * np.matmul(g,np.matmul(mInv ,np.matmul( h,np.matmul(mInv,np.transpose(g))))) - gam1*np.matmul(mInv,np.transpose(g))**3 )\
                / np.matmul( g,np.matmul(mInv,np.transpose(g)) )**2 - 1

        res_poly = [oneD_poly.params['a0'].value, -1*oneD_poly.params['k'].value, B0[0][0]*1.60219*10**(-19 + 30 - 9)*nat, Bp[0][0]]
        #res_poly = [oneD_poly.params['a0'].value, -1*oneD_poly.params['k'].value/(nat*nat_type), B0[0][0]*1.60219*10**(-19 + 30 - 9)/nat, Bp[0][0]]
        #res_poly = [oneD_poly.params['a0'].value, -1*oneD_poly.params['k'].value/nat, B0[0][0]*1.60219*10**(-19 + 30 - 9), Bp[0][0]]
    else:
        e_min=min(e_list)
        e_list=[ (e-e_min) for e in e_list ]    
        oneD_poly = one_var_fit(param_list,e_list)

        res_poly = [oneD_poly.params['a0'].value, -1*oneD_poly.params['k'].value/nat, 0., 0.]

    result_dict={"name":name,"functional":functional,"BM_poly":BM_poly,"poly":oneD_poly,"res_poly":res_poly,"e_list":e_list,"param_list":param_list,"eRy_list":eRy_list,"v_list":v_list,"geo":geo,"vrat":vrat,"nat":nat}
    return result_dict

def two_dim_bulkmod(name,functional,param_list,v_list,e_list,eRy_list,geo,vrat,nat):
    """TODO: Docstring for two_dim_bulkmod.

    :param_list: TODO
    :v_list: TODO
    :e_list: TODO
    :returns: TODO

    """
    BM_poly=[[0],[1]]

    if sym!="3L-ab":
        twoD_poly = two_var_fit(param_list,e_list)
        #twoD_poly = two_var_Refit(twoD_poly,param_list,e_list)

        # TODO: BM FIT NOT YET IMPLEMENTED, this is simply placeholder
        #BM_poly=bmfit_1D(v_list,e_list)
        mInv = inv(np.array([[twoD_poly.params['m11'].value, twoD_poly.params['m12'].value], [twoD_poly.params['m12'].value, twoD_poly.params['m22'].value]]))

        [V0,g,h] = ana_expr([twoD_poly.params['a0'].value, twoD_poly.params['c0'].value],geo,vrat)

        B0 = V0 / np.matmul( g,np.matmul(mInv,np.transpose(g)) )

        F = np.matmul(mInv,np.transpose(g))
        gam111 = twoD_poly.params['gam111'].value
        gam222 = twoD_poly.params['gam222'].value
        gam112 = twoD_poly.params['gam112'].value
        gam122 = twoD_poly.params['gam122'].value

        B = gam111*F[0]**3 + gam222*F[1]**3 + 3*gam112*F[0]**2*F[1]**1 + 3*gam122*F[0]**1*F[1]**2

        Bp = V0 * (3 * np.matmul(g,np.matmul(mInv ,np.matmul( h,np.matmul(mInv,np.transpose(g))))) - B )\
                / np.matmul( g,np.matmul(mInv,np.transpose(g)) )**2 - 1

        res_poly = [twoD_poly.params['a0'].value, twoD_poly.params['c0'].value, -1*twoD_poly.params['k'].value/nat, B0[0][0]*1.60219*10**(-19 + 30 - 9)/nat, Bp[0][0]]
    else:
        e_min=min(e_list)
        e_list=[ (e-e_min) for e in e_list ]    
        twoD_poly = two_var_fit(param_list,e_list)

        #mInv = inv(np.array([[twoD_poly.params['m11'].value, twoD_poly.params['m12'].value], [twoD_poly.params['m12'].value, twoD_poly.params['m22'].value]]))

        #[V0,g,h] = ana_expr([twoD_poly.params['a0'].value, twoD_poly.params['c0'].value],sym)


        #B0 = V0 / np.matmul( g,np.matmul(mInv,np.transpose(g)) )

        #F = np.matmul(mInv,np.transpose(g))
        #gam111 = twoD_poly.params['gam111'].value
        #gam222 = twoD_poly.params['gam222'].value
        #gam112 = twoD_poly.params['gam112'].value
        #gam122 = twoD_poly.params['gam122'].value

        #B = gam111*F[0]**3 + gam222*F[1]**3 + 3*gam112*F[0]**2*F[1]**1 + 3*gam122*F[0]**1*F[1]**2

        #Bp = V0 * (3 * np.matmul(g,np.matmul(mInv ,np.matmul( h,np.matmul(mInv,np.transpose(g))))) - B )\
        #        / np.matmul( g,np.matmul(mInv,np.transpose(g)) )**2 - 1

        res_poly = [twoD_poly.params['a0'].value, twoD_poly.params['c0'].value, 0, 0, 0]

    result_dict={"name":name,"functional":functional,"BM_poly":BM_poly,"poly":twoD_poly,"res_poly":res_poly,"e_list":e_list,"param_list":param_list,"eRy_list":eRy_list,"v_list":v_list,"geo":geo,"vrat":vrat,"nat":nat}
    return result_dict

#def pressure_calc(E_poly,param_list):
def pressure_calc(E_poly,param_list):
    """Calculates the hydrostatic (?) pressure along the parameters specified.

    :E_poly: TODO
    :param_list: TODO
    :returns: TODO

    """
    
    E_der = E_der_2D(E_poly,param_list)
    #print(E_der)
    #[V_0,g,H]=ana_expr([E_poly.params['a0'].value, E_poly.params['c0'].value],'hcp')
    #pressure=[-1*E_der[0]/g[0],-1*E_der[1]/g[1]]
    #print(pressure)

    V_der=[]
    pressure=[[],[]]
    for i in range(len(param_list[0])):
        [V_i,g_i,h_e] = ana_expr([param_list[0][i], param_list[1][i]],'hcp')
        V_der.append([g_i])
        #print(g_i[0][0])
        pressure[0].append(-1*E_der[0][i]/g_i[0][0])
        pressure[1].append(-1*E_der[1][i]/g_i[0][1])

    #print(V_der)
    print(pressure[0][-1]*1.60219*10**(-19 + 30 - 9))
    print(pressure[1][-1]*1.60219*10**(-19 + 30 - 9))
    
    return pressure


# ------------------ ANALYTICAL EXPRESSIONS --------------------
def ana_expr(p_list,geo,vrat):
    """
    Calculates the analytically derived expressions for V_0,g and H defined in paper by Ziambaras 
    and Schröder Physical Review B 68. 064112

    Parameters
    ----------
    param_list : List of parameters as, where sym=hcp,bcc or hcp
       [a_0,c_0,sym] 

    Returns
    -------
    out : list
        [V_0,g,H]
    """
    if geo=="2L":
        q=1./1.*np.sin(np.pi/3)
        return [q*p_list[0]**3*p_list[1] , np.array([[q*3*p_list[0]**2*p_list[1],q*p_list[0]**3]]) , np.array([ [ q*6*p_list[0]*p_list[1], q*3*p_list[0]**2 ] , [ q*3*p_list[0]**2 , 0 ] ])]
    elif geo=="1L":
        q=1./vrat
        print("QQQQ",q)
        return [q*p_list**3, np.array([[q*3*p_list**2]]), np.array([[ q*6*p_list ]])]
    #TODO fix the expressions below
    elif "3L" in geo:
        return [float('nan'),float('nan'),float('nan')]
    else:
        return [float('nan'),float('nan'),float('nan')]
#------------------------ MAIN --------------------------------------------------
def main(argv):

    parser = argparse.ArgumentParser(description='Calculates the bulk modulus.')
    parser.add_argument('-dir', type=str, help='Identifier for the output directories - most commonly the functional used.')
    parser.add_argument('-i'  , type=str, action="store",  help='Identifier for the output files.')
    parser.add_argument('-ref', type=str, action="store",  help='Identifier for the atomic reference calculation.')
    parser.add_argument('-geo', type=str, action="store",  help='Identifier for the symmetry.')

    parser.add_argument('-list',type=str, action="store",  help='Specify a file containing a list of calculations to process.')
    parser.add_argument('-file',type=str, action="store",  help='Specify a file containing ALL RESULTS of calculations.')
    parser.add_argument('-exp',type=str, action="store" ,  help='Supply experimental values of the cell parameters.')

    parser.add_argument('-sys_name', type=str, action="store",  help='Name of system; to be used in table.')
    parser.add_argument('-func_name',type=str, action="store",  help='Name of functional; to be used in table.')
    parser.add_argument('-show_plot' , action="store_true", help='Show interactive plot.')

    parser.add_argument('-results'   , action="store_true", help='Print ASCII table of data and fitting.')
    parser.add_argument('-fit_report', action="store_true", help='Print results of poly. fit method.')
    parser.add_argument('-show_input', action="store_true", help='Print results of poly. fit method.')

    parser.add_argument('-table', type=str, help='OPTIONAL - Name of .tex table containing the bulkmod fit.')
    parser.add_argument('-plot',  type=str, help='OPTIONAL - Name_of_the_violin_plot.pdf plot.')
    parser.add_argument('-violin_plot',  type=str, help='OPTIONAL - Name_of_the_violin_plot.pdf plot.')
    parser.add_argument('-bar_plot',  type=str, help='OPTIONAL - Name_of_the_violin_plot.pdf plot.')

    args = parser.parse_args()
    calc_list=[]
    result_dict_list=[]

    print()

    # Read list of inputs, that is calculations to be processed
    if args.list:
        inp_list = open(args.list,"r")
        for line in inp_list:
            line_list = line.split(",")
            #TESTING
            #!print("Input line:\n",line_list,"\n")
            in_geo    = line_list[0].strip()
            #in_exp    = float(line_list[1].strip())
            in_dir    = line_list[1].strip()
            in_i      = line_list[2].strip()
            in_ref    = line_list[3].strip()
            if "+" in in_ref:
                in_ref=in_ref.split("+")
            else:
                in_ref=[in_ref]
            in_sys    = line_list[4].strip()
            in_fun    = line_list[5].strip()

            calc_list.append({"in_dir":in_dir,"in_i":in_i,"in_ref":in_ref,"in_sys":in_sys,"in_fun":in_fun,"in_exp":exp_val(in_sys,0),"in_geo":in_geo})
    # Read inputs from std.in/command line 
    elif (args.dir) and (args.i) and (args.ref) and (args.geo):
        if (args.sys_name) and (args.func_name) and (args.exp):
            calc_list.append([args.dir,args.i,[args.ref],args.sys_name,args.func_name,args.geo,args.exp])
        elif args.exp:
            calc_list.append([args.dir,args.i,[args.ref],"N/A","N/A",args.geo,args.exp])
        else:
            calc_list.append([args.dir,args.i,[args.ref],"N/A","N/A",args.geo])
    elif (args.file) and (args.geo):
        calc_list.append(["N/A","N/A","N/A","N/A","N/A",args.geo])
    else:
        print("Wrong inputs..")
        sys.exit()

    # Make analysis of calculations
    for inputs in calc_list:
        #print("CALC LIST",calc_list)
        #!print("Currently processing",inputs,"\n")
        in_dir = inputs["in_dir"];in_i=inputs["in_i"];in_ref=inputs["in_ref"];in_sys=inputs["in_sys"];in_fun=inputs["in_fun"];in_exp=inputs["in_exp"];in_geo=inputs["in_geo"]

        # TODO: FIX THIS TO BE GENERAL?
        if args.file:
            if in_geo == "2L":
                param_list = [[],[]]
                e_list     = []
                v_list     = []
                eRy_list   = []
                with open(args.file) as inputFile:
                    for line in inputFile:
                        param_list[0].append( float(re.findall("\d+\.\d+",line)[0]))
                        param_list[1].append(float(re.findall("\d+\.\d+",line)[1])*float(re.findall("\d+\.\d+",line)[0]))
                        e_list.append(float( re.findall("-\d+\.\d+",line)[-1]))
        else:
            e_list,param_list,eRy_list,v_list,vrat,nat,nat_type=read_input(in_dir,in_i,in_ref,in_geo)
            if (len(v_list)!=len(e_list)) and (in_geo!="3L-c"):
                print('len(v_list) neq len(e_list)')
                print(v_list)

        if args.show_input:
            if (in_geo == "2L"):
                print('a\t\t','c/a\t\t','E_list')
                for i in range(len(param_list[0])):
                    print(round(param_list[0][i],5),round(param_list[1][i],5),round(e_list[i],4),sep="\t\t")
            elif (in_geo=="3L-ab"):
                print('a\t\t','b\t\t','E_list')
                for i in range(len(param_list[0])):
                    print(round(param_list[0][i],5),round(param_list[1][i],5),round(e_list[i],5),sep="\t\t")
                print()
                print("MINIMUM ENERGY VALUES")
                print(param_list[0][e_list.index(min(e_list))],param_list[1][e_list.index(min(e_list))],min(e_list))
            else:
                print('a\t\t','E_list')
                #print('V_list\t\t\t\t\t','e_list\t\t\t\t','eRy\t\t\t\t')
                for i in range(len(param_list)):
                    print(round(param_list[i],5),round(e_list[i],5),sep="\t\t")
        #sys.exit()

        # Fitting
        if (in_geo=="1L") or (in_geo=="3L-c"):
            #!print()
            result_dict=one_dim_bulkmod(in_sys,in_fun,param_list,v_list,e_list,eRy_list,in_geo,vrat,nat,nat_type)
            result_dict.update({"exp":exp_val(in_sys,0)})
        # Print fit_report
            if args.fit_report:
                #!print() 
                print(result_dict["poly"].fit_report())
                #!print()
        elif (in_geo=="2L") or (in_geo=="3L-ab"):
            #!print()
            result_dict=two_dim_bulkmod(in_sys,in_ref,param_list,v_list,e_list,eRy_list,in_geo,vrat,nat,nat_type)
            result_dict.update({"exp":exp_val(in_sys,0)})
        # Print fit_report
            if args.fit_report:
                #!print() 
                print(result_dict["poly"].fit_report())
                #!print()
        else:
            print("Code not yet implemented for this structure")

        result_dict_list.append(result_dict)
    #TESTING
    #sys.exit()

#------------------------------------------GENERATE OUTPUTS---------------------------------------------------------
        # Generate ASCII table of 
    if args.results:
        for res_dict in result_dict_list:
            #print(res_dict)
            if (in_geo=="1L") or (in_geo=="3L-c"):
                print()
                print("   "+res_dict["name"]+"\t\t"+res_dict["functional"])
                print("\t\t Murnaghan Fit")
                print_ascii_output_table(res_dict["BM_poly"][0],res_dict["geo"],"bm",res_dict["vrat"],res_dict["nat"])
            print()
            print("\t\t Fourth Order Poly. Fit")
            print_ascii_output_table(res_dict["res_poly"],res_dict["geo"],"poly",res_dict["vrat"],res_dict["nat"])
            print()

    # Generate TEX table
    if args.table:
        #generate_tex_table_light(args.table, result_dict_list, in_geo)
        generate_tex_table_light_transposed(args.table, result_dict_list, in_geo)
        #generate_tex_table(args.table, result_dict_list, in_geo)
        print('Tex table ',args.table,'.tex was created.\n')
        #print()

    # PDF Plot
    if args.plot:
        if (in_geo=="1L") or (in_geo=="3L-c"):
            pdf_plot_oneD(args.plot,result_dict_list,'ko','b-','ko','r-')
        elif (in_geo=="2L") or (in_geo=="3L-ab"):
            pdf_plot_twoD(args.plot,result_dict_list[0])
        else:
            sys.exit("Not yet implemented PDF plot for this structure")

    # PDF VIOLIN Plot
    if args.violin_plot:
        if (in_geo=="1L") or (in_geo=="3L-c"):
            pdf_violin_plot_oneD(args.violin_plot,result_dict_list,'ko','b-','ko','r-')
        elif (in_geo=="2L") or (in_geo=="3L-ab"):
            sys.exit("Not yet implemented PDF plot for 2d struct")

    # PDF BAR Plot
    if args.bar_plot:
        if (in_geo=="1L") or (in_geo=="3L-c"):
            pdf_bar_plot_oneD_OLD(args.bar_plot,result_dict_list,'ko','b-','ko','r-')
        elif (in_geo=="2L") or (in_geo=="3L-ab"):
            sys.exit("Not yet implemented PDF plot for 2d struct")

    # Heads-up plot
    if args.show_plot:
        print()
        if (in_geo=="1L") or (in_geo=="3L-c"):
            nbr_uniq_systems=len(np.unique([ res_dict["name"] for res_dict  in result_dict_list]))
            if(nbr_uniq_systems>1):
                fig = pdf_plot_oneD("NO PDF",result_dict_list)
            else:
                fig = data_plot_oneD(result_dict_list)
            #fig = data_plot_oneD(param_list,v_list,e_list,BM_poly,oneD_poly,'k--o','ro','kx','-')
            plt.show()
        elif ((in_geo=="2L") or (in_geo=="3L-ab")) and (len(result_dict_list)==1) :
            #CAN SO FAR ONLY PLOT ONE SET OF DATA
            #res_list=result_dict_list[0]
            #v_list=res_list[-1]; e_list=res_list[-4]; param_list=res_list[-3]; lm_poly=res_list[2];BM_poly=res_list[1];
            #fig = data_plot_surface(param_list,e_list,v_list,'k--o','ro','kx','-')
            fig= data_plot_twoD(result_dict_list[0])
            plt.show()
        else:
            sys.exit("Not plottable...")

#-------------------------------------------------------------------------------

if __name__ == "__main__": 
    main(sys.argv[1:])
