import os, getopt, sys, math, numpy as np, glob, argparse, re

def get_vol_factor(sym_index):
    """Gets the symmetry in string-form

    :sym_index: TODO
    :returns: TODO

    """
    if sym_index==0:
        return 1
    elif sym_index==1:
        return 1
    elif sym_index==2:
        return 4
    elif sym_index==3:
        return 2
    elif sym_index==4:
        return 1

def extract_data(output_file_path,geo):
    """Returns the total energy from the molecule [mol] 

    Parameters
    ----------
    :mol: string containing name of molecule
    :ecut: energy to use as cutof in *.in (ecutwfc)
    :kpoints: energy to use as cutof in *.in (ecutwfc)

    Returns
    -------
    :returns: TODO

    """
    fin_tot_ene=0.0
    l_par1=0.0
    l_par2=0.0
    l_par3=0.0
    sym_index=0
    nat=0

    try:
        out_file=open(output_file_path,"r")
        #print(output_file_path)
        #for num_line,line in enumerate(out_file):
        lines = iter(out_file)
        for line in lines:
            if ("bravais-lattice index" in line):
                sym_index=float(re.findall("\\d+", line)[0])
            elif ("number of atoms/cell" in line):
                nat=float(re.findall("\\d+", line)[0])
            elif ("number of atomic types" in line):
                nat_type=float(re.findall("\\d+", line)[0])
                #if nat_type==nat:
                #    nat
            elif (geo == "1L"):
                if (("total energy" in line) or ("Final energy" in line)) and ("=" in line) and ("!" in line): 
                    fin_tot_ene=float(re.findall("-\d+\.\d+", line)[0])
                if "celldm(1)" in line:
                    l_par1=float(re.findall("\d+\.\d+", line)[0])
            elif geo == "2L":
                if (("total energy" in line) or ("Final energy" in line)) and ("=" in line) and ("!" in line): 
                    fin_tot_ene=float(re.findall("-\d+\.\d+", line)[0])
                if "celldm(1)" in line:
                    #a
                    l_par1=float(re.findall("\d+\.\d+", line)[0])
                    #c/a
                    l_par2=float(re.findall("\d+\.\d+", line)[2])
            elif geo == "3L-c":
                nat=1#?
                if (("total energy" in line) or ("Final energy" in line)) and ("=" in line) and ("!" in line): 
                    fin_tot_ene=float(re.findall("-\d+\.\d+", line)[0])
                if "CELL_PARAMETERS" in line:
                    cell_param_lines = [next(lines) for _ in range(3)]
                    l_par3=float(re.findall("\d+\.\d+", cell_param_lines[-1])[0])
                    l_par2=float(re.findall("\d+\.\d+", cell_param_lines[-1])[1])
                    l_par1=float(re.findall("\d+\.\d+", cell_param_lines[-1])[2])
                    #print(l_par1)
            elif geo == "3L-ab":
                nat=1#?
                if (("total energy" in line) or ("Final energy" in line)) and ("=" in line) and ("!" in line): 
                    fin_tot_ene=float(re.findall("-\d+\.\d+", line)[0])
                if "celldm(1)" in line:
                    scaling_factor=float(re.findall("\d+\.\d+", line)[0])
                    #print('scaling factor',scaling_factor)
                    cell_param_lines = [next(lines) for _ in range(6)]
                    #print(cell_param_lines)
                    l_par1=float(re.findall("\d+\.\d+", cell_param_lines[-3])[0])*scaling_factor
                    l_par2=float(re.findall("\d+\.\d+", cell_param_lines[-2])[1])*scaling_factor
                    l_par3=float(re.findall("\d+\.\d+", cell_param_lines[-1])[2])*scaling_factor
        out_file.close()
        return fin_tot_ene,[l_par1,l_par2,l_par3],nat,sym_index,nat_type
    except IOError:
        print("There is no %s file to use."%(output_file_path))	
        return

#def extract_data(output_file_path,sym):
#    """Returns the total energy from the molecule [mol] 
#
#    Parameters
#    ----------
#    :mol: string containing name of molecule
#    :ecut: energy to use as cutof in *.in (ecutwfc)
#    :kpoints: energy to use as cutof in *.in (ecutwfc)
#
#    Returns
#    -------
#    :returns: TODO
#
#    """
#    fin_tot_ene=0.0
#    l_par1=0.0
#    l_par2=0.0
#    l_par3=0.0
#    try:
#        out_file=open(output_file_path,"r")
#        #print(output_file_path)
#        #for num_line,line in enumerate(out_file):
#        lines = iter(out_file)
#        for line in lines:
#            if (sym == "bcc") or (sym == "fcc") or (sym=="dia") or (sym=="halite"):
#                if (("total energy" in line) or ("Final energy" in line)) and ("=" in line) and ("!" in line): 
#                    #print(line)
#                    fin_tot_ene=float(re.findall("-\d+\.\d+", line)[0])
#                if "celldm(1)" in line:
#                    #print(line)
#                    l_par1=float(re.findall("\d+\.\d+", line)[0])
#                    #print(l_par1)
#            elif sym == "hcp":
#                if (("total energy" in line) or ("Final energy" in line)) and ("=" in line) and ("!" in line): 
#                    fin_tot_ene=float(re.findall("-\d+\.\d+", line)[0])
#                if "celldm(1)" in line:
#                    #a
#                    l_par1=float(re.findall("\d+\.\d+", line)[0])
#                    #c/a
#                    l_par2=float(re.findall("\d+\.\d+", line)[2])
#            elif sym == "ortho-c":
#                if (("total energy" in line) or ("Final energy" in line)) and ("=" in line) and ("!" in line): 
#                    fin_tot_ene=float(re.findall("-\d+\.\d+", line)[0])
#                if "CELL_PARAMETERS" in line:
#                    cell_param_lines = [next(lines) for _ in range(3)]
#                    l_par3=float(re.findall("\d+\.\d+", cell_param_lines[-1])[0])
#                    l_par2=float(re.findall("\d+\.\d+", cell_param_lines[-1])[1])
#                    l_par1=float(re.findall("\d+\.\d+", cell_param_lines[-1])[2])
#                    #print(l_par1)
#            elif sym == "ortho-ab":
#                if (("total energy" in line) or ("Final energy" in line)) and ("=" in line) and ("!" in line): 
#                    fin_tot_ene=float(re.findall("-\d+\.\d+", line)[0])
#                if "celldm(1)" in line:
#                    scaling_factor=float(re.findall("\d+\.\d+", line)[0])
#                    #print('scaling factor',scaling_factor)
#                    cell_param_lines = [next(lines) for _ in range(6)]
#                    #print(cell_param_lines)
#                    l_par1=float(re.findall("\d+\.\d+", cell_param_lines[-3])[0])*scaling_factor
#                    l_par2=float(re.findall("\d+\.\d+", cell_param_lines[-2])[1])*scaling_factor
#                    l_par3=float(re.findall("\d+\.\d+", cell_param_lines[-1])[2])*scaling_factor
#                    #print(l_par1)
#                    #print(l_par2)
#                    #print(l_par3)
#        out_file.close()
#        return fin_tot_ene,[l_par1,l_par2,l_par3]
#    except IOError:
#        print("There is no %s file to use."%(output_file_path))	
#        return

# TODO : Do I need to sort the lists? 
def read_input(output_dir,output_name,ref_list,geo):
    """Reads lattice parameters and total energy (subtracting the atomic reference energy) for all outputfiles in subdirectories containing 'output_dir'
    
    Parameters
    ----------
    :output_dir: String used to find the directories containing the results of interest. 
    Should be "." or "./" if code is run in directory containing only subdirectories of interest.

    :output_name: Unique string contained in all the results file names - often this is the name of the functional used in the calculation.
    :ref_list: Directory containing the results of the atomic reference calculation.

    Returns
    -------
    :returns: 
    param_list - listtor with all the lattice parameters, 
    e_list - listtor with all the energies subtracted by atomic reference, 
    eRy_list - listtor with all energies w/o subtracting atomic reference. 
        """
    bohrToAng=0.529177
    rydToev=13.6057

    print()
    print("output_dir",output_dir)
    print("output_name",output_name)

    eRy_list=[]
    e_list=[]
    v_list=[]
    param_list=[]
    ref_ene=0
    sym_index=0
    nat=0

    if output_dir[-1]=="*":
    #strip until first/ = output_dir
        dir_path=output_dir[:output_dir.rfind("/")]
        output_dir=output_dir[output_dir.rfind("/"):].strip("/").rstrip("*")
        print('output_dir',output_dir)
        print('dir_path',dir_path)
        dir_list = next(os.walk(dir_path))[1]
    else:
        dir_list = next(os.walk('./'))[1]
        dir_path="."
        print('dir_list',dir_list)

    # Find the atomic reference energy atom_ref_ene, if 0 is passed the atom_ref_ene=0, if file or directory is passed,
    # search for energy there.
    print( "Searching for atomic energy" )
    for atom_ref in ref_list:
        if atom_ref == "0":
            ref_ene=0
        elif os.path.isfile(atom_ref):
            #print(atom_ref)
            atom_ref_ene = extract_data(atom_ref,geo)
            ref_ene+= atom_ref_ene[0]
        elif os.path.isdir(atom_ref):
            atom_ref = atom_ref.rstrip("/")
            ref_path=glob.glob('%s/*%s*'%(atom_ref,output_name))[0]
            atom_ref_ene = extract_data(ref_path,geo)
            ref_ene+= atom_ref_ene[0]
        else:
            print('./*%s*/*%s*'%(atom_ref,output_name))
            ref_path=glob.glob('./*%s*/*%s*'%(atom_ref,output_name))[0]
            atom_ref_ene = extract_data(ref_path,geo)
            ref_ene+= atom_ref_ene[0]

    print('Atomic reference energy is ',ref_ene,'Ry.\n')

    # Find the energy and lattice parameters of the QE calculations
    # IF script is run in directory containing the output files for the calculations of interest
    # If output_dir is path to directory containing calculations of interest
    if os.path.isdir(output_dir):
        output_dir=output_dir.rstrip("/")
        for ending in ["out","log",""]:
            inputF_id='*%s*'%(output_name) if ending=="" else '*%s*.%s'%(output_name,ending) # width in inches
            #print("Checking:",'%s/%s'%(output_dir,inputF_id))
            in_file_path_list=glob.glob('%s/%s'%(output_dir,inputF_id))
            if (len(in_file_path_list) > 0):
                for output_file_path in in_file_path_list:
                    energy,a_par,nat,sym_index,nat_type = extract_data(output_file_path,geo)
                    vrat=get_vol_factor(sym_index)*nat
                    if (energy == 0 ):
                        print("Found no energy for file " + output_file_path)
                        print("Omitting this result!\n")
                    elif (a_par == 0):
                        print("Found no lattice parameter for file " + output_file_path)
                        print("Omitting this result!\n")
                    else:
                        eRy_list.append(energy)
                        if geo == "2L":
                            e_list.append(rydToev*(energy-nat*ref_ene))
                            param_list.append([bohrToAng*a_par[0],a_par[1]])
                        elif geo=="3L-c":
                            e_list.append(rydToev*(energy-nat*ref_ene))
                            param_list.append(a_par[0])
                            v_list.append(a_par[0],a_par[1],a_par[2])
                        elif geo=="3L-ab":
                            e_list.append(rydToev*(energy-nat*ref_ene))
                            param_list.append([bohrToAng*a_par[0],bohrToAng*a_par[1]])
                            v_list.append(a_par[0],a_par[1],a_par[2])
                        else:
                            if len(ref_list)==nat:
                                e_list.append(rydToev*(energy-ref_ene))
                            else:
                                e_list.append(rydToev*(energy-nat*ref_ene))
                            param_list.append(bohrToAng*a_par[0])
            else:
                print("No outputfiles ending with .%s\n"%(ending))
        param_list,e_list=zip(*sorted(zip(param_list,e_list)))
        param_list,eRy_list=zip(*sorted(zip(param_list,eRy_list)))
    #If scrip is exec at parent node and output_dir is specifier for *multiple subdirectories* containing ouputs 
    elif (len(dir_list)>0):
        #dir_list contains all subdirs of dir_path, now we filter out those result dirs which matches output_dir 
        #filtered_list= [string for string in dir_list if output_dir in string]
        filtered_list= [res_dir for res_dir in dir_list if res_dir.startswith(output_dir)]
        for directory in filtered_list:
            #print('%s/*%s*/*%s*'%(dir_path,directory,output_name))
            #print(glob.glob('%s/*%s*/%s'%(dir_path,directory,output_name)))
            #output_file_path=glob.glob('%s/*%s*/*%s*'%(dir_path,directory,output_name))[0]
            #print("Checking:",'%s/%s/%s'%(dir_path,directory,output_name))
            #print('directory',directory)
            output_file_path=glob.glob('%s/%s/%s'%(dir_path,directory,output_name))[0]
            #print("output_file_path",output_file_path)
            energy,a_par,nat,sym_index,nat_type = extract_data(output_file_path,geo)
            vrat=get_vol_factor(sym_index)*nat
            if (energy == 0 ):
                print("Found no energy for file " + output_file_path)
                print("Omitting this result!\n")
            elif (a_par == 0):
                print("Found no lattice parameter for file " + output_file_path)
                print("Omitting this result!\n")
            else:
                eRy_list.append(energy)
                if geo == "2L":
                    e_list.append(rydToev*(energy-nat*ref_ene))
                    param_list.append([bohrToAng*a_par[0],a_par[1]])
                elif geo=="3L-c":
                    e_list.append(rydToev*(energy-nat*ref_ene))
                    param_list.append(a_par[0])
                    v_list.append(a_par[0],a_par[1],a_par[2])
                elif geo=="3L-ab":
                    e_list.append(rydToev*(energy-nat*ref_ene))
                    param_list.append([bohrToAng*a_par[0],bohrToAng*a_par[1]])
                    v_list.append(a_par[0],a_par[1],a_par[2])
                else:
                    if len(ref_list)==nat:
                        #e_list.append((1/(vrat*nat_type))*rydToev*(energy-1*ref_ene))
                        e_list.append(1/nat_type*rydToev*(energy-1*ref_ene))
                    else:
                        #e_list.append((1/nat)*rydToev*(energy-nat*ref_ene))
                        e_list.append((1/1)*rydToev*(energy/2-ref_ene))
                    param_list.append(bohrToAng*a_par[0])
            #print()
        param_list,e_list=zip(*sorted(zip(param_list,e_list)))
        param_list,eRy_list=zip(*sorted(zip(param_list,eRy_list)))
    else:
        print("No output files were found")

    # Volume list construction
    #vrat=get_vol_factor(sym_index)*nat
    vrat=get_vol_factor(sym_index)
    if geo=="1L":
        for a in param_list:
            v_list.append((1/vrat)*a**3)
    elif geo == "2L":
        p1=[]
        p2=[]
        for a in param_list:
            v_list.append((1/vrat)*a[0]**2*np.sin(np.pi/3)*a[1])
            p1.append(a[0]) # a
            p2.append(a[1]) # c
        param_list=[p1,p2]
    #TODO Fix the below volume list!
    elif geo == "3L-ab":
        p1=[]
        p2=[]
        for a in param_list:
            p1.append(a[0]) # a
            p2.append(a[1]) # c
        param_list=[p1,p2]

    return e_list,param_list,eRy_list,v_list,vrat,nat,nat_type
