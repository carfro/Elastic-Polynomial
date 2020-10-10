import numpy as np
import sys
from   prettytable import PrettyTable

def generate_tex_table_light_transposed(output_name,result_dict_list,geo):
    """Generates a table of the opt lattice parameter, elastic energy and bulk modulus 

    Parameters
    ----------
    :output_name: name/path to outputfile
    :pfit: listtor of constants from third order fitting

    -------
    Returns
    :returns: NULL

    """


    if (output_name.split(".")[-1]=="tex"):
        data_table = open(output_name,'w+')
    else:
        data_table = open(output_name+".tex",'w+')

    species_list        = [ res_dict["name"] for res_dict in result_dict_list ]
    functional_list     = [ res_dict["functional"] for res_dict in result_dict_list ]
    unique_species      = np.unique(species_list)
    unique_functionals  = np.unique(functional_list).tolist()
    #unique_functionals.sort(key  = len)

    len_func  = len(unique_functionals)
    len_spec  = len(unique_species)
    nbr_columns = len_spec*(2 +len_func) 

    index_species       = [[i for i, e in enumerate(species_list) if e == u_s] for u_s in unique_species]
    index_functionals   = [[i for i, e in enumerate(functional_list) if e == u_s] for u_s in unique_functionals]
    plot_order          = [[i for j in range(len_func) for i in index_functionals[j] if i in S] for S in index_species]

    #Makes the tables in the same order as the input
    plot_order.sort()
    [sublist.sort() for sublist in plot_order]
    
    max_length          = 27
    table_length        = len(index_species) + len(index_functionals)

    def find_max_list(list):
        list_len = [len(i) for i in list]
        return max(list_len)

    col_tot='c'*(find_max_list(plot_order)+2)
    data_table.write("\\begin{tabular}{l|"+col_tot+"}\n")

    i=0
    for plot_list in plot_order:
        print(plot_list)
        data_table.write("\hline\\\\[-3.5mm]\n")
        row_2 = []
        row_3 = ["$a_0$ \hspace{2.2mm} [Å]"]
        row_4 = ["$E_{coh}$ [eV]"]
        row_5 = ["$B_0$ \hspace{1.9mm} [GPa]"]
        row_6 = ["$B'$ "]
        for cur_res in plot_list:
            res_dict = result_dict_list[cur_res];BM_poly=res_dict["BM_poly"];pfit=res_dict["res_poly"];geo=res_dict["geo"]
            print(res_dict["name"],res_dict["functional"])

            if geo=="1L":
                if i==0:
                    row_2.append("\\textbf{"+res_dict["name"]+"}")
                    row_2.append(res_dict["functional"])
                else:  
                    row_2.append(res_dict["functional"])
                for j in range(len(pfit)):
                    if j==0:
                        row_3.append("$"+'{:.3f}'.format(round(pfit[0],3))+"$")
                    elif j==1:
                        row_4.append("$"+'{:.3f}'.format(round(pfit[j],3))+"$")
                    elif j==2:
                        row_5.append("$"+'{:.1f}'.format(round(pfit[j],1))+"$")
                    else:
                        row_6.append("$"+'{:.2f}'.format(round(pfit[j],2))+"$")
            elif geo == "2L":
                if i==0:
                    row1.append("\\textbf{"+res_dict["name"]+"}")
                    row2.append(res_dict["functional"])
                else:  
                    row1.append("\\textbf{}")
                    row2.append(res_dict["functional"])
                for j in range(len(pfit)):
                    if j==0:
                        row_3.append("$"+'{:.3f}'.format(round(pfit[0],3))+"$")
                    elif j==1:
                        row_4.append("$"+'{:.3f}'.format(round(pfit[j],3))+"$")
                    elif j==2:
                        row_5.append("$"+'{:.1f}'.format(round(pfit[j],1))+"$")
                    else:
                        row_6.append("$"+'{:.2f}'.format(round(pfit[j],2))+"$")
            i+=1
        #Add experimental columns
        if len( res_dict )>10:
            # ONLY INCLUDES EXP*
            #for exp_values in [exp_table_row(res_dict["name"]),exp_corrected_table_row(res_dict["name"])]:
            for exp_values in [exp_corrected_table_row(res_dict["name"])]:
                j=0
                for item in exp_values:
                    if j==1:
                        row_2.append(item)
                    elif j==2:
                        row_3.append(item)
                    elif j==3:
                        row_4.append(item)
                    elif j==4:
                        row_5.append(item)
                    elif j==5:
                        row_6.append(item)
                    j+=1
        print("i",i)
        data_table.write(' & '.join((word.ljust(max_length) for word in row_2)) + " \\\\ \n")
        data_table.write("\hline\\\[-5pt]\n")
        #data_table.write("\hline\hline%\n")
        data_table.write(' & '.join((word.ljust(max_length) for word in row_3)) + " \\\\[+1.5mm] \n")
        data_table.write(' & '.join((word.ljust(max_length) for word in row_4)) + " \\\\[+1.5mm] \n")
        data_table.write(' & '.join((word.ljust(max_length) for word in row_5)) + " \\\\[+1.5mm] \n")
        data_table.write(' & '.join((word.ljust(max_length) for word in row_6)) + " \\\\[+1.5mm] \n")
        #data_table.write("\hline\hline%\n")
        data_table.write("\hline%\n")

        i=0

    data_table.write("\end{tabular}\n")
    data_table.close()
    pass

def generate_tex_table_light(output_name,result_dict_list,geo):
    """Generates a table of the opt lattice parameter, elastic energy and bulk modulus 

    Parameters
    ----------
    :output_name: name/path to outputfile
    :pfit: listtor of constants from third order fitting

    -------
    Returns
    :returns: NULL

    """


    if (output_name.split(".")[-1]=="tex"):
        data_table = open(output_name,'w+')
    else:
        data_table = open(output_name+".tex",'w+')

    species_list        = [ res_dict["name"] for res_dict in result_dict_list ]
    functional_list     = [ res_dict["functional"] for res_dict in result_dict_list ]
    unique_species      = np.unique(species_list)
    unique_functionals  = np.unique(functional_list).tolist()
    unique_functionals.sort(key  = len)

    len_func  = len(unique_functionals)
    len_spec  = len(unique_species)
    size_grid = len_func*len_spec

    index_species       = [[i for i, e in enumerate(species_list) if e == u_s] for u_s in unique_species]
    index_functionals   = [[i for i, e in enumerate(functional_list) if e == u_s] for u_s in unique_functionals]
    plot_order          = [[ i for j in range(len_func) for i in index_functionals[j] if i in S] for S in index_species]
    #Makes the tables in the same order as the input
    plot_order.sort()
    [sublist.sort() for sublist in plot_order]
    
    max_length           = 20


    data_table.write("\\begin{tabular}{lc|cccc}\n")
    data_table.write("\hline\\\\[-3.5mm]\n")
    header  = ["System","Functional","$a_0$ [Å]","$E_0$ [eV]","$B_0$ [GPa]","$B'$ [N/A] \\\\ \n"]
    data_table.write( " & ".join(word.ljust(max_length) for word in header ) )
    data_table.write("\hline\hline%\n")

    i=0
    for plot_list in plot_order:
        for cur_res in plot_list:
            res_dict = result_dict_list[cur_res];BM_poly=res_dict["BM_poly"];pfit=res_dict["res_poly"];geo=res_dict["geo"];vrat=res_dict["vrat"]

            #print("SYMM IN TABEL",geo,"VRAT",vrat)
            if (geo == "1L"):
                if i==0:
                    row = ["\\textbf{"+res_dict["name"]+"}",res_dict["functional"]]
                    print(res_dict["name"])
                else:  
                    row = ["\\textbf{}",res_dict["functional"]]
                for j in range(len(pfit)):
                    if j==0:
                        row.append("$"+'{:.3f}'.format(round(pfit[0],3))+"$")
                        print(res_dict["functional"],'{:.3f}'.format(round(pfit[0]*1.8897259886,3)),'Bhor')
                    elif j==1:
                        row.append("$"+'{:.3f}'.format(round(pfit[j],3))+"$")
                    elif j==2:
                        row.append("$"+'{:.1f}'.format(vrat*round(pfit[j]*1.60219*10**(-19 + 30 - 9),1))+"$")
                    else:
                        row.append("$"+'{:.2f}'.format(round(pfit[j],2))+"$")
                data_table.write(' & '.join(word.ljust(max_length) for word in row)  + " \\\\ \n")
            elif geo == "2L":
                #vrat = 1./2.
                if i==0:
                    row = ["\\textbf{"+res_dict["name"]+"}",res_dict["functional"]]
                else:  
                    row = ["\\textbf{}",res_dict["functional"]]
                for j in range(len(pfit)):
                    if j==0:
                        row_3.append("$"+'{:.3f}'.format(round(pfit[0],3))+"$")
                    elif j==1:
                        row_4.append("$"+'{:.3f}'.format(round(pfit[j],3))+"$")
                    elif j==2:
                        row_5.append("$"+'{:.1f}'.format(vrat*round(pfit[j]*1.60219*10**(-19 + 30 - 9),1))+"$")
                    else:
                        row_6.append("$"+'{:.2f}'.format(round(pfit[j],2))+"$")
                data_table.write(' & '.join(row) + " \\\\ \n")
            i+=1
        if len( res_dict )>10:
            #data_table.write("\\textbf{}  & \\textbf{Experiment} & $\\mathbf{%s}$ & $\\mathbf{X}$ & $\\mathbf{X}$ & $\\mathbf{X}$ \\\\\n"%(str(res_dict["exp"])))
            #data_table.write(' & '.join(word.ljust(max_length) for word in exp_table_row(res_dict["name"]))  + " \\\\ \n")
            data_table.write(' & '.join(word.ljust(max_length) for word in exp_corrected_table_row(res_dict["name"]))  + " \\\\ \n")
        print("i",i)
        data_table.write("\hline\n")
        i=0
    data_table.write("\end{tabular}\n")
    data_table.close()
    pass

def generate_tex_table(output_name,result_dict_list,sym):
    """Generates a table of the opt lattice parameter, elastic energy and bulk modulus 

    Parameters
    ----------
    :output_name: name/path to outputfile
    :pfit: listtor of constants from third order fitting

    -------
    Returns
    :returns: NULL

    """


    if (output_name.split(".")[-1]=="tex"):
        data_table = open(output_name,'w+')
    else:
        data_table = open(output_name+".tex",'w+')

    species_list        = [ res_dict["name"] for res_dict in result_dict_list ]
    functional_list     = [ res_dict["functional"] for res_dict in result_dict_list ]
    unique_species      = np.unique(species_list)
    unique_functionals  = np.unique(functional_list).tolist()
    unique_functionals.sort(key  = len)

    len_func  = len(unique_functionals)
    len_spec  = len(unique_species)
    size_grid = len_func*len_spec

    index_species       = [[i for i, e in enumerate(species_list) if e == u_s] for u_s in unique_species]
    index_functionals   = [[i for i, e in enumerate(functional_list) if e == u_s] for u_s in unique_functionals]
    plot_order          = [[ i for j in range(len_func) for i in index_functionals[j] if i in S] for S in index_species]

    data_table.write("\\begin{tabular}{lc|cccc|cccc}\n")
    data_table.write("\hline\\\\[-3.5mm]\n")
    #data_table.write("%s & & & Bulkdata & & \\\\ \n")
    #data_table.write("& $a_0$ [Å] & $E_0$ [eV] & $B_0$ [GPa] & $B'$ [N/A] & $a_0$ [Å] & $E_0$ [eV] & $B_0$ [GPa] & $B'$ [N/A] \\\\ \n")
    data_table.write("System  & Method & $a_0 $ [Å] & $E_0$ [eV] & $B_0$ [GPa] & $B'$ [N/A] & $a_0^{\\text{Mur}}$ & $E_0^{\\text{Mur}}$ [eV] & $B_0^{\\text{Mur}}$ [GPa] & $B'^{\\text{Mur}}$ [N/A]\\\\ \n")
    data_table.write("\hline\hline%\n")

    i=0
    for plot_list in plot_order:
        for cur_res in plot_list:
            res_dict = result_dict_list[cur_res];BM_poly=res_dict["BM_poly"];pfit=res_dict["res_poly"];geo=res_dict["geo"];vrat=res_dict["vrat"]

            if geo=="1L":
                if i==0:
                    row = ["\\textbf{%s}  & %s \t\t"%(res_dict["name"],res_dict["functional"])]
                else:  
                    row = ["\\textbf{}  & %s \t\t"%(res_dict["functional"])]
                #row2 = ["Murnaghan fit \t\t& \\textbf{%s %s}"%(res_dict["name"],res_dict["functional"])]
                for j in range(len(pfit)):
                    if j==0:
                        row.append("$"+str(round(pfit[0],3))+"$")
                    elif j==2:
                        row.append("$"+'{:.1f}'.format(round(pfit[j],1))+"$")
                        #row.append("$"+str(round(pfit[j]*1.60219*10**(-19 + 30 - 9),3))+"$")
                    else:
                        row.append("$"+str(round(pfit[j],3))+"$")
                for j in range(len(BM_poly[0])):
                    if j==0:
                        row.append("$"+str(round((BM_poly[0][j]/vrat)**(1./3),3))+"$")
                    elif j==2:
                        row.append("$"+str(round(BM_poly[0][j]*1.60219*10**(-19 + 30 - 9),1))+"$")
                    else:
                        row.append("$"+str(round(BM_poly[0][j],3))+"$")
                #data_table.write(' & '.join(row2) + " \\\\ \n")
                data_table.write(' & '.join(row)  + " \\\\ \n")
            elif geo=="2L":
                # TODO : MURNAGHAN FIT 
                vrat = 1./2.
                row = ["\\textbf{%s}  &  %s "%(res_dict["name"],res_dict["functional"])]
                for j in range(len(pfit)):
                    if j==0:
                        row.append(str(pfit[0]))
                    elif j==2:
                        row.append(str(round(pfit[j]*1.60219*10**(-19 + 30 - 9),3)))
                    else:
                        row.append(str(round(pfit[j],2)))
                data_table.write(' & '.join(row) + " \\\\ \n")
            i+=1
        if len( res_dict )>10:
            data_table.write("\\textbf{}  & \\textbf{Experiment} & $\\mathbf{%s}$ & $\\mathbf{X}$ & $\\mathbf{X}$ & $\\mathbf{X}$ & $-$ & $-$ & $-$ & $-$ \\\\\n"%(str(res_dict["exp"])))
        print("i",i)
        data_table.write("\hline\n")
        i=0
    data_table.write("\end{tabular}\n")
    data_table.close()
    pass

def print_ascii_output_table(pfit,geo,sort,vrat,nat):
    """Prints an ascii table of the opt lattice parameter, elastic energy and bulk modulus 

    Parameters
    ----------
    :pfit: listtor of constants from third order fitting

    Returns
    -------
    :returns: NULL
    """
    if geo=="1L":
        #vrat = 1./2. if sym=="bcc" else 1./4. # width in inches
        ascii_table = PrettyTable()
        ascii_table.field_names=["a_0 [Å]"," E_0 [eV] "," B_0 [GPa] "," B' [N/A]"]
        row = []
        for i in range(len(pfit)):
            if (i==0):
                if sort=="bm":
                    row.append('{:.3f}'.format(round((pfit[i]*vrat)**(1./3),3)))
                else:
                    row.append('{:.3f}'.format(round(pfit[i],3)))
            elif (i==1):
                if sort=="bm":
                    row.append('{:.3f}'.format(round((-1*pfit[i]/nat),3)))
                else:
                    row.append('{:.3f}'.format(round(pfit[i],3)))
            elif (i==2):
                if sort=="bm":
                    row.append('{:.1f}'.format(round(pfit[i]*1.60219*10**(-19 + 30 - 9)/nat,2)))
                else:
                    row.append('{:.1f}'.format(round(pfit[i],1)))
            elif (i==3):
                row.append('{:.2f}'.format(round(pfit[i],2)))
        ascii_table.add_row(row)
        print(ascii_table)
    elif (geo == "2L"):
        #vrat = 1./2.
        ascii_table = PrettyTable()
        ascii_table.field_names=["a_0 [Å]","c/a" ," E_0 [eV] "," B_0 [GPa] "," B' [N/A]"]
        row = []
        for i in range(len(pfit)):
            if (i==0):
                if sort=="bm":
                    row.append('{:.3f}'.format(round((pfit[i]*vrat)**(1./3),3)))
                else:
                    row.append('{:.3f}'.format(round(pfit[i],3)))
            elif (i==1):
                if sort=="bm":
                    row.append('{:.3f}'.format(round((-1*pfit[i]/nat),3)))
                else:
                    row.append('{:.3f}'.format(round(pfit[i],3)))
            elif (i==2):
                if sort=="bm":
                    row.append('{:.1f}'.format(round(pfit[i]*1.60219*10**(-19 + 30 - 9)/nat,2)))
                else:
                    row.append('{:.1f}'.format(round(pfit[i],1)))
            elif (i==3):
                row.append('{:.2f}'.format(round(pfit[i],2)))
        ascii_table.add_row(row)
    elif (sym=="3L-ab"):
        ascii_table = PrettyTable()
        ascii_table.field_names=["a_0 [Å]","b_0 [Å]" ," E_0 [eV] "," B_0 [GPa] "," B' [N/A]"]
        row = []
        for i in range(len(pfit)):
            if (i==0):
                if sort=="bm":
                    row.append('{:.3f}'.format(round((pfit[i]*vrat)**(1./3),3)))
                else:
                    row.append('{:.3f}'.format(round(pfit[i],3)))
            elif (i==1):
                if sort=="bm":
                    row.append('{:.3f}'.format(round((-1*pfit[i]/nat),3)))
                else:
                    row.append('{:.3f}'.format(round(pfit[i],3)))
            elif (i==2):
                if sort=="bm":
                    row.append('{:.1f}'.format(round(pfit[i]*1.60219*10**(-19 + 30 - 9)/nat,2)))
                else:
                    row.append('{:.1f}'.format(round(pfit[i],1)))
            elif (i==3):
                row.append('{:.2f}'.format(round(pfit[i],2)))
        ascii_table.add_row(row)
        print(ascii_table)
    pass

def print_ascii_energy_table(param_list,e_list,eRy_list):
    """Prints an ascii table of the opt lattice parameter, elastic energy and bulk modulus 

    Parameters
    ----------
    :pfit: listtor of constants from third order fitting

    Returns
    -------
    :returns: NULL

    """

    ascii_table = PrettyTable()

    ascii_table.field_names=["a_0 [Å]"," E-E_atom [eV] "," E [Ry] "]

    for i in range(len(eRy_list)):
        row=[]
        row.append(str(param_list[i]))
        row.append(str(e_list[i]))
        row.append(str(eRy_list[i]))
        ascii_table.add_row(row)
    print(ascii_table)

    pass

def exp_table_row(system):
    """Finds the tex row of experimental values for the system specified.

    :system: System your looking for
    :returns: tex table row values

    """
    exp_values={
            "C":["\\textbf{}","\\textbf{Experiment}"," $\\mathbf{3.566}$"," $\\mathbf{7.799}$","$\\mathbf{454.7}$","$\\mathbf{X}$"],
            "Si":["\\textbf{}","\\textbf{Experiment}"," $\\mathbf{5.430}$"," $\\mathbf{4.743}$ ","$\\mathbf{100.8}$","$\\mathbf{X}$"],
            "SiC":["\\textbf{}","\\textbf{Experiment}","$\\mathbf{4.360}$"," $\\mathbf{6.687}$ ","$\\mathbf{229.1}$","$\\mathbf{X}$"],
            "Ge":["\\textbf{}","\\textbf{Experiment}","$\\mathbf{5.652}$"," $\\mathbf{3.899}$","$\\mathbf{77.3}$","$\\mathbf{X}$"],
            "GaAs":["\\textbf{}","\\textbf{Experiment}","$\\mathbf{5.648}$"," $\\mathbf{3.436}$","$\\mathbf{76.7}$","$\\mathbf{X}$"], 
            "NaCl":["\\textbf{}","\\textbf{Experiment}","$\\mathbf{5.594}$","$\\mathbf{3.372}$","$\\mathbf{27.6}$","$\\mathbf{X}$"],
            "NaF":["\\textbf{}","\\textbf{Experiment}","$\\mathbf{4.609}$","$\\mathbf{4.026}$","$\\mathbf{53.1}$","$\\mathbf{X}$"],
            "LiCl":["\\textbf{}	","\\textbf{Experiment}","$\\mathbf{5.088}$","$\\mathbf{3.632}$","$\\mathbf{38.7}$","$\\mathbf{X}$"],
            "LiF":["\\textbf{}","\\textbf{Experiment}","$\\mathbf{4.010}$","$\\mathbf{4.542}$","$\\mathbf{76.3}$","$\\mathbf{X}$"],
            "MgO":["\\textbf{}","\\textbf{Experiment}","$\\mathbf{4.203}$","$\\mathbf{5.363}$","$\\mathbf{165}$","$\\mathbf{X}$"],
            "Ag":["\\textbf{}","\\textbf{Experiment}","$\\mathbf{4.086}$","$\\mathbf{2.487}$","$\\mathbf{100.7}$","$\\mathbf{4.73}$"],
            "Au":["\\textbf{}","\\textbf{Experiment}","$\\mathbf{4.079}$","$\\mathbf{3.109}$","$\\mathbf{173.2}$","$\\mathbf{6.40}$"],
            "Cu":["\\textbf{}","\\textbf{Experiment}","$\\mathbf{3.615}$","$\\mathbf{3.513}$","$\\mathbf{137.0}$","$\\mathbf{4.88}$"],
            "Pd":["\\textbf{}","\\textbf{Experiment}","$\\mathbf{3.888}$","$\\mathbf{3.928}$","$\\mathbf{180.8}$","$\\mathbf{5}$"],
            "Rh":["\\textbf{}","\\textbf{Experiment}","$\\mathbf{3.796}$","$\\mathbf{5.873}$","$\\mathbf{270.4}$","$\\mathbf{4.5}$"],
            "Al":["\\textbf{}","\\textbf{Experiment}","$\\mathbf{4.049}$","$\\mathbf{3.389}$","$\\mathbf{77.1}$","$\\mathbf{4.45}$"],
            "Pt":["\\textbf{}","\\textbf{Experiment}","$\\mathbf{3.924}$","$\\mathbf{5.845}$","$\\mathbf{278.3}$","$\\mathbf{5.18}$"],
            "Ni":["\\textbf{}","\\textbf{Experiment}","$\\mathbf{3.523}$","$\\mathbf{4.436}$","$\\mathbf{186.0}$","$\\mathbf{4}$"],
            "Na":["\\textbf{}","\\textbf{Experiment}","$\\mathbf{4.226}$","$\\mathbf{1.109}$","$\\mathbf{6.8}$","$\\mathbf{4.13}$"],
            "Fe":["\\textbf{}","\\textbf{Experiment}","$\\mathbf{2.867}$","$\\mathbf{4.280}$","$\\mathbf{175.1}$","$\\mathbf{4.6}$"],
            "Li":["\\textbf{}","\\textbf{Experiment}","$\\mathbf{x}$","$\\mathbf{x$","$\\mathbf{x}$","$\\mathbf{x}$"],
            }
    
    return exp_values[system]


def exp_corrected_table_row(system):
    """Finds the tex row of back-corrected experimental values for the system specified.

    :system: System your looking for
    :returns: tex table row values

    """
    exp_values_cor={
            "C":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{3.543}$","$\\mathbf{7.583}$","$\\mathbf{443}$","$\\mathbf{X}$"],
            "Si":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{5.416}$","$\\mathbf{4.681}$","$\\mathbf{99.2}$","$\\mathbf{X}$"],
            "SiC":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{4.342}$","$\\mathbf{6.488}$","$\\mathbf{225}$","$\\mathbf{X}$"],
            "Ge":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{5.640}$","$\\mathbf{3.863}$","$\\mathbf{75.8}$","$\\mathbf{X}$"],
            "GaAs":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{5.638}$","$\\mathbf{3.393}$","$\\mathbf{75.6}$","$\\mathbf{X}$"],
            "NaCl":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{5.565}$","$\\mathbf{3.341}$","$\\mathbf{26.6}$","$\\mathbf{X}$"],
            "NaF":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{4.579}$","$\\mathbf{3.978}$","$\\mathbf{51.4}$","$\\mathbf{X}$"],
            "LiCl":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{5.056}$","$\\mathbf{3.591}$","$\\mathbf{35.4}$","$\\mathbf{X}$"],
            "LiF":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{3.964}$","$\\mathbf{4.471}$","$\\mathbf{69.8}$","$\\mathbf{X}$"],
            "MgO":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{4.184}$","$\\mathbf{5.271}$","$\\mathbf{169.8}$","$\\mathbf{X}$"],
            "Ag":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{4.070}$","$\\mathbf{2.964}$","$\\mathbf{105.7}$","$\\mathbf{4.73}$"],
            "Au":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{4.067}$","$\\mathbf{3.835}$","$\\mathbf{182.0}$","$\\mathbf{6.40}$"],
            "Cu":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{3.595}$","$\\mathbf{3.472}$","$\\mathbf{144.3}$","$\\mathbf{4.88}$"],
            "Pd":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{3.876}$","$\\mathbf{3.700}$","$\\mathbf{187.2}$","$\\mathbf{5}$"],
            "Rh":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{3.786}$","$\\mathbf{5.876}$","$\\mathbf{277.1}$","$\\mathbf{4.5}$"],
            "Al":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{4.022}$","$\\mathbf{3.431}$","$\\mathbf{72.2}$","$\\mathbf{4.45}$"],
            "Pt":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{3.917}$","$\\mathbf{5.866}$","$\\mathbf{285.5}$","$\\mathbf{5.18}$"],
            "Ni":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{3.510}$","$\\mathbf{4.477}$","$\\mathbf{192.5}$","$\\mathbf{4}$"],
            "Fe":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{2.855}$","$\\mathbf{4.322}$","$\\mathbf{168.3}$","$\\mathbf{4.6}$"],
            "Na":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{4.205}$","$\\mathbf{1.119}$","$\\mathbf{7.9}$","$\\mathbf{4.13}$"],
            "Li":["\\textbf{}","\\textbf{Experiment*}","$\\mathbf{3.443}$","$\\mathbf{1.669}$","$\\mathbf{13.1}$","$\\mathbf{5.51}$"],
            }
    
    return exp_values_cor[system]

def exp_val(system,val):
    """Finds the tex row of back-corrected experimental values for the system specified.

    :system: System your looking for
    :returns: tex table row values

    """
    exp_values_cor={
            "C"   :[3.543,7.583,443  ,"X"],
            "Si"  :[5.416,4.681,99.2 ,"X"],
            "SiC" :[4.342,6.488,225  ,"X"],
            "Ge"  :[5.640,3.863,75.8 ,"X"],
            "GaAs":[5.638,3.393,75.6 ,"X"],
            "NaCl":[5.565,3.341,26.6 ,"X"],
            "NaF" :[4.579,3.978,51.4 ,"X"],
            "LiCl":[5.056,3.591,35.4 ,"X"],
            "LiF" :[3.964,4.471,69.8 ,"X"],
            "MgO" :[4.184,5.271,169.8,"X"],
            "Ag"  :[4.070,2.964,105.7,4.73],
            "Au"  :[4.067,3.835,182.0,6.40],
            "Cu"  :[3.595,3.472,144.3,4.88],
            "Pd"  :[3.876,3.700,187.2,5],
            "Rh"  :[3.786,5.876,277.1,4.5],
            "Al"  :[4.022,3.431,72.2 ,4.45],
            "Pt"  :[3.917,5.866,285.5,5.18],
            "Ni"  :[3.510,4.477,192.5,4],
            "Fe"  :[4.226,1.109,6.8,4.13],
            "Na"  :[4.205,1.119,7.9,4.13],
            "Li"  :[3.443,1.669,13.1,3.51],
            }
    
    return exp_values_cor[system][val]

