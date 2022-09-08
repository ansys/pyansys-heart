# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 15:31:51 2021

@author: cshao
"""

import os


def check_int(s):
    # function take as input a string, and check if the string is an integer
    if s[0] in ("-", "+"):
        return s[1:].isdigit()
    return s.isdigit()


def extract_components_variable_from_cellml_file(cellml_file):
    # take as input the cell_ml txt file and return the component as a dictionnary,
    # The dictionnary associates each component to its variables ([var_name,unit,initial_value,pub_in_out,priv_in_out])
    with open(cellml_file, "r") as f:
        cellml_lines = f.readlines()
    comp_cellml = []
    comp_line = False
    for line in cellml_lines:

        if "def comp" in line:
            comp_line = True
        elif "enddef;" in line:
            comp_line = False
        if len(line.split()) > 0:

            if comp_line and line.split()[0][0:2] != "//":

                comp_cellml.append(line)
    comp_variable = {}

    for i in comp_cellml:
        if "def comp" in i:
            comp_name = i.split()[2]
            comp_variable[comp_name] = []

        else:
            if comp_name == "cell_geometry":
                print(i)

                # unit_attribut.append(i.split()[1])
            if "var" in i:

                var = [i.split()[1].replace(":", "").replace("};", "")]
                if var[0] != "time":
                    var.append(i.split()[2].replace(";", ""))
                    if "init:" in i:

                        if check_int(i.split()[4].replace(",", "").replace("};", "")):

                            var.append(i.split()[4].replace(",", "").replace("};", "") + ".0")
                        else:
                            var.append(i.split()[4].replace(",", "").replace("};", ""))
                    else:
                        var.append("No init")
                    if "pub:" in i:
                        var.append(
                            i[i.index("pub:") :].split()[1].replace(",", "").replace("};", "")
                        )
                    else:
                        var.append("No pub")
                    if "priv:" in i:
                        var.append(
                            i[i.index("priv:") :].split()[1].replace(",", "").replace("};", "")
                        )
                    else:
                        var.append("No priv")
                    if comp_name == "cell_geometry":
                        print(i)

                    comp_variable[comp_name].append(var)

    return comp_variable


def extract_components_equation_from_cellml_file(cellml_file):
    # take as input the cell_ml txt file and return the component as a dictionnary,
    # The dictionnary associates each component to a list where each element contains the lines corresponding to one equation
    with open(cellml_file, "r") as f:
        cellml_lines = f.readlines()
    comp_cellml = []
    comp_line = False
    for line in cellml_lines:

        if "def comp" in line:
            comp_line = True
        elif "enddef;" in line:
            comp_line = False
        if len(line.split()) > 0:
            if comp_line and line.split()[0][0:2] != "//":
                comp_cellml.append(line)

    comp_equation = {}
    sel_line = False
    equation = []

    for i in comp_cellml:

        if "def comp" in i:
            comp_name = i.split()[2]
            comp_equation[comp_name] = []

        else:
            if "= sel" in i:
                sel_line = True
            elif "endsel;" in i:
                sel_line = False
                comp_equation[comp_name].append(equation)
                equation = []

            if sel_line:
                equation.append(i)

            elif "=" in i:
                equation = [i]
                comp_equation[comp_name].append(equation)
                equation = []

    return comp_equation


def extract_unit_from_cellml_file(cellml_file):
    # take as input the cell_ml txt file and return the component as a dictionnary,
    # The dictionnary associates each component to a list where each element contains the lines corresponding to one equation

    with open(cellml_file, "r") as f:
        cellml_lines = f.readlines()
    unit_cellml = []
    unit_line = False
    for line in cellml_lines:

        if "def unit" in line:
            unit_line = True
        elif "enddef;" in line:
            unit_line = False

        if unit_line:
            unit_cellml.append(line)
    unit = {}
    for i in unit_cellml:
        if i == "\n":
            A = 1
        elif "def unit" in i:
            unit_name = i.split()[2].replace(";", "")
            unit[unit_name] = []

        else:
            # unit_attribut.append(i.split()[1])
            if "pref" in i:
                pref = i[i.index("pref") :].split()[1].replace(",", "").replace("};", "")
            else:
                pref = ""
            if "expo" in i:

                expo = i[i.index("expo") :].split()[1].replace("};", "")
            else:
                expo = ""
            unit_attribute = [i.split()[1].replace(";", ""), pref, expo]
            unit[unit_name].append(unit_attribute)
    return unit


def extract_mapping_from_cellml_file(cellml_file):
    # take as input the cellml txt file
    # return a dictionnary specifying for each coupled components the shared variables
    # map_cellml([component1,component2])=[[shared_variable_component1,shared_variable_component2],...]
    with open(cellml_file, "r") as f:
        cellml_lines = f.readlines()
    comp_cellml = []
    comp_line = False
    for line in cellml_lines:

        if "def map" in line:
            comp_line = True
        elif "enddef;" in line:
            comp_line = False
        if len(line.split()) > 0:
            if comp_line and line.split()[0][0:2] != "//":
                comp_cellml.append(line)

    map_cellml = {}

    for i in comp_cellml:

        if "def map" in i:
            comp_name_1 = i.split()[3]
            comp_name_2 = i.split()[5]
            map_cellml[(comp_name_1, comp_name_2)] = []

        else:

            if "time" not in i:
                var_1 = i.split()[1].replace(";", "")
                var_2 = i.split()[3].replace(";", "")
                map_cellml[(comp_name_1, comp_name_2)].append([var_1, var_2])

    return map_cellml


def replace_lower_case_by_capital(chaine):
    # take as input a string, remove each '_' and capitalize the first letter of the string in between the '_'
    # ex: replace_lower_case_by_capital(ion_channel_gate)=ionChannelGate
    liste = chaine.split("_")
    new_chaine = ""
    for i in liste:
        new_chaine += i.capitalize()
    return new_chaine


def convert_symbole_cellml_to_modelica(quantity, pref, expo):
    # convert the unit of cellml into the unit of modelica to have the abrevation of the units instead of the full name

    q = ""
    p = ""
    e = ""
    if quantity == "metre":
        q = "m"
    elif quantity == "coulomb":
        q = "C"
    elif quantity == "joule":
        q = "J"
    elif quantity == "litre":
        q = "l"
    elif quantity == "ampere":
        q = "A"
    elif quantity == "farad":
        q = "F"
    elif quantity == "microF":
        q = "uF"
    elif quantity == "mole":
        q = "Mol"
    elif quantity == "millimolar":
        q = "mMol"
    elif quantity == "second":
        q = "s"
    elif quantity == "volt":
        q = "V"
    elif quantity == "siemens":
        q = "S"
    elif quantity == "micrometre3":
        q = "um3"
    elif quantity == "picoA":
        q = "pA"
    elif quantity == "cm2":
        q = "cm2"
    elif quantity == "kelvin":
        q = "K"
    elif quantity == "liter":
        q = "l"
    elif quantity == "hour":
        q = "h"
    else:
        q = quantity

    if pref == "centi":
        p = "c"
    elif pref == "milli":
        p = "m"
    elif pref == "micro":
        p = "u"
    elif pref == "nano":
        p = "n"
    elif pref == "pico":
        p = "p"
    else:
        p = pref

    if expo == "-1":
        e = "/"
        return e + p + q
    else:
        e = expo
        return p + q + e


def convert_equation_cellml_to_modelica(eq, var_in_out, comp):
    # convert cellml equation to modelica format
    # replace ln() function from cellml
    if "ln(" in eq:
        eq = eq.replace("ln(", "log(")
    if "{" in eq:
        eq_split = eq.replace("{", "}").split("}")
        eq_modelica = ""
        for i in range(len(eq_split)):
            if i % 2 == 0:
                eq_modelica += eq_split[i]
        eq = eq_modelica

    # replace sqr() function from cellml
    while "sqr(" in eq:
        a = eq.index("sqr(") + 4
        compt = 0
        while eq[a] != ")" or compt != 0:

            if eq[a] == "(":
                compt += 1
            if eq[a] == ")":
                compt += -1
            a += 1
        eq = eq.replace(eq[eq.index("sqr(") : a + 1], "(" + eq[eq.index("sqr(") + 4 : a] + ")^2")

    # replace pow() function from cellml
    while "pow(" in eq:
        exp_index = eq.index("pow(") + 4
        compt = 0
        end_index = eq.index("pow(") + 4
        while eq[exp_index] != ",":
            exp_index += 1
        while eq[end_index] != ")" or compt != 0:
            if eq[end_index] == "(":
                compt += 1
            if eq[end_index] == ")":
                compt += -1
            end_index += 1

        eq = eq.replace(
            eq[eq.index("pow(") : end_index + 1],
            "("
            + eq[eq.index("pow(") + 4 : exp_index]
            + ")^("
            + eq[exp_index + 1 : end_index].replace(" ", "")
            + ")",
        )
    # replace ode() from cellml
    while "ode(" in eq:
        exp_index = eq.index("ode(") + 4
        end_index = eq.index("ode(") + 4
        while eq[exp_index] != ",":
            exp_index += 1
        while eq[end_index] != ")":
            end_index += 1
        eq = eq.replace(
            eq[eq.index("ode(") : end_index + 1],
            "der(" + eq[eq.index("ode(") + 4 : exp_index] + ")",
        )
    # replace in and out from private compartemnt
    for j in var_in_out:
        if eq not in comp[len(comp) - len(var_in_out) :] and j in eq:
            eq_split = eq.split(j)
            eq_temp = ""
            for i in range(len(eq_split) - 1):

                if (
                    eq_split[i + 1].isalpha()
                    or eq_split[i + 1][0] == "_"
                    or eq_split[i][-1].isalpha()
                    or eq_split[i][-1] == "_"
                ):
                    eq_temp += eq_split[i] + j
                else:
                    eq_temp += eq_split[i] + j + "_in"
            eq_temp += eq_split[-1]
            eq = eq_temp
    eq_temp = ""

    return eq


##### Variable of the model #####
in_path = "D:/PKPD/Model/Yang_2006/"
cellml_file = "D:/PKPD/Model/Yang_2006/Yang_2006.txt"
library_name = "Yang_2006"
modelica_comp_name = library_name + ".Cell"
# real_input = [['membrane','stim_period']]
real_input = []
os.chdir(in_path)


# extract units, variables, equations and mapping of each cellml components
units = extract_unit_from_cellml_file(cellml_file)
comp_var = extract_components_variable_from_cellml_file(cellml_file)
comp_equation = extract_components_equation_from_cellml_file(cellml_file)
map_cellml = extract_mapping_from_cellml_file(cellml_file)
for i in real_input:
    for j in range(len(comp_var[i[0]])):
        if comp_var[i[0]][j][0] == i[1]:
            comp_var[i[0]][j][2] = "No init"
            comp_var[i[0]][j][3] = "in"
            map_cellml[(i[0], "")] = [[i[1], i[1]]]

#####create folder and other .mo file of the library ####
if not os.path.exists("./" + library_name):
    os.mkdir("./" + library_name)

os.chdir("./" + library_name)
with open("package.mo", "w") as f:
    f.write("within ;\npackage " + library_name + ' "CellsModel extension for modelica."\n')
    f.write(
        'extends Modelica.Icons.Package;\n\nannotation (uses(Modelica(version="3.2.2")), Documentation(info="<html>\n<p>This is the Modelica library group to provide the utilities for using Modelica models in Twin Builder environment.</p>\n</html>"));'
    )
    f.write("\nend " + library_name + ";\n")
with open("package.order", "w") as f:
    f.write("Cell")

if not os.path.exists("./Cell"):
    os.mkdir("./Cell")

os.chdir("./Cell")
with open("package.mo", "w") as f:
    f.write("within " + library_name + ";\npackage Cell\nend " + library_name + ";")
with open("package.order", "w") as f:
    # f.write('Components\nUnits\nInterfaces')
    f.write("Components\nUnits\n" + library_name + "_cell\n")

#######      Creation of the Units file in modelica     ######
if len(units) != 0:
    with open("Units.mo", "w") as f:
        f.write("within " + modelica_comp_name + ";\npackage Units\n")
        for unit in units:

            u = replace_lower_case_by_capital(unit)
            s = ""
            if units[unit] != []:
                for i in units[unit]:
                    s += convert_symbole_cellml_to_modelica(i[0], i[1], i[2])
                if s[0] == "/":
                    s = "1" + s
                if s.count("/") > 1:
                    s = ".".join(s.rsplit("/", 1))
                f.write(
                    "type "
                    + u
                    + ' = Real ( final quantity="'
                    + u
                    + '" , final unit = "'
                    + s
                    + '");\n'
                )
        f.write("end Units;")


#######      Creation of the Component file in modelica     ######

# write the component file
comp_var_in_out = {}
with open("Components.mo", "w") as f:
    f.write("within " + modelica_comp_name + ";\npackage Components\n")
    for comp in comp_var:
        comp_var_in_out[comp] = []
        f.write(
            " model "
            + replace_lower_case_by_capital(comp)
            + ' "'
            + replace_lower_case_by_capital(comp)
            + '"\n'
        )
        compute_var = []
        var_in_out = []
        for var in comp_equation[comp]:
            if var != []:
                if "ode" not in var[0]:
                    compute_var.append(var[0].split()[0])
                else:
                    compute_var.append(var[0].split()[0].replace("ode(", "").replace(",", ""))
        # declare variable
        for var in comp_var[comp]:
            # declare inputs
            if var[3] == "in" or var[4] == "in":
                if var[4] != "out" and var[3] != "out":
                    f.write("    Modelica.Blocks.Interfaces.RealInput " + var[0] + ";\n")
                else:
                    f.write("    Modelica.Blocks.Interfaces.RealInput " + var[0] + "_in;\n")
                    f.write("    Modelica.Blocks.Interfaces.RealOutput " + var[0] + "_out;\n")
                    comp_equation[comp].append(
                        ["        " + var[0] + "_out = " + var[0] + "_in;\n"]
                    )
                    var_in_out.append(var[0])
                    comp_var_in_out[comp].append(var[0])

            # declare parameters and intern variables
            elif var[3] == "No pub" and var[4] == "No priv":
                if var[0] not in compute_var:
                    if var[1] != "dimensionless":
                        f.write("    parameter Real " + var[0] + "=" + var[2] + ";\n")
                        # f.write('    parameter Units.'+replace_lower_case_by_capital(var[1])+' '+var[0]+'='+var[2]+';\n')
                    else:
                        f.write("    parameter Real " + var[0] + "=" + var[2] + ";\n")

                else:
                    if var[1] != "dimensionless":
                        if var[2] != "No init":
                            f.write("    Real " + var[0] + "(start=" + var[2] + ");\n")
                        else:
                            f.write("    Real " + var[0] + ";\n")
                        # f.write('    Units.'+replace_lower_case_by_capital(var[1])+' '+var[0]+';\n')
                    else:
                        if var[2] != "No init":
                            f.write("    Real " + var[0] + "(start=" + var[2] + ");\n")
                        else:
                            f.write("    Real " + var[0] + ";\n")
            # declare outputs
            elif var[3] == "out" and var[0] not in compute_var:

                if var[1] != "dimensionless":
                    f.write("    Modelica.Blocks.Interfaces.RealOutput " + var[0] + ";\n")
                    comp_equation[comp].append(["        " + var[0] + " = " + var[2] + ";\n"])
                    # f.write('    parameter Units.'+replace_lower_case_by_capital(var[1])+' '+var[0]+'='+var[2]+';\n')

                else:
                    f.write("    Modelica.Blocks.Interfaces.RealOutput " + var[0] + ";\n")
                    comp_equation[comp].append(["        " + var[0] + " = " + var[2] + ";\n"])

            elif var[3] == "No pub" and var[4] == "out":
                if var[2] == "No init":
                    f.write("    Modelica.Blocks.Interfaces.RealInput " + var[0] + "_in;\n")
                    f.write("    Modelica.Blocks.Interfaces.RealInput " + var[0] + "_out;\n")
                    comp_equation[comp].append(
                        ["        " + var[0] + "_out = " + var[0] + "_in;\n"]
                    )
                    var_in_out.append(var[0])
                    comp_var_in_out[comp].append(var[0])
                else:
                    f.write("    Modelica.Blocks.Interfaces.RealInput " + var[0] + "_out;\n")
                    comp_equation[comp].append(["        " + var[0] + "_out = " + var[2] + ";\n"])
                    var_in_out.append(var[0])
                    comp_var_in_out[comp].append(var[0])
            else:
                if var[2] == "No init":
                    f.write("    Modelica.Blocks.Interfaces.RealOutput " + var[0] + ";\n")
                else:
                    f.write(
                        "    Modelica.Blocks.Interfaces.RealOutput "
                        + var[0]
                        + "(start="
                        + var[2]
                        + ");\n"
                    )

        # write equation
        f.write("\nequation\n\n")

        for eq in comp_equation[comp]:
            # write equation with no 'if' condition
            if "sel" not in eq[0]:

                eq = convert_equation_cellml_to_modelica(eq[0], var_in_out, comp_equation[comp])

                f.write(eq)
            # write equation with 'if' condition
            else:
                var_sel = eq[0].split()[0]
                condition_sel = convert_equation_cellml_to_modelica(
                    eq[1].replace("case", "").replace(":", "").replace("\n", ""),
                    var_in_out,
                    comp_equation[comp],
                )
                temp = ""
                for k in condition_sel.split():
                    if check_int(k):
                        k = k + ".0"
                    temp += k + " "

                condition_sel = temp
                # case of equation containing '==' condition, the '=' sign in modelica compare only string, for real need to be a double condition with <= and >=
                if "==" in condition_sel:
                    condition_sel = (
                        condition_sel.split()[0]
                        + " >= "
                        + condition_sel.split()[2]
                        + " and "
                        + condition_sel.split()[0]
                        + " <= "
                        + condition_sel.split()[2]
                    )
                equation_1 = convert_equation_cellml_to_modelica(
                    eq[2], var_in_out, comp_equation[comp]
                )
                equation_2 = convert_equation_cellml_to_modelica(
                    eq[-1], var_in_out, comp_equation[comp]
                )
                f.write("	if " + condition_sel + " then \n")
                f.write("		" + var_sel + " = " + equation_1.replace(" ", ""))
                f.write("	else\n")
                f.write("		" + var_sel + " = " + equation_2.replace(" ", ""))
                f.write("	end if;\n")

        f.write("end " + replace_lower_case_by_capital(comp) + ";\n\n")
    f.write("end Components;\n")

# os.chdir('../../')
# write sim model:
comp_var_in_out[""] = []
with open(library_name + "_cell.mo", "w") as f:
    f.write("model " + library_name + "_cell\n")
    # declare the different compartment of the model
    for comp in comp_var:

        f.write(
            "    "
            + modelica_comp_name
            + ".Components."
            + replace_lower_case_by_capital(comp)
            + " "
            + replace_lower_case_by_capital(comp)
            + ";\n"
        )
    for i in real_input:
        f.write("    Modelica.Blocks.Interfaces.RealInput " + i[1] + ";\n")
    f.write("equation\n")
    # write the shared variables between the components
    for mapp in map_cellml:
        for var in map_cellml[mapp]:
            if mapp[1] == "":
                f.write(
                    "   connect("
                    + replace_lower_case_by_capital(mapp[0])
                    + "."
                    + var[0]
                    + ", "
                    + var[1]
                    + ");\n"
                )
            elif var[0] not in comp_var_in_out[mapp[0]] and var[1] not in comp_var_in_out[mapp[1]]:
                f.write(
                    "   connect("
                    + replace_lower_case_by_capital(mapp[0])
                    + "."
                    + var[0]
                    + ", "
                    + replace_lower_case_by_capital(mapp[1])
                    + "."
                    + var[1]
                    + ");\n"
                )
            elif var[0] in comp_var_in_out[mapp[0]]:
                f.write(
                    "   connect("
                    + replace_lower_case_by_capital(mapp[0])
                    + "."
                    + var[0]
                    + "_out, "
                    + replace_lower_case_by_capital(mapp[1])
                    + "."
                    + var[1]
                    + ");\n"
                )
            elif var[1] in comp_var_in_out[mapp[1]]:
                f.write(
                    "   connect("
                    + replace_lower_case_by_capital(mapp[0])
                    + "."
                    + var[0]
                    + ", "
                    + replace_lower_case_by_capital(mapp[1])
                    + "."
                    + var[1]
                    + "_out);\n"
                )

    f.write("end  " + library_name + "_cell;")
