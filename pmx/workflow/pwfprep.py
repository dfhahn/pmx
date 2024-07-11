#!/usr/bin/env python

"""
This script is used to prepare the protein and ligand files, directories, and the 
necessary yaml files of target, ligands and edges for free energy calculations
following the pmx workflow.
It is possible to activate file conversion from a schrodinger project (convert) and 
to create the corresponding directories and yaml files (prep). Required input is 
the target name and the directory in which to create the output directory. Optionally,
you can provide names of properties that are stored for target and ligands to be
extracted from the input files. By default, all properties are extracted.
"""

from argparse import ArgumentParser
import csv
from datetime import date
import os
import shlex
import shutil
import subprocess
import sys

import yaml
from rdkit import Chem

__author__ = "Kenneth Goossens"
__version__ = "1.0.0"
__email__ = "goossens_kenny@hotmail.com"

def convert_schrod(target_name, directory, flag2):
    if not flag2:
        return
    print("file conversion function invoked.")
    schrodinger_path=(os.environ['SCHRODINGER'])
    command1 = shlex.split(f'{schrodinger_path}/run -FROM scisol extract_pv.py\
        -include-pred-dg {target_name}.fmp {target_name}.maegz')
    command2 = shlex.split(f'{schrodinger_path}/utilities/sdconvert\
        -imae {target_name}.maegz -osd {target_name}.sdf')
    command3 = shlex.split(f'{schrodinger_path}/run -FROM scisol\
        fmp2excel.py {target_name}.fmp -csv -o {target_name}')
    command4 = shlex.split(f'{schrodinger_path}/utilities/pdbconvert\
        -imae {target_name}.maegz -opdb {target_name}.pdb -n 1')
    for command in [command1, command2, command3, command4]:
        process = subprocess.Popen(command, stdout=subprocess.PIPE).stdout.read()
    return

def ligand_files_setup(target_name, directory, flag3):
    if not flag3:
        return
    print("ligand file setup invoked.")
    dict_file = os.path.join(directory, "FEP_workflow_pmx/targets.yml")
    with open(dict_file, "r") as f:
        target_dict = yaml.safe_load(f)
    target_dir = target_dict[target_name]["dir"]
    path_to_ligands = os.path.join(directory, "FEP_workflow_pmx", target_dir, "02_ligands/")
    mols = Chem.SDMolSupplier(f"{target_name}.sdf")
    for count, mol in enumerate(mols):
        if count == 0:
            continue 
        name = mol.GetProp("_Name")
        writer = Chem.SDWriter(f"{name}.sdf")
        writer.write(mol)
        mol_path = os.path.join(path_to_ligands, name, "crd")
        if os.path.isfile(mol_path):
            print(mol_path, "is a file")
        os.makedirs(mol_path, exist_ok=True)
        os.replace(f"{name}.sdf", os.path.join(mol_path, f"{name}.sdf"))
    return


def main(target_name, directory, flag1, *props):
    """Optional args are the names of the SDF file headers that you would like to extract.
    Extracts all headers by default. Directory is the absolute path to the directory in 
    which the input files are located.
    """

    """Directory creation"""

    if not flag1:
        return
    print("main function invoked.\n")
    print(f"Creating directories in '{dir_path}'.")
    for directory in ["00_data", "01_protein", "02_ligands", "03_hybrid"]:
        os.makedirs(os.path.join(path_target, directory), exist_ok=True)
    os.makedirs(os.path.join(path_target, "01_protein/crd"), exist_ok=True)

    """Append target to targets.yml"""

    target_dict = {target_name:{"date":str(today),
                                "dir":target_name+"_"+today,
                                "name":target_name}}
    tgts_yml_path = os.path.join(output_dir, "targets.yml")
    if not os.path.exists(tgts_yml_path):
        print("NOTE: targets.yml did not exist yet. New file created.")
        with open(tgts_yml_path, 'w') as file:
            pass        
    with open(tgts_yml_path, "r+") as file:
        check = file.read()
        if target_name + ":" in check:
            print(f"{target_name} already exists in the dataset. No new instance was created.")
        else:
            print("writing...")
            file.write(yaml.dump(target_dict)) #prints date with quotation marks

    """Edges.yml"""

    data_path = os.path.join(path_target, "00_data")
    edge_dict = {}
    with open(t_edges, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            new_edge = {"edge_" + row["Ligand1"] + "_" + row["Ligand2"]:\
                        {"ligand_a":row["Ligand1"],\
                        "ligand_b":row["Ligand2"]}}
            edge_dict.update(new_edge)
    with open(os.path.join(data_path, "edges_tmp.yml") ,"w+") as file: 
        file.write(yaml.dump(edge_dict))
#removing quotation marks from ligand names
    with open(os.path.join(data_path, "edges_tmp.yml"), "r") as file,\
        open(os.path.join(data_path, "edges.yml"), "w") as file_out:
        for line in file:
            file_out.write(line.replace('\'', ''))
    os.remove(os.path.join(data_path, "edges_tmp.yml"))

    """Target.yml"""

    tgt_yaml_path = os.path.join(data_path, "target.yml")
    tgt_dict = {"name":target_name, "date":today}
    dict_name = None
    dict_value = None
    if not os.path.exists(tgt_yaml_path):
        with open(t_structures, "r") as file:
            linecount = 0
            for line in file:
                if "$$$$" in line:
                    break
                if "> <" in line:
                    dict_name = line.strip("<> \n").replace('>','')
                    linecount += 1
                    continue                
                elif linecount == 1:
                    dict_value = line.strip('\n')
                    linecount = 0
                if not args.properties_target == None:
                    if dict_name in args.properties_target:        
                        tgt_dict[dict_name] = dict_value
                else:
                    if not dict_name == None:
                        tgt_dict[dict_name] = dict_value
            with open(tgt_yaml_path, 'w') as file:
                file.write(yaml.dump(tgt_dict))
    else:
        print("WARNING: target.yml file already exists. Proceeding with workflow,"\
              " but no new target.yml file was created.")

    """Ligands"""

    lig_yaml_path = os.path.join(data_path,"ligands.yml")
    if not os.path.exists(lig_yaml_path):
        ref_dict = []
        lig_dict = {}
        ligcount = 0
        linecount = 0
        with open(t_ligands, newline = '') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                ref_dict.append(row["Ligand"])
        with open(t_structures, "r") as file:
            for line in file:
                if "$$$$" in line:
                    if ligcount != 0:
                        lig_dict.update(ligand)
                    ligcount += 1
                    try:
                        ligand = {ref_dict[ligcount-1]:{'name':ref_dict[ligcount-1]}}
                    except IndexError:
                        break
                    linecount = 0    
                if "> <" in line and ligcount != 0:
                    linecount = 0
                    dict_name = line.strip("<> \n").replace('>','').replace('[V] ','')
                    linecount += 1
                    continue                
                elif linecount == 1:
                    dict_value = line.strip('\n')
                    linecount = 0
                    if not args.properties_ligands == None:
                        if dict_name in args.properties_ligands:        
                            ligand[ref_dict[ligcount - 1]][dict_name] = dict_value
                    else:
                        ligand[ref_dict[ligcount - 1]][dict_name] = dict_value
 
        with open(os.path.join(data_path,"ligands_tmp.yml"), 'w+') as file:
            file.write(yaml.dump(lig_dict))
        with open(os.path.join(data_path,"ligands_tmp.yml"),"r") as file,\
            open(os.path.join(lig_yaml_path),"w") as file_out:
            for line in file:
                file_out.write(line.replace('\'', ''))
        os.remove(os.path.join(data_path,"ligands_tmp.yml"))
        print('Metadata of',len(lig_dict), 'ligands written to ligands.yml')
    else:
        print("WARNING: ligands.yml file already exists. Proceeding with workflow,"\
              " but no new ligands.yml file was created.")
    return

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("-t", "--target",
                        dest="target",
                        type=str,
                        help="The name of the TARGET without extension.",
                        metavar="TARGET")
    parser.add_argument("-d", "--directory",
                        dest="directory",
                        type=str,
                        help="The full path to the DIRECTORY with input files.",
                        metavar="DIRECTORY")
    parser.add_argument("-pt", "--properties_target",
                        dest="properties_target",
                        type=str,
                        nargs="+",
                        help="Properties to extract from the target SDF file.",
                        metavar="TARGET_PROPERTIES")
    parser.add_argument("-pl", "--properties_ligands",
                        dest="properties_ligands",
                        type=str,
                        nargs="+",
                        help="Properties to extract from the ligands in the SDF file.",
                        metavar="LIGAND_PROPERTIES")
    parser.add_argument("-prep",
                        dest="flag1",
                        default=False,
                        action='store_true',
                        help="If true, folder preparation function will run.",)
    parser.add_argument("-convert",
                        dest="flag2",
                        default=False,
                        action='store_true',
                        help="If true, schrodinger conversion function will run.")
    parser.add_argument("-lig_setup",
                        dest="flag3",
                        default=False,
                        action='store_true',
                        help=("If true, function to extract and move ligands "
                              "from sdf file will run."))

    args = parser.parse_args()
    
    dt = date.today()
    today = dt.strftime("%Y-%m-%d")

    dir_path = args.directory
    output_dir = os.path.join(dir_path, "FEP_workflow_pmx")
    path_target = os.path.join(output_dir, args.target + "_" + today) 

    t_structures = os.path.join(dir_path, args.target + ".sdf")
    t_edges = os.path.join(dir_path, args.target + "_ddG.csv")
    t_ligands = os.path.join(dir_path, args.target + "_dG.csv")

    convert_schrod(args.target, args.directory, args.flag2)
    main(args.target, args.directory, args.flag1, args.properties_target, args.properties_ligands)
    ligand_files_setup(args.target, args.directory, args.flag3)
