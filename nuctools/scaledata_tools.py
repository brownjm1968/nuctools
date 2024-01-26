
import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET
import h5py as h5
import os,sys,json
import matplotlib.pyplot as plt
import pkgutil

from . import chem_tools as ct
from . import physics_tools as pt

__all__ = ['plot_h5scale_xs',"get_cross_section","get_std_comp","get_std_comp_nat_abund","get_zlist",
           "get_zaidlist","calc_num_densities","write_tsl_table","get_scaleza_name_thermal",
           "get_scaleza_name_metastable","get_scaleza_name_specialNuclei","append_xml_to_table",
           "get_xml_root","get_fastmat_thermal","get_single_mat"]

def plot_h5scale_xs(filename,scaleid,temp,emin=2.1e7,mt=None):
    """
    Plot xs of every MT (reaction number) in SCALE data H5 format

    Parameters
    ----------
    filename : str
        The path to the H5 library file
    scaleid : int
        The integer representing a single isotope or compound, e.g. 1001 
        for hydrogen, 5008016 for the oxygen in BeO 
    temp : float
        The temperature, assumed to only need one digit after the decimal
    emin : float, optional
        Minimum energy of XS must be below this value to be plotted. 
    mt : int, optional
        The MT ENDF/SCALE reaction number to be plotted. Default is all within 
        energy limit defined

    Returns
    -------
    None

    Examples
    --------
    >>> import nuctools as nuc
    >>> 
    >>> scaleid = 5008016
    >>> temperature = 293.6
    >>> nuc.plot_all_h5scale_xs('n_008016.h5',temperature,scaleid)
    """
    with h5.File(filename,'r+') as f:
        plt.figure()
        if mt is None:
            for key in f['n_{:0>7d}_{:0>6.1f}'.format(scaleid,temp)]:
                try:
                    e  = f['n_{:0>7d}_{:0>6.1f}/{}/energy'.format(scaleid,temp,key)][()]
                    xs = f['n_{:0>7d}_{:0>6.1f}/{}/xs'.format(scaleid,temp,key)][()]
                    if e.min() < emin:
                        plt.plot(e,xs,label='{}'.format(key))
                except:
                    print("key = '{}' cannot be plotted".format(key))
        else:
            try:
                e  = f['n_{:0>7d}_{:0>6.1f}/mt_{:0>4d}/energy'.format(scaleid,temp,mt)][()]
                xs = f['n_{:0>7d}_{:0>6.1f}/mt_{:0>4d}/xs'.format(scaleid,temp,mt)][()]
                if e.min() < emin:
                    plt.plot(e,xs,label='{}'.format("mt={}".format(mt)))
            except:
                print("key = 'mt_{:0>4d}' cannot be plotted".format(mt))
        plt.xlabel('Energy [eV]')
        plt.ylabel("Cross section [b]")
        plt.xscale('log')
        plt.yscale('log')
        plt.legend(ncol=3)

def get_cross_section(filename,scaleid,temp,mt):
    """
    Plot xs of every MT (reaction number) in SCALE data H5 format

    Parameters
    ----------
    filename : str
        The path to the H5 library file
    scaleid : int
        The integer representing a single isotope or compound, e.g. 1001 
        for hydrogen, 5008016 for the oxygen in BeO 
    temp : float
        The temperature, assumed to only need one digit after the decimal
    emin : float, optional
        Minimum energy of XS must be below this value to be plotted. 
    mt : int, optional
        The MT ENDF/SCALE reaction number to be plotted. Default is all within 
        energy limit defined

    Returns
    -------
    data : DataFrame
        A Pandas DataFrame with keys: "e" and "cs" for energy and cross section

    Examples
    --------
    >>> import numpy as np
    >>> 
    >>> scaleid = 5008016
    >>> temperature = 293.6
    >>> data = nuc.get_cross_section('n_008016.h5',temperature,scaleid)

    """
    with h5.File(filename,'r') as f:
        try:
            data = pd.DataFrame({
                "e" : f['n_{:0>7d}_{:0>6.1f}/mt_{:0>4d}/energy'.format(scaleid,temp,mt)][()],
                "cs" : f['n_{:0>7d}_{:0>6.1f}/mt_{:0>4d}/xs'.format(scaleid,temp,mt)][()]
            })
        except:
            print("key = 'n_{:0>7d}_{:0>6.1f}/mt_{:0>4d}' cannot be found".format(scaleid,temp,mt))
            df = pd.DataFrame(f['nuclide_md'][:])
            print("Temperatures available for SCALEID = {}:".format(scaleid))
            print(df[df['zaid']==scaleid]['temperature'])
            raise ValueError("Temperature was not found.")
    return data



def get_std_comp(datadir):
    """
    Return a std composition table pandas DataFrame

    assumes the directory and file names are static in the 
    SCALE data directory

    Parameters
    ----------
    datadir : str
        The directory where the "scale data" exists. This is a 
        specific Git repo

    Returns
    -------
    stdcomp : DataFrame
        A table of the standard composition from SCALE data
    """

    std_composition = datadir+'library_inputs/compoz.inp'
    begin = 0
    max_rows = 3527 # num rows in table not including mixtures
    with open(std_composition,'r+') as f:
        for i,line in enumerate(f):
            if "SCID          ROTH.    ICP      NCZA     ATPM     END" in line:
                begin_natural = i+1
    stdcomp = pd.read_fwf(std_composition,skiprows=begin_natural,nrows=max_rows,
                          colspecs=( (0,14), (15,22), (23,27), (29,39), (40,48), (49,55)),
                          names=['name','rho_theor','compound_flag','scaleid','wtpct','end'])
    stdcomp.fillna(0,inplace=True)
    stdcomp.scaleid = stdcomp.scaleid.astype('int')

    return stdcomp

def get_std_comp_nat_abund(datadir):
    """
    Return a std composition table pandas DataFrame for the 
    natural abundance table

    assumes the directory and file names are static in the 
    SCALE data directory

    Parameters
    ----------
    datadir : str
        The directory where the "scale data" exists. This is a 
        specific Git repo

    Returns
    -------
    stdcomp : DataFrame
        A table of the standard composition of natural abundance
        from SCALE data
    """

    std_composition = datadir+'library_inputs/compoz.inp'
    begin_natural = 0
    max_rows = 361 # num rows with natural abundance
    with open(std_composition,'r+') as f:
        for i,line in enumerate(f):
            if "NZN         ISZA   ABWP      END" in line:
                begin_natural = i+1
    stdcomp = pd.read_fwf(std_composition,skiprows=begin_natural,nrows=max_rows,names=['natid','isoid','frac','end'])
    stdcomp.fillna(0,inplace=True)
    stdcomp.isoid = stdcomp.isoid.astype('int')

    return stdcomp

def get_zlist(symbol_list):
    """
    Get list of proton numbers matching the list of symbols given
    in the input

    Parameters
    ----------
    symbol_list : list
        A list of strings with elemental symbols (properly caps), 
        e.g. Am, H, Ca, etc.
    
    Returns
    -------
    zlist : list
        A list of ints for the number of protons in each element
    """
    zlist = []
    f_protons = os.path.join(os.path.dirname(__file__),'data','protons.json')
    with open(f_protons,'r') as f:
        protons = json.load(f)
    for symbol in symbol_list:
        zlist.append(protons[symbol])
    return zlist

def get_zaidlist(zlist,stdcomp):
    """
    Get the list of ZAIDs, relevant to ENDF/SCALE ID system, using the
    standard composition to get the natural composition of isotopes

    Parameters
    ----------
    zlist : list
        List of atomic numbers for unique elements in the salt mixture
    stdcomp : DataFrame
        DataFrame object made exactly like the output of `get_std_comp_nat_abund`

    Returns
    -------
    zaidlist : list
        The list of ZAIDs for the mixture of salts
    """
    zaidlist = []
    natid = 0.0
    for z in zlist:
        for i,za in enumerate(stdcomp.isoid):
            # determine if TSL or metastable (skipping metas!)
            if stdcomp.natid.iloc[i] != 0:
                natid = stdcomp.natid.iloc[i]
            if( z == za//1000 and natid < 1000000):
                zaidlist.append(za)
    return zaidlist

def calc_num_densities(mol_df,mass,area,scaledatadir):
    """
    Calculate number densities for every SCALE ZAID based on 
    mass fractions of molecules input, mass of sample, and 
    area of sample

    Parameters
    ----------
    mol_df : DataFrame
        A Pandas DataFrame with column names: "molecule", and "mass_frac"
        which describe the molecules (with proper capitalization) in the
        sample and the mass fractions of each molecule. 
    mass : float
        The mass of the sample in grams
    area : float
        The area of the sample perpendicular to the beam (neutrons) in cm
    scaledatadir : str
        The path to the SCALE data directory. (This dependency could be 
        removed by including a natural abundance file in nuctools in 
        future releases)

    Returns
    -------
    numdensity_iso : dict
        A dictionary with keys of ZAIDs and values of number density for 
        every ZAID
    
    Example
    -------
    >>> import nuctools as nuc
    >>> import pandas as pd
    >>> mol = pd.DataFrame({"CaO":0.5,"Fe2O3":0.5})
    >>> mass = 40 # grams
    >>> area = 4.5 # cm
    >>> n = nuc.calc_num_densities(mol,mass,area)
    """
    f_protons = os.path.join(os.path.dirname(__file__),'data','protons.json')
    f_unary_list = os.path.join(os.path.dirname(__file__),'data','unary-list.json')
    with open(f_unary_list,'r') as f:
        unary_dict = json.load(f)
        unary_stoich = unary_dict['unary_stoich']
        unary_MM = unary_dict['unary_mol_mass']
    with open(f_protons,'r') as f:
        protons = json.load(f)
    stdcomp = get_std_comp_nat_abund(scaledatadir)

    # get unique elements in sample
    full_elem_list = ct.get_single_elem_list(mol_df.molecule)
    # get proton numbers for all elements in sample
    zlist = [protons[el] for el in full_elem_list]
    # get unique naturally abundant ZAIDs for the elements in sample
    zaidlist = get_zaidlist(zlist,stdcomp)
    # initialize number density dict with zeros
    numdensity_iso = {}
    for zaid in zaidlist:
        numdensity_iso[zaid] = 0.0

    # --- molecule loop
    for i,molec in enumerate(mol_df.molecule):
        MF_mol = mol_df.mass_frac.iloc[i]
        MM_mol = unary_MM[molec]
        elem_list = ct.get_single_elem_list([molec])
        molec_stoich = unary_stoich[molec]
        # --- elements in molecule loop
        for j,elem in enumerate(elem_list):
            z = protons[elem]
            s_elem = molec_stoich[j]
            elem_zaidlist = get_zaidlist([z],stdcomp)
            # --- isotopes in element loop
            for zaid in elem_zaidlist:
                gamma_iso = 0.01 * stdcomp.loc[stdcomp.isoid==zaid].frac.iloc[0]
                numdensity_iso[zaid] += gamma_iso * s_elem * MF_mol * mass * pt.Na / area / MM_mol
    return numdensity_iso

def write_tsl_table(library_master_file,output_file,stdcomp):
    """
    Read library master file to find TSLs and temperatures
    available and write to file

    Parameters
    ----------
    library_master_file : str
        The path for the library master file in HDF5
    output_file : str
        The path where to write a markdown table
    stdcomp : DataFrame
        The SCALE standard composition information which 
        can be placed in a DataFrame with `get_std_comp`

    Returns
    -------
    None


    """
    tsl_table = [['scaleid','name','temps','**notes**']]
    with h5.File(library_master_file,'r') as f:
        last_scaleid = None
        last_scaleid_pos = 0
        for i,nuclide in enumerate(f):
            if nuclide == 'metadata' or nuclide == "nuclide_md" or "p_" in nuclide:
                continue
            found_tsl = False
            scaleid = nuclide.split('_')[1]
            temp = nuclide.split('_')[2]
            is_water = (int(scaleid) == 1001) or (int(scaleid) == 1002)
            if (int(scaleid) < 100*1000) and not is_water:
                continue
            for mtlist in f[nuclide]['reaction_md']:
                mt = mtlist[0]
                if mt == 1007 or mt == 1008:
                    found_tsl = True
                    break
            if found_tsl:
                if last_scaleid == scaleid:
                    temp = str(float(temp))
                    tsl_table[last_scaleid_pos][2] += ", {}".format(temp)
                else:
                    name = stdcomp.name[stdcomp.scaleid==int(scaleid)].to_numpy()[0]
                    temp = str(float(temp))
                    tsl_table.append([str(int(scaleid)),name,temp,''])
                    last_scaleid_pos += 1
                last_scaleid = scaleid
    print("Number of TSLs found = ",last_scaleid_pos)
    tsl_df = pd.DataFrame(tsl_table)
    header = tsl_df.iloc[0]
    tsl_df = tsl_df.iloc[1:len(tsl_df)]
    tsl_df.columns = header
    tsl_df.to_markdown(output_file,tablefmt="grid",index=False,maxcolwidths=[None,None, 45,None])

def get_single_mat(root,endfmat):
    """
    Get the XML root object using python xml

    Parameters
    ----------
    root : Element
        The path to the XML file
    endfmat : str
        The ENDF MAT number identifying the element we want

    Returns
    -------
    child : Element
        A child Element object from the root Element
    """
    for child in root:
        endf = child.attrib['endf']
        if endf == endfmat:
            return child

def get_fastmat_thermal(root,realza,endfmat):
    """
    Read XML config file from SCALE data repo and get the SCALEID
    and ENDF MAT number for nuclei listed under "thermal" to find
    the MAT number for fast evaluation

    Parameters
    ----------
    root : Element object
        Element object from an ElementTree object from the python xml 
        package
    realza : str
        string of an int for the ZA of a material in ENDF format
    endfmat : str
        string of int for the ENDF MAT number

    Returns
    -------
    endf : str
        The MAT number for the fast evaluation
    """
    for child in root:
        # print(child.tag)
        if child.tag == 'thermal':
            for grandchild in child:
                # print(grandchild.keys())
                temp_realza = grandchild.get('realza')
                endf = grandchild.get('endf')
                # print(temp_realza,endf)
                if temp_realza == realza and endf == endfmat:
                    for greatgrandchild in grandchild:
                        endf = greatgrandchild.get('endf')
                        return endf

def get_xml_root(filename):
    """
    Get the XML root object using python xml

    Parameters
    ----------
    filename : str
        The path to the XML file

    Returns
    -------
    root : Element
        The root Element object from ElementTree
    """
    tree = ET.parse(filename)
    root = tree.getroot()
    return root

def get_scaleza_name_specialNuclei(root,realza,endfmat):
    """
    Read XML config file from SCALE data repo and get the SCALEID
    and ENDF MAT number for nuclei listed under "specialNuclei"

    Parameters
    ----------
    root : Element object
        Element object from an ElementTree object from the python xml 
        package
    realza : str
        string of an int for the ZA of a material in ENDF format
    endfmat : str
        string of int for the ENDF MAT number

    Returns
    -------
    scaleza, name : tuple
        A tuple of strings, the SCALEID and the name for SCALE

    Examples
    --------
    >>> import nuctools as nuc
    >>> e80root = nuc.get_xml_root('ENDF-8.0.xml_config')
    >>> nuc.get_scaleza_name_specialNuclei(e80root,'1001','125')
    """
    for child in root:
        #print(child.tag)
        if child.tag == 'specialNuclei':
            for grandchild in child:
                #print(grandchild.keys())
                temp_realza = grandchild.get('realza')
                endf = grandchild.get('endf')
                #print(temp_realza,endf)
                if temp_realza == realza and endf == endfmat:
                    scaleza = grandchild.get('scaleza')
                    name = grandchild.get('name')
                    return scaleza,name
def get_scaleza_name_thermal(root,realza,endfmat):
    """
    Read XML config file from SCALE data repo and get the SCALEID
    and ENDF MAT number for nuclei listed under "thermal"

    Parameters
    ----------
    root : Element object
        Element object from an ElementTree object from the python xml 
        package
    realza : str
        string of an int for the ZA of a material in ENDF format
    endfmat : str
        string of int for the ENDF MAT number

    Returns
    -------
    scaleza, name : tuple
        A tuple of strings, the SCALEID and the name for SCALE

    Examples
    --------
    >>> import nuctools as nuc
    >>> e80root = nuc.get_xml_root('ENDF-8.0.xml_config')
    >>> nuc.get_scaleza_name_thermal(e80root,'126','425')
    """
    names,scalezas = [],[]
    for child in root:
        # print(child.tag)
        if child.tag == 'thermal':
            for grandchild in child:
                # print(grandchild.keys())
                temp_realza = grandchild.get('realza')
                endf = grandchild.get('endf')
                # print(temp_realza,endf)
                if temp_realza == realza and endf == endfmat:
                    for greatgrandchild in grandchild:
                        scaleza = greatgrandchild.get('scaleza')
                        name = greatgrandchild.get('name')
                        scalezas.append(scaleza)
                        names.append(name)
                    return scalezas,names
def get_scaleza_name_metastable(root,realza,endfmat):
    """
    Read XML config file from SCALE data repo and get the SCALEID
    and ENDF MAT number for nuclei listed under "metastable"

    Parameters
    ----------
    root : Element object
        Element object from an ElementTree object from the python xml 
        package
    realza : str
        string of an int for the ZA of a material in ENDF format
    endfmat : str
        string of int for the ENDF MAT number

    Returns
    -------
    scaleza, name : tuple
        A tuple of strings, the SCALEID and the name for SCALE

    Examples
    --------
    >>> import nuctools as nuc
    >>> e80root = nuc.get_xml_root('ENDF-8.0.xml_config')
    >>> nuc.get_scaleza_name_metastable(e80root,'61148','6153')
    """
    for child in root:
        # print(child.tag)
        if child.tag == 'metastable':
            for grandchild in child:
                temp_realza = grandchild.get('realza')
                endf = grandchild.get('endf')
                # print(temp_realza,endf)
                if temp_realza == realza and endf == endfmat:
                    scaleza = grandchild.get('scaleza')
                    return scaleza
def append_xml_to_table(table,prot_dict,root,configroot,elibrary,additional_lib=False):
    """
    Append row of values of interest for table based on XML files

    Parameters
    ----------
    table : DataFrame
        DataFrame with members: "scaleid","name","gamprod","bondfac","gamxs","elibrary"
    prot_dict : dict
        python dictionary associating chemical symbols (e.g. Li) with number of protons
        (e.g. 3)
    root : Element
        root of main XML file listing for SCALE data, make with function `get_xml_root` above
    configroot : Element
        root of XML config file, make with function `get_xml_root` above
    elibrary : str
        Which endf library name you want to put in the table
    additional_lib : bool, optional
        is this the first lib appending to table
    
    Returns
    -------
    None : None
        It modifies the table in the input parameter

    """
    for child in root:
        scaleids,names=None,None # for thermal
        za = child.attrib['za']
        meta = child.attrib['metaStable']
        endf = child.attrib['endf']
        awi = child.attrib['awi']
        if za=='1':
            continue # skip neutron eval
        if awi=='0.0':
            continue # skip photo-atomic
        
        # gamma production
        try:
            awp0 = child.attrib['AWP0']
        except:
            awp0 = "no"
        # is TSL
        try:
            file7 = child.attrib['file7']
        except:
            file7 = "no"
        
        tag = child.attrib['tag']
        if file7 == "no":
            name = tag
            scaleid = za
        else:
            scaleids,names = get_scaleza_name_thermal(configroot,za,endf)

        # is metastable
        if meta == "true" and file7 == "no":
            scaleid = get_scaleza_name_metastable(configroot,za,endf)
        
        # special nuclei
        if za == '1001' and endf == '125':
            scaleid,name = get_scaleza_name_specialNuclei(configroot,za,endf)
        if za == '1002' and endf == '128':
            scaleid,name = get_scaleza_name_specialNuclei(configroot,za,endf)
    
        # ---- if TSL
        if file7 == "yes":
            for i in range(len(scaleids)):
                scaleid = scaleids[i]
                name = names[i]
                fastendf = get_fastmat_thermal(configroot,za,endf)
                fastelem = get_single_mat(root,fastendf)
                prot = int(fastelem.attrib['za'])//1000
                gamxs = prot_dict[str(prot)].lower()
                gamprod = "yes"
                bondfac = "yes"
                if additional_lib:
                    if scaleid in table['scaleid'].unique():
                        table.loc[table.scaleid==scaleid,'elibrary'] += ", " + elibrary
                    else:
                        # graphite's can still be called using same name
                        if '6000' in scaleid: 
                            elib = elibrary + "*"
                            print(elib)
                        else:
                            elib = elibrary
                        table.loc[len(table)] = [scaleid,name,gamprod,bondfac,gamxs,elib]
                else:        
                    table.loc[len(table)] = [scaleid,name,gamprod,bondfac,gamxs,elibrary]
        # ---- if *not* TSL
        else:
            prot = int(za)//1000
            gamxs = prot_dict[str(prot)].lower()
            gamprod = awp0
            bondfac = "yes"
            if additional_lib:
                if scaleid in table['scaleid'].unique():
                    table.loc[table.scaleid==scaleid,'elibrary'] += ", " + elibrary
                else:
                    # graphite's can still be called using same name
                    if '6000' in scaleid: 
                        elib = elibrary + "*"
                        print(elib)
                    else:
                        elib = elibrary
                    table.loc[len(table)] = [scaleid,name,gamprod,bondfac,gamxs,elib]
            else:        
                table.loc[len(table)] = [scaleid,name,gamprod,bondfac,gamxs,elibrary]







