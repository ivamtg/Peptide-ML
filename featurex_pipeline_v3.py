from pydpi.pydrug import Chem
from pydpi.drug import *
from pydpi.protein import PseudoAAC

#template: {"A":,"C":,"D":,"E":,"F":,"G":,"H":,"I":,"K":,"L":,"M":,"N":,"P":,"Q":,"R":,"S":,"T":,"V":,"W":,"Y":}

_Hydrophobicity={"A":0.62,"R":-2.53,"N":-0.78,"D":-0.90,"C":0.29,"Q":-0.85,"E":-0.74,"G":0.48,"H":-0.40,"I":1.38,"L":1.06,"K":-1.50,"M":0.64,"F":1.19,"P":0.12,"S":-0.18,"T":-0.05,"W":0.81,"Y":0.26,"V":1.08}

_hydrophilicity={"A":-0.5,"R":3.0,"N":0.2,"D":3.0,"C":-1.0,"Q":0.2,"E":3.0,"G":0.0,"H":-0.5,"I":-1.8,"L":-1.8,"K":3.0,"M":-1.3,"F":-2.5,"P":0.0,"S":0.3,"T":-0.4,"W":-3.4,"Y":-2.3,"V":-1.5}

_residuemass={"A":15.0,"R":101.0,"N":58.0,"D":59.0,"C":47.0,"Q":72.0,"E":73.0,"G":1.000,"H":82.0,"I":57.0,"L":57.0,"K":73.0,"M":75.0,"F":91.0,"P":42.0,"S":31.0,"T":45.0,"W":130.0,"Y":107.0,"V":43.0}

_pK1={"A":2.35,"C":1.71,"D":1.88,"E":2.19,"F":2.58,"G":2.34,"H":1.78,"I":2.32,"K":2.20,"L":2.36,"M":2.28,"N":2.18,"P":1.99,"Q":2.17,"R":2.18,"S":2.21,"T":2.15,"V":2.29,"W":2.38,"Y":2.20}

_pK2={"A":9.87,"C":10.78,"D":9.60,"E":9.67,"F":9.24,"G":9.60,"H":8.97,"I":9.76,"K":8.90,"L":9.60,"M":9.21,"N":9.09,"P":10.6,"Q":9.13,"R":9.09,"S":9.15,"T":9.12,"V":9.74,"W":9.39,"Y":9.11}

_pI={"A":6.11,"C":5.02,"D":2.98,"E":3.08,"F":5.91,"G":6.06,"H":7.64,"I":6.04,"K":9.47,"L":6.04,"M":5.74,"N":10.76,"P":6.30,"Q":5.65,"R":10.76,"S":5.68,"T":5.60,"V":6.02,"W":5.88,"Y":5.63}

"NEW FEATURES BELOW"
#hydrogen bonds: 0 = does not participate in H-bonds, 1=hydrogen-bond donor, -1=hydrogen-bond acceptor, 2=can be either donor or acceptor
_hydrogenbonds= {"A":0,"C":2,"D":0,"E":0,"F":0,"G":0,"H":2,"I":0,"K":0,"L":0,"M":0,"N":-1,"P":0,"Q":-1,"R":0,"S":1,"T":1,"V":0,"W":0,"Y":1}

#disulfide bond only formed with C, but may be worth looking into?
#A: would need to come late since we can't quantify that rn


#aromatic: 0 = not aromatic, 1-3 organized from least to most hydrophobic
_aromaticity={"A":0,"C":0,"D":0,"E":0,"F":3,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":2,"Y":1}

#van der Waals volume in A^3: space enclosed by the van der Waals spheres of the constituent atoms
#related to the packing density of the interior of a macromolecule:
    #ratio of vdw vol to vol actually occupied when in solvent
_vanderwaalsvolume = {"A":67,"C":86,"D":91,"E":109,"F":135,"G":48,"H":118,"I":124,"K":135,"L":124,"M":124,"N":96,"P":90,"Q":114,"R":148,"S":90,"T":93,"V":105,"W":163,"Y":141}

#pi-pi interactions

#"transmembrane tendency" https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2242586/
#check out above link w/ raul and cesar

#KAPPA DESCRIPTORS VERY GENERAL, bc it doesn't take into account composition of protein

#helical penalty kJ/mol (important bc when in an hydrophobic enviro, are generally alpha helix)
_helicalpenalty = {"A":0,"C":2.85,"D":2.89,"E":1.67,"F":2.26,"G":4.18,"H":2.55,"I":1.72,"K":1.09,"L":0.88,"M":1.00,"N":2.72,"P":13.22,"Q":1.63,"R":0.88,"S":2.09,"T":2.76,"V":2.55,"W":2.05,"Y":2.22}

#how to indentify sequential motifs that correlate w/ fxn

#polar angle

#protein concentration and impact on activity (can't take into account, but interesting and look into)

#salt composition of solution (ionic force of the media)

#Zimm Bragg Model (these values determined experimentally, but unknown why they are the way they are)
#s: propagation paramter for the equilib constant for transforming a coil residue to a helical residue
    #at the end of a helical sequence

#sigma*s: equilib constant for initiation a helical unit in a sequence of coil residues

##sigma: the nucleation parameter

def extract_named_descriptors_of_seq(sequence):
    '''
    Returns a map ("descriptor" -> value) of descriptors when given a sequence of aminoacids (string)
    :param sequence:
    :return:
    '''
    #mol = Chem.MolFromSequence(str(sequence))
    res = {}
    sequence=str(sequence)
    if len(sequence) > 3:
        print(sequence)
        #print(Autocorrelation.CalculateGearyAutoTotal(sequence))
        res = PseudoAAC.GetPseudoAAC(sequence, lamda=3, weight=0.05, AAP=[_Hydrophobicity, _hydrophilicity, _residuemass, _pK1, _pK2, _pI])
        mol = Chem.MolFromSequence(str(sequence))
        #res = geary.GetGearyAuto(mol)
        res.update(kappa.GetKappa(mol))
        res.update(charge.GetCharge(mol))
        #res.update(moran.GetMoranAuto(mol))
        #res.update(moreaubroto.GetMoreauBrotoAuto(mol))
        res.update(molproperty.GetMolecularProperty(mol))
        #res.update(moe.GetMOE(mol))
        res.update(basak.Getbasak(mol))

        #print(res)
        #input()
    return res


def num_vector_from_descriptor_vector(descriptor_vector):
    '''
    Given a sequence (map) of ("descriptor_name" -> value), returns a vector of values
    :param descriptor_vector:
    :return:
    '''
    x = []
    for k, v in descriptor_vector.items():
        x.append(v)
    return x

if __name__ == "__main__":
    print("TODO")
