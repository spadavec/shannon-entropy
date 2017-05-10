import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import csv
import numpy
from rdkit.Chem import DataStructs
from random import randint
from operator import itemgetter
import math
import csv
import argparse



def get_train_actives(mol_list):
    """ Split the training dataset so that 20\% most active compounds are
    considered active """

    # Take only 20% most active
    mol_list.sort(key=lambda x: float(x[1]))
    num_to_keep = int(math.ceil(0.2*len(mol_list)))
    actives = mol_list[0:num_to_keep]

    return actives, mol_list

def get_index_totals(fps):
    """ Counts the number of '1' bits at every FP index of a given mol """
    se_i = []

    indexes = [x for x in range(1024)]
    idx_counts = {key: 0 for key in indexes}
    ct = 0
    for fp in fps:
        for i,x in enumerate(fp):
            if x is '1':
                ct += 1
                idx_counts[i] += 1

    return idx_counts


def check_addition(actives_np_fps,m):
     """ Checks the change in shannon entropy of a given molecule 'm' against
        a set of active molecules (actives_np_fps) """

     frequencies = []
     se_i = []
     mol = Chem.MolFromSmiles(m[0])

     fp = AllChem.GetMorganFingerprintAsBitVect(mol,6,nBits=1024)
     np_fps = convert_fps(fp)
     actives_np_fps.append(np_fps)

     total_len = float(len(actives_np_fps))

     flattened = [item for sublist in actives_np_fps for item in sublist]


     # Get the number of 1s
     ones_dict = get_index_totals(flattened)
     ones = ones_dict.values()

     # Calculate frequency of 1s in every position
     for x in ones:
         frequencies.append(float(x/total_len))

     for x in frequencies:
         if x in {0,1}:
             se = 0
         else:
             # Use Shannon Entropy definition as detailed in Publication
             se = -x*math.log(x,2) - (1-x)*math.log(1-x,2)

         se_i.append(se)

     se_t = sum(se_i)

     # Remove the appended 'test molecule'
     del actives_np_fps[-1]

     return se_t




def read_csv(filename):
    smiles = []
    ic50 = []

    with open(filename) as ifile:
        data = csv.reader(ifile, delimiter=',')

        for line in data:
            smiles.append(line[0])
            ic50.append(float(line[1]))

    return smiles, ic50

def convert_fps(fp):
    """ Converts RDKit Fingerprints to numpy array """
    np_fps = []
    array = numpy.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, array)
    np_fps.append(''.join([str(int(x)) for x in array]))

    return np_fps


def write_csv(iter, arr):
    l = iter+1
    rounds = [x for x in range(l)]

    final = zip(rounds,arr)

    with open('output.csv', 'a') as ofile:
        writer = csv.writer(ofile)
        writer.writerows(final)



def end_condition(mol_list, ic50_cutoff, total_num_actives, percents):
    """ Checks for end condition; whether all 'active' molecules have been
    recovered (top five percent most potent) """

    ct = 0
    for x in mol_list:
        if x[1] <= ic50_cutoff:
            ct += 1

    fraction_found = float(float(ct)/float(total_num_actives))*100
    percents.append(fraction_found)

    if ct < total_num_actives:
        return False
    else:
        print ""
        print "All actives found"
        return True


def grab_lowest_entropy_idx(se, num_to_grab):
    """ Finds index of test compounds in test list which have lowest shannon
    entropy contribution to the training set"""

    idx = []
    sorted_se = sorted(se)

    cutoff_se = sorted_se[num_to_grab]

    ct = 0

    if cutoff_se == 0:
        idx = [x for x in range(num_to_grab)]
    else:
        for i,x in enumerate(se):
            if float(x) <= cutoff_se:
                if ct < num_to_grab:
                    idx.append(i)
                    ct += 1

    return idx

def main():
    percents = []
    test = []
    train = []
    iteration_num = 1

    # Set run type; 'fixed' maxes selection pool to 96 (to mimic cell plates)
    # If set to anything else, will select as many as it sees fit 
    run_type = 'fixed'

    # Argument parser
    parser = argparse.ArgumentParser(description='argument parser for shannon entropy search')
    parser.add_argument('--input', type=str, dest="input_file")
    args = parser.parse_args()
    input_csv_name = args.input_file



    print ""
    print "Recovering actives from dataset..."

    # Grab data
    smiles, ic50 = read_csv(input_csv_name)

    # Get total number of actives to be retrieved
    # Set "actives" to be 5% most potent compounds
    sorted_ic50 = sorted(ic50)
    total_num_actives = int(math.ceil(len(sorted_ic50)*float(0.05)))
    ic50_cutoff = sorted_ic50[total_num_actives]
    print "The potency cutoff for this dataset is {}".format(ic50_cutoff)
    print "There are a total of {} actives in this dataset ".format(total_num_actives)
    print "Picking 2 random molecules to seed the search..."

    ceiling = len(ic50)

    #Select 2 random molecules to start the process
    rand_indexes = []
    for x in range(0,2):
        rand_indexes.append(randint(0,ceiling))

    # Create training set
    for x in rand_indexes:
        train.append([smiles[x], ic50[x]])

    # Create the test set
    for i,x in enumerate(smiles):
        if i not in rand_indexes:
            test.append([smiles[i], ic50[i]])

    print "Begining search based on shannon entropy..."
    while not end_condition(train, ic50_cutoff, total_num_actives, percents):
        se_idxs = []
        se = []
        print ""
        print "-"*25
        print "Iteration #{}".format(iteration_num)
        print "-"*25

        # Determine all the actives in the training set
        actives, train_mols = get_train_actives(train)

        # Generate all necessary training data
        train_mols = [Chem.MolFromSmiles(x[0]) for x in train_mols]
        train_fps = [AllChem.GetMorganFingerprintAsBitVect(x,6,nBits=1024)
            for x in train_mols]
        train_np_fps = [convert_fps(x) for x in train_fps]

        # Generate intial actives
        actives_mols = [Chem.MolFromSmiles(x[0]) for x in actives]
        actives_fps = [AllChem.GetMorganFingerprintAsBitVect(x,6,nBits=1024)
            for x in actives_mols]
        actives_np_fps = [convert_fps(x) for x in actives_fps]


        # Get all of the entropies of the test molecules
        for x in test:
            se.append(check_addition(actives_np_fps, x))


        if run_type == 'free':
            num_to_grab = int(math.ceil(int(len(test)*0.2)))
        else:
            num_to_grab = 96

        # Get the indexes of the lowest entropy additions
        se_idxs = grab_lowest_entropy_idx(se, num_to_grab)

        # Append the new compounds to train
        for i,x in enumerate(test):
            if i in se_idxs:
                train.append(x)

        # Delete the new compounds from test
        for index in sorted(se_idxs, reverse=True):
            del test[index]


        len_train = float(len(train))
        len_test = float(len(test))
        len_total = len_test+ len_train
        fraction = (len_train/len_total)*100
        print "Total fraction searched: {}% ".format(int(fraction))
        print "Total fraciton of actives found: {}% ".format(int(percents[-1]))


        iteration_num += 1

    write_csv(iteration_num, percents)


if __name__ == "__main__":
    main()
