# shannon-entropy

Sandbox for attempt to create an iterative method that recovers active molecules from a dataset using the shannon entropy contribution of test molecules to the already existing active molecules in the training set. This is losely based off off Wang et. al (2009)[1].


Requirements 
---------------

[rdkit](http://www.rdkit.org/)


How to Run
----------------

```
python shannon.py --input inputfile.csv 
``` 

where __inputfile.csv__ has the format of: 

```
SMILES_1, potency_value1
SMILES_2, potency_value2
...
SMILES_n, potency_valuen
```

Review
--------------

The script will take 2 random molcules from the input file, and take the 
most potent molcule and designate it 'active'. It will then look at all the remaining molecules in the file and look for compounds which minimally increase the shannon entropy of the dataset, were it to be included, using the equation: 

```
shannon entropy = sum(-x*log(x) - (1-x)*log(1-x)) 
```

where x is the frequency of 'on' bits at every index of the 1024-bit fingerprint of the molcule set (Note that the entropy is calculated before and after the addition of a test molecule, and those with lost dEntropy are added in batches of 96). Here, the Log is base 2, and values of x of 0 or 1 result in a dEntropy of 0.  


References 
----------------
[1] J. Chem. Inf. Model., 2009, 49 (7), pp 1687–1691



