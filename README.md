# Transposon Detection in *C. elegans*
Scripts for simulation and detection of Transposons

#### Simulations
Simulate transposon insertions and exsisions with BAMSurgeon and RSVSim and evaluate results
```
bash run_simulations.sh
bash analyze_simulations.sh
```


#### Set Up Directory Structure
```
mkdir data
mkdir results
mkdir piRNA
mkdir tables
```

#### Run Transposon Detection Programs
Run TEMP to detect insertion events relative to the reference and TE-Locate to detect absences or maintenance of transposons relative to the reference. In data directory:
```
bash run_TransposonCaller.sh
```

#### Process Transposon Caller Outputs
* Filter out redundant calls
* Apply read support and population frequency filters
* Convert transposon names
* Collapse transposons
* Produce consistently formatted output files
```
bash run_ProcessTransposonCallerOutput.sh
```
#### Define Traits and Classify
* Check for and resolve contradictory calls
* Merge transposons
* Generate matrix of calls
* Filter monomorphic sites
* Account for sites with no coverage
* Calculate trait metrics
* Check transposon classification and prevalence
* Find intersection with genes
```
bash CleanTransposons.sh
```

#### GWAS Results
Obtain mapping results file and analyze GWAS results in context of genes of interest and piRNAs
```
bash PostMappings.sh
```

