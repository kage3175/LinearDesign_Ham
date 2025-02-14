<img src="/pic/baidu_research_logo.jpg"  width="40%" alt="Baidu Research Logo">

#  Algorithm for Optimized mRNA Design Improves Stability and Immunogenicity (LinearDesign)
![GitHub all releases](https://img.shields.io/github/downloads/LinearDesignSoftware/LinearDesign/total)


This repository contains the source code for the LinearDesign project.

He Zhang†, Liang Zhang†, Ang Lin†, Congcong Xu†, Ziyu Li, Kaibo Liu, Boxiang Liu, Xiaopin Ma, Fanfan Zhao, Huiling Jiang, Chunxiu Chen, Haifa Shen, Hangwen Li*, David H. Mathews*, Yujian Zhang*, Liang Huang†*<sup>#</sup>. Algorithm for Optimized mRNA Design Improves Stability and Immunogenicity. Nature [https://doi.org/10.1038/s41586-023-06127-z](https://doi.org/10.1038/s41586-023-06127-z) (2023)

† contributed equally, 
\* corresponding authors, 
<sup>#</sup> lead corresponding author

For questions, please contact the lead corresponding author at <liang.huang.sh@gmail.com>.

## Dependencies
Clang 11.0.0 (or above) or GCC 4.8.5 (or above)

python2.7

LinearPartition should be installed in same directory with LinearDesign, and compiled.
https://github.com/LinearFold/LinearPartition

## To Compile
```
make
```

## To Run
The LinearDesign program can be run with:
```
echo SEQUENCE | ./lineardesign [OPTIONS]

OR

cat FASTA_FILE | ./lineardesign [OPTIONS]
```

OPTIONS:
```
--lambda LAMBDA or -l LAMBDA
```
Set LAMBDA, a hyperparameter balancing MFE and CAI. (default 0.0)
```
--codonusage FILE_NAME or -c FILE_NAME
```
Import a Codon Usage Frequency Table. See "codon_usage_freq_table_human.csv" for the format.
(default: using human codon usage frequency table)
```
--verbose or -v
```
Print out more details. (default False)

For Macbook, users may encounter a pop-up message at the first run.
For Mac-M1 system, the message is:
```
"LinearDesign_Mac_M1.so" can't be opened because Apple cannot check it for malicious software.
```
For Mac-Intel system, the message is:
```
"LinearDesign_Mac_Intel.so" cannot be opened because it is from an unidentified developer.
```
If so, please go to "System Preferences -> Security & Privacy -> General" to allow LinearDesign-Mac-M1.so (or LinearDesign-Mac-Intel.so) to open.

## Example: Single Sequence Design
```
echo MNDTEAI | ./lineardesign
Penalties: None
mRNA sequence:  AUGAACGAUACGGAGGCGAUC
mRNA structure: ......(((.((....)))))
mRNA folding free energy: -1.10 kcal/mol; mRNA CAI: 0.695
```

## Example: Multiple Sequences Design with Option --lambda (-l)
```
cat testseq | ./lineardesign --lambda 3
>seq1
mRNA sequence:  AUGCCAAACACCCUGGCAUGCCCC
mRNA structure: ((((((.......)))))).....
mRNA folding free energy: -6.00 kcal/mol; mRNA CAI: 0.910

>seq2
mRNA sequence:  AUGCUGGAUCAGGUGAACAAGCUGAAGUACCCAGAGGUGAGCCUGACCUGA
mRNA structure: .....((.((((((..((...(((.......)))..))..))))))))...
mRNA folding free energy: -13.50 kcal/mol; mRNA CAI: 0.979
```

## Example: Option --codonusage (-c)
```
Penalties: None
echo MNDTEAI | ./lineardesign -l 0.3 --codonusage codon_usage_freq_table_yeast.csv
mRNA sequence:  AUGAAUGAUACGGAAGCGAUC
mRNA structure: ......(((.((....)))))
mRNA folding free energy: -1.10 kcal/mol; mRNA CAI: 0.670
```

## Example: Option --verbose (-v)
```
echo MNDTEAI | ./lineardesign --verbose
Penalties: None
Input protein: MNDTEAI
Using lambda = 0; Using codon frequency table = codon_usage_freq_table_human.csv
mRNA sequence:  AUGAACGAUACGGAGGCGAUC
mRNA structure: ......(((.((....)))))
mRNA folding free energy: -1.10 kcal/mol; mRNA CAI: 0.695
Runtime: 0.002 seconds
```

## Example: Option --penalty (-p), for lineardesign
```
echo MNDTEAI | ./lineardesign -p 2~4,6~8
Penalties: 2~4,6~8
mRNA sequence:  AUGAACGACACCGAGGCCAUC
mRNA structure: (((..(........)..))).
mRNA folding free energy: -0.30 kcal/mol
mRNA CAI: 1.000
```

## Example: Option --penalty (-p), for lineardesign_legacy
```
echo MNDTEAI | ./lineardesign_legacy -p 2~4
Penalties: 2~4
mRNA sequence:  AUGAACGACACCGAGGCCAUC
mRNA structure: (((..(........)..))).
mRNA folding free energy: -0.30 kcal/mol
mRNA CAI: 1.000
```

## Example: Option --score (-s)
```
echo MNDTEAI | ./lineardesign -s
Penalties: None
mRNA sequence:  AUGAACGAUACGGAGGCGAUC
mRNA structure: ......(((.((....)))))
mRNA folding free energy: -1.10 kcal/mol
mRNA CAI: 0.695
score: 5.337136760587025
```

## Example: Option --gapsearch (-g) for lineardesign_legacy
```
echo MNDTEAIVVMDY | ./lineardesign_legacy -g 2
Penalties: 10~12
mRNA sequence:  AUGAAUGAUACGGAGGCCAUCGUGGUCAUGGAUUAU
mRNA structure: ....(((((.(.(.(((((...))))).).))))))
mRNA folding free energy: -6.90 kcal/mol
mRNA CAI: 0.810

Penalties: 11~13
mRNA sequence:  AUGAAUGAUACGGAGGCCAUCGUGGUCAUGGAUUAU
mRNA structure: ....(((((.(.(.(((((...))))).).))))))
mRNA folding free energy: -6.90 kcal/mol
mRNA CAI: 0.810

...

Penalties: 23~25
mRNA sequence:  AUGAAUGAUACCGAGGCCAUCGUGGUCAUGGAUUAU
mRNA structure: ....(((((.(((.(((((...))))).))))))))
mRNA folding free energy: -10.20 kcal/mol
mRNA CAI: 0.887
```

## Example: Option --radial (-r) for lineardesign_legacy
```
cat rnapol | ./lineardesign_legacy -r rnapol
## rnapol.sorted saved in ./result directory
```

## Optional: Option --max_core_num (-m) with --radial option
```
cat rnapol | ./lineardesign_legacy -r rnapol -m 10
## max cpu core number to use for -r option
```


## Declarations
Baidu Research has filed a patent for the LinearDesign algorithm that lists He Zhang, Liang Zhang, Ziyu Li, Kaibo Liu, Boxiang Liu, and Liang Huang as inventors.
