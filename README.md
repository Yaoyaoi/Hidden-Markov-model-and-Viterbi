Implementation of Viterbi using python

### command line to run the program 
1. run the program with: 
    ` python3 ViterbiHmm.py `
2. then follow the tips to enter the files' name, like:
    `hmm file name:example.hmm`
    `fasta file name:example.fa`
* fasta file format:
contain multiple lines or just two lines.
    ```
    >test1
    AACCGGA
    Or
    >test2
    AACGG
    AGCT
    AAACGTA
    ```

### environment 
Only use the math package
python 3.6.8

### output 
My program outputs both in the command line and the file named 'output.txt'.If the command line can not show all the outputs, you can check the file 'output.txt'