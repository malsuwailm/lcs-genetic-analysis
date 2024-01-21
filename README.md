## Running the LCS algorithm for DNA sequences

This is a Python program that generates and tests the LCS algorithm for DNA sequence pairs. The program reads input sequences from a file specified by the user and writes the output to another file. To evaluate the algorithm's robustness, the program also generates extreme test cases with predetermined lengths and errors.
Prerequisites

Make sure you have Python 3 installed on your machine, with libraries random, time and matplotlib.

## Running the program:

    1. Download the project files and extract them to a folder.
    2. Open a terminal window and navigate to the folder where the project files are located.
    3. Run the following command to execute the program:

```
python lcs.py
```

    1. The program will then execute the LCS algorithm for all pairs of DNA sequences and print the results to the console. 2. The total number of comparisons made during the execution of the algorithm will also be displayed.
    3. The program will write the input sequences, pairwise LCS results, total number of comparisons, a asymptotic complexity table to the output file and showcase a graph.

## Generating additional test cases

The program automatically generates three additional test cases: an extreme case with a sequence of one thousand single nucleotides, a complexity case with an arbitrarily generated sequence of 500 characters in length, and a case with an error character to test error handling. These automated test cases are stored in a separate input file named inputExtremeCases.txt.

To add user custom test cases, you can append them to the input.txt file in the following format:

Test Case Name = DNA sequence

Replace "Test Case Name" with a unique name for the test case, and "DNA sequence" with the DNA sequence you want to test. The program will automatically read these test cases from the file and include them in the execution.
