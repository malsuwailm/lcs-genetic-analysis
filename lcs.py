'''

This project is a Python program that generates and tests the LCS algorithm for DNA sequence pairs. The program 
reads input sequences from a file specified by the user and writes the output to another file. To evaluate the algorithm's 
robustness, the program also generates extreme test cases with predetermined lengths and errors.

The program initially scans the input sequences from a file and verifies that they are DNA sequences by ensuring that 
they contain only the characters 'A', 'C', 'G', and 'T'. It then generates three additional test cases: an extreme case 
with a sequence of one thousand single nucleotides, a complexity case with an arbitrarily generated sequence of 500 characters in length, 
and a case with an error character to test error handling.

The program calculates the LCS for each pair of DNA sequences and publishes the result along with the number of comparisons 
made during the algorithm's execution. If any of the input sequences contain invalid characters, the program detects the 
exception and instead of executing the LCS algorithm, it displays an error message.

The program also evaluates the asymptotic complexity of the LCS algorithm by executing it on pairs of increasing-length 
sequences and documenting the execution time and number of comparisons made during each iteration. The output is a table 
of the results.

The program then writes the input sequences, pairwise LCS results, total number of comparisons, and asymptotic complexity 
table to an output file.
'''

import random
import time
import matplotlib.pyplot as plt


def generate_random_sequence(length):
    """
    Generates a random DNA sequence of the specified length.
    :param length: The length of the DNA sequence to generate.
    :return: A randomly generated DNA sequence.
    """
    return ''.join(random.choice('ACTG') for _ in range(length))


def test_asymptotic_complexity(min_len, max_len, step):
    """
    Test the asymptotic complexity of the LCS algorithm for sequences of increasing lengths, and return the results.

    Parameters:
    min_len (int): the minimum length of the sequences to test
    max_len (int): the maximum length of the sequences to test
    step (int): the step size to use when increasing the sequence length

    Returns:
    list: a list of lists, where each inner list contains the length of the sequences, the execution time of the LCS
          algorithm, and the number of comparisons made
    """

    results = []

    try:
        for length in range(min_len, max_len + 1, step):
            seq1 = generate_random_sequence(length)
            seq2 = generate_random_sequence(length)

            start_time = time.time()
            _, comparisons = longest_common_subsequence(seq1, seq2)
            execution_time = time.time() - start_time

            results.append([length, execution_time, comparisons])
    except Exception as e:
        print(f"Error occurred while generating random sequence: {e}")
        return []

    return results

def plot_asymptotic_cost(results):
    """
    Plot the asymptotic costs observed (execution time and number of comparisons) for the LCS algorithm
    as a function of input sequence length.

    Parameters:
    results (list): a list of lists, where each inner list contains the length of the sequences, the execution time
                    of the LCS algorithm, and the number of comparisons made

    Returns:
    None
    """
    # Extract the data
    lengths = [result[0] for result in results]
    execution_times = [result[1] for result in results]
    comparisons = [result[2] for result in results]

    # Create the plot
    fig, ax1 = plt.subplots()

    # Plot execution times
    ax1.set_xlabel("Length")
    ax1.set_ylabel("Execution Time", color="tab:blue")
    ax1.plot(lengths, execution_times, color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")

    # Create a second y-axis for comparisons
    ax2 = ax1.twinx()
    ax2.set_ylabel("Comparisons", color="tab:red")
    ax2.plot(lengths, comparisons, color="tab:red")
    ax2.tick_params(axis="y", labelcolor="tab:red")

    # Set the title
    plt.title("Asymptotic Costs Observed")

    # Show the plot
    plt.show()


def print_results_table(results):
    """
    Print the results of the asymptotic complexity test in a tabular format.

    Parameters:
    results (list): a list of lists, where each inner list contains the length of the sequences, the execution time
                    of the LCS algorithm, and the number of comparisons made.

    Returns:
    None
    """
    print("Length | Execution Time | Comparisons")
    print("-------------------------------------")
    for row in results:
        print(f"{row[0]:<6} | {row[1]:<14.6f} | {row[2]}")


def generate_extreme_case():
    """
    Generate an extreme test case with a sequence of 1000 nucleotides, with one randomly selected nucleotide.

    Parameters:
    None

    Returns:
    str: a string of 1000 nucleotides
    """
    nucleotides = ['A', 'T', 'C', 'G']
    selected_nucleotide = random.choice(nucleotides)
    seq = selected_nucleotide * 1000
    return seq


def generate_complexity_case(n):
    """
    Generate a test case with a random sequence of the specified length.

    Parameters:
    n (int): the length of the sequence to generate

    Returns:
    str: a string of length n containing random characters from the 'ACTG' sequence
    """
    alphabet = ['A', 'C', 'G', 'T']
    seq = ''.join(random.choices(alphabet, k=n))
    return seq


def generate_input_with_errors():
    """
    Generate a test case with an error character to test error handling.

    Parameters:
    None

    Returns:
    str: a string of 10 random characters from the 'ACTG' sequence, followed by a single random character that is not
         from the 'ACTG' sequence
    """
    valid_chars = ['A', 'C', 'G', 'T']
    error_chars = ['B', 'D', 'E', 'F', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'U', 'V', 'W', 'X', 'Y', 'Z']
    seq = ''.join(random.choices(valid_chars, k=10)) + ''.join(random.choices(error_chars, k=1))
    return seq


def validate_sequence(seq):
    """
    Validates that a given DNA sequence only contains valid characters.
    :param seq: The DNA sequence to validate.
    :raises ValueError: If the sequence contains an invalid character.
    """
    for char in seq:
        if char not in ['A', 'C', 'G', 'T']:
            raise ValueError(f"Invalid character '{char}' in sequence '{seq}'")
        
   
def append_test_case_to_file(filename, test_case_label, sequence):
    """
    Appends a test case label and its corresponding sequence to a given file.
    :param filename: The name of the file to append to.
    :param test_case_label: The label for the test case.
    :param sequence: The DNA sequence for the test case.
    """
    with open(filename, 'a') as file:
        file.write(f"{test_case_label} = {sequence}\n")


def longest_common_subsequence(X, Y):
    """
    Finds the longest common subsequence between two DNA sequences.
    :param X: The first DNA sequence.
    :param Y: The second DNA sequence.
    :return: A tuple containing the longest common subsequence and the number of character comparisons made.
    """
    try:
        m, n = len(X), len(Y)
        c = [[0] * (n + 1) for _ in range(m + 1)]
        comparisons = 0

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if X[i - 1] == Y[j - 1]:
                    c[i][j] = c[i - 1][j - 1] + 1
                else:
                    c[i][j] = max(c[i - 1][j], c[i][j - 1])
                comparisons += 1

        return c, comparisons
    except Exception as e:
        print(f"Error in function 'longest_common_subsequence': {e}")
        return None, 0


def print_lcs(c, X, Y, i, j):
    """
    Recursively constructs the longest common subsequence based on the input subsequence matrix.
    :param c: The subsequence matrix generated by the longest_common_subsequence function.
    :param X: The first DNA sequence.
    :param Y: The second DNA sequence.
    :param i: The current index of the first sequence.
    :param j: The current index of the second sequence.
    :return: The longest common subsequence of the two DNA sequences.
    """
    try:
        if i == 0 or j == 0:
            return ""
        if X[i - 1] == Y[j - 1]:
            return print_lcs(c, X, Y, i - 1, j - 1) + X[i - 1]
        if c[i - 1][j] > c[i][j - 1]:
            return print_lcs(c, X, Y, i - 1, j)
        return print_lcs(c, X, Y, i, j - 1)
    except Exception as e:
        print(f"Error in function 'print_lcs': {e}")
        return "" 


def read_sequences_from_user_file(filename):
    """
    Reads sequences from a file in the format 'name = sequence', and returns a list of the sequences.

    Args:
    - filename: A string representing the name of the file to read.

    Returns:
    - A list of strings representing the sequences.
    """
    sequences = []
    try:
        with open(filename, 'r') as file:
            for line in file:
                if '=' in line:
                    _, sequence = line.strip().split('=')
                    sequences.append(sequence.strip())
    except FileNotFoundError as e:
        print(f"Error: {e}")
    return sequences


def write_output_to_file(filename, strings, results, total_comparisons, table):
    """
    Writes the input sequences, pairwise LCS results, total comparisons, and asymptotic complexity table to a file.

    Args:
    - filename: A string representing the name of the output file.
    - strings: A list of strings representing the input sequences.
    - results: A list of strings representing the pairwise LCS results.
    - total_comparisons: An integer representing the total number of comparisons performed.
    - table: A list of lists representing the asymptotic complexity table.
    """
    try:
        with open(filename, 'w') as file:
            file.write("Input sequences:\n")
            for s in strings:
                file.write(f"{s}\n")

            file.write("\nPairwise LCS results:\n")
            for r in results:
                file.write(f"{r}\n")

            file.write(f"\nTotal comparisons: {total_comparisons}\n")

            # Add table to the output
            file.write("\nAsymptotic Complexity Table:\n")
            file.write("Length | Execution Time | Comparisons\n")
            file.write("-------------------------------------\n")
            for row in table:
                file.write(f"{row[0]:<6} | {row[1]:<14.6f} | {row[2]}\n")
    except Exception as e:
        print(f"Error in function 'write_output_to_file': {e}")


filename = 'input.txt'
filename_extreme_cases = 'inputExtremeCases.txt'

# Read sequences from the user-specified input file
strings = read_sequences_from_user_file(filename)

# Add pre-defined test cases for extreme cases and complexity
extreme_case = generate_extreme_case()
append_test_case_to_file(filename_extreme_cases, 'S_extreme', extreme_case)
strings.append(extreme_case)

complexity_case = generate_complexity_case(500)
append_test_case_to_file(filename_extreme_cases, 'S_complex', complexity_case)
strings.append(complexity_case)

input_with_errors = generate_input_with_errors()
append_test_case_to_file(filename_extreme_cases, 'S_error', input_with_errors)
strings.append(input_with_errors)

total_comparisons = 0
results = []

# Find the LCS for all pairs of sequences, and print the results
for i in range(len(strings)):
    for j in range(i + 1, len(strings)):
        try:
            validate_sequence(strings[i])
            validate_sequence(strings[j])

            c, comparisons = longest_common_subsequence(strings[i], strings[j])
            lcs = print_lcs(c, strings[i], strings[j], len(strings[i]), len(strings[j]))
            total_comparisons += comparisons
            result = f"LCS of {strings[i]} and {strings[j]}: {lcs} (Comparisons: {comparisons})"
            results.append(result)
            print(result)
        except ValueError as e:
            print(f"Error: {e}")

print(f"Total comparisons: {total_comparisons}")

# Write the output to a file, including the input sequences, pairwise LCS results, total comparisons, and complexity table
output_filename = 'output.txt'
results_table = test_asymptotic_complexity(100, 500, 100)
print_results_table(results_table)
plot_asymptotic_cost(results_table)
write_output_to_file(output_filename, strings, results, total_comparisons, results_table)