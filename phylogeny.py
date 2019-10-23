import argparse
import os
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--inputFile", dest="input_file", required=True,
        help="input fasta file"
    )

    args = parser.parse_args()

    # Get the path to the directory of this script
    script_dir = os.path.dirname(__file__)
    script_dir = script_dir + "/"

    # Get the path to the fasta file
    input_file_path = os.path.abspath(args.input_file)

    # Open the input fasta file
    fasta_file = open(input_file_path, "r")

    # Extract the identifiers and sequences from the fasta file
    identifiers = []
    sequences = []

    for line in fasta_file:
        if line.startswith(">"):
            line = line.rstrip()
            # Remove the ">" from the string and add it to the identifiers list
            identifiers.append(line[1:])
        else:
            line = line.rstrip()
            sequences.append(line)

    # Calculate the genetic distance between every pair of sequences and store in a dictionary
    genetic_distance_dict = {}

    for i in range(0, len(identifiers)):
        genetic_distance_list = []
        for j in range(0, len(sequences)):
            genetic_distance = calculate_genetic_distance(sequences[i], sequences[j])
            genetic_distance_list.append(genetic_distance)
        genetic_distance_dict[identifiers[i]] = genetic_distance_list

    # Create a tab-delimited file of genetic distance in the directory where this script lives
    genetic_distance_file = open(script_dir + "genetic_distances.txt", "w")

    # Make the header for the file, start with a blank "cell"
    genetic_distance_file.write("\t")
    genetic_distance_file.write("\t".join(identifiers))
    genetic_distance_file.write("\n")
    # Iterate through all the items in genetic_distance_dict and write it to the result file
    for key, val in genetic_distance_dict.items():
        genetic_distance_file.write(key + "\t")
        # Convert the list of float into string so it can be written to a file
        val = [str(i) for i in val]
        genetic_distance_file.write("\t".join(val))
        genetic_distance_file.write("\n")

    # Calculate the Q matrix and store the Q distances in a dictionary
    Q_matrix_dict = {}
    for i in range(0, len(identifiers)):
        Q_matrix_dict[identifiers[i]] = []
        # Calculate the sum of distances to the ith sequence
        sum_distance_i_k = sum(genetic_distance_dict[identifiers[i]])
        for j in range(0, len(identifiers)):
            # No need to calculate Q for the same sequences
            # Have the Q_distance be 0 as a placeholder
            if i == j:
                Q_matrix_dict[identifiers[i]].append(0)
            else:
                # Calculate the sum of distances to the jth sequence
                sum_distance_j_k = sum(genetic_distance_dict[identifiers[j]])
                # Calculate the Q distance
                Q_distance = ((len(identifiers) - 2) * genetic_distance_dict[identifiers[i]][j] -
                              sum_distance_i_k - sum_distance_j_k)
                # Append the Q distance as a list, which will be the values for the ith sequence key
                Q_matrix_dict[identifiers[i]].append(Q_distance)

    # Find the minimum distance in the Q matrix
    # Initialize the overall minimum Q distance as 0
    min_q_distance = 0
    for key in Q_matrix_dict.keys():
        # Check the minimum Q distance for a given sequence
        min_q_distance_for_sequence = min(Q_matrix_dict[key])
        # If the minimum Q distance for the sequence is smaller than the overall min Q distance,
        # replace the overall minimum Q distance
        if min_q_distance_for_sequence < min_q_distance:
            min_q_distance = min_q_distance_for_sequence
            # Store the sequence (key) where this happens
            min_q_distance_key = key
            # Find and store the index of the min Q distance in the list of values for the key
            min_q_distance_index = Q_matrix_dict[key].index(min_q_distance)
            # See against which sequence the min Q value was calculated relative to the key
            min_q_distance_partner = identifiers[min_q_distance_index]
    print(min_q_distance)
    print(min_q_distance_key)
    print(min_q_distance_partner)

    fasta_file.close()
    genetic_distance_file.close()


def calculate_genetic_distance(first_sequence, second_sequence):
    """
    Calculates the genetic distance between two sequences. If the two sequences is not the same length,
    the program will exit with an error
    :param first_sequence:
    :param second_sequence:
    :return: genetic_distance
    """
    similarity_score = 0
    if len(first_sequence) != len(second_sequence):
        sys.exit("Length of two sequences needs to be the same!")
    elif len(first_sequence) == len(second_sequence):
        for position in range(0, len(first_sequence)):
            if first_sequence[position] == second_sequence[position]:
                similarity_score += 1
        genetic_distance = 1 - (similarity_score / len(first_sequence))
        return genetic_distance


if __name__ == "__main__":
    main()
