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

    print(genetic_distance_dict)

    fasta_file.close()


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
