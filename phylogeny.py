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
    genetic_distance_dict, genetic_distance_file = create_genetic_distance_dict_and_file(identifiers, script_dir,
                                                                                         sequences)

    # Initialize by starting at root 120 being the first node to be joined
    root = 120
    edges_dict = {}

    # Loop until all of the neighbors are joined
    while root > len(identifiers):
        # Calculate the Q matrix and store the Q distances in a dictionary
        q_matrix_dict = calculate_q_distance(genetic_distance_dict)

        # Find the minimum distance in the Q matrix
        min_q_distance, min_q_distance_key, min_q_distance_partner = find_neighbors(q_matrix_dict)

        print(min_q_distance)
        print(min_q_distance_key)
        print(min_q_distance_partner)

        # Join the neighbors and calculate the distance to the newly formed root
        edges_dict = join_neighbors(edges_dict, genetic_distance_dict, min_q_distance_key, min_q_distance_partner, root)

        # Update genetic_distance_dict
        genetic_distance_dict = update_genetic_distance(genetic_distance_dict, min_q_distance_key,
                                                        min_q_distance_partner, root)
        # Subtract 1 from the root number and repeat
        root -= 1

    print(edges_dict)

    fasta_file.close()
    genetic_distance_file.close()


def update_genetic_distance(genetic_distance_dict, min_q_distance_key, min_q_distance_partner, root):
    """
    Update the genetic_distance_dict with the newly formed root
    :param genetic_distance_dict:
    :param min_q_distance_key:
    :param min_q_distance_partner:
    :param root:
    :return: Updated genetic_distance_dict
    """
    # Add the new root to the dictionary
    genetic_distance_dict[str(root)] = {}
    # Calculate the distance to the new node
    for key in genetic_distance_dict.keys():
        # No need to update the tip/root that is being joined
        if key == min_q_distance_key or key == min_q_distance_partner:
            continue
        # No need to do this for the newly formed root also because we will be copying the distances as we go
        if key == str(root):
            continue
        else:
            genetic_distance_dict[key][str(root)] = (0.5 * (genetic_distance_dict[min_q_distance_key][key] +
                                                            genetic_distance_dict[min_q_distance_partner][key] -
                                                            genetic_distance_dict[min_q_distance_key][
                                                                min_q_distance_partner]))
            # Copy the same distance to the nested dictionary of the newly formed root
            genetic_distance_dict[str(root)][key] = genetic_distance_dict[key][str(root)]
            # Delete the joined tip/root from the nested dictionary
            del genetic_distance_dict[key][min_q_distance_key]
            del genetic_distance_dict[key][min_q_distance_partner]
    # Delete the joined tip/root from the genetic distance dictionary
    del genetic_distance_dict[min_q_distance_key]
    del genetic_distance_dict[min_q_distance_partner]
    return genetic_distance_dict


def join_neighbors(edges_dict, genetic_distance_dict, min_q_distance_key, min_q_distance_partner, root):
    """
    Join the neighbors, calculate the distance from the neighbors to the newly formed root
    and store the information in the edges_dict dictionary
    :param edges_dict: Updates the already existing edges_dict
    :param genetic_distance_dict:
    :param min_q_distance_key:
    :param min_q_distance_partner:
    :param root:
    :return: Updated edges_dict
    """
    # Calculate the distance to the new node from the first branch
    first_branch_distance = (0.5 * genetic_distance_dict[min_q_distance_key][min_q_distance_partner] +
                             (1 / (2 * (len(list(genetic_distance_dict.keys())) - 2))) *
                             (sum(genetic_distance_dict[min_q_distance_key].values()) -
                              sum(genetic_distance_dict[min_q_distance_partner].values())))
    # Calculate the second branch distance
    second_branch_distance = genetic_distance_dict[min_q_distance_key][min_q_distance_partner] - first_branch_distance
    # Store the joined node information, start at number 120 and decrease every cycle
    edges_dict[str(root)] = {min_q_distance_key: first_branch_distance,
                             min_q_distance_partner: second_branch_distance}
    return edges_dict


def find_neighbors(Q_matrix_dict):
    """
    Find the tip/root with the smallest Q distance (neighbors)
    :param Q_matrix_dict:
    :return:
    """
    # Initialize the overall minimum Q distance as 0
    min_q_distance = 0
    for key in Q_matrix_dict.keys():
        # Check the minimum Q distance for a given sequence
        min_q_distance_for_sequence = min(Q_matrix_dict[key].values())
        # If the minimum Q distance for the sequence is smaller than the overall min Q distance,
        # replace the overall minimum Q distance
        if min_q_distance_for_sequence < min_q_distance:
            min_q_distance = min_q_distance_for_sequence
            # Store the sequence (key) where this happens
            min_q_distance_key = key
            # See against which tip or root the min Q value was calculated relative to the key
            for nested_key, nested_value in Q_matrix_dict[key].items():
                if nested_value == min_q_distance:
                    min_q_distance_partner = nested_key
    return min_q_distance, min_q_distance_key, min_q_distance_partner


def calculate_q_distance(genetic_distance_dict):
    """
    Calculates the Q distance between all tip/root and store it in a nested dictionary.
    Uses the same tip/root number as the genetic_distance_dict
    :param genetic_distance_dict:
    :return:
    """
    q_matrix_dict = {}
    for first_key in genetic_distance_dict.keys():
        node_i = first_key
        # Created a nested dictionary inside the Q_matrix_dict, the key will be the tip or root number
        q_matrix_dict[node_i] = {}
        # Calculate the sum of distances to the ith sequence
        sum_distance_i_k = sum(genetic_distance_dict[node_i].values())
        for second_key in genetic_distance_dict.keys():
            node_j = second_key
            # No need to calculate Q for the same sequences
            if first_key == second_key:
                continue
            else:
                # Calculate the sum of distances to the jth sequence
                sum_distance_j_k = sum(genetic_distance_dict[node_j].values())
                # Calculate the Q distance
                Q_distance = ((len(list(genetic_distance_dict.keys())) - 2) * genetic_distance_dict[node_i][node_j] -
                              sum_distance_i_k - sum_distance_j_k)
                # Append the root or tip number as the key and Q distance as the value
                q_matrix_dict[node_i][node_j] = Q_distance
    return q_matrix_dict


def create_genetic_distance_dict_and_file(identifiers, script_dir, sequences):
    """
    The key for the genetic_distance_dict will represent the tip or root number and will start with 1,
    where 1 is a tip and represents the first sequence in the fasta file.
    Until the dictionary gets updated, it will start out with the 61 sequences in the fasta file.
    The internal root node will start with 62 since there are 61 sequences in the fasta file.
    There are 59 (61 - 2) nodes possible in this phylogenetic tree, so the first pair to be joined will be
    root 120 (61 + 59).
    :param identifiers:
    :param script_dir:
    :param sequences:
    :return:
    """
    genetic_distance_dict = {}
    for i in range(0, len(identifiers)):
        # Create a nested dictionary that will store the genetic distances with the tip number as key
        nested_genetic_distance_dict = {}
        for j in range(0, len(sequences)):
            genetic_distance = calculate_genetic_distance(sequences[i], sequences[j])
            nested_genetic_distance_dict[str(j + 1)] = genetic_distance
        genetic_distance_dict[str(i + 1)] = nested_genetic_distance_dict
    # Create a tab-delimited file of genetic distance in the directory where this script lives
    genetic_distance_file = open(script_dir + "genetic_distances.txt", "w")
    # Make the header for the file, start with a blank "cell"
    genetic_distance_file.write("\t")
    genetic_distance_file.write("\t".join(identifiers))
    genetic_distance_file.write("\n")
    # Iterate through all the items in genetic_distance_dict and write it to the result file
    for key, val in genetic_distance_dict.items():
        sequence = identifiers[int(key) - 1]
        genetic_distance_file.write(sequence + "\t")
        # Iterate through the nested dictionary to write the genetic distances to a file
        for nested_value in genetic_distance_dict[key].values():
            genetic_distance_file.write(str(nested_value) + "\t")
        genetic_distance_file.write("\n")
    return genetic_distance_dict, genetic_distance_file


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
