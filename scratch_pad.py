identifiers = ['152801', '3855581', '4345542', '1009440', '996519', '990269', '954340', '947627', '816087', '809189',
               '806723', '720729', '664158', '593544', '586332', '565897', '555079', '553272', '536135', '523025',
               '511006', '348167', '319214', '286686', '273179', '260976', '208595', '197222', '188674', '166838',
               '164216', '146154', '144260', '106543', '100149', '11136', '5364', '698164', '1920715', '2034243',
               '2263183', '2678325', '847142', '779409', '2963287', '3148512', '3183891', '1136782', '646961',
               '3906666', '822879', '4164932', '4314491', '4378332', '4387271', '4397763', '4431405', '4439229',
               '4448694', '4468812', '4472218']

edges_file = open("/Users/m006703/Class/CSCI5481/Homework3/test_edges.txt", "r")

newick_string = ""
for line in edges_file:
    line = line.rstrip()
    line_item = line.split("\t")
    ancestral_node = line_item[0]
    descendant_node = line_item[1]
    node_distance = line_item[2]
    # Check to see if the ancestral node is the root
    if ancestral_node == str(len(identifiers) + 1):
        try:
            # If there are items already in the lists, do...
            #FIXME
            if len(ancestral_node_list) > 0:
                # Reverse the order of the lists for postorder since the edges.txt was preorder
                ancestral_node_list.reverse()
                descendant_node_list.reverse()
                node_distance_list.reverse()
                ancestral_node_count = []
                newick_string += "("
                for index in range(0, len(ancestral_node_list)):
                    # If the descendant node is a tip, do...
                    if int(descendant_node_list[index]) <= len(identifiers):
                        descendant_node_key = (int(descendant_node_list[index]) - 1)
                        descendant_node_identifier = identifiers[descendant_node_key]
                        if ancestral_node_list[index] not in ancestral_node_count:
                            newick_string += "(" + descendant_node_identifier + ":" + node_distance_list[index] + ","
                            ancestral_node_count.append(ancestral_node_list[index])
                        else:
                            newick_string += descendant_node_identifier + ":" + node_distance_list[index] + ")"
                    # If the descendant node is a node
                    elif int(descendant_node_list[index]) > len(identifiers):
                        if ancestral_node_list[index] not in ancestral_node_count:
                            newick_string += ":" + node_distance_list[index] + ","
                            try:
                                if ancestral_node_list[index] == ancestral_node_list[index + 1]:
                                    index_for_parenthesis = newick_string.rfind("(")
                                    newick_string = (newick_string[:(index_for_parenthesis + 1)] + "(" +
                                                     newick_string[(index_for_parenthesis + 1):])
                                else:
                                    newick_string = "(" + newick_string
                            except IndexError:
                                continue
                            ancestral_node_count.append(ancestral_node_list[index])
                        else:
                            newick_string += ":" + node_distance_list[index] + ")"
        except NameError:
            # If the list doesn't exist, create the list and add the root information
            ancestral_node_list = []
            descendant_node_list = []
            node_distance_list = []
            ancestral_node_list.append(ancestral_node)
            descendant_node_list.append(descendant_node)
            node_distance_list.append(node_distance)
    # If the ancestral node is not a root, add the node information to the lists
    else:
        ancestral_node_list.append(ancestral_node)
        descendant_node_list.append(descendant_node)
        node_distance_list.append(node_distance)

print(newick_string)