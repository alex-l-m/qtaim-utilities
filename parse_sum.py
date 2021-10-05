import re
import sys
import csv

variable_description_re = r"^ *([^ ^=][^=]*[^ ^=]) *= *([^ ^=][^=]*[^ ^=])"

def parse_variables(lines):
    variable_description_dict = dict()
    for line in lines:
        variable, description = \
                re.match(variable_description_re, line).groups()
        variable_description_dict.update([(variable.strip(), description.strip())])
    return variable_description_dict

def parse_atomic_properties(title, header, rows, variable_description_dict):
    # Some special cases: some tables which include variables that are never
    # given descriptions
    if title == "Nuclear Charges and Cartesian Coordinates":
        variable_description_dict.update(\
                [("Charge", ""), ("X", ""), ("Y", ""), ("Z", "")])
    #elif title == "Eigenvalues and Eigenvectors of Atomic Traceless Quadrupole Moment Tensors":
    #    variable_description_dict.update([("n", "")])

    # Names of atom name properties (such as Atom A, etc).  Also includes "n",
    # as a hack, because really this should be row key, not atom names
    atom_colnames = []
    potential_atom_colnames = ["Atom", "Atom A", "Atom B", "Atom C", "n"]
    for atom_colname in potential_atom_colnames:
        # Need boundaries, otherwise any variable with an "n" in it will get
        # included. Kind of weird that "Atom" gets included even if it's really
        # "Atom A"
        if re.search(r"\b{}\b".format(atom_colname), header) is not None:
            atom_colnames.append(atom_colname)
    colnames = atom_colnames + list(variable_description_dict.keys())
    colname_starts = [header.find(variable) for variable in colnames]
    colname_ends = [start + len(variable) for start, variable in zip(colname_starts, colnames)]
    colname_middles = [(start + end)/2 for start, end in zip(colname_starts, colname_ends)]
    threed_property_table = []
    one_atom_property_table = []
    two_atom_property_table = []
    three_atom_property_table = []
    for row in rows:
        atom_table = []
        for match in re.finditer(r"[^ ]+\b", row):
            start, end = match.span()
            value_middle = (start + end) / 2
            variable, distance = min([\
                (colname, abs(colname_middle - value_middle)) \
                for colname, colname_middle in zip(colnames, colname_middles)\
                ], key=lambda x: x[1])
            if "Atom" in variable or variable == "n":
                atom_table.append((variable, match.group(0)))
            else:
                # Values of the atom names (like C1, H1, H2, etc)
                atom_names = [i[1] for i in sorted(atom_table, key=lambda x: x[0])]
                line = [title] + atom_names + [variable, match.group(0)]
                if "n" in atom_colnames:
                    threed_property_table.append(line)
                elif len(atom_names) == 1:
                    one_atom_property_table.append(line)
                elif len(atom_names) == 2:
                    two_atom_property_table.append(line)
                elif len(atom_names) == 3:
                    three_atom_property_table.append(line)
    return one_atom_property_table, two_atom_property_table, three_atom_property_table, threed_property_table

state = "other"

title_buffer= None
possible_title_buffer = None
variable_buffer = []
header_buffer = None
data_row_buffer = []
one_atom_property_table = []
two_atom_property_table = []
three_atom_property_table = []
threed_property_table = []

for line in open(sys.argv[1]):
    if re.match("^[^ ].*: *$", line) is not None and \
            (state == "other" or state == "variable_description_line" or \
            state == "variable_description_interruption" or state == "possible_title" \
            or state == "blank_after_total"):
        state = "possible_title"
    elif re.match("^-+ *$", line) is not None:
        if state == "possible_title":
            state = "dashes_after_title"
        elif state == "variable_description_line" or state == "variable_description_interruption":
            state = "dashes_after_variable_descriptions"
        elif state == "header":
            state = "dashes_after_header"
        elif state == "data_row":
            state = "dashes_after_data"
    elif re.match(variable_description_re, line) is not None and \
            (state == "dashes_after_title" or state == "variable_description_line" or \
            state == "variable_description_interruption" or state == "possible_title"):
        state = "variable_description_line"
    elif re.match("^Atom", line) is not None:
        state = "header"
    elif state == "variable_description_line" or state == "variable_description_interruption":
        state = "variable_description_interruption"
    elif re.search(r"\w+", line) is not None and \
            (state == "dashes_after_header" or state == "data_row"):
        state = "data_row"
    elif re.search("Total", line) is not None and \
            (state == "dashes_after_data" or state == "blank_after_total"):
        state = "total_line"
    elif re.match("^ *$", line) is not None and \
            state == "total_line":
        state = "blank_after_total"
    elif "=" not in line and ":" not in line and \
            re.match("^ *$", line) is None and state == "blank_after_total":
        state = "data_row"
    else:
        state = "other"

    # Parse the table if we are no longer in a table
    # Can't do this in a possible title because those can occur during variable
    # descriptions
    if state == "other" or state == "dashes_after_title":
        if len(data_row_buffer) > 0:
            this_one_atom_property_table, this_two_atom_property_table, this_three_atom_property_table, this_threed_property_table = \
                    parse_atomic_properties(title_buffer, header_buffer, data_row_buffer, \
                    parse_variables(variable_buffer))
            one_atom_property_table += this_one_atom_property_table
            two_atom_property_table += this_two_atom_property_table
            three_atom_property_table += this_three_atom_property_table
            threed_property_table += this_threed_property_table
        title_buffer= None
        variable_buffer = []
        header_buffer = None
        data_row_buffer = []

    if state == "possible_title":
        possible_title_buffer = re.match("[^:]*", line).group(0).strip()
    if state == "dashes_after_title":
        title_buffer = possible_title_buffer
    elif state == "variable_description_line":
        variable_buffer.append(line)
    elif state == "header":
        header_buffer = line
    elif state == "data_row":
        data_row_buffer.append(line)

prefix = sys.argv[2]

with open(prefix + "_oneatom.csv", "w") as f:
    csv.writer(f).writerows(one_atom_property_table)

with open(prefix + "_twoatom.csv", "w") as f:
    csv.writer(f).writerows(two_atom_property_table)

with open(prefix + "_threeatom.csv", "w") as f:
    csv.writer(f).writerows(three_atom_property_table)

with open(prefix + "_3d.csv", "w") as f:
    csv.writer(f).writerows(threed_property_table)
