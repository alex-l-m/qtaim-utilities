import re
import sys

def parse_variables(lines):
    variable_description_dict = dict()
    for line in lines:
        variable, description = \
                re.match(r"^ *([^ ^=][^=]*[^ ^=]) *= *([^ ^=][^=]*[^ ^=])", line).groups()
        variable_description_dict.update([(variable.strip(), description.strip())])
    return variable_description_dict

def parse_atomic_properties(header, rows, variable_description_dict):
    atom_colnames = []
    for atom_colname in ["Atom A", "Atom B", "Atom C"]:
        if atom_colname in header:
            atom_colnames.append(atom_colname)
    colnames = atom_colnames + list(variable_description_dict.keys())
    colname_starts = [header.find(variable) for variable in colnames]
    colname_ends = [start + len(variable) for start, variable in zip(colname_starts, colnames)]
    colname_middles = [(start + end)/2 for start, end in zip(colname_starts, colname_ends)]
    property_table = []
    for row in rows:
        atom_table = []
        for match in re.finditer(r"[^ ]+\b", row):
            start, end = match.span()
            value_middle = (start + end) / 2
            variable, distance = min([\
                (colname, abs(colname_middle - value_middle)) \
                for colname, colname_middle in zip(colnames, colname_middles)\
                ], key=lambda x: x[1])
            if "Atom" in variable:
                atom_table.append((variable, match.group(0)))
            else:
                atom_names = [i[1] for i in sorted(atom_table, key=lambda x: x[0])]
                property_table.append(atom_names + [variable, match.group(0)])
    return property_table

state = "other"

for line in open(sys.argv[1]):
    print(line, end="")
    if re.search(": *$", line) is not None and \
            (state == "other" or state == "variable_description_line" or \
            state == "variable_description_interruption" or state == "possible_title"):
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
    elif re.search("=", line) is not None and \
            (state == "dashes_after_title" or state == "variable_description_line" or \
            state == "variable_description_interruption" or state == "possible_title"):
        state = "variable_description_line"
    elif re.match("^Atom A", line) is not None:
        state = "header"
    elif state == "variable_description_line" or state == "variable_description_interruption":
        state = "variable_description_interruption"
    elif re.search(r"\w+", line) is not None and \
            (state == "dashes_after_header" or state == "data_row"):
        state = "data_row"
    else:
        state = "other"
    print(state, end="\n\n")
