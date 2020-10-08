import re

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
