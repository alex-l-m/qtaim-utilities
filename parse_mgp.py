import sys
import re
import csv

inpath = sys.argv[1]

# Initialize table
table = [("critical_point", "property", "value")]

critical_point = None
for line in open(inpath):
    # Critical point descriptions are ended by empty lines
    if line.strip() == "":
        critical_point = None

    # Critical point descriptions begin by specifying the critical point number
    cp_match = re.match(r"CP# (\d+)", line)
    if cp_match is not None:
        critical_point = int(cp_match.group(1))
        # Remove the critical point marker from the line so that it doesn't get
        # considered as part of a variable name
        line = line[cp_match.end():]

    if critical_point is not None:

        # Detect three-dimensional vectors
        vector_variables = set()
        for variable, x, y, z in re.findall(\
                r"([a-zA-Z-_][^=]+[a-zA-Z-_]) *= *([+\-0-9.Ee]+) +([+\-0-9.Ee]+) +([+\-0-9.Ee]+)",\
                line):
            # Record that this is a vector variable so it is not added later as a scalar
            vector_variables.add(variable)
            # Add each coordinate individually to the table
            table.append((critical_point, variable + "_x", float(x)))
            table.append((critical_point, variable + "_y", float(y)))
            table.append((critical_point, variable + "_z", float(z)))

        # Detect scalars
        for variable, value in re.findall(r"([a-zA-Z-_][^=]+[a-zA-Z-_]) *= *([+\-0-9.Ee]+)", line):
            # Only add to the table if we didn't already add this as a vector
            # variable
            if variable not in vector_variables:
                table.append((critical_point, variable, float(value)))

        # Detect type
        type_match = re.search(r"Type += +\(3, *([\-+0-9]+)\) +(\w+)(( +\w+)+)", line)
        if type_match is not None:
            total_sign, type_string, atoms = type_match.groups()[:3]
            table.append((critical_point, "total_sign", int(total_sign)))
            table.append((critical_point, "type_string", type_string))
            for i, atom in enumerate(atoms.split()):
                table.append((critical_point, "atom_{}".format(i+1), atom.strip()))

# Send CSV table to stdout
csv.writer(sys.stdout).writerows(table)
