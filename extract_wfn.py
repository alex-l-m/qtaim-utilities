import sys

# Check for a "--functional" argument
if len(sys.argv) > 2:
    if sys.argv[2] == "--functional":
        functional = sys.argv[3]
    else:
        raise ValueError("Unknown argument: " + sys.argv[2])
else:
    functional = None

start_line = "----- TOP OF INPUT FILE FOR BADER'S AIMPAC PROGRAM -----\n"
end_line = "----- END OF INPUT FILE FOR BADER'S AIMPAC PROGRAM -----\n"

in_section = False
for line in open(sys.argv[1]):
    if line == end_line:
        break
    if in_section:
        if functional is None:
            print(line, end="")
        else:
            print(line.replace("NUCLEI", f"NUCLEI {functional}"), end="")
    if line == start_line:
        in_section = True
