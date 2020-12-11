import sys

start_line = "----- TOP OF INPUT FILE FOR BADER'S AIMPAC PROGRAM -----\n"
end_line = "----- END OF INPUT FILE FOR BADER'S AIMPAC PROGRAM -----\n"

in_section = False
for line in open(sys.argv[1]):
    if line == end_line:
        break
    if in_section:
        print(line, end="")
    if line == start_line:
        in_section = True
