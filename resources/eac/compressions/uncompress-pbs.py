uncompress

# read header bytes
file_header = read(2)
output_size = read(1) << 16 | read(2)
num_entries = read(1)

# 1. read the tree and build the tables
build_tree()

