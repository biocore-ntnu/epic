from sys import argv

bed_file = argv[1]

for line in open(bed_file):
    chromosome, start, end, name, score, strand = line.split()
    start, end = int(start), int(end)
    out = "\t".join(str(s)
                    for s in [chromosome, start, start + 1, chromosome, end,
                              end + 1, name, score, strand, strand])
    print(out)
