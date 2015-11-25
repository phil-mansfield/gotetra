with open("shell.dat") as fp: s = fp.read()
lines = s.split("\n")
tokens = [[t for t in line.split(" ") if len(t) > 0] for line in lines]
snaps = map(int, [line[-1] for line in tokens[:-1]])
ids = map(int, [line[-2].split(")")[-1] for line in tokens[:-1]])
tokens = tokens[1:]
tokens = [line[1:] for (i, line) in enumerate(tokens)]
for i in xrange(len(tokens) - 1):
    tokens[i] = tokens[i][:-1]
for i in xrange(len(tokens)):
    tokens[i][-1] = tokens[i][-1].split(")")[0]
print tokens
