import numpy as np

def get_col(file_name, *cols):
    """ Function `get_col` returns a tuple containing the conents of each of the
    specified columns in a textfile containing space-separated floats.
    """
    with open(file_name) as fp: s = fp.read()
    lines = s.split("\n")
    def get_word(line):
        words = [word for word in line.split(" ") if word != ""]
        return [words[col] for col in cols]

    res =  map(np.array, zip(*[map(float, get_word(line))
                               for line in lines if line != "" and
                               line[0] != "#"]))
    if len(res) == 1: return res[0]
    return res
