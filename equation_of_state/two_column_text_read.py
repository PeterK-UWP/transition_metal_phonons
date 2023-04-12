import numpy as np


def two_column_text_read(file_name):
    try:
        file = open(file_name)
    except OSError as error:
        print(error)
        return

    content = file.readlines()
    data = np.zeros([2, (len(content))], float)
    n = 0

    for line in content:
        elements = line.split()
        data[0, n] = float(elements[0])
        data[1, n] = float(elements[1])
        n += 1

    return data