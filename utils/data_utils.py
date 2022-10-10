import cooler


def read_loops(loop_path):
    """
    read the loop annotations from a .bedpe file
    :param loop_path: path, path to a loop annotation file
    :return: list of list, list of loop positions
    """
    loops = []
    with open(loop_path, 'r') as loop_file:
        for line in loop_file:
            loops.append(line.split('\t')[:6])
            # loops[i] = ['chr1', 'x_start', 'x_end', 'chr1', 'y_start', 'y_end']
            # Note that all loop[i][j] are strings instead of integers

    return loops


def load_loops(loop_path, chr_num, start=None, end=None):
    """
    select a list of loop positions from a chromosome
    :param loop_path: path, path to a loop annotation file
    :param chr_num: string, number of the chromosomes, e.g., '1'
    :param start: int, mark the start of the range to find loops, None by default
    :param end: int, mark the end of the range to find loop, None by default
    :return: list of list, list of loop positions, e.g., [['chr1', 21000, 21200],]
    """
    raw_loops = read_loops(loop_path)
    # skip header line
    if raw_loops[0][1] == 'x1':
        raw_loops = raw_loops[1:]
    # check for empty file
    if len(raw_loops) == 0:
        raise Exception('Input File Error: the provided file contains no loop calls')
    # determine the format of the chromosome representation
    if len(raw_loops[0][0]) >= 4:  # format: 'chr1', 'chr12', etc.
        chr_num = f'chr{chr_num}'
    # else: the format is '1', '12', etc.

    loops = []

    if start is not None and end is not None:
        # only select the loops within the range specified by the input
        if start > end:
            raise Exception('Incorrect Arguments: the starting coordinate must be smaller than the ending one')

        for loop in raw_loops:
            if loop[0] != chr_num:
                continue

            x = int((int(loop[1]) + int(loop[2])) * 0.5)
            y = int((int(loop[4]) + int(loop[5])) * 0.5)

            if start < x < end and start < y < end:
                loops.append([loop[0], x, y])
                loops.append([loop[0], y, x])

        else:
            # load all the loops on this chromosome
            for loop in raw_loops:
                if loop[0] != chr_num:
                    continue

                x = int((int(loop[1]) + int(loop[2])) * 0.5)
                y = int((int(loop[4]) + int(loop[5])) * 0.5)

                loops.append([loop[0], x, y])
                loops.append([loop[0], y, x])

    print(f'{int(len(loops))} loops loaded from {chr_num}')
    return loops


def load_contact_map(map_path, chr_num, start=None, end=None, resolution=10000):
    """
    load a contact map from a .cool or .mcool file
    :param map_path: path, path to a .cool/.mcool file
    :param chr_num: string, number of the chromosome, e.g., '1'
    :param start: int, starting coordinate of the contact map to be loaded, None by default
    :param end: int, ending coordinate of the contact map to be loaded, None by default
    :param resolution: int, resolution of the contact map to be loaded, 10kb by default
    :return: array, the contact map as an array
    """

    if '.mcool' in map_path:
        raw_mat = cooler.Cooler(f'{map_path}::/resolutions/{resolution}')
        matrix = raw_mat.matrix(balance=False).fetch('{}'.format(chr_num))
    elif '.cool' in map_path:
        print('Warning: Note that the program will load the contact map stored in the user-specified .cool file,'
              'which may or may not fit the specified the resolution.')
        raw_mat = cooler.Cooler(map_path)
        matrix = raw_mat.matrix(balance=False).fetch(chr_num)  # chr_num is a string, e.g., '1'
    else:
        raise Exception('Incorrect Format: the format of the contact map must be .cool or .mcool')

    if start is not None and end is not None:
        # load part of the contact map
        if start > end:
            raise Exception('Incorrect Arguments: the starting coordinate must be smaller than the ending one')
        return matrix[start:end, start:end]

    elif start is not None and end is not None:
        # load the entire contact map
        return matrix
