"""
@File   : util_data.py
@Author : Siyuan Chen
@Date   : 2025/3/10
@Desc   : 
"""
import os
import time
import cooler
# import hicstraw
from tqdm.auto import tqdm


def write_bedpe(loop_list,file_name,resolution):
    with open(file_name, 'w') as f:
        for loop in loop_list:
            score = loop[3]
            chr_num = loop[0]
            if loop[1] > loop[2]:
                start_1 = int(loop[2]) * resolution
                end_1 =  int(loop[2] + 1) * resolution
                start_2 =  int(loop[1]) * resolution
                end_2 = int(loop[1] + 1) * resolution
            else:
                start_1 =  int(loop[1]) * resolution
                end_1 =  int(loop[1] + 1) * resolution
                start_2 =  int(loop[2]) * resolution
                end_2 =  int(loop[2] + 1) * resolution

            line = "{}\t{}\t{}\tchr{}\t{}\t{}\t{}\n".format(chr_num, start_1, end_1, chr_num, start_2, end_2, score)
            f.write(line)

def read_bedpe(bedpe_path):
    """
    read the loop annotations from a .bedpe file

    Args:
        bedpe_path (string): path to a .bedpe loop annotation file

    Returns:
         list of list, list of loop positions
    """
    lines = []

    with open(bedpe_path, 'r') as bedpe_file:
        # select the first 6 entries: ['chr1', 'x_start', 'x_end', 'chr1', 'y_start', 'y_end']
        # all entries are strings instead of integers
        counter = 0
        for line in bedpe_file:
            counter +=1
            line_info = line.rstrip().split()[:7]
            # Threshold Score for IMR90 and GM12878
            lines.append(line_info[0:6])

    print("Before {}".format(counter))
    # sort the list by the chromosome name
    lines = sorted(lines, key=lambda x: x[0])
    print("After {}".format(len(lines)))
    print(f'[info] {len(lines)} loop annotations read from {bedpe_path}')
    return lines


def load_loops(loop_path, chr_name):
    """
    select a list of loop coordinates from a chromosome

    Args:
        loop_path (string):, path to a loop annotation file
        chr_name (string): name of the chromosomes, e.g., 'chr1'

    Returns:
        list of list, list of loop positions
    """
    all_loops = read_bedpe(loop_path)
    loops = []

    for loop in all_loops:
        if loop[0] != chr_name:
            continue

        x = int((int(loop[1]) + int(loop[2])) * 0.5)
        y = int((int(loop[4]) + int(loop[5])) * 0.5)

        loops.append([loop[0], x, y])
        loops.append([loop[0], y, x])

    print(f'[info] {int(len(loops))} loop annotations loaded from {chr_name}')

    return loops


def read_hic(hic_path, chr_name, load_all=True,
             x_start=-1, x_end=-1, y_start=-1, y_end=-1, resolution=10000):
    """
        load the contact map from a .hic file

        Args:
            hic_path (string): path to a  .hic/.mcool/.cool file
            chr_name (string): name of the chromosome, e.g., 'chr1'
            load_all (bool): whether to load the entire contact matrix or not, True by default
            x_start (int): starting x-coordinate (row index) of the chromosome
            x_end (int): ending x-coordinate (row index) of the contact map
            y_start (int): starting y-coordinate (col index) of the chromosome
            y_end (int): ending y-coordinate (col index) of the contact map
            resolution (int): resolution of the contact map, 100kb by default

        Returns:
             numpy array: contact map as a matrix
        """

    chr_num = chr_name[3:]

    raise Exception(f'[error] .hic file loading is not implemented at this moment')

    # hic_obj = hicstraw.HiCFile(hic_path)
    #
    # mzd = hic_obj.getMatrixZoomData(f'{chr_num}', f'{chr_num}',
    #                                 'observed', 'NONE', 'BP', resolution)
    #
    # print(f'[info] loading contact map from .hic file')
    #
    # chr_dict = {c.name: c for c in hic_obj.getChromosomes()}  # chr number as a str: chr object
    # print(chr_dict)
    #
    # if load_all:
    #     # load the entire contact map
    #     print(f'[info] loading the entire contact map')
    #     print(chr_dict[f'{chr_num}'].length)
    #     # chr_len = (chr_dict[f'{chr_num}'].length // resolution + 1) * resolution
    #     # print(f'range up to: {chr_len}')
    #     contact_list = mzd.getRecords(0, chr_dict[f'{chr_num}'].length, 0, chr_dict[f'{chr_num}'].length)
    #     print('check point 1')
    #     print(len(contact_list))
    #     print(contact_list[0])
    #
    #     data = np.array([float(contact.counts) for contact in contact_list])
    #     row = np.array([contact.binX // resolution for contact in contact_list])
    #     col = np.array([contact.binY // resolution for contact in contact_list])
    #     print(row[:5])
    #     print(col[:5])
    #
    #     matrix = csr_matrix(data, (row, col))

        # max_row_idx = max(contact.binX // resolution for contact in contact_list)
        # max_col_idx = max(contact.binY // resolution for contact in contact_list)

        # matrix = np.zeros((max_row_idx + 1, max_col_idx + 1))
        # print('check point 2')
        #
        # bin_x = np.array([contact.binX for contact in contact_list]) // resolution
        # bin_y = np.array([contact.binY for contact in contact_list]) // resolution
        # counts = np.array([contact.counts for contact in contact_list])
        # print('check point 3')
        #
        # # since we are loading the entire contact map, we can assume it to be symmetric
        # matrix[bin_x, bin_y] = counts
        # matrix[bin_y, bin_x] = counts

    # else:
    #     # load a part of the contact map
    #     print(f'[info] retrieving data within range:\n'
    #           f'row ({x_start}, {x_end}), column: ({y_start}, {y_end})')
    #
    #     matrix = mzd.getRecordsAsMatrix(x_start, x_end, y_start, y_end)
    #
    # print(f'[debug] contact map data retrieved')

    # return matrix


def read_cool(cool_path, chr_name, load_all=True,
              x_start=-1, x_end=-1, y_start=-1, y_end=-1, resolution=10000, balance=False):
    """
    load the contact map from a .mcool/.cool file

    Args:
        cool_path (string): path to a  .mcool/.cool file
        chr_name (string): name of the chromosome, e.g., 'chr1'
        load_all (bool): whether to load the entire contact matrix or not, True by default
        x_start (int): starting x-coordinate (row index) of the chromosome
        x_end (int): ending x-coordinate (row index) of the contact map
        y_start (int): starting y-coordinate (col index) of the chromosome
        y_end (int): ending y-coordinate (col index) of the contact map
        resolution (int): resolution of the contact map, 10kb by default
        balance (bool): whether to use a balanced matrix, False by default

    Returns:
         numpy array: contact map as a matrix
    """

    # read as .cool or .mcool
    if '.cool' in cool_path:
        raw_mat = cooler.Cooler(cool_path)
    elif '.mcool' in cool_path:
        print(cool_path)
        raw_mat = cooler.Cooler(f'{cool_path}::/resolutions/{resolution}')
    else:
        raise Exception(f'[error] unrecognized file extension in: {cool_path}')

    # load all the contact map or only a part
    print(f'[info] using balanced matrix: {bool(balance)}')

    if load_all:
        try:
            matrix = raw_mat.matrix(balance=balance).fetch(chr_name)
        except ValueError:
            matrix = raw_mat.matrix(balance=balance).fetch(chr_name[3:])
    else:
        # cooler can only load a symmetric matrix
        # hence, we load a larger matrix and return only a part of it
        start = min(x_start, y_start)
        end = max(x_end, y_end)

        try:
            large_matrix = raw_mat.matrix(balance=balance).fetch((chr_name, start, end))
        except ValueError:
            large_matrix = raw_mat.matrix(balance=balance).fetch((chr_name[3:], start, end))

        matrix = large_matrix[(x_start - start) // resolution:(x_end - start) // resolution,
                              (y_start - start) // resolution:(y_end - start) // resolution]

    return matrix


def load_contact_map(cm_path, chr_name, load_all=True,
                     x_start=-1, x_end=-1, y_start=-1, y_end=-1, resolution=10000, balance=False):
    """
    wrapper function around read_hic and read_cool
    to load the square contact map from a .hic/.mcool/.cool file

    Args:
        cm_path (string): path to a  .hic/.mcool/.cool file
        chr_name (string): name of the chromosome, e.g., 'chr1'
        load_all (bool): whether to load the entire contact matrix or not, True by default
        x_start (int): starting x-coordinate (row index) of the chromosome
        x_end (int): ending x-coordinate (row index) of the contact map
        y_start (int): starting y-coordinate (col index) of the chromosome
        y_end (int): ending y-coordinate (col index) of the contact map
        resolution (int): resolution of the contact map, 10kb by default
        balance (bool): whether to use a balanced matrix, False by default

    Returns:
         numpy array: contact map as a matrix
    """

    if '.hic' in cm_path:
        matrix = read_hic(cm_path, chr_name, load_all, x_start, x_end, y_start, y_end, resolution)

    elif 'cool' in cm_path:
        matrix = read_cool(cm_path, chr_name, load_all, x_start, x_end, y_start, y_end, resolution, balance)

    else:
        raise Exception(f'[error] unrecognized file extension in: {cm_path}')

    return matrix


if __name__ == '__main__':
    start_time = time.perf_counter()
    # mat = read_hic('../data/GSE63525_K562_combined.hic', 'chr1')
    mat = load_contact_map('../data/k562-hic.cool', 'chr1', load_all=False,
                           x_start=24080000, x_end=29200000, y_start=24050000, y_end=29170000)
    a = mat[10:100, 10:100]
    # mat = load_contact_map('../data/k562-hic.mcool', 'chr1')
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.6f} seconds")
    print(mat.shape)
    print(a.shape)

