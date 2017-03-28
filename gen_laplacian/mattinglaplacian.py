import argparse, glob, os
import numpy as np
import scipy.sparse as sps
import scipy.ndimage as spi
import scipy.misc as spm
import scipy.io as spio


def getlaplacian1(i_arr: np.ndarray, consts: np.ndarray, epsilon: float = 0.0000001, win_size: int = 1):
    neb_size = (win_size * 2 + 1) ** 2
    h, w, c = i_arr.shape
    n = h
    m = w
    img_size = w * h
    consts = spi.morphology.grey_erosion(consts, footprint=np.ones(shape=(win_size * 2 + 1, win_size * 2 + 1)))

    indsM = np.reshape(np.array(range(img_size)), newshape=(h, w), order='F')
    tlen = (-consts[win_size:-win_size, win_size:-win_size] + 1).sum() * (neb_size ** 2)

    row_inds = np.zeros(tlen)
    col_inds = np.zeros(tlen)
    vals = np.zeros(tlen)
    l = 0
    for j in range(win_size, w - win_size):
        for i in range(win_size, h - win_size):
            if consts[i, j]:
                continue
            win_inds = indsM[i - win_size:i + win_size + 1, j - win_size: j + win_size + 1]
            win_inds = win_inds.ravel(order='F')
            win_i = i_arr[i - win_size:i + win_size + 1, j - win_size: j + win_size + 1, :]
            win_i = win_i.reshape((neb_size, c), order='F')
            win_mu = np.mean(win_i, axis=0).reshape(1, win_size * 2 + 1)
            win_var = np.linalg.inv(
                np.matmul(win_i.T, win_i) / neb_size - np.matmul(win_mu.T, win_mu) + epsilon / neb_size * np.identity(
                    c))

            win_i2 = win_i - win_mu
            tvals = (1 + np.matmul(np.matmul(win_i2, win_var), win_i2.T)) / neb_size

            ind_mat = np.broadcast_to(win_inds, (neb_size, neb_size))
            row_inds[l: (neb_size ** 2 + l)] = ind_mat.ravel(order='C')
            col_inds[l: neb_size ** 2 + l] = ind_mat.ravel(order='F')
            vals[l: neb_size ** 2 + l] = tvals.ravel(order='F')
            l += neb_size ** 2

    vals = vals.ravel(order='F')
    row_inds = row_inds.ravel(order='F')
    col_inds = col_inds.ravel(order='F')
    a_sparse = sps.csr_matrix((vals, (row_inds, col_inds)), shape=(img_size, img_size))

    sum_a = a_sparse.sum(axis=1).T.tolist()[0]
    a_sparse = sps.diags([sum_a], [0], shape=(img_size, img_size)) - a_sparse

    return a_sparse


def im2double(im):
    min_val = np.min(im.ravel())
    max_val = np.max(im.ravel())
    return (im.astype('float') - min_val) / (max_val - min_val)


def reshape_img(in_img, l=512):
    in_h, in_w, _ = img.shape
    if in_h > in_w:
        h2 = l
        w2 = int(in_w * h2 / in_h)
    else:
        w2 = l
        h2 = int(in_h * w2 / in_w)

    return spm.imresize(in_img, (h2, w2))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-dir", "--directory", help="Path to examples")
    parser.add_argument("-out", "--out_directory", help="Path to examples")
    args = parser.parse_args()

    for filename in glob.iglob(os.path.join(args.directory, '*.png')):
        img = spi.imread(filename, mode="RGB")
        img = im2double(reshape_img(img, 700))
        h, w, c = img.shape
        out_filename = 'Input_Laplacian_3x3_1e-7_CSR' + os.path.basename(filename).replace("in", "").replace(".png", "") + ".csv"
        CSR = getlaplacian1(img, np.zeros(shape=(h, w)), 1e-7, 1)
        COO = CSR.tocoo()
        zipped = zip(COO.row + 1, COO.col + 1, COO.data)
        with open(os.path.join(args.out_directory, out_filename), 'w') as out_file:
            for row, col, val in zipped:
                out_file.write("%d,%d,%.15f\n" % (row, col, val))
