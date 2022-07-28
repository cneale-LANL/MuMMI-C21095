import numpy as np
import scipy
import scipy.spatial

def nearest_neighbor_inds_kd_tree(AA, XX0, kk):
    tree = scipy.spatial.cKDTree(AA)
    ii = tree.query(XX0, k=kk, workers=-1)[1]
    # R indexing starts from 1
    return ii + 1
