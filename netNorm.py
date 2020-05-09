"""Main function of netNorm for the paper: Estimation of Connectional Brain Templates using Selective Multi
View Network Normalization
    Details can be found in:
    (1) the original paper https://www.ncbi.nlm.nih.gov/pubmed/31622839
        Salma Dhifallah, and Islem Rekik.
    ---------------------------------------------------------------------
    This file contains the implementation of three key steps of our netNorm framework:
        netNorm(sourceGraph, number_of_subjects, number_of_regions)
                Inputs:
                        sourceGraph: (n × m x t x t) matrix stacking the source graphs of all subjects
                                     n the total number of views
                                     m the number of subjects
                                     t number of regions
                Output:
                        CBT:         (t x t) matrix representing the cortical brain template
    (2) Dependencies: please install the following libraries:
        - matplotlib
        - numpy
        - snfpy (https://github.com/rmarkello/snfpy)
        - scikitlearn
    ---------------------------------------------------------------------
    Copyright 2019 Salma Dhifallah, Sousse University.
    Please cite the above paper if you use this code.
    All rights reserved.
    """

import numpy as np
import snf
from sklearn import preprocessing
import matplotlib.pyplot as plt

def netNorm(v, nbr_of_sub, nbr_of_regions):
    nbr_of_feat = int((np.square(nbr_of_regions) - nbr_of_regions) / 2)

    def minmax_sc(x):
        min_max_scaler = preprocessing.MinMaxScaler()
        x = min_max_scaler.fit_transform(x)
        return x

    def upper_triangular():
        All_subj = np.zeros((nbr_of_sub, len(v), nbr_of_feat))
        for i in range(len(v)):
            for j in range (nbr_of_sub):

                subj_x = v[i, j, :, :]
                subj_x = np.reshape(subj_x, (nbr_of_regions, nbr_of_regions))
                subj_x = minmax_sc(subj_x)
                subj_x = subj_x[np.triu_indices(nbr_of_regions, k = 1)]
                subj_x = np.reshape(subj_x, (1, 1, nbr_of_feat))
                All_subj[j, i, :] = subj_x

        return All_subj

    def distances_inter(All_subj):
        theta = 0
        distance_vector = np.zeros(1)
        distance_vector_final = np.zeros(1)
        x = All_subj
        for i in range(nbr_of_feat): #par rapport ll number of ROIs
            ROI_i = x[:, :, i]
            ROI_i = np.reshape(ROI_i, (nbr_of_sub, nbr_of_views)) #1,3
            for j in range(nbr_of_sub):
                subj_j = ROI_i[j:j+1, :]
                subj_j = np.reshape(subj_j, (1, nbr_of_views))
                for k in range(nbr_of_sub):
                    if k != j:
                        subj_k = ROI_i[k:k+1, :]
                        subj_k = np.reshape(subj_k, (1, nbr_of_views))

                        for l in range(nbr_of_views):
                            if l ==0:
                                distance_euclidienne_sub_j_sub_k = np.square(subj_k[:, l:l+1] - subj_j[:, l:l+1])
                            else:
                                distance_euclidienne_sub_j_sub_k = distance_euclidienne_sub_j_sub_k + np.square(subj_k[:, l:l+1] - subj_j[:, l:l+1])

                            theta +=1
                if j ==0:
                    distance_vector = np.sqrt(distance_euclidienne_sub_j_sub_k)
                else:
                    distance_vector = np.concatenate((distance_vector, np.sqrt(distance_euclidienne_sub_j_sub_k)), axis=0)

            if i ==0:
                distance_vector_final = distance_vector
            else:
                distance_vector_final = np.concatenate((distance_vector_final, distance_vector), axis=1)

        print(theta)
        return distance_vector_final


    def minimum_distances(distance_vector_final):
        x = distance_vector_final

        for i in range(nbr_of_feat):
            for j in range(nbr_of_sub):
                minimum_sub = x[j:j+1, i:i+1]
                minimum_sub = float(minimum_sub)
                for k in range(nbr_of_sub):
                    if k != j:
                        local_sub = x[k:k+1, i:i+1]
                        local_sub = float(local_sub)
                        if local_sub < minimum_sub:
                            general_minimum = k
                            general_minimum = np.array(general_minimum)
            if i == 0:
                final_general_minimum = np.array(general_minimum)
            else:
                final_general_minimum = np.vstack((final_general_minimum, general_minimum))

        final_general_minimum = np.transpose(final_general_minimum)

        return final_general_minimum

    def new_tensor(final_general_minimum, All_subj):
        y = All_subj
        x = final_general_minimum
        for i in range(nbr_of_feat):
            optimal_subj = x[:, i:i+1]
            optimal_subj = np.reshape(optimal_subj, (1))
            optimal_subj = int(optimal_subj)
            if i ==0:
                final_new_tensor = y[optimal_subj: optimal_subj+1, :, i:i+1]
            else:
                final_new_tensor = np.concatenate((final_new_tensor, y[optimal_subj: optimal_subj+1, :, i:i+1]), axis=2)

        return final_new_tensor


    def make_sym_matrix(nbr_of_regions, feature_vector):

        nbr_of_regions = nbr_of_regions
        feature_vector = feature_vector
        my_matrix = np.zeros([nbr_of_regions,nbr_of_regions], dtype=np.double)

        my_matrix[np.triu_indices(nbr_of_regions, k=1)] = feature_vector
        my_matrix = my_matrix + my_matrix.T
        my_matrix[np.diag_indices(nbr_of_regions)] = 0

        return my_matrix

    def re_make_tensor(final_new_tensor, nbr_of_regions):
        x =final_new_tensor
        x = np.reshape(x, (nbr_of_views, nbr_of_feat))
        for i in range (nbr_of_views):
            view_x = x[i, :]
            view_x = np.reshape(view_x, (1, nbr_of_feat))
            view_x = make_sym_matrix(nbr_of_regions, view_x)
            view_x = np.reshape(view_x, (1, nbr_of_regions, nbr_of_regions))
            if i ==0:
                tensor_for_snf = view_x
            else:
                tensor_for_snf = np.concatenate((tensor_for_snf, view_x), axis=0)
        return tensor_for_snf

    def create_list(tensor_for_snf):
        x =tensor_for_snf
        for i in range(nbr_of_views):
            view = x[i, :, :]
            view = np.reshape(view, (nbr_of_regions, nbr_of_regions))
            list = [view]
            if i ==0:
                list_final = list
            else:
                list_final = list_final +list
        return list_final

    def cross_subjects_cbt(fused_network, nbr_of_exemples):
        final_cbt = np.zeros((nbr_of_exemples, nbr_of_feat))
        x = fused_network
        x = x[np.triu_indices(nbr_of_regions, k=1)]
        x = np.reshape(x, (1, nbr_of_feat))
        for i in range(nbr_of_exemples):
            final_cbt[i, :] = x


        return final_cbt

    Upp_trig = upper_triangular()
    Dis_int = distances_inter(Upp_trig)
    Min_dis = minimum_distances(Dis_int)
    New_ten = new_tensor(Min_dis, Upp_trig)
    Re_ten = re_make_tensor(New_ten, nbr_of_regions)
    Cre_lis = create_list(Re_ten)
    fused_network = snf.snf((Cre_lis), K=20)
    fused_network = minmax_sc(fused_network)
    np.fill_diagonal(fused_network, 0)
    fused_network = np.array(fused_network)
    return fused_network

nbr_of_sub = int(input('Please select the number of subjects: '))
nbr_of_regions = int(input('Please select the number of regions: '))
nbr_of_views= int(input('Please select the number of views: '))

v = np.random.rand(nbr_of_views, nbr_of_sub,nbr_of_regions, nbr_of_regions)
A = netNorm(v, nbr_of_sub, nbr_of_regions)
print(A)

mx = A.max()
mn = A.min()
print(mx)
print(mn)
plt.pcolor(A, vmin=mn, vmax=mx)
plt.imshow(A)
plt.show()
