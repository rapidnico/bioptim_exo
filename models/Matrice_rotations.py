import bioviz
from matplotlib import pyplot as plt
from enums import Models
import numpy as np
import biorbd
from scipy import stats


def matrice_1_rotation(nbSolide=int,
                       n_frames=int,
                       biorbd_model_path=str,
                       amplitude_finale = int
                       ):
    biorbd_model = biorbd.Model(biorbd_model_path)
    Q0 = np.zeros((biorbd_model.nbQ()))

    all_q = np.zeros((biorbd_model.nbQ(), n_frames))
    all_q[11, :] = np.linspace(0, amplitude_finale, n_frames)
    q_elevation = all_q[11, :]

    for i in range(0, len(q_elevation)):
        all_q[0, i] = q_elevation[i] * (-0.63355 / 2.61799)  # same as unrot_scap_R2
        all_q[1, i] = q_elevation[i] * (0.322013 / 3.14159)  # same as unrot_scap_R3
        all_q[2, i] = q_elevation[i] * (-0.322013 / 3.14159)
        all_q[3, i] = q_elevation[i] * (0.633555 / 2.61799)
        all_q[4, i] = q_elevation[i] * (-0.128282 / 2.61799)  # same as unrot_hum_R2
        all_q[5, i] = q_elevation[i] * (1.03673 / 2.61799)  # same as unrot_hum_R3
        all_q[6, i] = q_elevation[i] * (0.46603 / 2.61799)  # same as unrot_hum_R1
        all_q[7, i] = q_elevation[i] * (-0.46603 / 2.61799)
        all_q[8, i] = q_elevation[i] * (-1.03673 / 2.61799)
        all_q[9, i] = q_elevation[i] * (0.128282 / 2.61799)

    Table_euler_angle = np.zeros((n_frames, 3))
    for j in range(0, n_frames):
        Q = all_q[:, j]

        T0 = biorbd_model.globalJCS(Q0, nbSolide).rot().transpose().to_array()

        T = biorbd_model.globalJCS(Q, nbSolide).rot().to_array()

        Matrice_rotation = np.matmul(T0, T)
        M = biorbd.Rotation(Matrice_rotation[0, 0], Matrice_rotation[0, 1], Matrice_rotation[0, 2],
                            Matrice_rotation[1, 0], Matrice_rotation[1, 1], Matrice_rotation[1, 2],
                            Matrice_rotation[2, 0], Matrice_rotation[2, 1], Matrice_rotation[2, 2])
        Euler_angles = biorbd.Rotation.toEulerAngles(M, 'xyz').to_array()
        Table_euler_angle[j, :] = Euler_angles

    return Euler_angles, Table_euler_angle, Matrice_rotation




def matrice_2_rotation(nbSolide1=int,
                       nbSolide2 = int,
                       n_frames=int,
                       biorbd_model_path=str,
                       amplitude_finale=int
                       ):


    biorbd_model = biorbd.Model(biorbd_model_path)
    sequence = 'xyz'
    Q0 = np.zeros((biorbd_model.nbQ()))

    all_q = np.zeros((biorbd_model.nbQ(), n_frames))
    all_q[11, :] = np.linspace(0, amplitude_finale, n_frames)
    q_elevation = all_q[11, :]

    for i in range(0, len(q_elevation)):
        all_q[0, i] = q_elevation[i] * (-0.63355 / 2.61799)  # same as unrot_scap_R2
        all_q[1, i] = q_elevation[i] * (0.322013 / 3.14159)  # same as unrot_scap_R3
        all_q[2, i] = q_elevation[i] * (-0.322013 / 3.14159)
        all_q[3, i] = q_elevation[i] * (0.633555 / 2.61799)
        all_q[4, i] = q_elevation[i] * (-0.128282 / 2.61799)  # same as unrot_hum_R2
        all_q[5, i] = q_elevation[i] * (1.03673 / 2.61799)  # same as unrot_hum_R3
        all_q[6, i] = q_elevation[i] * (0.46603 / 2.61799)  # same as unrot_hum_R1
        all_q[7, i] = q_elevation[i] * (-0.46603 / 2.61799)
        all_q[8, i] = q_elevation[i] * (-1.03673 / 2.61799)
        all_q[9, i] = q_elevation[i] * (0.128282 / 2.61799)

    Table_euler_angle = np.zeros((n_frames, 3))
    Table_euler_angle1 = np.zeros((n_frames, 3))
    Table_euler_angle_2 = np.zeros((n_frames, 3))

    for j in range(0, n_frames):
        Q = all_q[:, j]
        T0_1 = biorbd_model.globalJCS(Q0, nbSolide1).rot().transpose().to_array()

        T_1 = biorbd_model.globalJCS(Q, nbSolide1).rot().to_array()

        Matrice_rotation1 = np.matmul(T0_1, T_1)
        M1 = biorbd.Rotation(Matrice_rotation1[0, 0], Matrice_rotation1[0, 1], Matrice_rotation1[0, 2],
                            Matrice_rotation1[1, 0], Matrice_rotation1[1, 1], Matrice_rotation1[1, 2],
                            Matrice_rotation1[2, 0], Matrice_rotation1[2, 1], Matrice_rotation1[2, 2])
        Euler_angles1 = biorbd.Rotation.toEulerAngles(M1, sequence).to_array()
        Table_euler_angle1[j, :] = Euler_angles1



        T0_2 = biorbd_model.globalJCS(Q0, nbSolide2).rot().to_array()

        # on recupere la matrice rotation inverse de R2 pour l,appliquer a la matrice de R3


        T0_2_base = np.matmul(T0_1, T0_2)  # on annule la rotation du prenier solide sur le 2eme et on sort la matrice de rotation originale du bioMod
        T0_2_base_inverse = T0_2_base.transpose()

        # On effectue les meme operations sur la matrice
        T_2 = biorbd_model.globalJCS(Q, nbSolide2).rot().to_array()

        T_1_undo = T_1.transpose()

        T_1_base = np.matmul(T_1_undo, T_2)

        Matrice_rotation_2 = np.matmul(T0_2_base_inverse, T_1_base)
        M2 = biorbd.Rotation(Matrice_rotation_2[0, 0], Matrice_rotation_2[0, 1], Matrice_rotation_2[0, 2],
                            Matrice_rotation_2[1, 0], Matrice_rotation_2[1, 1], Matrice_rotation_2[1, 2],
                            Matrice_rotation_2[2, 0], Matrice_rotation_2[2, 1], Matrice_rotation_2[2, 2])
        Euler_angles_2 = biorbd.Rotation.toEulerAngles(M2, sequence).to_array()
        Table_euler_angle_2[j, :] = Euler_angles_2
    return Table_euler_angle1, Table_euler_angle_2





def matrice_3_rotation(nbSolide1=int,
                       nbSolide2 = int,
                       nbSolide3 = int,
                       n_frames=int,
                       biorbd_model_path=str,
                       amplitude_finale=int
                       ):


    biorbd_model = biorbd.Model(biorbd_model_path)
    sequence = 'xyz'
    Q0 = np.zeros((biorbd_model.nbQ()))

    all_q = np.zeros((biorbd_model.nbQ(), n_frames))
    all_q[11, :] = np.linspace(0, amplitude_finale, n_frames)
    q_elevation = all_q[11, :]

    for i in range(0, len(q_elevation)):
        all_q[0, i] = q_elevation[i] * (-0.63355 / 2.61799)  # same as unrot_scap_R2
        all_q[1, i] = q_elevation[i] * (0.322013 / 3.14159)  # same as unrot_scap_R3
        all_q[2, i] = q_elevation[i] * (-0.322013 / 3.14159)
        all_q[3, i] = q_elevation[i] * (0.633555 / 2.61799)
        all_q[4, i] = q_elevation[i] * (-0.128282 / 2.61799)  # same as unrot_hum_R2
        all_q[5, i] = q_elevation[i] * (1.03673 / 2.61799)  # same as unrot_hum_R3
        all_q[6, i] = q_elevation[i] * (0.46603 / 2.61799)  # same as unrot_hum_R1
        all_q[7, i] = q_elevation[i] * (-0.46603 / 2.61799)
        all_q[8, i] = q_elevation[i] * (-1.03673 / 2.61799)
        all_q[9, i] = q_elevation[i] * (0.128282 / 2.61799)

    Table_euler_angle = np.zeros((n_frames, 3))
    Table_euler_angle1 = np.zeros((n_frames, 3))
    Table_euler_angle_2 = np.zeros((n_frames, 3))
    Table_euler_angle_3 = np.zeros((n_frames, 3))

    for j in range(0, n_frames):
        Q = all_q[:, j]
        T0_1 = biorbd_model.globalJCS(Q0, nbSolide1).rot().transpose().to_array()

        T_1 = biorbd_model.globalJCS(Q, nbSolide1).rot().to_array()

        Matrice_rotation1 = np.matmul(T0_1, T_1)
        M1 = biorbd.Rotation(Matrice_rotation1[0, 0], Matrice_rotation1[0, 1], Matrice_rotation1[0, 2],
                            Matrice_rotation1[1, 0], Matrice_rotation1[1, 1], Matrice_rotation1[1, 2],
                            Matrice_rotation1[2, 0], Matrice_rotation1[2, 1], Matrice_rotation1[2, 2])
        Euler_angles1 = biorbd.Rotation.toEulerAngles(M1, sequence).to_array()
        Table_euler_angle1[j, :] = Euler_angles1

        ##### Solide 2

        T0_2 = biorbd_model.globalJCS(Q0, nbSolide2).rot().to_array()

        # on recupere la matrice rotation inverse du 1er solide pour l'appliquer a la matrice de R3

        T0_2_base = np.matmul(T0_1, T0_2)  # on annule la rotation du prenier solide sur le 2eme et on sort la matrice de rotation originale du bioMod
        T0_2_base_inverse = T0_2_base.transpose()

        # On effectue les meme operations sur la matrice
        T_2 = biorbd_model.globalJCS(Q, nbSolide2).rot().to_array()

        T_1_undo = T_1.transpose()

        T_2_base = np.matmul(T_1_undo, T_2)

        Matrice_rotation_2 = np.matmul(T0_2_base_inverse, T_2_base)
        M2 = biorbd.Rotation(Matrice_rotation_2[0, 0], Matrice_rotation_2[0, 1], Matrice_rotation_2[0, 2],
                            Matrice_rotation_2[1, 0], Matrice_rotation_2[1, 1], Matrice_rotation_2[1, 2],
                            Matrice_rotation_2[2, 0], Matrice_rotation_2[2, 1], Matrice_rotation_2[2, 2])
        Euler_angles_2 = biorbd.Rotation.toEulerAngles(M2, sequence).to_array()
        Table_euler_angle_2[j, :] = Euler_angles_2


        ##### Solide 3

        T0_3 = biorbd_model.globalJCS(Q0, nbSolide3).rot().to_array()

        # On se sert des matrices inverses dans le biomode des bases du 1er et 2eme solide : T0_2_base_inverse et T0_1

        T0_3_base_1 = np.matmul(T0_1, T0_3)
        T0_3_base_2 = np.matmul(T0_2_base_inverse, T0_3_base_1)

        # On fait de meme pour la matrice qui subit la rotation avec T_1 et T_2_base qu'on doit inverser

        T_3 = biorbd_model.globalJCS(Q0, nbSolide3).rot().to_array()
        T_1_undo = T_1.transpose()
        T_2_base_undo = T_2_base.transpose()

        T_3_base_1 = np.matmul(T_1_undo, T_3)
        T_3_base_2 = np.matmul(T_2_base_undo, T_3_base_1)

        Matrice_rotation_3 = np.matmul(T0_2_base_inverse, T_2_base)
        M3 = biorbd.Rotation(Matrice_rotation_3[0, 0], Matrice_rotation_3[0, 1], Matrice_rotation_3[0, 2],
                            Matrice_rotation_3[1, 0], Matrice_rotation_3[1, 1], Matrice_rotation_3[1, 2],
                            Matrice_rotation_3[2, 0], Matrice_rotation_3[2, 1], Matrice_rotation_3[2, 2])
        Euler_angles_3 = biorbd.Rotation.toEulerAngles(M3, sequence).to_array()
        Table_euler_angle_3[j, :] = Euler_angles_3




    return Table_euler_angle1, Table_euler_angle_2, Table_euler_angle_3