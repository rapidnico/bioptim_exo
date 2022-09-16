import bioviz
from matplotlib import pyplot as plt
from enums import Models
import numpy as np
import biorbd
from scipy import stats
from Matrice_rotations import (
matrice_1_rotation,
matrice_2_rotation,
matrice_3_rotation
)
manually_animate = False

# Load the model
model_name = Models.Stanford_VA_upper_limb_model_0_40.value
biorbd_model = bioviz.Viz(model_name, show_floor=False, show_global_ref_frame=False)


# Create a movement
n_frames = 20
all_q = np.zeros((biorbd_model.nQ, n_frames))
all_q[11, :] = np.linspace(0, 3.14, n_frames)
q_elevation = all_q[11, :]

# Apply the movement to the DoF of the Scapula and Clavicula
q_sterno_clav_R2 = np.zeros(n_frames)
q_sterno_clav_R3 = np.zeros(n_frames)
q_accromioclaviculaire_R1 = np.zeros(n_frames)
q_accromioclaviculaire_R2 = np.zeros(n_frames)
q_accromioclaviculaire_R3 = np.zeros(n_frames)
unrot_scap_R2 = np.zeros(n_frames)
unrot_scap_R3 = np.zeros(n_frames)
unrot_hum_R1 = np.zeros(n_frames)
unrot_hum_R2 = np.zeros(n_frames)
unrot_hum_R3 = np.zeros(n_frames)

for i in range(0, len(q_elevation)):
    all_q[0, i] = q_elevation[i]* (-0.63355/2.61799) # same as unrot_scap_R2
    all_q[1, i] = q_elevation[i] * (0.322013/3.14159) # same as unrot_scap_R3
    all_q[2, i] = q_elevation[i] * (-0.322013 / 3.14159)
    all_q[3, i] = q_elevation[i] * (0.633555 / 2.61799)
    all_q[4, i] = q_elevation[i] * (-0.128282 / 2.61799)  # same as unrot_hum_R2
    all_q[5, i] = q_elevation[i] * (1.03673 / 2.61799)  # same as unrot_hum_R3
    all_q[6, i] = q_elevation[i] * (0.46603 / 2.61799)  # same as unrot_hum_R1
    all_q[7, i] = q_elevation[i] * (-0.46603 / 2.61799)
    all_q[8, i] = q_elevation[i] * (-1.03673 / 2.61799)
    all_q[9, i] = q_elevation[i] * (0.128282 / 2.61799)


# Animate
# if manually_animate:
#     i = 0
#     while biorbd_model.vtk_window.is_active:
#         biorbd_model.set_q(all_q[i, :])
#         i = (i+1) % n_frames
# else:
#     biorbd_model.load_movement(all_q)
#     biorbd_model.exec()
#
# np.set_printoptions(precision=10)
#
model_name = Models.Stanford_VA_upper_limb_model_0_40.value
biorbd_model = biorbd.Model(model_name)
Q0 = np.zeros((biorbd_model.nbQ()))
####################################################################


print(biorbd_model.nbQ())
ALL_Q_humerus = np.zeros((n_frames, 3))


Angles_table_sternoclav_R2 = np.zeros((n_frames, 3))
Angles_table_sternoclav_R3 = np.zeros((n_frames, 3))
Angles_table_accromioclav_R2 = np.zeros((n_frames, 3))
Angles_table_accromioclav_R3 = np.zeros((n_frames, 3))
Angles_table_accromioclav_R1 = np.zeros((n_frames, 3))

#A_sternoclav_R2, Angles_table_sternoclav_R2, Angles_sternoclav_R2 = matrice_1_rotation(5, n_frames, model_name, 3.14)
# A_accromio_R2, Angles_table_accromioclav_R2, angles_accrmio_R2 = matrice_1_rotation(19, n_frames, model_name, 0, 3.14)
# Angles_table_sternoclav_R2, Angles_table_sternoclav_R3 = matrice_2_rotation(5, 6, n_frames, model_name, 0 ,3.14)
# Angles_table_accromioclav_R2, Angles_table_accromioclav_R3, Angles_table_accromioclav_R1 = matrice_3_rotation(19, 20, 21, n_frames, model_name, 0 ,3.14)


for j in range(0, n_frames):
    Q = all_q[:, j]

    T0_clav = biorbd_model.globalJCS(Q, 9).rot().to_array()
    T0_clav_inv = T0_clav.transpose()

    T_clavphant = biorbd_model.globalJCS(Q, 16).rot().to_array()
    T_clavphant_inv = T_clavphant.transpose()

    T_scap = biorbd_model.globalJCS(Q, 20).rot().to_array()

    clav_T_clavphant = np.matmul(T0_clav_inv, T_clavphant)
    clavphant_T_scap = np.matmul(T_clavphant_inv, T_scap)

    clav_T_scap = np.matmul(clav_T_clavphant, clavphant_T_scap)

    Matrice_rotation = np.matmul(T0_clav, clav_T_scap)
    M = biorbd.Rotation(Matrice_rotation[0, 0], Matrice_rotation[0, 1], Matrice_rotation[0, 2],
                        Matrice_rotation[1, 0], Matrice_rotation[1, 1], Matrice_rotation[1, 2],
                        Matrice_rotation[2, 0], Matrice_rotation[2, 1], Matrice_rotation[2, 2])
    Euler_angles = biorbd.Rotation.toEulerAngles(M, 'xyz').to_array()
    Angles_table_sternoclav_R2[j, :] = Euler_angles

    #   elevation hum
    Rototrans_matrix_world_humerus = biorbd_model.globalJCS(Q, 52)
    Rototrans_matrix_world_humerus_angles = biorbd.RotoTrans_toEulerAngles(Rototrans_matrix_world_humerus, 'xyz')
    r_humerus = Rototrans_matrix_world_humerus_angles.to_array()
    ALL_Q_humerus[j,:]=r_humerus # put the matrix of the rotation of the humerus regarding the world in a table (3 column, j lines)
#
# ##########################################################################################
pente_sternoclav_R2 = np.zeros((n_frames-1, 3))
# pente_sternoclav_R3 = np.zeros((n_frames-1, 1))
# pente_accromioclav_R2 = np.zeros((n_frames-1, 1))
# pente_accromioclav_R3 = np.zeros((n_frames-1, 1))
# pente_accromioclav_R1 = np.zeros((n_frames-1, 1))
#
# # calcul de pentes :
for j in range(1, n_frames):
    pente_sternoclav_R2[j-1,0] = Angles_table_sternoclav_R2[j, 0] / ALL_Q_humerus[j, 0]
    pente_sternoclav_R2[j - 1, 1] = Angles_table_sternoclav_R2[j, 1] / ALL_Q_humerus[j, 0]
    pente_sternoclav_R2[j - 1, 2] = Angles_table_sternoclav_R2[j, 2] / ALL_Q_humerus[j, 0]

#     pente_sternoclav_R3[j-1] = Angles_table_sternoclav_R3[j, 1] /  -ALL_Q_humerus[j, 0]
#     pente_accromioclav_R2[j-1] = Angles_table_accromioclav_R2[j, 0] / -ALL_Q_humerus[j, 0]
#     pente_accromioclav_R3[j-1] = Angles_table_accromioclav_R3[j, 1] / -ALL_Q_humerus[j, 0]
#     pente_accromioclav_R1[j-1] = Angles_table_accromioclav_R1[j, 2] / -ALL_Q_humerus[j,0]
#
pente_moy_sternoclav_R2 = np.zeros((n_frames-1, 3))
pente_moy_sternoclav_R2[:, 0] = (sum((pente_sternoclav_R2[:, 0]))) / (n_frames-1)
pente_moy_sternoclav_R2[:, 1] = (sum((pente_sternoclav_R2[:, 1]))) / (n_frames-1)
pente_moy_sternoclav_R2[:, 2] = (sum((pente_sternoclav_R2[:, 2]))) / (n_frames-1)
# pente_moy_sternoclav_R3 = (sum((pente_sternoclav_R3))) / (n_frames-1)
# pente_moy_accromioclav_R2 = (sum((pente_accromioclav_R2))) / (n_frames-1)
# pente_moy_accromioclav_R3 = (sum((pente_accromioclav_R3))) / (n_frames-1)
# pente_moy_accromioclav_R1 = (sum((pente_accromioclav_R1))) / (n_frames-1)
#
print("La pente moyenne de la rotation de sternoclav R2 selon l'elevation de l'humerus est :")
print(pente_moy_sternoclav_R2)
print(pente_sternoclav_R2)
# print("La pente moyenne de la rotation de sternoclav R3 selon l'elevation de l'humerus est :")
# print(pente_moy_sternoclav_R3)
# print("La pente moyenne de la rotation de accromioclav R2 selon l'elevation de l'humerus est :")
# print(pente_moy_accromioclav_R2)
# print("La pente moyenne de la rotation de accromioclav R3 selon l'elevation de l'humerus est :")
# print(pente_moy_accromioclav_R3)
# print("La pente moyenne de la rotation de accromioclav R1 selon l'elevation de l'humerus est :")
# print(pente_moy_accromioclav_R1)
#
#
fig, axsY = plt.subplots(2, 3)

axsY[0, 0].plot(ALL_Q_humerus[:,0], Angles_table_sternoclav_R2[:,0], 'tab:red')
axsY[0, 0].set_title("rot sterno clav x par elev hum y")
axsY[0, 1].plot(ALL_Q_humerus[:,0], Angles_table_sternoclav_R2[:,1], 'tab:red')
axsY[0, 1].set_title("rot sterno clav y par elev hum y")
axsY[0, 2].plot(ALL_Q_humerus[:,0], Angles_table_sternoclav_R2[:,2], 'tab:red')
axsY[0, 2].set_title("rot accromio clav x par elev hum x")
# axsY[1, 0].plot(ALL_Q_humerus[:,0], Angles_table_accromioclav_R3[:,1], 'tab:red')
# axsY[1, 0].set_title("rot accromio clav y par elev hum x")
# axsY[1, 1].plot(ALL_Q_humerus[:,0], Angles_table_accromioclav_R1[:,2], 'tab:red')
# axsY[1, 1].set_title("rot accromio clav Z par elev hum x")
#
#
plt.show()
