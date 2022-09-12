import bioviz
from matplotlib import pyplot as plt
from enums import Models
import numpy as np
import biorbd
manually_animate = False

# Load the model
biorbd_viz = Models.Stanford_VA_upper_limb_model_0_40.value
b = bioviz.Viz(biorbd_viz, show_floor=False, show_global_ref_frame=False)


# Create a movement
n_frames = 20
all_q = np.zeros((b.nQ, n_frames))
all_q[11, :] = np.linspace(0, 3.14, n_frames)
q_elevation = all_q[11, :]

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
    q_sterno_clav_R2[i] = q_elevation[i]* (-0.63355/2.61799) #same as unrot_scap_R2
    q_sterno_clav_R3[i] = q_elevation[i] * (0.322013/3.14159)#same as unrot_scap_R3
    unrot_scap_R3[i] = q_elevation[i] * (-0.322013 / 3.14159)
    unrot_scap_R2[i] = q_elevation[i]* (0.633555/2.61799)
    q_accromioclaviculaire_R2[i] = q_elevation[i] * (-0.128282/2.61799)#same as unrot_hum_R2
    q_accromioclaviculaire_R3[i] = q_elevation[i] * (1.03673 / 2.61799)  # same as unrot_hum_R3
    q_accromioclaviculaire_R1[i] = q_elevation[i] * (0.46603 / 2.61799)  # same as unrot_hum_R1
    unrot_hum_R1[i] = q_elevation[i] * (-0.46603 / 2.61799)
    unrot_hum_R3[i] = q_elevation[i] * (-1.03673 / 2.61799)
    unrot_hum_R2[i] = q_elevation[i] * (0.128282 / 2.61799)

#
all_q[0, :] = q_sterno_clav_R2
all_q[1, :] = q_sterno_clav_R3
all_q[2, :] = unrot_scap_R3
all_q[3, :] = unrot_scap_R2
all_q[4, :] = q_accromioclaviculaire_R2
all_q[5, :] = q_accromioclaviculaire_R3
all_q[6, :] = q_accromioclaviculaire_R1
all_q[7, :] = unrot_hum_R1
all_q[8, :] = unrot_hum_R3
all_q[9, :] = unrot_hum_R2


if manually_animate:
    i = 0
    while b.vtk_window.is_active:
        b.set_q(all_q[i, :])
        i = (i+1) % n_frames
else:
    b.load_movement(all_q)
    b.exec()

np.set_printoptions(precision=10)

model_name = Models.Stanford_VA_upper_limb_model_0_40.value

biorbd_model = biorbd.Model(model_name)

# Now the model is loaded as a biorbd object
print(biorbd_model.nbQ())
ALL_Q_humerus = np.zeros((n_frames, 3))
ALL_Q_sterno_clav_R2 = np.zeros((n_frames, 3))
ALL_Q_sterno_clav_R3 = np.zeros((n_frames, 3))
ALL_Q_unrot_scap_R3 = np.zeros((n_frames, 3))
ALL_Q_unrot_scap_R2 = np.zeros((n_frames, 3))
ALL_Q_accromioclaviculaire_R2 = np.zeros((n_frames, 3))
ALL_Q_accromioclaviculaire_R3 = np.zeros((n_frames, 3))
ALL_Q_accromioclaviculaire_R1 = np.zeros((n_frames, 3))
ALL_Q_unrot_hum_R1 = np.zeros((n_frames, 3))
ALL_Q_unrot_hum_R3 = np.zeros((n_frames, 3))
ALL_Q_unrot_hum_R2 = np.zeros((n_frames, 3))


for j in range(0, n_frames):
    Q = all_q[:, j]
    Rototrans_matrix_world_humerus = biorbd_model.globalJCS(Q, 11)
    Rototrans_matrix_world_humerus_angles = biorbd.RotoTrans_toEulerAngles(Rototrans_matrix_world_humerus, 'xyz')
    r_humerus = Rototrans_matrix_world_humerus_angles.to_array()
    ALL_Q_humerus[j,:]=r_humerus

    Rototrans_matrix_world_sterno_clav_R2 = biorbd_model.globalJCS(Q, 0)
    Rototrans_matrix_world_sterno_clav_R2_angles = biorbd.RotoTrans_toEulerAngles(Rototrans_matrix_world_sterno_clav_R2, 'xyz')
    r_sterno_clav_R2 = Rototrans_matrix_world_sterno_clav_R2_angles.to_array()
    ALL_Q_sterno_clav_R2[j,:]=r_sterno_clav_R2

    Rototrans_matrix_world_sterno_clav_R3 = biorbd_model.globalJCS(Q, 1)
    Rototrans_matrix_world_sterno_clav_R3_angles = biorbd.RotoTrans_toEulerAngles(Rototrans_matrix_world_sterno_clav_R3, 'xyz')
    r_sterno_clav_R3 = Rototrans_matrix_world_sterno_clav_R3_angles.to_array()
    ALL_Q_sterno_clav_R3[j, :] = r_sterno_clav_R3

    Rototrans_matrix_world_unrot_scap_R3 = biorbd_model.globalJCS(Q, 2)
    Rototrans_matrix_world_unrot_scap_R3_angles = biorbd.RotoTrans_toEulerAngles(Rototrans_matrix_world_unrot_scap_R3, 'xyz')
    r_unrot_scap_R3 = Rototrans_matrix_world_unrot_scap_R3_angles.to_array()
    ALL_Q_unrot_scap_R3[j, :] = r_unrot_scap_R3

    Rototrans_matrix_world_unrot_scap_R2 = biorbd_model.globalJCS(Q, 3)
    Rototrans_matrix_world_unrot_scap_R2_angles = biorbd.RotoTrans_toEulerAngles(Rototrans_matrix_world_unrot_scap_R2, 'xyz')
    r_unrot_scap_R2 = Rototrans_matrix_world_unrot_scap_R2_angles.to_array()
    ALL_Q_unrot_scap_R2[j, :] = r_unrot_scap_R2

    Rototrans_matrix_world_accromioclaviculaire_R2 = biorbd_model.globalJCS(Q, 4)
    Rototrans_matrix_world_accromioclaviculaire_R2_angles = biorbd.RotoTrans_toEulerAngles(Rototrans_matrix_world_accromioclaviculaire_R2, 'xyz')
    r_accromioclaviculaire_R2 = Rototrans_matrix_world_accromioclaviculaire_R2_angles.to_array()
    ALL_Q_accromioclaviculaire_R2[j, :] = r_accromioclaviculaire_R2

    Rototrans_matrix_world_accromioclaviculaire_R3 = biorbd_model.globalJCS(Q, 5)
    Rototrans_matrix_world_accromioclaviculaire_R3_angles = biorbd.RotoTrans_toEulerAngles(Rototrans_matrix_world_accromioclaviculaire_R3, 'xyz')
    r_accromioclaviculaire_R3 = Rototrans_matrix_world_accromioclaviculaire_R3_angles.to_array()
    ALL_Q_accromioclaviculaire_R3[j, :] = r_accromioclaviculaire_R3

    Rototrans_matrix_world_accromioclaviculaire_R1 = biorbd_model.globalJCS(Q, 6)
    Rototrans_matrix_world_accromioclaviculaire_R1_angles = biorbd.RotoTrans_toEulerAngles(Rototrans_matrix_world_accromioclaviculaire_R1, 'xyz')
    r_accromioclaviculaire_R1 = Rototrans_matrix_world_accromioclaviculaire_R1_angles.to_array()
    ALL_Q_accromioclaviculaire_R1[j, :] = r_accromioclaviculaire_R1

    Q = all_q[:, j]
    Rototrans_matrix_world_unrot_hum_R1 = biorbd_model.globalJCS(Q, 7)
    Rototrans_matrix_world_unrot_hum_R1_angles = biorbd.RotoTrans_toEulerAngles(Rototrans_matrix_world_unrot_hum_R1, 'xyz')
    r_unrot_hum_R1 = Rototrans_matrix_world_unrot_hum_R1_angles.to_array()
    ALL_Q_unrot_hum_R1[j, :] = r_unrot_hum_R1

    Q = all_q[:, j]
    Rototrans_matrix_world_unrot_hum_R3 = biorbd_model.globalJCS(Q, 8)
    Rototrans_matrix_world_unrot_hum_R3_angles = biorbd.RotoTrans_toEulerAngles(Rototrans_matrix_world_unrot_hum_R3, 'xyz')
    r_unrot_hum_R3 = Rototrans_matrix_world_unrot_hum_R3_angles.to_array()
    ALL_Q_unrot_hum_R3[j, :] = r_unrot_hum_R3

    Q = all_q[:, j]
    Rototrans_matrix_world_unrot_hum_R2 = biorbd_model.globalJCS(Q, 9)
    Rototrans_matrix_world_unrot_hum_R2_angles = biorbd.RotoTrans_toEulerAngles(Rototrans_matrix_world_unrot_hum_R2, 'xyz')
    r_unrot_hum_R2 = Rototrans_matrix_world_unrot_hum_R2_angles.to_array()
    ALL_Q_unrot_hum_R2[j, :] = r_unrot_hum_R2

print(ALL_Q_humerus)

# Q_X_Membres = [ALL_Q_sterno_clav_R2[:, 0],
#                 ALL_Q_sterno_clav_R3[:, 0],
#                 ALL_Q_unrot_scap_R3[:, 0],
#                 ALL_Q_unrot_scap_R2[:, 0],
#                 ALL_Q_accromioclaviculaire_R2[:, 0],
#                 ALL_Q_accromioclaviculaire_R3[:, 0],
#                 ALL_Q_accromioclaviculaire_R1[:, 0],
#                 ALL_Q_unrot_hum_R1[:, 0],
#                 ALL_Q_unrot_hum_R3[:, 0],
#                 ALL_Q_unrot_hum_R2[:, 0]
#                ]
# Q_Y_Membres = [ALL_Q_sterno_clav_R2[:, 1],
#                 ALL_Q_sterno_clav_R3[:, 1],
#                 ALL_Q_unrot_scap_R3[:, 1],
#                 ALL_Q_unrot_scap_R2[:, 1],
#                 ALL_Q_accromioclaviculaire_R2[:, 1],
#                 ALL_Q_accromioclaviculaire_R3[:, 1],
#                 ALL_Q_accromioclaviculaire_R1[:, 1],
#                 ALL_Q_unrot_hum_R1[:, 1],
#                 ALL_Q_unrot_hum_R3[:, 1],
#                 ALL_Q_unrot_hum_R2[:, 1]
#                ]
# Q_Z_Membres = [ALL_Q_sterno_clav_R2[:, 2],
#                 ALL_Q_sterno_clav_R3[:, 2],
#                 ALL_Q_unrot_scap_R3[:, 2],
#                 ALL_Q_unrot_scap_R2[:, 2],
#                 ALL_Q_accromioclaviculaire_R2[:, 2],
#                 ALL_Q_accromioclaviculaire_R3[:, 2],
#                 ALL_Q_accromioclaviculaire_R1[:, 2],
#                 ALL_Q_unrot_hum_R1[:, 2],
#                 ALL_Q_unrot_hum_R3[:, 2],
#                 ALL_Q_unrot_hum_R2[:, 2]
#                ]

for i in range(0,2):
    for j in range(0,2):

        fig, axs = plt.subplots(2, 5)

        axs[0, 0].plot(ALL_Q_humerus[:,i], ALL_Q_sterno_clav_R2[:,2], 'tab:red')
        axs[0, 0].set_title("sterno_clav_R2 " + j + " par elev hum " + i)

        axs[0, 1].plot(ALL_Q_humerus[:,i], ALL_Q_sterno_clav_R3[:,2], 'tab:red')
        axs[0, 1].set_title("sterno_clav_R3 " + j + " par elev hum " + i)

        axs[0, 2].plot(ALL_Q_humerus[:,i], ALL_Q_unrot_scap_R3[:,2], 'tab:red')
        axs[0, 2].set_title("unrot_scap_R3 " + j + " par elev hum " + i)

        axs[0, 3].plot(ALL_Q_humerus[:,i], ALL_Q_unrot_scap_R2[:,2], 'tab:red')
        axs[0, 3].set_title("unrot_scap_R2 " + j + " par elev hum " + i)

        axs[0, 4].plot(ALL_Q_humerus[:,i], ALL_Q_accromioclaviculaire_R2[:,2], 'tab:red')
        axs[0, 4].set_title("accromioclaviculaire_R2 " + j + " par elev hum " + i)

        axs[1, 0].plot(ALL_Q_humerus[:,i], ALL_Q_accromioclaviculaire_R3[:,2], 'tab:red')
        axs[1, 0].set_title("accromioclaviculaire_R3 " + j + " par elev hum " + i)

        axs[1, 1].plot(ALL_Q_humerus[:,i], ALL_Q_accromioclaviculaire_R1[:,2], 'tab:red')
        axs[1, 1].set_title("accromioclaviculaire_R1 " + j + " par elev hum " + i)

        axs[1, 2].plot(ALL_Q_humerus[:,i], ALL_Q_unrot_hum_R1[:,2], 'tab:red')
        axs[1, 2].set_title("unrot_hum_R1 " + j + " par elev hum " + i)

        axs[1, 3].plot(ALL_Q_humerus[:,i], ALL_Q_unrot_hum_R3[:,2], 'tab:red')
        axs[1, 3].set_title("unrot_hum_R3 " + j + " par elev hum " + i)

        axs[1, 4].plot(ALL_Q_humerus[:,i], ALL_Q_unrot_hum_R2[:,2], 'tab:red')
        axs[1, 4].set_title("unrot_hum_R2 " + j + " par elev hum " + i)




plt.show()
# We have Rototrans_matrix_ulna_support, so, we are now able to transform the coordinates of any point in the support
# frame to the ulna frame with the following formula:
# Position_in_ulna = Rototrans_matrix_ulna_support * Position_in_support