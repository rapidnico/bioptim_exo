import bioviz
from matplotlib import pyplot as plt
from enums import Models
import numpy as np
import biorbd
from scipy import stats
manually_animate = False

# Load the model
biorbd_viz = Models.Stanford_VA_upper_limb_model_0_40.value
b = bioviz.Viz(biorbd_viz, show_floor=False, show_global_ref_frame=False)


# Create a movement
n_frames = 20
all_q = np.zeros((b.nQ, n_frames))
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
    q_sterno_clav_R2[i] = q_elevation[i]* (-0.63355/2.61799) #same as unrot_scap_R2
    q_sterno_clav_R3[i] = q_elevation[i] * (0.322013/3.14159)#same as unrot_scap_R3
    unrot_scap_R3[i] = q_elevation[i] * (-0.322013 / 3.14159)
    unrot_scap_R2[i] = q_elevation[i]* (0.633555/2.61799)
    q_accromioclaviculaire_R2[i] = q_elevation[i] * (-0.128282/2.61799) #same as unrot_hum_R2
    q_accromioclaviculaire_R3[i] = q_elevation[i] * (1.03673 / 2.61799)  # same as unrot_hum_R3
    q_accromioclaviculaire_R1[i] = q_elevation[i] * (0.46603 / 2.61799)  # same as unrot_hum_R1
    unrot_hum_R1[i] = q_elevation[i] * (-0.46603 / 2.61799)
    unrot_hum_R3[i] = q_elevation[i] * (-1.03673 / 2.61799)
    unrot_hum_R2[i] = q_elevation[i] * (0.128282 / 2.61799)

# Send the applied movements in the correct DoF
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

#Animate
# if manually_animate:
#     i = 0
#     while b.vtk_window.is_active:
#         b.set_q(all_q[i, :])
#         i = (i+1) % n_frames
# else:
#     b.load_movement(all_q)
#     b.exec()

# np.set_printoptions(precision=10)
#
model_name = Models.Stanford_VA_upper_limb_model_0_40.value
biorbd_model = biorbd.Model(model_name)
####################################################################


print(biorbd_model.nbQ())
ALL_Q_humerus = np.zeros((n_frames, 3))
# ALL_Q_clav_ = np.zeros((n_frames, 3))
# ALL_Q_scapula = np.zeros((n_frames, 3))
# ALL_Q_thorax= np.zeros((n_frames, 3))

Q0 = np.zeros((b.nQ))
Angles_table_sternoclav_R2 = np.zeros((n_frames, 3))
Angles_table_sternoclav_R3 = np.zeros((n_frames, 3))


for j in range(0, n_frames):
    Q = all_q[:, j]

#   Rot sternoclav R2
    T0_sternoclav_R2 = biorbd_model.globalJCS(Q0, 5)
    rot_T0_inverse_sternoclav_R2 = T0_sternoclav_R2.rot().transpose()
    T0_inverse_sternoclav_R2 = rot_T0_inverse_sternoclav_R2.to_array()

    rot_T_sternoclav_R2 = biorbd_model.globalJCS(Q, 5).rot()
    T_sternoclav_R2 = rot_T_sternoclav_R2.to_array()

    A_sternoclav_R2 = np.matmul(T0_inverse_sternoclav_R2, T_sternoclav_R2)
    M_sternoclav_R2 = biorbd.Rotation(A_sternoclav_R2[0, 0], A_sternoclav_R2[0, 1], A_sternoclav_R2[0, 2], A_sternoclav_R2[1, 0], A_sternoclav_R2[1, 1], A_sternoclav_R2[1, 2], A_sternoclav_R2[2, 0], A_sternoclav_R2[2, 1], A_sternoclav_R2[2, 2])
    Angles_sternoclav_R2 = biorbd.Rotation.toEulerAngles(M_sternoclav_R2, 'xyz').to_array()
    Angles_table_sternoclav_R2[j,:] = Angles_sternoclav_R2


    #   Rot sternoclav R3
    T0_sternoclav_R3 = biorbd_model.globalJCS(Q0, 6)
    rot_T0_inverse_sternoclav_R3 = T0_sternoclav_R3.rot().transpose()
    T0_inverse_sternoclav_R3 = rot_T0_inverse_sternoclav_R3.to_array()

    rot_T_sternoclav_R3 = biorbd_model.globalJCS(Q, 6).rot()
    T_sternoclav_R3 = rot_T_sternoclav_R3.to_array()

    A_sternoclav_R3 = np.matmul(T0_inverse_sternoclav_R3, T_sternoclav_R3)
    M_sternoclav_R3 = biorbd.Rotation(A_sternoclav_R3[0, 0], A_sternoclav_R3[0, 1], A_sternoclav_R3[0, 2], A_sternoclav_R3[1, 0], A_sternoclav_R3[1, 1], A_sternoclav_R3[1, 2], A_sternoclav_R3[2, 0], A_sternoclav_R3[2, 1], A_sternoclav_R3[2, 2])
    Angles_sternoclav_R3 = biorbd.Rotation.toEulerAngles(M_sternoclav_R3, 'xyz').to_array()
    Angles_table_sternoclav_R3[j, :] = Angles_sternoclav_R3


    #   elevation hum
    Rototrans_matrix_world_humerus = biorbd_model.globalJCS(Q, 52)
    Rototrans_matrix_world_humerus_angles = biorbd.RotoTrans_toEulerAngles(Rototrans_matrix_world_humerus, 'xyz')
    r_humerus = Rototrans_matrix_world_humerus_angles.to_array()
    ALL_Q_humerus[j,:]=r_humerus # put the matrix of the rotation of the humerus regarding the world in a table (3 column, j lines)

##########################################################################################



print(ALL_Q_humerus)
# print(ALL_Q_clav_)
# print(ALL_Q_thorax)
# print("La pente de la rotation X de la clavicule en fonction de l'elevation de l'humerus est :")
# print(slope_clav_X)
# print("La pente de la rotation X de la scapula en fonction de l'elevation de l'humerus est :")
# print(slope_scapula_X)

fig, axsY = plt.subplots(2, 3)

axsY[0, 0].plot(ALL_Q_humerus[:,0], Angles_table_sternoclav_R2[:,0], 'tab:red')
axsY[0, 0].set_title("rot clav x par elev hum x")
axsY[0, 1].plot(ALL_Q_humerus[:,0], Angles_table_sternoclav_R2[:,1], 'tab:red')
axsY[0, 1].set_title("rot clav y par elev hum x")
axsY[0, 2].plot(ALL_Q_humerus[:,0], Angles_table_sternoclav_R2[:,2], 'tab:red')
axsY[0, 2].set_title("rot clav z par elev hum x")
axsY[1, 0].plot(ALL_Q_humerus[:,0], Angles_table_sternoclav_R3[:,0], 'tab:red')
axsY[1, 0].set_title("rot clav x par elev hum x")
axsY[1, 1].plot(ALL_Q_humerus[:,0], Angles_table_sternoclav_R3[:,1], 'tab:red')
axsY[1, 1].set_title("rot clav y par elev hum x")
axsY[1, 2].plot(ALL_Q_humerus[:,0], Angles_table_sternoclav_R3[:,2], 'tab:red')
axsY[1, 2].set_title("rot clav z par elev hum x")
#
# fig, axsX = plt.subplots(2, 3)
#
# axsX[0, 0].plot(ALL_Q_humerus[:,1], ALL_Q_clav_[:,0], 'tab:red')
# axsX[0, 0].set_title("rot clav x par elev hum y")
# axsX[1, 0].plot(ALL_Q_humerus[:,1], ALL_Q_scapula[:,0], 'tab:red')
# axsX[1, 0].set_title("rot scapula x par elev hum y")
# axsX[0, 1].plot(ALL_Q_humerus[:,1], ALL_Q_clav_[:,1], 'tab:red')
# axsX[0, 1].set_title("rot clav y par elev hum y")
# axsX[1, 1].plot(ALL_Q_humerus[:,1], ALL_Q_scapula[:,1], 'tab:red')
# axsX[1, 1].set_title("rot scapula y par elev hum y")
# axsX[0, 2].plot(ALL_Q_humerus[:,1], ALL_Q_clav_[:,2], 'tab:red')
# axsX[0, 2].set_title("rot clav z par elev hum y")
# axsX[1, 2].plot(ALL_Q_humerus[:,1], ALL_Q_scapula[:,2], 'tab:red')
# axsX[1, 2].set_title("rot scapula z par elev hum y")
#
# fig, axsZ = plt.subplots(2, 3)
#
# axsZ[0, 0].plot(ALL_Q_humerus[:,2], ALL_Q_clav_[:,0], 'tab:red')
# axsZ[0, 0].set_title("rot clav x par elev hum z")
# axsZ[1, 0].plot(ALL_Q_humerus[:,2], ALL_Q_scapula[:,0], 'tab:red')
# axsZ[1, 0].set_title("rot scapula x par elev hum z")
# axsZ[0, 1].plot(ALL_Q_humerus[:,2], ALL_Q_clav_[:,1], 'tab:red')
# axsZ[0, 1].set_title("rot clav y par elev hum z")
# axsZ[1, 1].plot(ALL_Q_humerus[:,2], ALL_Q_scapula[:,1], 'tab:red')
# axsZ[1, 1].set_title("rot scapula y par elev hum z")
# axsZ[0, 2].plot(ALL_Q_humerus[:,2], ALL_Q_clav_[:,2], 'tab:red')
# axsZ[0, 2].set_title("rot clav z par elev hum z")
# axsZ[1, 2].plot(ALL_Q_humerus[:,2], ALL_Q_scapula[:,2], 'tab:red')
# axsZ[1, 2].set_title("rot scapula z par elev hum z")

# fig, axsX = plt.subplots(2, 2)
#
# axsX[0, 0].plot(ALL_Q_humerus[:,1], ALL_Q_clav_[:,0], 'tab:red')
# axsX[0, 0].plot(ALL_Q_humerus[:,1], fitLine_clav_X, 'tab:blue')
# axsX[0, 0].set_title("sterno_clav_R2 x par elev hum y et droite de regression ")
#
# axsX[0, 1].plot(ALL_Q_humerus[:,1], ALL_Q_scapula[:,0], 'tab:red')
# axsX[0, 1].plot(ALL_Q_humerus[:,1], fitLine_scap_X, 'tab:blue')
# axsX[0, 1].set_title("sterno_clav_R3 x par elev hum y et droite de regression")

# ############################################################################
# test regression
# slope_clav_X, intercept_clav_X, r_value_clav_X, p_value_clav_X, std_err_clav_X = stats.linregress(ALL_Q_clav_[:,0], ALL_Q_humerus[:,1])
# slope_scapula_X, intercept_scapula_X, r_value_scapula_X, p_value_scapula_X, std_err_scapula_X = stats.linregress(ALL_Q_scapula[:,0], ALL_Q_humerus[:,1])
#
# fitLine_clav_X = slope_clav_X*ALL_Q_humerus[:,1]
# fitLine_scap_X = slope_scapula_X*ALL_Q_humerus[:,1]



# slope_clav_Y, intercept_clav_Y, r_value_clav_Y, p_value_clav_Y, std_err_clav_Y = stats.linregress(ALL_Q_clav_[:,1],ALL_Q_humerus [:,1])
# slope_scapula_Y, intercept_scapula_Y, r_value_scapula_Y, p_value_scapula_Y, std_err_scapula_Y = stats.linregress(ALL_Q_scapula[:,1], ALL_Q_humerus[:,1])
# #
# fitLine_clav_Y = slope_clav_Y*ALL_Q_humerus[:,1] + intercept_clav_Y
# fitLine_scap_Y = slope_scapula_Y*ALL_Q_humerus[:,1] + intercept_scapula_Y
# #
# print("La pente de la rotation Y de la clavicule en fonction de l'elevation de l'humerus est :")
# print(slope_clav_Y)
# print("La pente de la rotation Y de la scapula en fonction de l'elevation de l'humerus est :")
# print(slope_scapula_Y)
# #
# fig, axsY = plt.subplots(2, 2)
#
# axsY[0, 0].plot(ALL_Q_humerus[:,1], ALL_Q_clav_[:,1], 'tab:red')
# axsY[0, 0].plot(ALL_Q_humerus[:,1], fitLine_clav_Y, 'tab:blue')
# axsY[0, 0].set_title("sterno_clav_R2 y par elev hum y et droite de regression")
#
# axsY[0, 1].plot(ALL_Q_humerus[:,1], ALL_Q_scapula[:,1], 'tab:red')
# axsY[0, 1].plot(ALL_Q_humerus[:,1], fitLine_scap_Y, 'tab:blue')
# axsY[0, 1].set_title("sterno_clav_R3 y par elev hum y et droite de regression")

# ##################################################################################
#
# slope_clav_Z, intercept_clav_Z, r_value_clav_Z, p_value_clav_Z, std_err_clav_Z = stats.linregress(ALL_Q_clav_[:,2], ALL_Q_humerus[:,1])
# slope_scapula_Z, intercept_scapula_Z, r_value_scapula_Z, p_value_scapula_Z, std_err_scapula_Z = stats.linregress(ALL_Q_scapula[:,2], ALL_Q_humerus[:,1])
# #
# # fitLine_clav_Z = slope_clav_Z*ALL_Q_clav_[:,2] + intercept_clav_Z
# # fitLine_scap_Z = slope_scapula_Z*ALL_Q_scapula[:,2] + intercept_scapula_Z
# #
# print("La pente de la rotation Z de la clavicule en fonction de l'elevation de l'humerus est :")
# print(slope_clav_Z)
# print("La pente de la rotation Z de la scapula en fonction de l'elevation de l'humerus est :")
# print(slope_scapula_Z)
# #
#
# fig, axsZ = plt.subplots(2, 2)
#
# axsZ[0, 0].plot(ALL_Q_humerus[:,1], ALL_Q_clav_[:,2], 'tab:red')
# axsZ[0, 0].plot(ALL_Q_humerus[:,1], fitLine_clav_Z, 'tab:blue')
# axsZ[0, 0].set_title("sterno_clav_R2 z par elev hum y")
#
# axsZ[0, 1].plot(ALL_Q_humerus[:,1], ALL_Q_scapula[:,2], 'tab:red')
# axsZ[0, 1].plot(ALL_Q_humerus[:,1], fitLine_scap_Z, 'tab:blue')
# axsZ[0, 1].set_title("sterno_clav_R3 z par elev hum y")


plt.show()
