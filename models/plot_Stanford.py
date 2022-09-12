# from matplotlib import pyplot as plt
#
# def plot_stanford(Q_hum, Q_autres):
#
#     figX, axs = plt.subplots(2, 5)
#
#     axs[0, 0].plot(Q_hum[0, :], ALL_Q_sterno_clav_R2[0, :], 'tab:red')
#     axs[0, 0].set_title("sterno_clav_R2 par elev hum x \n")
#
#     axs[0, 1].plot(Q_hum[0, :], ALL_Q_sterno_clav_R3[0, :], 'tab:red')
#     axs[0, 1].set_title("sterno_clav_R3 par elev hum x \n")
#
#     axs[0, 2].plot(Q_hum[0, :], ALL_Q_unrot_scap_R3[0, :], 'tab:red')
#     axs[0, 2].set_title("unrot_scap_R3 par elev hum x \n")
#
#     axs[0, 3].plot(Q_hum[0, :], ALL_Q_unrot_scap_R2[0, :], 'tab:red')
#     axs[0, 3].set_title("unrot_scap_R2 par elev hum x \n")
#
#     axs[0, 4].plot(Q_hum[0, :], ALL_Q_accromioclaviculaire_R2[0, :], 'tab:red')
#     axs[0, 4].set_title("accromioclaviculaire_R2 par elev hum x \n")
#
#     axs[1, 0].plot(Q_hum[0, :], ALL_Q_accromioclaviculaire_R3[0, :], 'tab:red')
#     axs[1, 0].set_title("accromioclaviculaire_R3 par elev hum x \n")
#
#     axs[1, 1].plot(Q_hum[0, :], ALL_Q_accromioclaviculaire_R1[0, :], 'tab:red')
#     axs[1, 1].set_title("accromioclaviculaire_R1 par elev hum x \n")
#
#     axs[1, 2].plot(Q_hum[0, :], ALL_Q_unrot_hum_R1[0, :], 'tab:red')
#     axs[1, 2].set_title("unrot_hum_R1 par elev hum x \n")
#
#     axs[1, 3].plot(Q_hum[0, :], ALL_Q_unrot_hum_R3[0, :], 'tab:red')
#     axs[1, 3].set_title("unrot_hum_R3 par elev hum x \n")
#
#     axs[1, 4].plot(Q_hum[0, :], ALL_Q_unrot_hum_R2[0, :], 'tab:red')
#     axs[1, 4].set_title("unrot_hum_R2 par elev hum x \n")
#
#     figY, axs = plt.subplots(2, 5)
#
#     axs[0, 0].plot(Q_hum[1, :], ALL_Q_sterno_clav_R2[1, :], 'tab:red')
#     axs[0, 0].set_title("sterno_clav_R2 par elev hum x \n")
#
#     axs[0, 1].plot(Q_hum[1, :], ALL_Q_sterno_clav_R3[1, :], 'tab:red')
#     axs[0, 1].set_title("sterno_clav_R3 par elev hum x \n")
#
#     axs[0, 2].plot(Q_hum[1, :], ALL_Q_unrot_scap_R3[1, :], 'tab:red')
#     axs[0, 2].set_title("unrot_scap_R3 par elev hum x \n")
#
#     axs[0, 3].plot(Q_hum[1, :], ALL_Q_unrot_scap_R2[1, :], 'tab:red')
#     axs[0, 3].set_title("unrot_scap_R2 par elev hum x \n")
#
#     axs[0, 4].plot(Q_hum[1, :], ALL_Q_accromioclaviculaire_R2[1, :], 'tab:red')
#     axs[0, 4].set_title("accromioclaviculaire_R2 par elev hum x \n")
#
#     axs[1, 0].plot(Q_hum[1, :], ALL_Q_accromioclaviculaire_R3[1, :], 'tab:red')
#     axs[1, 0].set_title("accromioclaviculaire_R3 par elev hum x \n")
#
#     axs[1, 1].plot(Q_hum[1, :], ALL_Q_accromioclaviculaire_R1[1, :], 'tab:red')
#     axs[1, 1].set_title("accromioclaviculaire_R1 par elev hum x \n")
#
#     axs[1, 2].plot(Q_hum[1, :], ALL_Q_unrot_hum_R1[1, :], 'tab:red')
#     axs[1, 2].set_title("unrot_hum_R1 par elev hum x \n")
#
#     axs[1, 3].plot(Q_hum[1, :], ALL_Q_unrot_hum_R3[1, :], 'tab:red')
#     axs[1, 3].set_title("unrot_hum_R3 par elev hum x \n")
#
#     axs[1, 4].plot(Q_hum[1, :], ALL_Q_unrot_hum_R2[1, :], 'tab:red')
#     axs[1, 4].set_title("unrot_hum_R2 par elev hum x \n")
#
#     figZ, axs = plt.subplots(2, 5)
#
#     axs[0, 0].plot(Q_hum[2, :], ALL_Q_sterno_clav_R2[2, :], 'tab:red')
#     axs[0, 0].set_title("sterno_clav_R2 par elev hum x \n")
#
#     axs[0, 1].plot(Q_hum[2, :], ALL_Q_sterno_clav_R3[2, :], 'tab:red')
#     axs[0, 1].set_title("sterno_clav_R3 par elev hum x \n")
#
#     axs[0, 2].plot(Q_hum[2, :], ALL_Q_unrot_scap_R3[2, :], 'tab:red')
#     axs[0, 2].set_title("unrot_scap_R3 par elev hum x \n")
#
#     axs[0, 3].plot(Q_hum[2, :], ALL_Q_unrot_scap_R2[2, :], 'tab:red')
#     axs[0, 3].set_title("unrot_scap_R2 par elev hum x \n")
#
#     axs[0, 4].plot(Q_hum[2, :], ALL_Q_accromioclaviculaire_R2[2, :], 'tab:red')
#     axs[0, 4].set_title("accromioclaviculaire_R2 par elev hum x \n")
#
#     axs[1, 0].plot(Q_hum[2, :], ALL_Q_accromioclaviculaire_R3[2, :], 'tab:red')
#     axs[1, 0].set_title("accromioclaviculaire_R3 par elev hum x \n")
#
#     axs[1, 1].plot(Q_hum[2, :], ALL_Q_accromioclaviculaire_R1[2, :], 'tab:red')
#     axs[1, 1].set_title("accromioclaviculaire_R1 par elev hum x \n")
#
#     axs[1, 2].plot(Q_hum[2, :], ALL_Q_unrot_hum_R1[2, :], 'tab:red')
#     axs[1, 2].set_title("unrot_hum_R1 par elev hum x \n")
#
#     axs[1, 3].plot(Q_hum[2, :], ALL_Q_unrot_hum_R3[2, :], 'tab:red')
#     axs[1, 3].set_title("unrot_hum_R3 par elev hum x \n")
#
#     axs[1, 4].plot(Q_hum[2, :], ALL_Q_unrot_hum_R2[2, :], 'tab:red')
#     axs[1, 4].set_title("unrot_hum_R2 par elev hum x \n")