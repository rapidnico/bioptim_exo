import bioviz
from enums import Models
import numpy as np
manually_animate = False

# Load the model
biorbd_viz = Models.Stanford_VA_upper_limb_model_0_40.value
b = bioviz.Viz(biorbd_viz, show_floor=False, show_global_ref_frame=False)




# Create a movement
n_frames = 100
all_q = np.zeros((b.nQ, n_frames))
all_q[4, :] = np.linspace(0, np.pi / 2, n_frames)

q_elevation = all_q[4, :]
q_sterno_clav_R2 = np.zeros(n_frames)
q_sterno_clav_R3 = np.zeros(n_frames)
q_accromioclaviculaire_R1 = np.zeros(n_frames)
q_accromioclaviculaire_R2 = np.zeros(n_frames)
q_accromioclaviculaire_R3 = np.zeros(n_frames)

for i in range(0, len(q_elevation)-1):
    q_sterno_clav_R2[i] = q_elevation[i]* 2.61799
    q_sterno_clav_R3[i] = q_elevation[i] * np.pi
    q_accromioclaviculaire_R1[i] = q_elevation[i] * 2.61799
    q_accromioclaviculaire_R2[i] = q_elevation[i] * 2.61799
    q_accromioclaviculaire_R3[i] = q_elevation[i] * 2.61799

if manually_animate:
    i = 0
    while b.vtk_window.is_active:
        b.set_q(all_q[i, :])
        i = (i+1) % n_frames
else:
    b.load_movement(all_q)
    b.exec()