from bioptim import OptimalControlProgram
import numpy as np
import pickle
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
import biorbd
from quat import quat2eul
from data.enums import Tasks
from models.enums import Models
from bioptim import (
    OptimalControlProgram,
    DynamicsFcn,
    Dynamics,
    Bounds,
    QAndQDotBounds,
    InitialGuess,
    ObjectiveFcn,
    ObjectiveList,
    ConstraintList,
    ConstraintFcn,
    Constraint,
    Objective,
    OdeSolver,
    CostType,
    Solver,
    Node,
)

import data.load_events as load_events
from lib import custom_load

def prepare_ocp(
    biorbd_model_path: Models.Stanford_VA_upper_limb_model_0_40.value,
    final_time: 5,
    n_shooting: 100,
    ode_solver: OdeSolver = OdeSolver.COLLOCATION(),
    n_threads: int = 4,
) -> OptimalControlProgram:
    """
    The initialization of an ocp

    Parameters
    ----------
    biorbd_model_path: str
        The path to the biorbd model
    final_time: float
        The time in second required to perform the task
    n_shooting: int
        The number of shooting points to define int the direct multiple shooting program
    ode_solver: OdeSolver = OdeSolver.RK4()
        Which type of OdeSolver to use
    use_sx: bool
        If the SX variable should be used instead of MX (can be extensive on RAM)
    n_threads: int
        The number of threads to use in the paralleling (1 = no parallel computing)

    Returns
    -------
    The OptimalControlProgram ready to be solved
    """

    biorbd_model = biorbd.Model(biorbd_model_path)

    # Add objective functions

    objective_functions = ObjectiveList()
    # objective_functions.add(ObjectiveFcn.Lagrange.PROPORTIONAL_CONTROL, key="tau", first_dof=11, second_dof=0, coef=(-0.63355/2.61799))
    # objective_functions.add(ObjectiveFcn.Lagrange.PROPORTIONAL_CONTROL, key="tau", first_dof=11, second_dof=1, coef=(0.322013/3.14159))
    # objective_functions.add(ObjectiveFcn.Lagrange.PROPORTIONAL_CONTROL, key="tau", first_dof=11, second_dof=2, coef=(-0.322013/3.14159))
    # objective_functions.add(ObjectiveFcn.Lagrange.PROPORTIONAL_CONTROL, key="tau", first_dof=11, second_dof=3, coef=(0.633555/2.61799))
    # objective_functions.add(ObjectiveFcn.Lagrange.PROPORTIONAL_CONTROL, key="tau", first_dof=11, second_dof=4, coef=(-0.128282/2.61799))
    # objective_functions.add(ObjectiveFcn.Lagrange.PROPORTIONAL_CONTROL, key="tau", first_dof=11, second_dof=5, coef=(1.03673/2.61799))
    # objective_functions.add(ObjectiveFcn.Lagrange.PROPORTIONAL_CONTROL, key="tau", first_dof=11, second_dof=6, coef=(0.46603/2.61799))
    # objective_functions.add(ObjectiveFcn.Lagrange.PROPORTIONAL_CONTROL, key="tau", first_dof=11, second_dof=7, coef=(-0.46603/2.61799))
    # objective_functions.add(ObjectiveFcn.Lagrange.PROPORTIONAL_CONTROL, key="tau", first_dof=11, second_dof=8, coef=(-1.03673/2.61799))
    # objective_functions.add(ObjectiveFcn.Lagrange.PROPORTIONAL_CONTROL, key="tau", first_dof=11, second_dof=9, coef=(0.128282/2.61799))
    # objective_functions.add(ObjectiveFcn.Lagrange.MINIMIZE_CONTROL, key="tau")
    # objective_functions.add(ObjectiveFcn.Mayer.SUPERIMPOSE_MARKERS, node=Node.END, first_marker="EPI_lat", second_marker="m1")


    # Dynamics
    dynamics = Dynamics(DynamicsFcn.TORQUE_DRIVEN)

    # Constraints
    constraints = ConstraintList()
    # constraints = Constraint(ConstraintFcn.SUPERIMPOSE_MARKERS, node=Node.END, first_marker="EPI_lat", second_marker="m1",index = 0)

    # Path constraint
    x_bounds = QAndQDotBounds(biorbd_model)
    # x_bounds[:, [0, -1]] = 0
    # x_bounds[0, -1] = 3.14
    # x_bounds[1, -1] = 0

    # Initial guess
    n_q = biorbd_model.nbQ()
    n_qdot = biorbd_model.nbQdot()
    x_init = InitialGuess([11] * (n_q + n_qdot))

    # Define control path constraint
    n_u = 14
    qddot_joints_min, qddot_joints_max, qddot_joints_init = -100, 100, 0
    u_bounds = Bounds([qddot_joints_min] * n_u, [qddot_joints_max] * n_u)

    u_init = InitialGuess([qddot_joints_init] * n_u)


    return OptimalControlProgram(
        biorbd_model_path,
        dynamics,
        n_shooting,
        final_time,
        x_init,
        u_init,
        x_bounds,
        u_bounds,
        objective_functions,
        constraints,

        ode_solver=ode_solver,
    )

def main():


        biorbd_model_path = Models.Stanford_VA_upper_limb_model_0_40.value
        # --- Prepare the ocp --- #
        ocp = prepare_ocp(biorbd_model_path, final_time=0.2, n_shooting=50)

        # Custom plots
        ocp.add_plot_penalty(CostType.ALL)
        solver = Solver.IPOPT(show_online_optim=False)
        sol = ocp.solve(solver)

        # --- Print ocp structure --- #
        #ocp.print(to_console=False, to_graph=False)

        # --- Solve the ocp --- #
        sol.animate(n_frames=100)

if __name__ == "__main__":
    main()
