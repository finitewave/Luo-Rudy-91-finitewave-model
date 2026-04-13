"""
This module provides a simple interface to run the Luo-Rudy 91 model in a 0D setting,
i.e., without spatial dimensions. It includes classes for defining stimulation protocols
and for managing the simulation of the Luo-Rudy 91 equations over time.
"""

from luo_rudy_91 import ops


class Stimulation:
    """
    Stimulus description for a 0D simulation.

    Parameters
    ----------
    t_start : float
        Start time (ms) of the first stimulus window.
    duration : float
        Duration (ms) of a single pulse.
    amplitude : float
        Pulse amplitude in the same units as du/dt contribution (typically "units/ms").

    Method
    ------
    stim(t: float) -> float
        Returns the instantaneous stimulus value at time t.

    """

    def __init__(self, t_start: float, duration: float, amplitude: float):
        self.t_start = t_start
        self.duration = duration
        self.amplitude = amplitude

    def stim(self, t: float) -> float:
        return self.amplitude if self.t_start <= t < self.t_start + self.duration else 0.0


class LuoRudy910D:
    """
    Luo-Rudy 91 0D model wrapper.
    """
    def __init__(self, dt: float, stimulations: list[Stimulation]):
        self.dt = dt
        self.stimulations = stimulations
        self.variables = ops.get_variables()
        self.parameters = ops.get_parameters()
        self.history = {s: [] for s in self.variables}
        self.stim_history = []
        self.times = []

    def step(self, i: int):
        """
        Perform a single time step update.

        Parameters
        ----------
        i : int
            Current time step index.
        """
        u_old = self.variables["u"]
        m_old = self.variables["m"]
        h_old = self.variables["h"]
        j_old = self.variables["j"]
        d_old = self.variables["d"]
        f_old = self.variables["f"]
        x_old = self.variables["x"]
        cai_old = self.variables["cai"]

        rhs, m_new, h_new, j_new, d_new, f_new, x_new, cai_new = ops.ionic_step(
            self.dt, u_old, m_old, h_old, j_old, d_old, f_old, x_old, cai_old,
            self.parameters["gna"], self.parameters["gsi"], self.parameters["gk"],
            self.parameters["gk1"], self.parameters["gkp"], self.parameters["gb"],
            self.parameters["ko"], self.parameters["ki"], self.parameters["nao"],
            self.parameters["nai"], self.parameters["R"], self.parameters["T"],
            self.parameters["F"], self.parameters["PR_NaK"],
            self.parameters["E_Na"], self.parameters["E_K1"]
        )
        stim_curr = self.dt * sum(stim.stim(t=self.dt*i) for stim in self.stimulations)
        self.stim_history.append(stim_curr)

        u_new = u_old + self.dt * rhs + stim_curr

        self.variables["u"] = u_new
        self.variables["m"] = m_new
        self.variables["h"] = h_new
        self.variables["j"] = j_new
        self.variables["d"] = d_new
        self.variables["f"] = f_new
        self.variables["x"] = x_new
        self.variables["cai"] = cai_new

    def run(self, t_max: float):
        """
        Run the simulation up to time t_max.
        
        Parameters
        ----------
        t_max : float
            Maximum simulation time.
        """
        n_steps = int(round(t_max/self.dt))
        for i in range(n_steps):
            self.step(i)
            self.times.append(self.dt * i)
            for s in self.variables:
                self.history[s].append(self.variables[s])