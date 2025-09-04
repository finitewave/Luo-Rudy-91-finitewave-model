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

    def step(self, i: int):
        """
        Perform a single time step update.

        Parameters
        ----------
        i : int
            Current time step index.
        """
        # Fast sodium current:
        self.variables["m"] += self.dt*ops.calc_dm(self.variables["u"], self.variables["m"])
        self.variables["h"] += self.dt*ops.calc_dh(self.variables["u"], self.variables["h"])
        self.variables["j"] += self.dt*ops.calc_dj(self.variables["u"], self.variables["j"])

        ina = ops.calc_ina(self.variables["u"], self.variables["m"], self.variables["h"], self.variables["j"], self.parameters["E_Na"], 
                           self.parameters["gna"])

        # Slow inward current:
        self.variables["d"] += self.dt*ops.calc_dd(self.variables["u"], self.variables["d"])
        self.variables["f"] += self.dt*ops.calc_df(self.variables["u"], self.variables["f"])
        
        isi = ops.calc_isk(self.variables["u"], self.variables["d"], self.variables["f"], self.variables["cai"], 
                           self.parameters["gsi"])

        self.variables["cai"] += self.dt*ops.calc_dcai(self.variables["cai"], isi)

        # Time-dependent potassium current:
        self.variables["x"] += self.dt*ops.calc_dx(self.variables["u"], self.variables["x"])

        ik = ops.calc_ik(self.variables["u"], self.variables["x"], self.parameters["ko"], self.parameters["ki"], 
                         self.parameters["nao"], self.parameters["nai"], self.parameters["PR_NaK"], 
                         self.parameters["R"], self.parameters["T"], self.parameters["F"], self.parameters["gk"])



        # Time-independent potassium current:
        ik1 = ops.calc_ik1(self.variables["u"], self.parameters["ko"], self.parameters["E_K1"], self.parameters["gk1"])

        # Plateau potassium current:
        ikp = ops.calc_ikp(self.variables["u"], self.parameters["E_K1"], self.parameters["gkp"])

        # Background current:
        ib = ops.calc_ib(self.variables["u"], self.parameters["gb"])

        # Total time-independent potassium current:
        self.variables["u"] += self.dt*(ops.calc_rhs(ina, 
                                                     isi, 
                                                     ik, 
                                                     ik1, 
                                                     ikp, 
                                                     ib)
                                        + sum(stim.stim(t=self.dt*i) for stim in self.stimulations))

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
            for s in self.variables:
                self.history[s].append(self.variables[s])