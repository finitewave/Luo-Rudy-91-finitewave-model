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
        u_old = self.variables["u"]
        m_old = self.variables["m"]
        h_old = self.variables["h"]
        j_old = self.variables["j"]
        d_old = self.variables["d"]
        f_old = self.variables["f"]
        x_old = self.variables["x"]
        cai_old = self.variables["cai"]

        # Gating derivatives from old state
        dm = ops.calc_dm(u_old, m_old)
        dh = ops.calc_dh(u_old, h_old)
        dj = ops.calc_dj(u_old, j_old)

        dd = ops.calc_dd(u_old, d_old)
        df = ops.calc_df(u_old, f_old)

        dx = ops.calc_dx(u_old, x_old)

        # Currents from old state
        ina = ops.calc_ina(
            u_old,
            m_old,
            h_old,
            j_old,
            self.parameters["E_Na"],
            self.parameters["gna"],
        )

        isi = ops.calc_isk(
            u_old,
            d_old,
            f_old,
            cai_old,
            self.parameters["gsi"],
        )

        ik = ops.calc_ik(
            u_old,
            x_old,
            self.parameters["ko"],
            self.parameters["ki"],
            self.parameters["nao"],
            self.parameters["nai"],
            self.parameters["PR_NaK"],
            self.parameters["R"],
            self.parameters["T"],
            self.parameters["F"],
            self.parameters["gk"],
        )

        ik1 = ops.calc_ik1(
            u_old,
            self.parameters["ko"],
            self.parameters["E_K1"],
            self.parameters["gk1"],
        )

        ikp = ops.calc_ikp(
            u_old,
            self.parameters["E_K1"],
            self.parameters["gkp"],
        )

        ib = ops.calc_ib(
            u_old,
            self.parameters["gb"],
        )

        dcai = ops.calc_dcai(cai_old, isi)

        stim_current = sum(stim.stim(t=self.dt * i) for stim in self.stimulations)

        du = ops.calc_rhs(
            ina,
            isi,
            ik,
            ik1,
            ikp,
            ib,
        ) + stim_current

        # Explicit update
        self.variables["m"] = m_old + self.dt * dm
        self.variables["h"] = h_old + self.dt * dh
        self.variables["j"] = j_old + self.dt * dj

        self.variables["d"] = d_old + self.dt * dd
        self.variables["f"] = f_old + self.dt * df

        self.variables["x"] = x_old + self.dt * dx
        self.variables["cai"] = cai_old + self.dt * dcai

        self.variables["u"] = u_old + self.dt * du

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