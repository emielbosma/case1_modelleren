import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
from random import gauss

def load_tumordata(filepath):
    """
    Reads and prepares the tumor growth dataset.
    """

    df = pd.read_csv(filepath)

    # Rename columns
    df = df.rename(columns={
        "ID": "Subject_ID",
        "DAYS": "Time_Days",
        "TRT": "Group",
        "VOL": "Tumor_Volume"
    })

    # Sort data by subject and measurement time
    df = df.sort_values(
        ["Subject_ID", "Time_Days"]
    ).reset_index(drop=True)

    return df

class DynamicModelling:
    """
    A class driven method to implement different
    techniques to visualize tumor growth activity data.

    Input: "generate_data/sir_dummy_dataset.csv" file containing tumor growth data.
    Output: Visualizations of tumor growth.
    """

    def __init__(self, data):
        self.data = data
        # Subject_ID (str) Group (str) Time_Days (int) Tumor_Volume_mm3 (float) Body_Weight_g (Float)
        # Or for SIR: Region (str), Time_days (int), Susceptible (int), Infected (int), Recovered (int), Population_Size (int)

    def plot_data(self, groupby="Group", x="Time_Days", y="Tumor_Volume_mm3"):
        """
        A function that plots the data that the class was initialized with.

        :param groupby: By which column the data will be divided into different groups and plots.
        :param x: Which column the x-as will be plotted against.
        :param y: Which column the y-as will be plotted against.
        :return: None, plots will be displayed.
        """

        for group, group_data in self.data.groupby(groupby):

            plt.figure(figsize=(8, 4))

            plt.scatter(
                group_data[x],
                group_data[y]
            )

            plt.xlabel(x.replace("_", " "))
            plt.ylabel(y.replace("_", " "))
            plt.title(f'{y.replace("_", " ")} over {x.replace("_", " ")} for {group.replace("_", " ")}')

            plt.grid()
            plt.show()

    """
    Workshop 1: differential equations.
    """

    # SIR-Model

    def SIR_model(self, beta=0.35, gamma=0.10, N=60, delta_t=1.0):
        """
        Parameters:
            beta(float) : Infection rate
            gamma(float) : Recovery rate
            N(int) : total steps integration
            delta_t(float) : time step size
        Returns:
            ts(list[float]) : time steps in days
            Ss (list[float]) : susceptible population
            Is (list[float]) : infected population
            Rs (list[float]) : recovered population
        """
        t = 0.0  # beginning time
        ts = [t]
        S, I, R = 0.99, 0.01, 0.00
        Ss, Is, Rs = [S], [I], [R]

        for _ in range(N):
            infections = beta * S * I * delta_t
            recoveries = gamma * I * delta_t
            t = t + delta_t
            S = S - infections
            I = I + infections - recoveries
            R = R + recoveries

            ts.append(t)
            Ss.append(S)
            Is.append(I)
            Rs.append(R)

        return ts, Ss, Is, Rs

    def plot_SIR(self, region="Region_North"):
        df = self.data[self.data['Region'] == region]
        ts, Ss, Is, Rs = self.SIR_model()

        plt.plot(ts, Ss, "b-", label="Susceptible")
        plt.plot(ts, Is, "r-", label="Infected")
        plt.plot(ts, Rs, "g-", label="Recovered")
        plt.plot(df['Time_Days'], df['Infected'], "ko", label="Observed Infected")
        plt.xlabel("Time (days)")
        plt.ylabel("Population")
        plt.title("SIR Model for North Region Tumor growth")
        plt.legend()
        plt.grid()
        plt.show()

    """
    Workshop 2: ODE-solvers
    """

    # Default parameters per growth model. 'c' (and for Von Bertalanffy also 'd')
    # are rate parameters; the others determine the shape of the curve.
    MODEL_DEFAULTS = {
        "lineair": {"c": 1.0},
        "exponentieel": {"c": 1.0},
        "mendelsohn": {"c": 1.0, "d": 0.5},
        "exponentieel_afvlakkend": {"c": 1.0, "Vmax": 1.0},
        "logistisch": {"c": 1.0, "Vmax": 1.0},
        "montroll": {"c": 1.0, "Vmax": 1.0, "d": 0.5},
        "allee": {"c": 1.0, "Vmin": 0.0, "Vmax": 1.0},
        "lineair_gelimiteerd": {"c": 1.0, "d": 1.0},
        "oppervlakte_gelimiteerd": {"c": 1.0, "d": 1.0},
        "von_bertalanffy": {"c": 1.0, "d": 1.0},
        "gompertz": {"c": 1.0, "Vmax": 1.0},
    }

    # Parameters that scale linearly with the growth rate (used when
    # estimating starting values from the data)
    RATE_PARAMS = {"von_bertalanffy": ["c", "d"]}

    def mathematical_methods(self, method, **kwargs):
        """
        Returns the differential equation dV/dt = f(t, V) of a growth model.

        Models (V = volume/size, t = time):
            lineair                  dV/dt = c
            exponentieel             dV/dt = c * V
            mendelsohn               dV/dt = c * V^d
            exponentieel_afvlakkend  dV/dt = c * (Vmax - V)
            logistisch               dV/dt = c * V * (Vmax - V)
            montroll                 dV/dt = c * V * (Vmax^d - V^d)
            allee                    dV/dt = c * (V - Vmin) * (Vmax - V)
            lineair_gelimiteerd      dV/dt = c * V / (V + d)
            oppervlakte_gelimiteerd  dV/dt = c * V / (V + d)^(1/3)
            von_bertalanffy          dV/dt = c * V^(2/3) - d * V
            gompertz                 dV/dt = c * V * ln(Vmax / V)

        :param method: Name of the model (see above, case-insensitive)
        :param kwargs: Model parameters (c, d, Vmax, Vmin); missing
                       parameters get their value from MODEL_DEFAULTS
        :return: Function f(t, V) that returns dV/dt
        """
        method = method.lower().replace(" ", "_").replace("-", "_")
        if method not in self.MODEL_DEFAULTS:
            raise ValueError(
                f"Unknown model '{method}'. Choose from: {list(self.MODEL_DEFAULTS)}"
            )

        p = {**self.MODEL_DEFAULTS[method], **kwargs}
        c = p["c"]
        d = p.get("d")
        Vmax = p.get("Vmax")
        Vmin = p.get("Vmin")
        pos = lambda V: max(V, 1e-12)  # prevents negative powers / log(0)

        models = {
            "lineair": lambda t, V: c,
            "exponentieel": lambda t, V: c * V,
            "mendelsohn": lambda t, V: c * pos(V) ** d,
            "exponentieel_afvlakkend": lambda t, V: c * (Vmax - V),
            "logistisch": lambda t, V: c * V * (Vmax - V),
            "montroll": lambda t, V: c * V * (Vmax ** d - pos(V) ** d),
            "allee": lambda t, V: c * (V - Vmin) * (Vmax - V),
            "lineair_gelimiteerd": lambda t, V: c * V / (V + d),
            "oppervlakte_gelimiteerd": lambda t, V: c * V / pos(V + d) ** (1 / 3),
            "von_bertalanffy": lambda t, V: c * pos(V) ** (2 / 3) - d * V,
            "gompertz": lambda t, V: c * V * np.log(Vmax / pos(V)),
        }
        return models[method]

    def estimate_params(self, method, t, y):
        """
        Estimates simple starting parameters for a growth model from the data.

        Shape parameters (Vmax, Vmin, d) are derived from the values in the data.
        The rate parameters are then scaled so that the slope of the model at the
        midpoint (V_mid) equals the average slope of the data.
        This is a rough estimate; optimization (Workshop 3) can improve it.

        :param method: Name of the model
        :param t: Array of time points
        :param y: Array of measured values
        :return: Dict with parameters for mathematical_methods
        """
        method = method.lower().replace(" ", "_").replace("-", "_")
        params = dict(self.MODEL_DEFAULTS[method])
        increasing = y[-1] >= y[0]

        # Shape parameters from the data
        if "Vmax" in params:
            params["Vmax"] = 1.1 * y.max() if increasing else 0.9 * y.min()
        if "Vmin" in params:
            params["Vmin"] = 0.5 * y.min()
        if method in ("lineair_gelimiteerd", "oppervlakte_gelimiteerd"):
            params["d"] = float(np.median(y))
        if method == "von_bertalanffy":
            # Equilibrium Vmax = (c/d)^3, so c = d * Vmax^(1/3)
            vmax = 1.1 * y.max() if increasing else 0.9 * y.min()
            params["d"] = 1.0
            params["c"] = max(vmax, 1e-12) ** (1 / 3)

        # Scale rate parameters to the average slope of the data
        slope = (y[-1] - y[0]) / (t[-1] - t[0])
        V_mid = (y[0] + y[-1]) / 2
        f_unit = self.mathematical_methods(method, **params)(t[0], V_mid)
        k = slope / f_unit if np.isfinite(f_unit) and f_unit != 0 else 1.0
        for key in self.RATE_PARAMS.get(method, ["c"]):
            params[key] *= k

        # Exponential: exact estimate r = ln(y_end / y_start) / T
        if method == "exponentieel" and y[0] > 0 and y[-1] > 0:
            params["c"] = np.log(y[-1] / y[0]) / (t[-1] - t[0])
        return params

    def ODE_methods(self, model, groupby="Region", ground="Time_Days",
                    references=["Susceptible", "Infected", "Recovered"], plot=True, show_data=True,
                    growth_model="exponentieel", params=None):
        """
        Calculates per group an ODE-solution (Euler, Heun or Runge-Kutta) based on
        the chosen growth model for every reference against the ground, and makes
        per reference a plot with a different line for every group.

        :param model: ODE-solver: 'Euler', 'Heun' or 'Runge-Kutta'
        :param groupby: Column on which will be grouped
        :param ground: Column on which Euler will work (mostly used as time)
        :param references: Columns where the Euler equation will be used on
        :param plot: Make the plots
        :param show_data: Put the raw-data on the plot
        :param growth_model: Growth model from mathematical_methods (e.g. 'logistisch')
        :param params: Optional dict with model parameters; missing parameters
                       are estimated from the data per group/reference
        :return: DataFrame with columns [groupby, ground, *references]
        """

        def prepare(data, y_name, h):
            data = data.sort_values(ground)
            t = data[ground].values.astype(float)
            y = data[y_name].values.astype(float)
            if y[0] == 0:  # Very bad fix for divide by zero error
                y[0] = 0.01
            # Estimate parameters from the data, overwrite with user-given ones
            model_params = self.estimate_params(growth_model, t, y)
            model_params.update(params or {})
            f = self.mathematical_methods(growth_model, **model_params)  # dy/dt = f(t, y)
            ts = np.arange(t.min(), t.max() + h, h)
            ys = np.zeros(len(ts))
            ys[0] = y[0]
            return f, ts, ys

        def euler_model(data, y_name, h=1):
            f, ts, ys = prepare(data, y_name, h)
            for i in range(len(ts) - 1):
                ys[i + 1] = ys[i] + h * f(ts[i], ys[i])
            return ts, ys

        def heun_model(data, y_name, h=1):
            f, ts, ys = prepare(data, y_name, h)
            for i in range(len(ts) - 1):
                k1 = f(ts[i], ys[i])
                k2 = f(ts[i] + h, ys[i] + h * k1)
                ys[i + 1] = ys[i] + h / 2 * (k1 + k2)
            return ts, ys

        def rungekutta_model(data, y_name, h=1):
            f, ts, ys = prepare(data, y_name, h)
            for i in range(len(ts) - 1):
                k1 = f(ts[i], ys[i])
                k2 = f(ts[i] + h / 2, ys[i] + h / 2 * k1)
                k3 = f(ts[i] + h / 2, ys[i] + h / 2 * k2)
                k4 = f(ts[i] + h, ys[i] + h * k3)
                ys[i + 1] = ys[i] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            return ts, ys

        results = []
        for group, group_data in self.data.groupby(groupby): # For every reference (column specified) it will run the specified model
            output = {}
            for reference in references:
                if model == 'Euler':
                    ts_model, output[reference] = euler_model(group_data, reference)
                elif model == 'Heun':
                    ts_model, output[reference] = heun_model(group_data, reference)
                elif model == 'Runge-Kutta':
                    ts_model, output[reference] = rungekutta_model(group_data, reference)

            group_df = pd.DataFrame(output)
            group_df.insert(0, ground, ts_model)
            group_df.insert(0, groupby, group)
            results.append(group_df)

        result = pd.concat(results, ignore_index=True)

        if plot:
            for reference in references:
                fig, ax = plt.subplots(figsize=(8, 5))
                for group, group_df in result.groupby(groupby):
                    line, = ax.plot(group_df[ground],
                                    group_df[reference],
                                    label=str(group))

                    if show_data: # To show the data points on the plot
                        obs = self.data[self.data[groupby] == group]
                        ax.scatter(obs[ground], obs[reference], s=12,
                                   color=line.get_color(), alpha=0.6)

                ax.set_title(f"{model} ({growth_model}): {reference} per {groupby}")
                ax.set_xlabel(ground)
                ax.set_ylabel(reference)
                ax.legend(title=groupby)
                ax.grid(alpha=0.3)
                fig.tight_layout()
            plt.show()

        return result

    """
    Workshop 3: Optimization
    """

    def _solve_ode(self, f, t_eval, y0, method="Runge-Kutta", h=1.0):
        """
        Solves dV/dt = f(t, V) numerically and returns the solution at the time points t_eval.

        :param f: Function f(t, V) from mathematical_methods
        :param t_eval: Time points at which the solution is needed
        :param y0: Starting value at t_eval.min()
        :param method: 'Euler', 'Heun' or 'Runge-Kutta'
        :param h: Step size of the solver
        :return: Array with model values at t_eval (NaN if the solution blows up)
        """

        t_eval = np.asarray(t_eval, dtype=float)
        ts = np.arange(t_eval.min(), t_eval.max() + h, h) # Start, stop, step
        ys = np.empty(len(ts))
        ys[0] = y0

        with np.errstate(over="ignore", invalid="ignore"):
            for i in range(len(ts) - 1):
                t, V = ts[i], ys[i]
                if method == "Euler":
                    ys[i + 1] = V + h * f(t, V)
                elif method == "Heun":
                    k1 = f(t, V)
                    k2 = f(t + h, V + h * k1)
                    ys[i + 1] = V + h / 2 * (k1 + k2)
                elif method == "Runge-Kutta":
                    k1 = f(t, V)
                    k2 = f(t + h / 2, V + h / 2 * k1)
                    k3 = f(t + h / 2, V + h / 2 * k2)
                    k4 = f(t + h, V + h * k3)
                    ys[i + 1] = V + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
                else:
                    raise ValueError("method must be 'Euler', 'Heun' or 'Runge-Kutta'")

                # Stop as soon as the solution becomes unstable
                if not np.isfinite(ys[i + 1]) or abs(ys[i + 1]) > 1e12:
                    return np.full(len(t_eval), np.nan)

        return np.interp(t_eval, ts, ys)

    def optimization(self, model="exponentieel", method="Runge-Kutta", groupby="TRT",
                     ground="DAYS", references=["VOL"], plot=False,
                     h=1.0, sigma=0.05, patience=500, max_iter=20000, params=None):
        """
        Fits a growth model from mathematical_methods per group to the data
        with Random Search (Workshop 3) and optionally plots the result.

        :param model: Growth model, a key of MODEL_DEFAULTS (e.g. 'gompertz')
        :param method: ODE-solver: 'Euler', 'Heun' or 'Runge-Kutta'
        :param groupby: Column on which will be grouped (e.g. treatment)
        :param ground: Time column
        :param references: List with the column to fit (only the first is used)
        :param plot: Make a plot of the data and the fitted curves
        :param h: Step size of the ODE-solver
        :param sigma: Relative step size of the random search
        :param patience: Stop after this many attempts without improvement
        :param max_iter: Maximum total number of attempts per group
        :param params: Optional dict with fixed starting values (overwrites the estimate)
        :return: dict with
                 'params' -> {group: optimized parameters}
                 'mse'    -> {group: MSE of that group}
                 'sse'    -> total sum of squared errors over all groups
                 'n'      -> total number of data points
                 'k'      -> total number of fitted parameters (params per group x groups)
        """
        model = model.lower().replace(" ", "_").replace("-", "_")
        y_name = references[0]

        def calc_mse(p):
            f = self.mathematical_methods(model, **p)
            y_model = self._solve_ode(f, x, y0, method=method, h=h)
            mse = np.mean((y - y_model) ** 2)
            return mse

        if plot:
            fig, ax = plt.subplots(figsize=(8, 5))

        fitted_params, fitted_mse = {}, {}
        sse_total, n_total, k_total = 0.0, 0, 0

        for group, data in self.data.groupby(groupby):
            data = data.sort_values(ground)
            x = data[ground].values
            y = data[y_name].values

            # Starting params
            means = data.groupby(ground)[y_name].mean()
            t_mean = means.index.values
            y_mean = means.values
            y0 = max(y_mean[0], 1e-2)  # Once again; a very bad fix for divide by zero
            y_mean[0] = y0

            best = self.estimate_params(model, t_mean, y_mean)
            best.update(params or {})
            best = {key: float(val) for key, val in best.items()}
            mse = calc_mse(best)

            # Simple Random Search adapted from Workshop 3-HTML
            tries, iterations = 0, 0
            while tries < patience and iterations < max_iter:
                # Because comparison runs this function as many times as the models we have, I added a max_iter
                # TODO: In the jupyter notebook we can maybe run the 3 best models with an bigger max_iter

                tries += 1
                iterations += 1

                new_params = {
                    key: val + gauss(0, sigma * max(abs(val), 1e-3))
                    for key, val in best.items()
                }
                new_mse = calc_mse(new_params)

                if new_mse < mse:
                    best, mse = new_params, new_mse
                    tries = 0

            fitted_params[group] = best
            fitted_mse[group] = mse
            sse_total = mse * len(y)
            n_total = len(y)
            k_total = len(best)

            if plot: # Added plot behind an if because I was tired of all those plots
                t_fine = np.linspace(x.min(), x.max(), 200)
                f = self.mathematical_methods(model, **best)
                y_fine = self._solve_ode(f, t_fine, y0, method=method, h=h)
                line, = ax.plot(t_fine, y_fine, label=f"{group} (MSE={mse:.3g})")
                ax.scatter(x, y, s=12, alpha=0.5, color=line.get_color())

        if plot:
            ax.set_xlabel(ground)
            ax.set_ylabel(y_name)
            ax.set_title(f"Optimized tumor growth per treatment: {model} ({method})")
            ax.legend(title=groupby)
            ax.grid(alpha=0.5)
            fig.tight_layout()
            plt.show()

        return {"params": fitted_params, "mse": fitted_mse,
                "sse": sse_total, "n": n_total, "k": k_total}

    """
    Workshop 4: model comparison
    """

    def model_comparison(self, groupby="TRT", ground="DAYS", references=["VOL"],
                         method="Runge-Kutta", models=None, plot=False, **kwargs):
        """
        Fits every growth model with plot_optimization and compares them with
        AIC, AICc and BIC (lower = better).

            AIC  = n * ln(SSE / n) + 2k
            AICc = AIC + 2k(k + 1) / (n - k - 1)
            BIC  = n * ln(SSE / n) + k * ln(n)

        with n = number of data points and k = number of fitted parameters.

        :param models: Optional list of models; default all models in MODEL_DEFAULTS
        :param plot: Plot the fit of every model
        :param kwargs: Passed on to plot_optimization (e.g. h, sigma, patience)
        :return: DataFrame per model with n, k, MSE, AIC, AICc and BIC, sorted on AICc
        """
        models = models or list(self.MODEL_DEFAULTS.keys())
        rows = {}

        for model in models:
            fit = self.optimization(model=model, method=method, groupby=groupby,
                                    ground=ground, references=references,
                                    plot=plot, **kwargs)
            n, k, sse = fit["n"], fit["k"], fit["sse"]
            mse = sse / n

            AIC = n * math.log(mse) + 2 * k
            AICc = AIC + (2 * k * (k + 1)) / (n - k - 1) if n - k - 1 > 0 else np.inf
            BIC = n * math.log(mse) + k * math.log(n)
            rows[model] = {"n": n, "k": k, "MSE": mse, "AIC": AIC, "AICc": AICc, "BIC": BIC}

        results = pd.DataFrame.from_dict(rows, orient="index").sort_values("AICc")
        return results


def main():
    # # Load tumor growth dataset
    # pd_df = load_tumordata(
    #     'generated_data/peerj-07-6983-s001.csv'
    # )
    # # SIR model
    # pd_df = pd.read_csv('generated_data/sir_dummy_dataset.csv')
    # # display(pd_df)
    # dm = DynamicModelling(pd_df)
    # dm.plot_data(groupby='Region', x='Time_Days', y='Infected')
    # # dm.plot_SIR()
    # dm.ODE_methods(model='Euler', groupby='Region', ground='Time_Days',
    #                references=['Susceptible', 'Infected', 'Recovered'])
    #
    # # All other models set-up
    # pd_df = pd.read_csv('generated_data/tumorgroei_dummy_dataset.csv')
    # # # display(pd_df)
    # dm = DynamicModelling(pd_df)
    # dm.ODE_methods(model='Euler', groupby='Group', references=['Tumor_Volume_mm3', 'Body_Weight_g'])
    #
    # # Different growth model with Runge-Kutta, optionally with custom parameters
    # dm.ODE_methods(model='Runge-Kutta', groupby='Group', references=['Tumor_Volume_mm3'],
    #                growth_model='gompertz')
    # # dm.ODE_methods(model='Heun', groupby='Group', references=['Tumor_Volume_mm3'],
    # #                growth_model='logistisch', params={'Vmax': 2000})

    # Workshop 3: PeerJ dataset
    pd_df = pd.read_csv('generated_data/peerj-07-6983-s001.csv')
    dm = DynamicModelling(pd_df)

    # dm.plot_optimization(model='gompertz', method='Runge-Kutta',
    #                      groupby='TRT', ground='DAYS', references=['VOL'], plot=True)
    # dm.plot_optimization_exp(groupby='TRT', ground='DAYS', references=['VOL'])  # Exponential model

    # Workshop 4: model comparison
    comparison = dm.model_comparison(groupby='TRT', ground='DAYS', references=['VOL'],
                                     method='Runge-Kutta')
    print(comparison.to_string())

if __name__ == "__main__":
    main()