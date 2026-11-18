import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import display
import numpy as np
import math

class DynamicModelling:
    def __init__(self, data):
        self.data = data
        # For generated data:
        # Subject_ID (str), Group (str), Time_Days (int), Tumor_Volume_mm3 (float), Body_Weight_g (Float)
        # Or for SIR: Region (str), Time_days (int), Suspectible (int), Infected (int), Recovered (int), Population_Size (int)
    """
    A class driven method to implement different 
    techniques to visualize tumor growth activity data.
    
    Input: "generate_data/sir_dummy_dataset.csv" file containing tumor growth data.
    Output: Visualizations of tumor growth. 
    """


    def plot_data(self, groupby = "Group", x = "Time_Days", y = "Tumor_Volume_mm3"):
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
    Workshop 1: differentuaalvergekijkingen.
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
        t = 0.0 # beginning time
        ts = [t]
        S, I, R = 0.99, 0.01, 0.00
        Ss, Is, Rs = [S], [I], [R]

        for _ in range(N):
            infections = beta * S * I * delta_t
            recoveries = gamma * I * delta_t
            t =  t + delta_t
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

    def Eulers_model(self, groupby="Region", ground="Time_Days", references=["Susceptible", "Infected", "Recovered"], plot=True, show_data=True):
        """
        Calculates per group an Euler-calculation based on exponential growth for
        every reference against the ground, and makes per reference a plot with a
        differente line for every group.

        :param groupby: Column on which will be grouped
        :param ground: Column on which Euler will work (mostly used as time)
        :param references: Columns where the Euler equation will be used on
        :param plot: Make the plots
        :param show_data: Put the raw-data on the plot
        :return: DataFrame with columns [groupby, ground, *references]
        """

        def model(data, y_name, h=1):
            data = data.sort_values(ground)
            t = data[ground].values
            V = data[y_name].values

            if V[0] == 0: # Very bad fix for divide by zero error
                V[0] = 0.01

            r = (np.log(V[-1] / V[0]) / # The average growth rate
                 (t[-1] - t[0])) # Over the (time)span of ground

            time_euler = np.arange(t.min(), t.max() + h, h)
            volume_euler = np.zeros(len(time_euler))
            volume_euler[0] = V[0]

            for i in range(len(time_euler) - 1):
                volume_euler[i + 1] = volume_euler[i] + h * r * volume_euler[i] # Calculates the derivative on every timepoint

            return time_euler, volume_euler

        results = []
        for group, group_data in self.data.groupby(groupby): # For every reference (column specified) it will run the Euler model
            output = {}
            for reference in references:
                time_euler, output[reference] = model(group_data, reference)

            group_df = pd.DataFrame(output)
            group_df.insert(0, ground, time_euler)
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

                ax.set_title(f"Euler-model: {reference} per {groupby}")
                ax.set_xlabel(ground)
                ax.set_ylabel(reference)
                ax.legend(title=groupby)
                ax.grid(alpha=0.3)
                fig.tight_layout()
            plt.show()

        return result


def main():
    # SIR model
    pd_df = pd.read_csv('generated_data/sir_dummy_dataset.csv')
    # display(pd_df)
    dm = DynamicModelling(pd_df)
    # dm.plot_data(groupby='Region', x='Time_Days', y='Infected')
    # dm.plot_SIR()
    dm.Eulers_model(groupby='Region', ground='Time_Days', references=['Susceptible', 'Infected', 'Recovered'])

    # All other models set-up
    pd_df = pd.read_csv('generated_data/tumorgroei_dummy_dataset.csv')
    # # display(pd_df)
    dm = DynamicModelling(pd_df)
    dm.Eulers_model(groupby='Group', references=['Tumor_Volume_mm3', 'Body_Weight_g'])


if __name__ == "__main__":
    main()