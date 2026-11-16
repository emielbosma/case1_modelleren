import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

"""
A class driven method to implement different 
techniques to visualize tumor growth activity data.

Input: "generate_data/sir_dummy_dataset.csv" file containing tumor growth data.
Output: Visualizations of tumor growth. 
"""


class DynamicModelling:
    def __init__(self, data):
        self.data = data
        # Subject_ID (str) Group (str) Time_Days (int) Tumor_Volume_mm3 (float) Body_Weight_g (Float)

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

def main():
    pd_df = pd.read_csv('case1_modelleren/generated_data/sir_dummy_dataset.csv')
    model = DynamicModelling(pd_df)
    model.plot_SIR()


if __name__ == "__main__":
    main()