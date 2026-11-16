import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import display

class DynamicModelling:
    def __init__(self, data):
        self.data = data
        # Subject_ID (str) Group (str) Time_Days (int) Tumor_Volume_mm3 (float) Body_Weight_g (Float)

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

def main():
    pd_df = pd.read_csv('generated_data/tumorgroei_dummy_dataset.csv')
    display(pd_df)
    dm = DynamicModelling(pd_df)
    dm.plot_data()

if __name__ == "__main__":
    main()