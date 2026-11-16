import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

class DynamicModelling:
    def __init__(self, data):
        self.data = data
        # Subject_ID (str) Group (str) Time_Days (int) Tumor_Volume_mm3 (float) Body_Weight_g (Float)

def main():
    pd_df = pd.read_csv('generated_data/tumorgroei_dummy_dataset.csv')
    DynamicModelling(pd_df)

if __name__ == "__main__":
    main()