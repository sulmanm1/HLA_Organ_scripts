import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, LogLogisticAFTFitter, WeibullAFTFitter, LogNormalAFTFitter, WeibullFitter, KaplanMeierFitter

def fit_coxph(df, covariates: list, left, right):
    from lifelines.plotting import plot_interval_censored_lifetimes
    ax=plot_interval_censored_lifetimes(df[left], df[right])
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Patient")
    ax.set_title("Grade 3 Interval Censored Data")
    plt.show()

def main():
    data=pd.read_csv("data/data_clean.csv")
    fit_coxph(data, [], "Time lower G3", "Time upper G3")

if __name__ == "__main__":
    main()