import pandas as pd
import numpy as np

def data_formatentry(entry):
    if not 'R' in str(entry):
        entry = '9'
    entry= str(entry).replace(" ", "")
    entry= str(entry).replace("Grade", "")
    entry= str(entry).replace("nan", "9")
    entry= str(entry).replace("R", "")
    return int(entry)
def data_modresults(data, columns):
    for column in columns:
        data[column] = data[column].apply(data_formatentry)
    return data
def data_modname(data, first, last):
    data[first] = data[first].apply(lambda x: x.lower().split(" ")[0])
    data[last] = data[last].apply(lambda x: x.lower().split(" ")[0])
    data["name"] = data[first] + " " + data[last]
    data = data.drop_duplicates(subset='name', keep="first")
    return data

def data_removeempty(data, columns):
    mask = (data[columns] == 9).all(axis=1)
    data_filtered = data[~mask]
    return data_filtered

def create_daydict(daylist):
    daydict={}
    for day in daylist:
        daydict[day]=int(day)*7
    return daydict

def time_to_grade(data, daydict, grade):
    colname_upper=f"Time upper G{grade}"
    data[colname_upper] = 999

    def funct_upper(**kwargs):
        for k in kwargs:
            if kwargs[k] >= grade and not kwargs[k] == 9:
                return daydict[k]
        return 999
    data[colname_upper] = data.apply(lambda row: funct_upper(**{k: row[k] for k in daydict}), axis=1)

    status_col=f"Status G{grade}"
    data[status_col] = data[colname_upper].apply(lambda x: 0 if x == 999 else 1)
    data[colname_upper] = data[colname_upper].replace(999, 364)

    ### writing lower bound
    colname_lower=f"Time lower G{grade}"
    data[colname_lower] = np.where(data[status_col] == 0, data[colname_upper], None)
    def funct_lower(row):
        if row[status_col]==0:
            return row[colname_upper]
        else:
            previous=0
            for k in daydict:
                if daydict[k]<row[colname_upper]:
                    previous=daydict[k]
                else:
                    return previous

    data[colname_lower]=data.apply(funct_lower, axis=1)
    return data

def max_grade(data, daydict):
    data["Max Grade"] = 0

    def funct(**kwargs):
        grade_arr = []
        for k in kwargs:
            if not kwargs[k] == 9:
                grade_arr.append(kwargs[k])
        return max(grade_arr)

    data["Max Grade"] = data.apply(lambda row: funct(**{k: row[k] for k in daydict}), axis=1)
    return data

def remempty(data, columns):
    print(data.iloc[data[columns].isnull().any(axis=1).to_string()])

def main():
    data=pd.read_csv("data/data.csv")
    data=data_modname(data,"First name","Last Name")

    result_columns=['1', '2', '3', '4', '6', '8', '10', '12', '26', '52']
    daydict=create_daydict(result_columns)

    data=data_modresults(data,result_columns)
    data=data_removeempty(data,result_columns)

    data=time_to_grade(data,daydict, 3)
    data = time_to_grade(data, daydict, 2)
    data = time_to_grade(data, daydict, 1)



    independent_vars = ['A Mismatch No.', 'B Mismatch No.', 'C Mismatch No.', 'DRB1 Mismatch No.', 'DRB345 Mismatch No.', 'DP Mismatch No.', 'DQ Mismatch No.']
    cumulative_cols = ['ClassI\nUniq# Mismatch\nAll', 'ClassI\n# Mismatch\nAll', 'ClassII\nUniq# Mismatch\nAll', 'ClassII\n# Mismatch\nAll']

    data=data.dropna(subset=independent_vars)
    data = max_grade(data, daydict)

    data.to_csv("data/data_clean.csv", index=False)




if __name__ == "__main__":
    main()