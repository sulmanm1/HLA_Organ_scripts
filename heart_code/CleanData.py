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

def data_add_allele(data):
    def format_HLAI(entry):
        HLA_line=entry.split("\n")[-1].replace(" ", "").split(",")
        HLA_dict={"A":[], "B":[], "C":[]}
        for entry in HLA_line:
            if entry.startswith('A'):
                HLA_dict['A'].append(entry.replace("A*", "").split(':')[0])
            elif entry.startswith('B'):
                HLA_dict['B'].append(entry.replace("B*", "").split(':')[0])
            elif entry.startswith('C'):
                HLA_dict['C'].append(entry.replace("C*", "").split(':')[0])
        for key in HLA_dict:
            if len(HLA_dict[key])==1:
                HLA_dict[key].append(HLA_dict[key][0])
        return HLA_dict

    def format_HLAII(entry):
        entry=entry.replace(" ", "")
        if "Not" in entry:
            entry=entry.split("\nNot")[0]
        HLA_line=entry.split("\n")[1].split(",")

        HLA_dict={"DRB1": [], "DRB3": [], "DRB4": [], "DRB5": [], "DPA1": [], "DPB1": [], "DQA1": [], "DQB1": []}
        for entry in HLA_line:
            if entry.startswith('DRB1'):
                HLA_dict['DRB1'].append(entry.replace("DRB1*", "").split(':')[0])
            elif entry.startswith('DRB3'):
                HLA_dict['DRB3'].append(entry.replace("DRB3*", "").split(':')[0])
            elif entry.startswith('DRB4'):
                HLA_dict['DRB4'].append(entry.replace("DRB4*", "").split(':')[0])
            elif entry.startswith('DRB5'):
                HLA_dict['DRB5'].append(entry.replace("DRB5*", "").split(':')[0])
            elif entry.startswith('DPA1'):
                HLA_dict['DPA1'].append(entry.replace("DPA1*", "").split(':')[0])
            elif entry.startswith('DPB1'):
                HLA_dict['DPB1'].append(entry.replace("DPB1*", "").split(':')[0])
            elif entry.startswith('DQA1'):
                HLA_dict['DQA1'].append(entry.replace("DQA1*", "").split(':')[0])
            elif entry.startswith('DQB1'):
                HLA_dict['DQB1'].append(entry.replace("DQB1*", "").split(':')[0])
        for key in HLA_dict:
            if len(HLA_dict[key])==1:
                HLA_dict[key].append(HLA_dict[key][0])
        return HLA_dict

    def matchI(list1, list2):
        list1.sort()
        list2.sort()
        count=0
        if not list1[0]==list2[0]:
            count+=1
        if not list1[1]==list2[1]:
            count+=1
        return count

    def matchII(list1, list2):
        if len(list1)==0 and len(list2)==0:
            return 0
        if len(list1)==0 or len(list2)==0:
            return 2
        else:
            list1.sort()
            list2.sort()
            count = 0
            if not list1[0] == list2[0]:
                count += 1
            if not list1[1] == list2[1]:
                count += 1
            return count

    def HLA_simplify(mm_dict):
        new_dict={}
        for key in mm_dict:
            if key=="DRB3" or key=="DRB4" or key=="DRB5":
                try:
                    new_dict["DRB345"]+=mm_dict[key]
                except:
                    new_dict["DRB345"]=mm_dict[key]
            elif key=="DPA1" or key=="DPB1":
                try:
                    new_dict["DP"]+=mm_dict[key]
                except:
                    new_dict["DP"]=mm_dict[key]
            elif key=="DQA1" or key=="DQB1":
                try:
                    new_dict["DQ"]+=mm_dict[key]
                except:
                    new_dict["DQ"]=mm_dict[key]
            else:
                new_dict[key]=mm_dict[key]
        return new_dict

    def funct_match(row):
        HLAIP_dict=format_HLAI(row['Donor class 1 HLA'])
        HLAID_dict=format_HLAI(row['Patient class 1 HLA'])
        HLAI_match={'A':0, "B":0, "C":0}
        for key in HLAIP_dict:
            HLAI_match[key]=matchI(HLAIP_dict[key], HLAID_dict[key])

        HLAIIP_dict=format_HLAII(row['Donor class 2 HLA'])
        HLAIID_dict=format_HLAII(row['Patient class 2 HLA'])
        HLAII_match={"DRB1": 0, "DRB3": 0, "DRB4": 0, "DRB5": 0, "DPA1": 0, "DPB1": 0, "DQA1": 0, "DQB1": 0}
        for key in HLAIIP_dict:
            HLAII_match[key]=matchII(HLAIIP_dict[key], HLAIID_dict[key])
        HLAII_match=HLA_simplify(HLAII_match)
        return [HLAI_match, HLAII_match]
    data[['HLAI', 'HLAII']] = pd.DataFrame(data.apply(lambda row: funct_match(row), axis=1).tolist(), index=data.index)
    data['A'] = data['HLAI'].apply(lambda x: x.get('A'))
    data['B'] = data['HLAI'].apply(lambda x: x.get('B'))
    data['C'] = data['HLAI'].apply(lambda x: x.get('C'))
    data['DRB1'] = data['HLAII'].apply(lambda x: x.get('DRB1'))
    data['DRB345'] = data['HLAII'].apply(lambda x: x.get('DRB345'))
    data['DP'] = data['HLAII'].apply(lambda x: x.get('DP'))
    data['DQ'] = data['HLAII'].apply(lambda x: x.get('DQ'))
    return data

def med_time(data):
    def averageG1(row):
        return (row["Time upper G1"] + row["Time lower G1"]) /2
    data["Med Time G1"] = data.apply(averageG1, axis=1)
    def averageG2(row):
        return (row["Time upper G2"] + row["Time lower G2"]) /2
    data["Med Time G2"] = data.apply(averageG2, axis=1)
    def averageG3(row):
        return (row["Time upper G3"] + row["Time lower G3"]) /2
    data["Med Time G3"] = data.apply(averageG3, axis=1)
    return data

def main():
    data=pd.read_csv("data2/data.csv")
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

    data=data_add_allele(data)
    data=med_time(data)
    print('creating csv')
    data.to_csv("data2/data_clean.csv", index=False)




print('running')
main()