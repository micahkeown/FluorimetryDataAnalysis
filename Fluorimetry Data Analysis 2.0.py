import numpy as np, pandas as pd, xml.etree.ElementTree as ET, re, matplotlib.pyplot as plt, seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import plotly.express as px
import dash
from dash import dcc, html
from dash.dependencies import Input, Output



filepath = r'C:\Users\micah\Downloads'
'''lipo_data = '013025 1 in 10 Lipo227 in DMSO MK Analyzed.xml'
FTPAL_data = '020425 1 in 10 F(ab)2 TPAL in DMSO MK Analyzed.xml'
linker_data = '013125 1 in 50 Fab2 Linker Sample 1 MK2 Analyzed.xml'

fluor_TPAL = 1600
fluor_0Lig = 600
y = [600, 600, 600, 800, 1400, 2800, 5700, 11000, 23000, 47000, 100000, 600]
#y = [600, 600, 600, 600, 800, 1400, 2800, 5700, 11000, 23000, 47000, 100000]
y.sort()
repeat = len(y)
diameter = 98.36

#FOR CHECKING DATA ONLY
#interpolated_TPAL = linker concentration of TPAL
#interpolated_0Lig = linker concentration of 0Lig
#lipo = BPD concentration of 0Lig
#FTPAL = BPD concentration of TPAL

#CHANGE BACK AFTER CHECKING
#interpolated_0Lig = 0.9727
#interpolated_TPAL = 3.243
#FTPAL = 367.65
#lipo = 3647.394

#DELETE WHEN ACTUALLY RUNNING'''

'''lipo_data = '021225 1 in 10 zero ligand lipo in DMSO RT3.xml'
FTPAL_data = '021225 1 in 10 mAb TPAL IRDye800 in DMSO RT4 (1).xml'
linker_data = '021225 1 in 50 mab linker IRDye800 in PBS 6.5 RT.xml'

fluor_TPAL = 570
fluor_0Lig = 390
#y = [2940000, 1730000.0, 1161030.0, 562540.0, 286610.0, 109360.0, 61650.0, 37240.0, 26790.0, 19330.0]
y = [234770,94560,36500,13490,6930,3310,1820,890,510,470,380,350]
y.sort()
diameter = 101.6
repeat = len(y)
#f'''


'''lipo_data = '021225 1 in 10 zero ligand lipo in DMSO RT3.xml'
FTPAL_data = '021225 1 in 10 mAb TPAL  AF in DMSO RT5.xml'
linker_data = '021225 1 in 50 mab linker AF in PBS 6.5 RT2.xml'

fluor_TPAL = 76330
fluor_0Lig = 25860
#y = [2940000, 1730000.0, 1161030.0, 562540.0, 286610.0, 109360.0, 61650.0, 37240.0, 26790.0, 19330.0]
y = [2940000, 1730000.0, 1161030.0, 562540.0, 286610.0, 109360.0, 61650.0, 37240.0, 26790.0, 19330.0]
y.sort()
diameter = 101.6
print(len(y))
repeat = len(y)
#'''

lipo_data = input('Enter Data File Name for 0-Ligand Liposome: ') + '.xml'
FTPAL_data = input('Enter Data File Name for TPAL: ') + '.xml'
linker_data = input('Enter Data File Name for Linker Group: ') + '.xml'


def xml_getter(filepath, file):
    data = filepath + r'/' + file
    
    # Parse the XML
    tree = ET.parse(data)
    root = tree.getroot()

    # Define namespace
    namespace = {"ss": "urn:schemas-microsoft-com:office:spreadsheet"}

    # Locate the table
    table = root.find(".//ss:Table", namespace)

    # Extract rows
    data = []
    for row in table.findall("ss:Row", namespace):
        cells = [cell.find("ss:Data", namespace).text if cell.find("ss:Data", namespace) is not None else "" 
             for cell in row.findall("ss:Cell", namespace)]
        data.append(cells)

    # Convert to Pandas DataFrame
    df = pd.DataFrame(data)
    return df

def preprocessing(df):
    # Set the first row as column headers
    df.columns = df.iloc[0]  # Use the first row as header
    df = df[1:].reset_index(drop=True)  # Remove header row from data
    df = df[['Wavelength (nm)', 'Absorbance']]  # Keep only the columns we need
    df = df[['Wavelength (nm)', 'Absorbance']].astype(float)
    
    df2 = pd.DataFrame(columns = df.columns)
    for i in range(0,799):
        df2.loc[i, :] = df.loc[800 - i, :]
    return df2

def max_abs_lipo(df2):
    abs_array = []
    for i in range(1, 20):
        abs_array.append(df2.loc[i + 480, 'Absorbance'])

    oLig_abs = np.max(abs_array)
    oLig_abs = oLig_abs.astype(float)
    print('Max Absorbance:', oLig_abs)
    return(oLig_abs)

def max_abs_linker(df2):
    abs_array = []
    for i in range(1, 20):
        abs_array.append(df2.loc[i + 70, 'Absorbance'])

    oLig_abs = np.max(abs_array)
    oLig_abs = oLig_abs.astype(float)
    print('Max Absorbance:', oLig_abs)
    return(oLig_abs)

def dil_factor(file):
    # Define the patterns to search for '1 in [number]'
    patterns = [r'1 in (\d+)', r'1in(\d+)']

    for pattern in patterns:
        match = re.search(pattern, file)
        if match:
            dil_factor = float(match.group(1))  # Extract the number and convert to float
            print('Dilution Factor:', dil_factor)
            return dil_factor  # Return immediately if a match is found

    print("Pattern not found")
    return None  # Return None if no match is found


def concentration(abs, file, dil_factor, df2):
    background = np.min(df2[df2['Absorbance'] > 0]['Absorbance'])
    print('Background Absorbance:', background)
    lipo = ['Lipo', '0Lig', '0-Lig', '0lig', '0-lig', '0 ligand', 
            'zero ligand', '0-ligand', '0 ligand', 'zero ligand', 'Zero Ligand']
    combined_lipo = '|'.join(lipo)
    TPAL = ['TPAL', 'tpal', 'TPAl', 'Tpal']
    combined_TPAL = '|'.join(TPAL)
    #combined_FTPAL = 'TPAL'
    Fab2 = ['Fab2 Linker', 'Fab2 Link', 'Fab2 Conjugate', 'Fab2Linker', 'F(ab)2 Linker']
    combined_Fab2 = '|'.join(Fab2)
    mAb = ['mAb Linker', 'MAB Linker', 'mAb linker', 'mab linker','mAb Conjugate']
    combined_mAb = '|'.join(mAb)
    match_lipo = re.search(combined_lipo, file)
    match_TPAL = re.search(combined_TPAL, file)
    match_Fab2 = re.search(combined_Fab2, file)
    match_mAb = re.search(combined_mAb, file)
    
    if match_lipo:
        BPD_exc = (34895)
        con = (abs - background) / BPD_exc  * 10**6 * dil_factor / 50 * 1000
        print('Concentration (nM):', con)
    elif match_TPAL:
        BPD_exc = (34895)
        con = (abs - background) / BPD_exc  * 10**6 * dil_factor / 50 * 1000
        print('Concentration (nM):', con)
    elif match_Fab2:
        protein_exc = (154000)
        
        con = (abs - background) / protein_exc  * 10**6 * dil_factor / 50 * 1000 /4
        print('Concentration (nM):', con)
        #print(abs, background, abs-background)
    elif match_mAb:
        protein_exc = (210000)
        con = (abs - background) / protein_exc  * 10**6 * dil_factor / 50 * 1000 / 4
        print('Concentration (nM):', con)
    else:
        print('Pattern not found')
    return(con)

def concentration_TPAL(abs, file, dil_factor, df2):
    background = np.min(df2[df2['Absorbance'] > 0]['Absorbance'])
    print('Background Absorbance:', background)
    
    BPD_exc = (34895)
    con = (abs - background) / BPD_exc  * 10**6 * dil_factor / 50 * 1000
    print('Concentration (nM):', con)
    return con
    
def concentration_linker(abs, file, dil_factor, df2):
    background = np.min(df2[df2['Absorbance'] > 0]['Absorbance'])
    print('Background Absorbance:', background)
    
    fab_vs_mab = input('Does the linker group include F(ab)2 or mAb? Put "F" for F(ab)2 and "M" for mAb: ')
    if fab_vs_mab == 'F':
        protein_exc = (154000)
        
        con = (abs - background) / protein_exc  * 10**6 * dil_factor / dil_factor * 1000 /4
        print('Concentration (nM):', con)
        
    elif fab_vs_mab == 'M':
        protein_exc = (210000)
        con = (abs - background) / protein_exc  * 10**6 * dil_factor / dil_factor * 1000 / 4
        print('Concentration (nM):', con)
    return con

def concentration_0Lig(abs, file, dil_factor, df2):
    background = np.min(df2[df2['Absorbance'] > 0]['Absorbance'])
    print('Background Absorbance:', background)
    
    BPD_exc = (34895)
    con = (abs - background) / BPD_exc  * 10**6 * dil_factor / 50 * 1000
    print('Concentration (nM):', con)
    return con
        
        
def linker_con(linker, repeat):
    linkers = []
    linker1 = linker
    linkers.append(linker1)


    linker1 = linker
    for i in range(0,repeat - 1):
        linker = linker/2
        linkers.append(linker)
    
    return(linkers) 

def get_y_features(y, fluor_TPAL, len_neg, len_pos):
    positive = []
    negative = []

    difference = (np.array(y) - fluor_TPAL)
    for i in range(len(difference)):
        if difference[i] > 0:
            positive.append(difference[i])
        
    for i in range(len(difference)):
        if difference[i] < 0:
            negative.append(difference[i])
    
    positive = sorted(positive)
    negative = sorted(negative)
    negative = negative[::-1]
    
    #print(f'Positive: {positive}\nNegative: {negative}')
    
    dif1 = positive[0:len_pos]
    dif2 = negative[0:len_neg]

    for _ in range(len(dif2)):
        dif1.append(dif2[_])
    

    dif1.sort()

    for i in range(len(dif1)):
        dif1[i] = dif1[i] + fluor_TPAL

    return dif1

def get_x_features(x, dif1, y, linkers): 
    x = np.array(x).reshape(-1,1)
    linkers = linkers[::-1]
    indnew = []
    for i in range(len(y)):
        if dif1[0] == y[i]:
            ind = i
            break
            #print(i)
    
        #indnew = [ind, ind + 1, ind+2, ind+3]
    
    xnew = []
    for o in range(len(dif1)):
        #xnew = [linkers[ind], linkers[ind+1], linkers[ind+2], linkers[ind+3]]
        if ind + o < len(linkers):
            xnew.append(linkers[ind + o])
    
    return xnew

def conjugation_calc(interpolated_0Lig, lipo, interpolated_TPAL, FTPAL, diameter):
    C613_lipo = interpolated_0Lig / lipo
    #print(C613_lipo)
    C613_TPAL = interpolated_TPAL / FTPAL
    #print(C613_TPAL)
    C619 = 17.69 * (((diameter / 2)**2) + (((diameter/2) - 5) ** 2))

    C614 = C613_TPAL - C613_lipo
    #print(C614)
    C615 = C614 * FTPAL
    #print(C615)

    C617 = C615 / FTPAL

    C621 = (0.6 / 100) * C619
    C622 = C617 * C621
    Conjugation_Ratio = C622
    

    C623 = C622 / 30 * 100
    Conjugation_Efficiency = C623
    return Conjugation_Ratio, Conjugation_Efficiency

def Regression_Analysis(xnew, ynew):
    lr = LinearRegression()
    lr.fit(xnew, ynew)
    yhat = lr.predict(xnew)
    #print(f'Best fit line: y = {lr.coef_}x + {lr.intercept_}')
    from sklearn.metrics import r2_score
    
    #print(f'R^2 value of best fit line: {r2_score(yhat, ynew)}')
    r2 = r2_score(yhat, ynew)
    '''plt.scatter(xnew, ynew, color = 'b')
    plt.plot(xnew, yhat, color = 'r')
    plt.show()'''
    return lr, r2

def stuff_for_best_regression_analysis(x11, y11):
    ynew = get_y_features(y, fluor_TPAL,x11, y11)
    dif = ynew
    xnew = get_x_features(x, dif, y, linkers)
    xnew = np.array(xnew).reshape(-1,1)
    lr, r2= Regression_Analysis(xnew, ynew)
    
    yhat22 = lr.predict(xnew)
    
    return xnew, ynew, lr, r2, yhat22

def best_regression_analysis():
    r2_list = []
    #For 2 above, 2 below
   
    
    
    xnew22, ynew22, lr22, r2_22, yhat22 = stuff_for_best_regression_analysis(2,2)
    r2_list.append(r2_22)
    xnew33, ynew33, lr33, r2_33, yhat33 = stuff_for_best_regression_analysis(3,3)
    r2_list.append(r2_33)
    
    xnew44, ynew44, lr44, r2_44, yhat44 = stuff_for_best_regression_analysis(4,4)
    r2_list.append(r2_44)
    xnew23, ynew23, lr23, r2_23, yhat23 = stuff_for_best_regression_analysis(2,3)
    r2_list.append(r2_23)
    xnew32, ynew32, lr32, r2_32, yhat32 = stuff_for_best_regression_analysis(3,2)
    r2_list.append(r2_32)
    xnew43, ynew43, lr43, r2_43, yhat43 = stuff_for_best_regression_analysis(4,3)
    r2_list.append(r2_43)
    xnew34, ynew34, lr34, r2_34, yhat34 = stuff_for_best_regression_analysis(3,4)
    r2_list.append(r2_34)
    xnew13, ynew13, lr13, r2_13, yhat13 = stuff_for_best_regression_analysis(1,3)
    r2_list.append(r2_13)
    xnew31, ynew31, lr31, r2_31, yhat31 = stuff_for_best_regression_analysis(3,1)
    r2_list.append(r2_31)
    
    
    if (np.max(r2_list) == r2_list[0]):
        xbest = xnew22
        ybest = ynew22
        yhatbest = yhat22
        best = r2_list[0]
        bestlr = lr22
    elif (np.max(r2_list) == r2_list[1]):
        xbest = xnew33
        ybest = ynew33
        yhatbest = yhat33
        best = r2_list[1]
        bestlr = lr33
        
    elif (np.max(r2_list) == r2_list[2]):
        xbest = xnew44
        ybest = ynew44
        yhatbest = yhat44
        best = r2_list[2]
        bestlr = lr44
    elif (np.max(r2_list) == r2_list[3]):
        xbest = xnew23
        ybest = ynew23
        yhatbest = yhat23
        best = r2_list[3]
        bestlr = lr23
    elif (np.max(r2_list) == r2_list[4]):
        xbest = xnew32
        ybest = ynew32
        yhatbest = yhat32
        best = r2_list[4]
        bestlr = lr32
    elif (np.max(r2_list) == r2_list[5]):
        xbest = xnew43
        ybest = ynew43
        yhatbest = yhat43
        best = r2_list[5]
        bestlr = lr43
    elif (np.max(r2_list) == r2_list[6]):
        xbest = xnew34
        ybest = ynew34
        yhatbest = yhat34
        best = r2_list[6]
        bestlr = lr34
    elif (np.max(r2_list) == r2_list[7]):
        xbest = xnew13
        ybest = ynew13
        yhatbest = yhat13
        best = r2_list[7]
        bestlr = lr13
    elif (np.max(r2_list) == r2_list[8]):
        xbest = xnew31
        ybest = ynew31
        yhatbest = yhat31
        best = r2_list[8]
        bestlr = lr31
    
    plt.scatter(xbest, ybest, color = 'b')
    plt.plot(xbest, yhatbest, color = 'r')
    plt.show()
    
    print(f'X values: {xbest}\nY values: {ybest}')
    print(f'\nr2 List: {r2_list}\n')
    print(f'Highest r2 value: {np.max(r2_list)}, Chosen r2 value: {best}, Position in list: {r2_list.index(max(r2_list))}')
    return bestlr

df_lipo = xml_getter(filepath, lipo_data)
df2_lipo = preprocessing(df_lipo)
oLig_abs = max_abs_lipo(df2_lipo)
dil_lipo = dil_factor(lipo_data)
lipo = concentration_0Lig(oLig_abs, lipo_data, dil_lipo, df2_lipo)
print(f'Liposome Concentration: {lipo}\n')

df_FTPAL = xml_getter(filepath, FTPAL_data)
df2_FTPAL = preprocessing(df_FTPAL)
FTPAL_abs = max_abs_lipo(df2_FTPAL)
dil_FTPAL = dil_factor(FTPAL_data)
FTPAL = concentration_TPAL(FTPAL_abs, FTPAL_data, dil_FTPAL, df2_FTPAL)
print(f'TPAL Concentration: {FTPAL}\n')

df_linker = xml_getter(filepath, linker_data)
df2_linker = preprocessing(df_linker)
linker_abs = max_abs_linker(df2_linker)
dil_linker = dil_factor(linker_data)
linker = concentration_linker(linker_abs, linker_data, dil_linker, df2_linker)
print(f'Linker Concentration for First Fluorimetry Dilution: {linker}\n')

repeat = int(input('How many data points do you have for your linker group? '))

# Create an empty list to store the fluorescence data points
y = []

highest = float(input('Enter highest fluorescence data point: '))
y.append(highest)
# Loop to collect the fluorescence data points
for i in range(repeat - 1):
    value = float(input(f'Enter fluorescence data point {i+2}: '))
    y.append(value)

# Print the collected data
print("Fluorescence data points:", y)

fluor_0Lig = float(input('Enter fluorescence value of 0-Ligand Liposome: '))
fluor_TPAL = float(input('Enter fluorescence value of TPAL: '))
#y = [fluor1, fluor2, fluor3, fluor4, fluor5, fluor6, fluor7, fluor8, fluor9, fluor10, fluor11, fluor12]
y.sort()
print(y)


linkers = linker_con(linker, repeat)
print(f'Linker concentrations: {linkers}')
#DELETE WHEN ACTUALLY RUNNING BECAUSE YOU
#INPUT THE VALUES ABOVE



#print(f'Positive difference: {positive}')
#print(f'Negative difference: {negative}')

x = linkers

lr = best_regression_analysis()


#interpolate values
interpolated_TPAL = (fluor_TPAL - lr.intercept_) / lr.coef_
interpolated_0Lig = (fluor_0Lig - lr.intercept_) / lr.coef_
print(f'Interpolated TPAL Concentration: {interpolated_TPAL}')
print(f'Interpolated 0-Ligand Liposome Concentration: {interpolated_0Lig}')

diameter = float(input('Enter diameter of 0-Ligand Liposome: '))





Conjugation_Ratio, Conjugation_Efficiency = conjugation_calc(interpolated_0Lig, lipo, interpolated_TPAL, FTPAL, diameter)
print(f'Conjugation Ratio: {Conjugation_Ratio}')
print(f'Conjugation Efficiency: {Conjugation_Efficiency}')

