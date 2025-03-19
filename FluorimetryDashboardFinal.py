import dash, numpy as np, matplotlib.pyplot as plt, re
from sklearn.linear_model import LinearRegression
from dash import dcc, html, dash_table
from dash.dependencies import Input, Output, State
import pandas as pd
import plotly.express as px
import xml.etree.ElementTree as ET
import base64
import io
import warnings
warnings.filterwarnings("ignore")


# Initialize Dash app
app = dash.Dash(__name__)

server = app.server

# Define app layout
app.layout = html.Div(
    style={
        'backgroundColor': '#3c5859',  # Dark background color
        'color': 'white',              # Light text color for contrast
        'fontFamily': 'Arial, sans-serif', 
        'height': '500vh', # Full viewport height
        'padding': '20px'  # Add some padding around the content
    },
    children = [
    html.H1("Fluorimetry Data Analysis", style={"textAlign": "center"}),

    # File Upload Component
    dcc.Upload(
        id="0ligandupload",
        children=html.Button("Upload 0-Ligand Liposome XML File", style={
                'fontSize': '20px', 'padding': '15px', 'borderRadius': '5px', 'backgroundColor': '#007BFF', 'color': 'white'}),
        multiple=False,
        style={"margin-bottom": "20px", 'fontFamily': 'Arial, sans-serif', 'fontSize': '18px', 'textAlign': 'center'}
    ),
    
      dcc.Upload(
        id="TPALupload",
        children=html.Button("Upload TPAL XML File", style={
                'fontSize': '20px', 'padding': '15px', 'borderRadius': '5px', 'backgroundColor': '#007BFF', 'color': 'white'}),
        multiple=False,  # Allow multiple files
        style={"margin-bottom": "20px", 'fontFamily': 'Arial, sans-serif', 'fontSize': '18px', 'textAlign': 'center'}
    ),
      
      dcc.Upload(
        id="linkerupload",
        children=html.Button("Upload Linker Group XML File", style={
                'fontSize': '20px', 'padding': '15px', 'borderRadius': '5px', 'backgroundColor': '#007BFF', 'color': 'white'}),
        multiple=False,  # Allow multiple files
        style={"margin-bottom": "20px", 'fontFamily': 'Arial, sans-serif', 'fontSize': '18px', 'textAlign': 'center'}
    ),
      
    html.Div([
        html.Label("Enter 'F' for F(ab)2 TPAL and 'M' for mAb TPAL:", style={'fontSize': '22px', 'color': 'white'}),
            dcc.Input(id="TPALtype", type="text", value='F', style={
                "margin-left": "10px","margin-right": "1px", "fontSize": "20px", 'width': '20px', 'padding': '5px', 'borderRadius': '5px','textAlign': 'center'})
        ], style={"padding": "5px",'width': '100%', 'display': 'flex', 'justifyContent': 'center'}),
    
      html.Div([
            html.Label("Enter highest fluorescence value: ", style={'fontSize': '22px', 'color': 'white'}),
            dcc.Input(id="fluor1", type="number", step="any", value=1000000, style={
            "margin-left": "10px", "margin-right": "10px", "fontSize": "20px", 'width': '100px', 'padding': '5px', 'borderRadius': '5px', 'textAlign': 'center'})
], style={"padding": "5px", 'width': '100%', 'display': 'flex', 'justifyContent': 'center', 'alignItems': 'center'}),

       
       html.Div([
    html.Label("Enter 2nd highest fluorescence value:", style={'fontSize': '22px', 'color': 'white'}),
    dcc.Input(id="fluor2", type="number", step="any", value=400000, style={
        "margin-left": "10px", "margin-right": "10px", "fontSize": "20px", 'width': '100px', 'padding': '5px', 
        'borderRadius': '5px', 'textAlign': 'center',
        '-webkit-appearance': 'none',  # Remove spinner in webkit browsers (Chrome, Safari)
        '-moz-appearance': 'textfield',  # Remove spinner in Firefox
        'appearance': 'none',  # Standard property to remove spinner
        'user-select': 'none',  # Disable text selection to prevent accidental changes
        'outline': 'none',  # Remove any outline
        'box-sizing': 'border-box'  # Ensure padding and border don't affect the input size
    }),
], style={"padding": "5px", 'width': '100%', 'display': 'flex', 'justifyContent': 'center', 'alignItems': 'center'}),

    
       html.Div([
        html.Label("Enter 3rd fluorescence value:", style={'fontSize': '22px', 'color': 'white'}),
        dcc.Input(id="fluor3", type="number", step="any", value=200000, style={
            "margin-left": "10px", "margin-right": "10px", "fontSize": "20px", 'width': '100px', 'padding': '5px', 'borderRadius': '5px', 'textAlign': 'center'})
], style={"padding": "5px", 'width': '100%', 'display': 'flex', 'justifyContent': 'center', 'alignItems': 'center'}),
    
       html.Div([
        html.Label("Enter 4th fluorescence value:", style={'fontSize': '22px', 'color': 'white'}),
        dcc.Input(id="fluor4", type="number", step="any", value=90000, style={
            "margin-left": "10px", "margin-right": "10px", "fontSize": "20px", 'width': '100px', 'padding': '5px', 'borderRadius': '5px', 'textAlign': 'center'})
], style={"padding": "5px", 'width': '100%', 'display': 'flex', 'justifyContent': 'center', 'alignItems': 'center'}),
    
    html.Div([
        html.Label("Enter 5th fluorescence value:", style={'fontSize': '22px', 'color': 'white'}),
        dcc.Input(id="fluor5", type="number", step="any", value=50000, style={
            "margin-left": "10px", "margin-right": "10px", "fontSize": "20px", 'width': '100px', 'padding': '5px', 'borderRadius': '5px', 'textAlign': 'center'})
], style={"padding": "5px", 'width': '100%', 'display': 'flex', 'justifyContent': 'center', 'alignItems': 'center'}),
    
    html.Div([
        html.Label("Enter 6th fluorescence value:", style={'fontSize': '22px', 'color': 'white'}),
        dcc.Input(id="fluor6", type="number", step="any", value=24000, style={
            "margin-left": "10px", "margin-right": "10px", "fontSize": "20px", 'width': '100px', 'padding': '5px', 'borderRadius': '5px', 'textAlign': 'center'})
], style={"padding": "5px", 'width': '100%', 'display': 'flex', 'justifyContent': 'center', 'alignItems': 'center'}),
    
    html.Div([
        html.Label("Enter 7th fluorescence value:", style={'fontSize': '22px', 'color': 'white'}),
        dcc.Input(id="fluor7", type="number", step="any", value=12500, style={
            "margin-left": "10px", "margin-right": "10px", "fontSize": "20px", 'width': '100px', 'padding': '5px', 'borderRadius': '5px', 'textAlign': 'center'})
], style={"padding": "5px", 'width': '100%', 'display': 'flex', 'justifyContent': 'center', 'alignItems': 'center'}),
    
    html.Div([
        html.Label("Enter 8th fluorescence value:", style={'fontSize': '22px', 'color': 'white'}),
        dcc.Input(id="fluor8", type="number", step="any", value=6000, style={
            "margin-left": "10px", "margin-right": "10px", "fontSize": "20px", 'width': '100px', 'padding': '5px', 'borderRadius': '5px', 'textAlign': 'center'})
], style={"padding": "5px", 'width': '100%', 'display': 'flex', 'justifyContent': 'center', 'alignItems': 'center'}),
    
    html.Div([
        html.Label("Enter 9th fluorescence value:", style={'fontSize': '22px', 'color': 'white'}),
        dcc.Input(id="fluor9", type="number", step="any", value=2700, style={
            "margin-left": "10px", "margin-right": "10px", "fontSize": "20px", 'width': '100px', 'padding': '5px', 'borderRadius': '5px', 'textAlign': 'center'})
], style={"padding": "5px", 'width': '100%', 'display': 'flex', 'justifyContent': 'center', 'alignItems': 'center'}),
    
    html.Div([
        html.Label("Enter 10th fluorescence value:",style={'fontSize': '22px', 'color': 'white'}),
        dcc.Input(id="fluor10", type="number", step="any", value=1400, style={
            "margin-left": "10px", "margin-right": "10px", "fontSize": "20px", 'width': '100px', 'padding': '5px', 'borderRadius': '5px', 'textAlign': 'center'})
], style={"padding": "5px", 'width': '100%', 'display': 'flex', 'justifyContent': 'center', 'alignItems': 'center'}),
    
    html.Div([
        html.Label("Enter 11th fluorescence value:",style={'fontSize': '22px', 'color': 'white'}),
        dcc.Input(id="fluor11", type="number", step="any", value=700, style={
            "margin-left": "10px", "margin-right": "10px", "fontSize": "20px", 'width': '100px', 'padding': '5px', 'borderRadius': '5px', 'textAlign': 'center'})
], style={"padding": "5px", 'width': '100%', 'display': 'flex', 'justifyContent': 'center', 'alignItems': 'center'}),
    
    html.Div([
        html.Label("Enter lowest fluorescence value:",style={'fontSize': '22px', 'color': 'white'}),
        dcc.Input(id="fluor12", type="number", step="any", value=700, style={
            "margin-left": "10px", "margin-right": "10px", "fontSize": "20px", 'width': '100px', 'padding': '5px', 'borderRadius': '5px', 'textAlign': 'center'})
], style={"padding": "5px", 'width': '100%', 'display': 'flex', 'justifyContent': 'center', 'alignItems': 'center'}),
    
    
    html.Div([
        html.Label("Enter TPAL fluorescence: ",style={'fontSize': '22px', 'color': 'white'}),
        dcc.Input(id="fluor_TPAL", type="number", step="any", value=1500, style={
            "margin-left": "10px", "margin-right": "10px", "fontSize": "20px", 'width': '100px', 'padding': '5px', 'borderRadius': '5px', 'textAlign': 'center'})
], style={"padding": "5px", 'width': '100%', 'display': 'flex', 'justifyContent': 'center', 'alignItems': 'center'}),
    
    html.Div([
        html.Label("Enter 0-Ligand Liposome Fluorescence: ",style={'fontSize': '22px', 'color': 'white'}),
        dcc.Input(id="fluor_0Lig", type="number", step="any", value=600, style={
            "margin-left": "10px", "margin-right": "10px", "fontSize": "20px", 'width': '100px', 'padding': '5px', 'borderRadius': '5px', 'textAlign': 'center'})
], style={"padding": "5px", 'width': '100%', 'display': 'flex', 'justifyContent': 'center', 'alignItems': 'center'}),
    
    html.Div([
        html.Label("Enter Liposome Diameter: ",style={'fontSize': '22px', 'color': 'white'}),
        dcc.Input(id="diameter", type="number", step="any", value=100, style={
            "margin-left": "10px", "margin-right": "10px", "fontSize": "20px", 'width': '100px', 'padding': '5px', 'borderRadius': '5px', 'textAlign': 'center'})
], style={"padding": "5px", 'width': '100%', 'display': 'flex', 'justifyContent': 'center', 'alignItems': 'center'}),

    

    # Graph output
    html.Div([html.P('Quick Results: ', style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'})]),
    
    html.Div([html.P('Files Chosen: ', style={'fontSize': '22px', 'color': 'white', 'marginTop': '60px'})]),
    html.Div(id="output_div", style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'}),
    html.Div(id="output_div1", style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'}),
    html.Div(id="output_div2", style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'}),
    
    html.Div([html.P('UV-Vis Concentrations: ', style={'fontSize': '22px', 'color': 'white', 'marginTop': '60px'})]),
    html.Div(id="concentration0lig", style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'}),
    html.Div(id="concentrationTPAL", style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'}),
    html.Div(id="concentrationlinker", style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'}),
    
    html.Div(id="graphs-container3"),
    
    html.Div([html.P('Interpolated Values: ', style={'fontSize': '22px', 'color': 'white', 'marginTop': '60px'})]),
    html.Div(id="inter0lig", style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'}),
    html.Div(id="interTPAL", style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'}),
    
    html.Div([html.P('Conjugation Info: ', style={'fontSize': '22px', 'color': 'white', 'marginTop': '60px'})]),
    html.Div(id="output_div3", style={'fontSize': '35px', 'color': 'white', 'marginTop': '30px'}),
    
    html.Div([html.P('Details: ', style={'fontSize': '30px', 'color': 'white', 'marginTop': '60px'})]),
    html.Div([html.P('This assumes that all fluorimetry dilutions are 1 in 50 for liposome groups. It also assumes that after UV-Vis of the linker group, a 1 in 4 dilution is used for fluorimetry. Artificial intelligence through the Scikit-Learn Library is used to determine line with the highest r^2 value from different combinations of data points for best results.', 
                     style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'})]),
    
    html.Div([html.P('Absorbance details: ', style={'fontSize': '22px', 'color': 'white', 'marginTop': '60px'})]),
    html.Div(id="abs0lig", style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'}),
    html.Div(id="absTPAL", style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'}),
    html.Div(id="abslinker", style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'}),
    html.Div(id="lig"),
    html.Div(id="TPAL"),
    html.Div(id="linker"),
    
    html.Div([html.P('Regression Details: ', style={'fontSize': '22px', 'color': 'white', 'marginTop': '60px'})]),
    html.Div(id="xlist", style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'}),
    html.Div(id="ylist", style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'}),
    html.Div(id="shortxlist", style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'}),
    html.Div(id="shortylist", style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'}),
    html.Div(id="r2list", style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'}),
    html.Div(id="shortr2list", style={'fontSize': '22px', 'color': 'white', 'marginTop': '20px'}),
    
    html.Div([html.P('Good luck!', style={'fontSize': '40px', 'color': 'white', 'marginTop': '60px'})]),
    
    
    
])

# Function to parse XML and extract data (Modified to match your `xml_getter` function)

def xml_getter(xml_content):
    """Extracts data from an XML file and converts it into a Pandas DataFrame."""
    try:
        # Parse XML content
        tree = ET.ElementTree(ET.fromstring(xml_content))
        root = tree.getroot()

        # Define namespace
        namespace = {"ss": "urn:schemas-microsoft-com:office:spreadsheet"}

        # Locate the table
        table = root.find(".//ss:Table", namespace)
        if table is None:
            return pd.DataFrame()  # Return an empty DataFrame if no table is found

        # Extract rows
        data = []
        for row in table.findall("ss:Row", namespace):
            cells = [cell.find("ss:Data", namespace).text if cell.find("ss:Data", namespace) is not None else "" 
                     for cell in row.findall("ss:Cell", namespace)]
            data.append(cells)

        # Convert to Pandas DataFrame
        df = pd.DataFrame(data)

        # Ensure column names are strings
        df.columns = df.columns.astype(str)

        return df
    except Exception as e:
        print(f"Error parsing XML: {e}")
        return pd.DataFrame()


    abs_array = []
    for i in range(1, 20):
        abs_array.append(df2.iloc[i + 70, 1])

    oLig_abs = np.max(abs_array)
    oLig_abs = oLig_abs.astype(float)
    
    return(oLig_abs)

def preprocessing(df):
    # Set the first row as column headers
    #df.columns = df.iloc[0]  # Use the first row as header
    df = df[1:].reset_index(drop=True)  # Remove header row from data
    #df = df[['Wavelength (nm)', 'Absorbance']]  # Keep only the columns we need
    #df = df[['Wavelength (nm)', 'Absorbance']].astype(float)
    
    df2 = pd.DataFrame(columns = df.columns)
    
    df2 = df.iloc[::-1].reset_index(drop = True)
    
    return df2

def max_abs_lipo(df2):
    abs_array = []
    for i in range(1, 20):
        abs_array.append(df2.iloc[i + 480, 1])

    oLig_abs = np.max(abs_array)
    oLig_abs = oLig_abs.astype(float)
    print('Max Absorbance Lipo:', oLig_abs)
    return(oLig_abs)

def max_abs_linker(df2):
    abs_array = []
    for i in range(1, 20):
        abs_array.append(df2.iloc[i + 70, 1])

    oLig_abs = np.max(abs_array)
    oLig_abs = oLig_abs.astype(float)
    print('Max Absorbance Linker:', oLig_abs)
    return(oLig_abs)

def dil_factor(file):
    # Define the patterns to search for '1 in [number]'
    patterns = [r'1 in (\d+)', r'1in(\d+)']

    for pattern in patterns:
        match = re.search(pattern, file)
        if match:
            dil_factor = float(match.group(1))  # Extract the number and convert to float
            print(f'Dilution Factor: {dil_factor}')
            return dil_factor  # Return immediately if a match is found

    
    return None  # Return None if no match is found

def concentration_TPAL(abs, file, dil_factor, df2):
    filtered_df = df2[df2.iloc[:, 1] > 0]  # Filter rows based on the second column
    background = np.min(filtered_df.iloc[:, 1])  # Select the second column and find the min value

    
    BPD_exc = (34895)
    con = (abs - background) / BPD_exc  * 10**6 * dil_factor / 50 * 1000
    print(f'Background: {background}')
    print(f'Concentration of TPAL: {con}')
    return con
    
def concentration_linker(abs, file, dil_factor, df2, fab_vs_mab):
    filtered_df = df2[df2.iloc[:, 1] > 0]  # Filter rows based on the second column
    background = np.min(filtered_df.iloc[:, 1])  # Select the second column and find the min value

    
    
    if fab_vs_mab == 'F':
        protein_exc = (154000)
        
        con = (abs - background) / protein_exc  * 10**6 * dil_factor / dil_factor * 1000 /4
        
        
    elif fab_vs_mab == 'M':
        protein_exc = (210000)
        con = (abs - background) / protein_exc  * 10**6 * dil_factor / dil_factor * 1000 / 4
        
    else:
        print('Invalid Input')
    print(f'Background: {background}')
    print(f'Concentration of linker: {con}')
    return con

def concentration_0Lig(abs, file, dil_factor, df2):
    filtered_df = df2[df2.iloc[:, 1] > 0]  # Filter rows based on the second column
    background = np.min(filtered_df.iloc[:, 1])  # Select the second column and find the min value

    
    
    BPD_exc = (34895)
    con = (abs - background) / BPD_exc  * 10**6 * dil_factor / 50 * 1000
    print(f'Background: {background}')
    print(f'Concentration of 0Ligand Liposome: {con}')
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
    
    return lr, r2

def stuff_for_best_regression_analysis(x11, y11, x, y, fluor_TPAL, linkers):
    ynew = get_y_features(y, fluor_TPAL,x11, y11)
    dif = ynew
    xnew = get_x_features(x, dif, y, linkers)
    xnew = np.array(xnew).reshape(-1,1)
    lr, r2= Regression_Analysis(xnew, ynew)
    
    yhat22 = lr.predict(xnew)
    
    return xnew, ynew, lr, r2, yhat22

def best_regression_analysis(graphs3, x, y, fluor_TPAL, linkers):
    r2_list = []
    #For 2 above, 2 below
   
    
    
    xnew22, ynew22, lr22, r2_22, yhat22 = stuff_for_best_regression_analysis(2,2, x, y, fluor_TPAL, linkers)
    r2_list.append(r2_22)
    xnew33, ynew33, lr33, r2_33, yhat33 = stuff_for_best_regression_analysis(3,3, x, y, fluor_TPAL, linkers)
    r2_list.append(r2_33)
    
    xnew44, ynew44, lr44, r2_44, yhat44 = stuff_for_best_regression_analysis(4,4, x, y, fluor_TPAL, linkers)
    r2_list.append(r2_44)
    xnew23, ynew23, lr23, r2_23, yhat23 = stuff_for_best_regression_analysis(2,3, x, y, fluor_TPAL, linkers)
    r2_list.append(r2_23)
    xnew32, ynew32, lr32, r2_32, yhat32 = stuff_for_best_regression_analysis(3,2, x, y, fluor_TPAL, linkers)
    r2_list.append(r2_32)
    xnew43, ynew43, lr43, r2_43, yhat43 = stuff_for_best_regression_analysis(4,3, x, y, fluor_TPAL, linkers)
    r2_list.append(r2_43)
    xnew34, ynew34, lr34, r2_34, yhat34 = stuff_for_best_regression_analysis(3,4, x, y, fluor_TPAL, linkers)
    r2_list.append(r2_34)
    xnew13, ynew13, lr13, r2_13, yhat13 = stuff_for_best_regression_analysis(1,3, x, y, fluor_TPAL, linkers)
    r2_list.append(r2_13)
    xnew31, ynew31, lr31, r2_31, yhat31 = stuff_for_best_regression_analysis(3,1, x, y, fluor_TPAL, linkers)
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
    
    
    
    '''xbest = np.array(xbest)
    ybest = np.array(ybest)
    
    xbest = xbest.flatten()
    ybest = ybest.flatten()
    
    print(f'Xbest: {xbest}, datatype: {type(xbest)}')
    print(f'Ybest: {ybest}, datatype: {type(ybest)}')
    
    m, b = np.polyfit(xbest, ybest, 1)
    print(f'm: {m}, b: {b}')
    print('Best Fit Line: y = ', m, 'x + ', b)
    
    y_fit = m * xbest + b
    
    df3 = pd.DataFrame({"X": xbest, "Y": ybest})
    print(df3.head())
    df_fit = pd.DataFrame({"X": xbest, "Y": yhatbest})  # Best-fit line points
    print(df_fit.head())'''
    
    # Create scatter plot using Plotly Express
    
    
    
    
    return bestlr, xbest, ybest, yhatbest


# Callback to process uploaded XML files
@app.callback(
    [Output("output_div", "children"),
     Output("output_div1", "children"),
     Output("output_div2", "children"),
     Output("concentration0lig", "children"),
     Output("concentrationTPAL", "children"),
     Output("concentrationlinker", "children"),
     Output("graphs-container3", 'children'),
     Output("inter0lig", "children"),
     Output("interTPAL", "children"),
     Output("output_div3", "children"),
     Output("abs0lig", "children"),
     Output("absTPAL", "children"),
     Output("abslinker", "children"),
     Output("lig", 'children'),
     Output("TPAL", 'children'),
     Output("linker", 'children'),
     Output("xlist", "children"),
     Output("ylist", "children"),
     Output("shortxlist", "children"),
     Output("shortylist", "children"),
     Output("r2list", "children"),
     Output("shortr2list", "children"),],
    [Input("0ligandupload", "contents"),
     Input("TPALupload", "contents"),
     Input("linkerupload", "contents"),
     Input("TPALtype", "value"),
     Input("fluor1", "value"),
     Input("fluor2", "value"),
     Input("fluor3", "value"),
     Input("fluor4", "value"),
     Input("fluor5", "value"),
     Input("fluor6", "value"),
     Input("fluor7", "value"),
     Input("fluor8", "value"),
     Input("fluor9", "value"),
     Input("fluor10", "value"),
     Input("fluor11", "value"),
     Input("fluor12", "value"),
     Input("fluor_TPAL", "value"),
     Input("fluor_0Lig", "value"),
     Input("diameter", "value")],
    [State("0ligandupload", "contents"),
     State("TPALupload", "contents"),
     State("linkerupload", "contents")]
)

def update_output(  # Accept all inputs and states
    contents_0ligand, contents_TPAL, contents_linker,
    TPALtype, fluor1, fluor2, fluor3, fluor4, fluor5, fluor6, fluor7, fluor8, fluor9, fluor10, fluor11, fluor12,
    fluor_TPAL, fluor_0Lig, diameter,
    filename_0ligand, filename_TPAL, filename_linker
):
    
    y = [fluor1, fluor2, fluor3, fluor4, fluor5, fluor6, fluor7, fluor8, fluor9, fluor10, fluor11, fluor12]
    y = [float(f) for f in y]
    y.sort()
    print(f'Fluor1: {fluor1}, datatype: {type(fluor1)}')
    print(f'Y: {y}')
    
    fluor_TPAL = float(fluor_TPAL)
    fluor_0Lig = float(fluor_0Lig)
    diameter = float(diameter)
    fab_vs_mab = TPALtype
    print(fab_vs_mab)
    
    
    
    # Initialize an empty list to store individual graphs
    
    graphs3 = []
    liggraph = []
    TPALgraph = []
    linkergraph = []
    
    
    df_0lig1 = pd.DataFrame()
    df_TPAL1 = pd.DataFrame()
    df_linker1 = pd.DataFrame()
    
    # Process the uploaded 0-ligand XML file
    if contents_0ligand:
        content_type, content_string = contents_0ligand.split(",")
        decoded = base64.b64decode(content_string)
        xml_content = decoded.decode("utf-8")

        # Parse XML and create DataFrame
        df_0lig = xml_getter(xml_content)
        
        if not df_0lig.empty:
            
            filename_0ligand = df_0lig.iloc[1, 2]
            print(filename_0ligand)
            df_0lig = df_0lig.iloc[:, :2]

            # Generate graph for this file
            if df_0lig.shape[1] >= 2:
                x_column, y_column = df_0lig.columns[0], df_0lig.columns[1]
                
                df_0lig[x_column] = pd.to_numeric(df_0lig[x_column], errors='coerce')
                df_0lig[y_column] = pd.to_numeric(df_0lig[y_column], errors='coerce')  
                
                df_combined = df_0lig[[x_column, y_column]]
                
                df_0lig1 = df_combined.copy()
                df_0lig1.columns = ['Wavelength (nm)', 'Absorbance']
                print(df_0lig1)
                fig = px.scatter(df_0lig1, x='Wavelength (nm)', y='Absorbance', 
                             title=f"Scatter Plot for {filename_0ligand}",
                             labels={"Wavelength (nm)": "Wavelength (nm)", "Absorbance": "Absorbance"})
                liggraph.append(dcc.Graph(figure=fig))
                
                
                
    
    if contents_TPAL:
        content_type, content_string = contents_TPAL.split(",")
        decoded = base64.b64decode(content_string)
        xml_content = decoded.decode("utf-8")

        # Parse XML and create DataFrame
        df_TPAL = xml_getter(xml_content)
        
        if not df_TPAL.empty:
            filename_TPAL = df_TPAL.iloc[1, 2]
            df_TPAL = df_TPAL.iloc[:, :2]

            # Generate graph for this file
            if df_TPAL.shape[1] >= 2:
                x_column, y_column = df_TPAL.columns[0], df_TPAL.columns[1]
                
                df_TPAL[x_column] = pd.to_numeric(df_TPAL[x_column], errors='coerce')
                df_TPAL[y_column] = pd.to_numeric(df_TPAL[y_column], errors='coerce')  
                
                df_combined1 = df_TPAL[[x_column, y_column]]
                
                df_TPAL1 = df_combined1.copy()
                df_TPAL1.columns = ['Wavelength (nm)', 'Absorbance']
                fig = px.scatter(df_TPAL1, x='Wavelength (nm)', y='Absorbance', 
                             title=f"Scatter Plot for {filename_TPAL}",
                             labels={"Wavelength (nm)": "Wavelength (nm)", "Absorbance": "Absorbance"})
                TPALgraph.append(dcc.Graph(figure=fig))
                
    # Return the table data and multiple graphs
    if contents_linker:
        content_type, content_string = contents_linker.split(",")
        decoded = base64.b64decode(content_string)
        xml_content = decoded.decode("utf-8")

        # Parse XML and create DataFrame
        df_linker = xml_getter(xml_content)
        
        if not df_linker.empty:
            filename_linker = df_linker.iloc[1, 2]
            df_linker = df_linker.iloc[:, :2]

            # Generate graph for this file
            if df_linker.shape[1] >= 2:
                x_column, y_column = df_linker.columns[0], df_linker.columns[1]
                
                df_linker[x_column] = pd.to_numeric(df_linker[x_column], errors='coerce')
                df_linker[y_column] = pd.to_numeric(df_linker[y_column], errors='coerce')  
                
                df_combined2 = df_linker[[x_column, y_column]]
                
                df_linker1 = df_combined2.copy()
                df_linker1.columns = ['Wavelength (nm)', 'Absorbance']
                fig = px.scatter(df_linker1, x='Wavelength (nm)', y='Absorbance', 
                             title=f"Scatter Plot for {filename_linker}",
                             labels={"Wavelength (nm)": "Wavelength (nm)", "Absorbance": "Absorbance"})
                linkergraph.append(dcc.Graph(figure=fig))
                
                
    
    
    lipo_data = filename_0ligand
    FTPAL_data = filename_TPAL
    linker_data = filename_linker
    
    df2_lipo = preprocessing(df_0lig1)
    oLig_abs = max_abs_lipo(df2_lipo)
    dil_lipo = dil_factor(lipo_data)
    lipo = concentration_0Lig(oLig_abs, lipo_data, dil_lipo, df2_lipo)
    

    
    df2_FTPAL = preprocessing(df_TPAL1)
    FTPAL_abs = max_abs_lipo(df2_FTPAL)
    dil_FTPAL = dil_factor(FTPAL_data)
    FTPAL = concentration_TPAL(FTPAL_abs, FTPAL_data, dil_FTPAL, df2_FTPAL)
    
    
    df2_linker = preprocessing(df_linker1)
    linker_abs = max_abs_linker(df2_linker)
    dil_linker = dil_factor(linker_data)
    linker = concentration_linker(linker_abs, linker_data, dil_linker, df2_linker, fab_vs_mab)
    
    
    
    
    
    
    
    repeat = 12
    linkers = linker_con(linker, repeat)
    x = linkers
    print(f'X: {x}')
    
    
    lr, xbest, ybest, yhatbest = best_regression_analysis(graphs3, x, y, fluor_TPAL, linkers)
    
    xbest = xbest.tolist()
    xbest = [item[0] for item in xbest]
    yhatbest = yhatbest.tolist()
    
    print(f'xbest: {xbest}, datatype: {type(xbest)}')
    print(f'ybest: {ybest}, datatype: {type(ybest)}')
    print(f'yhatbest: {yhatbest}, datatype: {type(yhatbest)}')
    print(f'coef: {lr.coef_}, intercept: {lr.intercept_}')
    print(f'xbest length: {len(xbest)}, ybest length: {len(ybest)}, yhatbest length: {len(yhatbest)}')
    
    fig = px.scatter(x=xbest, y=ybest, title="Scatter Plot of Fluorimetry Data with Best Fit Line from Regression Analysis", labels={"X": "X Values", "Y": "Y Values"})
    fig.update_layout(xaxis_title = 'Linker Concentration (nM)', yaxis_title = 'Fluorescence (RFUs)')
    
    print('scatter working')
    
    # Add best-fit line
    fig.add_scatter(x=xbest, y=yhatbest, mode="lines", name="Best Fit Line")
    print('fig working')
    
    graphs3.append(dcc.Graph(figure=fig))
    
    
    interpolated_TPAL = (fluor_TPAL - lr.intercept_) / lr.coef_
    interpolated_0Lig = (fluor_0Lig - lr.intercept_) / lr.coef_
    
    print(f'Interpolated TPAL: {interpolated_TPAL}')
    print(f'Interpolated 0Lig: {interpolated_0Lig}')
    

    Conjugation_Ratio, Conjugation_Efficiency = conjugation_calc(interpolated_0Lig, lipo, interpolated_TPAL, FTPAL, diameter)
    
    print(f'Conjugation Ratio: {Conjugation_Ratio}')
    
    print(type(Conjugation_Efficiency))
    print(str(Conjugation_Efficiency))
    
    x = np.array(x).reshape(-1, 1)
    y = np.array(y).reshape(-1, 1)
    xbest = np.array(xbest).reshape(-1, 1)
    ybest = np.array(ybest).reshape(-1, 1)
    
    
    return (f'0Ligand Liposome UV-Vis File Name: {lipo_data}', 
            f'TPAL UV-Vis File Name: {FTPAL_data}', 
            
            f'Linker UV-Vis File Name: {linker_data}', 
            f'0Ligand Concentration: {lipo.round(3)} nM',
            f'TPAL Concentration: {FTPAL.round(3)} nM',
            f'Linker Concentration: {linker.round(3)} nM',
            graphs3, 
            f'Interpolated 0Ligand: {interpolated_0Lig[0].round(3)} nM',
            f'Interpolated TPAL: {interpolated_TPAL[0].round(3)} nM',
            f'Conjugation Ratio: {Conjugation_Ratio[0].round(3)}, Conjugation Efficiency: {Conjugation_Efficiency[0].round(3)}',
            f'Max Absorbance 0Ligand Liposome: {oLig_abs.round(3)}, 0Ligand Liposome UV-Vis Dilution Factor: {dil_lipo}',
            f'Max Absorbance TPAL: {FTPAL_abs.round(3)}, TPAL UV-Vis Dilution Factor: {dil_FTPAL}',
            f'Max Absorbance Linker: {linker_abs.round(3)}, Linker UV-Vis Dilution Factor: {dil_linker}',
            liggraph, TPALgraph, linkergraph,
            f'X Values: {x}', f'Y Values: {y}', f'Best X Values for Regression: {xbest}', f'Best Y Values for Regression: {ybest}', f'R^2 Values List from Various Combinations: {lr.score(xbest, ybest)}', f'Best R^2 Value: {lr.score(xbest, ybest)}')
    
    
    '''lipo_data = filename_0ligand
    FTPAL_data = filename_TPAL
    linker_data = filename_linker
    
    df2_lipo = preprocessing(df_0lig1)
    oLig_abs = max_abs_lipo(df2_lipo)
    dil_lipo = dil_factor(lipo_data)
    lipo = concentration_0Lig(oLig_abs, lipo_data, dil_lipo, df2_lipo)
    

    
    df2_FTPAL = preprocessing(df_TPAL1)
    FTPAL_abs = max_abs_lipo(df2_FTPAL)
    dil_FTPAL = dil_factor(FTPAL_data)
    FTPAL = concentration_TPAL(FTPAL_abs, FTPAL_data, dil_FTPAL, df2_FTPAL)
    
    
    df2_linker = preprocessing(df_linker1)
    linker_abs = max_abs_linker(df2_linker)
    dil_linker = dil_factor(linker_data)
    linker = concentration_linker(linker_abs, linker_data, dil_linker, df2_linker, fab_vs_mab)'''
    
# Run app
if __name__ == "__main__":
    app.run(debug=True)
