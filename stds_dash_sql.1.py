#!/usr/bin/env python3

#works OK updates seem to be linked to the DB wen pressing refresh button. 


#Dependencies
import dash
import dash_auth
from dash.dependencies import Input, Output, State
from dash import dash_table
from dash import dcc
from dash import html
import pandas as pd
import plotly.graph_objs as go
import numpy as np
import plotly.express as px
import statistics as st
import mysql.connector
import time
from datetime import datetime
import csv




#read the data from the mysql database with the hd200 samples
def get_sql():
    #Use MySQL Connector for establishing link to DB in function, so that we may update live to the user
    mydb = mysql.connector.connect(host="localhost", database = 'HD200',user=" ", passwd=" ")
    query = "SELECT CallData.pass_filter,  CallData.afreq,  CallData.coverage,  CallData.norm_count,  CallData.sample, VarData.name, RunInfo.IonWF_version, RunInfo.name, RunInfo.filedate, Transcripts.name , HGVS.transcript, HGVS.HGVSc, HGVS.HGVSp, Genes.name FROM VarData LEFT JOIN HGVS ON HGVS.id = VarData.hgvs LEFT JOIN CallData ON VarData.id = CallData.variant LEFT JOIN RunInfo ON CallData.sample = RunInfo.id LEFT JOIN Transcripts ON Transcripts.id = HGVS.transcript LEFT JOIN Genes ON VarData.gene = Genes.id;"
    df = pd.read_sql(query,mydb)
    mydb.close()

    #read-in data and change duplicated column headers
    df.columns = ['pass_filter','afreq','coverage','norm_count','sample','variant','IonWF_version','samplename','filedate','trname','transcript','HGVSc', 'HGVSp','gene']

    #get only HD200 and seracare samples
    df = df[df['samplename'].str.contains("HD200_SSEQ_CONTROL", na=False)]

    #Get the variants of interest
    df = df[df['variant'].isin(["chr1_115256530_G_T_snp_1",\
    "chr3_41266101_C_A_snp_1",\
    "chr3_41266133_CCTT_C_del_3",\
    "chr3_178936091_G_A_snp_1",\
    "chr3_178952085_A_G_snp_1",\
    "chr4_55599321_A_T_snp_1",\
    "chr7_55241707_G_A_snp_1",\
    "chr7_55259515_T_G_snp_1",\
    "chr7_140453136_A_T_snp_1",\
    "chr12_25398281_C_T_snp_1",\
    "chr12_25398284_C_T_snp_1",\
    "chr1_154142876_C_C[CHR1:156844363[_Fusion_None",\
    "chr1_156100564_G_G[CHR1:156844697[_Fusion_None",\
    "chr1_205649522_C_C]CHR7:140494267]_Fusion_None",\
    "chr2_42522656_G_G]CHR2:29446394]_Fusion_None",\
    "chr2_113992971_C_C[CHR3:12421203[_Fusion_None",\
    "chr3_100451516_G_G[CHR1:156844362[_Fusion_None",\
    "chr4_1808623_G_G[CHR4:1739412[_Fusion_None",\
    "chr4_1808661_C_C]CHR7:97991744]_Fusion_None",\
    "chr4_1808661_C_C[CHR4:1741428[_Fusion_None",\
    "chr4_25665952_G_G]CHR6:117645578]_Fusion_None",\
    "chr5_149784243_C_C]CHR6:117645578]_Fusion_None",\
    "chr10_51582939_G_G[CHR10:43612031[_Fusion_None",\
    "chr10_61665880_G_G[CHR10:43612032[_Fusion_None",\
    "chr12_12022903_G_G]CHR15:88483984]_Fusion_None",\
    "chr21_42880008_C_C]CHR21:39956869]_Fusion_None",\
    "chr7_116411708_G_G[CHR7:116414934[_RNAExonVariant_None",\
    "chr7_55087058_G_G[CHR7:55223522[_RNAExonVariant_None",\
    "chr10_32306071_C_C[CHR10:43609928[_Fusion_None"])]
    #Change fusion gene names to something more interpretable
    df.replace(["chr1_154142876_C_C[CHR1:156844363[_Fusion_None",\
    "chr1_156100564_G_G[CHR1:156844697[_Fusion_None",\
    "chr1_205649522_C_C]CHR7:140494267]_Fusion_None",\
    "chr2_42522656_G_G]CHR2:29446394]_Fusion_None",\
    "chr2_113992971_C_C[CHR3:12421203[_Fusion_None",\
    "chr3_100451516_G_G[CHR1:156844362[_Fusion_None",\
    "chr4_1808623_G_G[CHR4:1739412[_Fusion_None",\
    "chr4_1808661_C_C]CHR7:97991744]_Fusion_None",\
    "chr4_1808661_C_C[CHR4:1741428[_Fusion_None",\
    "chr4_25665952_G_G]CHR6:117645578]_Fusion_None",\
    "chr5_149784243_C_C]CHR6:117645578]_Fusion_None",\
    "chr10_51582939_G_G[CHR10:43612031[_Fusion_None",\
    "chr10_61665880_G_G[CHR10:43612032[_Fusion_None",\
    "chr12_12022903_G_G]CHR15:88483984]_Fusion_None",\
    "chr21_42880008_C_C]CHR21:39956869]_Fusion_None",\
    "chr7_116411708_G_G[CHR7:116414934[_RNAExonVariant_None",\
    "chr7_55087058_G_G[CHR7:55223522[_RNAExonVariant_None",\
    "chr10_32306071_C_C[CHR10:43609928[_Fusion_None"],["TPM3(7) - NTRK1(10)",\
    "LMNA(2) - NTRK1(11)",\
    "SLC45A3(1) - BRAF(8)",\
    "EML4(13) - ALK(20)",\
    "PAX8(9) - PPARG(2)",\
    "TFG(5) - NTRK1(10)",\
    "FGFR3(17) - TACC3(10)",\
    "FGFR3(17) - BAIAP2L1(2)",\
    "FGFR3(17) - TACC3(11)",\
    "SLC34A2(4) - ROS1(34)",\
    "CD74(6) - ROS1(34)",\
    "NCOA4(7) - RET(12)",\
    "CCDC6(1) - RET(12)",\
    "ETV6(5) - NTRK3(15)",\
    "TMPRSS2(1) - ERG(2)",\
    "MET(13) - MET(15)",\
    "EGFR(1) - EGFR(8)",\
    "KIF5B(24) - RET(11)"], inplace = True)

    #Get only workflow 5.16
    df = df.loc[df['IonWF_version'] == 1]
    #Sort by date
    df = df.sort_values(["filedate","variant"])
    df = df.reset_index(drop=True)

    #Drop uninformative columns, add stdev for numerical values and group the entire table by variant. means and modes are used, except the samples names, which are the sum of unique names. 
    tabledf = df.drop(labels=['IonWF_version','pass_filter','transcript'], axis =1)
    tabledf['afreq'] = tabledf['afreq'].fillna(0)
    tabledf['afreq'] = tabledf['afreq'].astype(float)
    tabledf['norm_count'] = tabledf['norm_count'].fillna(0)
    tabledf['afreq_normcount'] = tabledf['afreq'] + tabledf['norm_count']
    tabledf['sd'] = tabledf.groupby('variant').afreq_normcount.transform('std')
    tabledf['upper_bound'] = tabledf['afreq_normcount'] + 2*(tabledf['sd'])
    tabledf['lower_bound'] = tabledf['afreq_normcount'] - 2*(tabledf['sd'])
    return(tabledf)

#initial call for the database query

def make_table():
    tabledf = get_sql()
    tabledf1 = tabledf.groupby('variant', as_index=False).agg({'samplename': 'nunique', 'coverage': ['mean'], 'afreq_normcount': ['mean'], 'trname': pd.Series.mode, 'HGVSc': pd.Series.mode, 'HGVSp':pd.Series.mode, 'filedate': ['max'], 'gene': pd.Series.mode, 'sd': ['mean'], 'upper_bound': ['mean'], 'lower_bound':['mean']})
    tabledf1.columns = tabledf1.columns.droplevel(1)
    neworder = ['variant','gene','afreq_normcount','sd', 'upper_bound', 'lower_bound','coverage','trname','HGVSc','HGVSp','samplename','filedate']
    tabledf = tabledf[neworder]
    tabledf1 = tabledf1[neworder]
    tabledf = tabledf.applymap(lambda x: round(x, 5) if isinstance(x, (int, float)) else x)
    tabledf1 = tabledf1.applymap(lambda x: round(x, 5) if isinstance(x, (int, float)) else x)
    # Needs to be jsonized to store in an empty div to gain speed, otherwise you are stuck building the dataframes 0 and 1 in each callback, which is expensive.
    return tabledf.to_json(double_precision = 5), tabledf1.to_json(double_precision = 5)

#create the app layout
app = dash.Dash(__name__)

VALID_USERNAME_PASSWORD_PAIRS = {
    '###': '###',
    '###': '###',
    '###': '###',
    '###': '###'
}

auth = dash_auth.BasicAuth(
    app,
    VALID_USERNAME_PASSWORD_PAIRS
)
#Avoid global variables, instead store in empty divs

#defining the layout
def serve_layout():
    tab0, tab1 = make_table()
    t0 = pd.read_json(tab0)
    t1 = pd.read_json(tab1)
    t0 = t0.applymap(lambda x: round(x, 5) if isinstance(x, (int, float)) else x)
    t1 = t1.applymap(lambda x: round(x, 5) if isinstance(x, (int, float)) else x)
    return html.Div(children=[
    #First is a title and a refresh button
    #Then a datatable with selectable rows for later graphs with callback
    html.Div([html.Button('Refresh', id='apply-button', n_clicks=0),
              html.Div(id='output-container-button', children='Click the button to update.')]),
    html.H1(children="HD200 / SeraCare controls"),
    dcc.Dropdown(id = 'drpdown', options = [{'label': i, 'value': i} for i in t0['samplename'].unique()[-50:]], value = t0['samplename'].tolist()[-1]),
    # add the conditionnal styling to allele frequencies that are out-of-bounds (2SD)
    dash_table.DataTable(id = 'table-2',
        columns=[{'name':i, 'id':i, 'deletable': False} for i in t0.columns if i not in ['id', 'afreqsd','normcountsd']],
        page_size=10,
        style_data_conditional=[
        {
            'if': {
                'filter_query':'{afreq_normcount} < {lower_bound}',
                'column_id':'afreq_normcount'
            },
            'backgroundColor': '#85144b',
            'color':'white'
        }, {
            'if': {
                'filter_query': '{afreq_normcount} > {upper_bound}',
                'column_id':'afreq_normcount'
            },
            'backgroundColor': '#85144b',
            'color':'white'
        }]),
    html.Br(),
    dcc.Graph(id = 'subplot-div', style={'width': '250vh'}),
    html.Div(dcc.Input(id='input-box', type='text')),
    html.Button('Submit', id='button'),
    html.Div(id='output-container-button2',children='No comments to display'),
    html.Br(),
    html.Button('Delete', id='button2'),
    html.Br(),
    dash_table.DataTable(
    id='table',
    columns=[
        {'name': i, 'id': i, 'deletable': False} for i in t1.columns
        # omit the id column
        if i != 'id'
    ],
    #This to_dict function adds a level to column headers, making them tuples if not dropped after using groupby. See above.
    data=t1.to_dict('records'),
    editable=False,
    filter_action="native",
    sort_action="native",
    sort_mode='multi',
    row_selectable='single',
    row_deletable=False,
    selected_rows=[],
    page_action='native',
    page_current= 0,
    page_size= 10,
    ),
    html.Br(),
    html.Div(id='LJ_graph'),
    dcc.Graph(id = 'LJ_graph2'),
    html.Br(),
    html.Div(id='dummy1'),
    html.Br(),
    html.Div(id='dummy2'),
    html.Div(id = 'store1',style={'display': 'none'}, children = tab0), #empty divs allow us to only compute / load the expensive database filtering upon initialization, rather than in each callback
    html.Div(id = 'store2',style={'display': 'none'}, children = tab1)
    ])

app.layout = serve_layout

#Callbacks and functions


#Function for printing the active variant on the screen
@app.callback(
    Output(component_id='LJ_graph', component_property='children'),
    Input('table', 'data'),
    Input('table', 'selected_rows'))

def print_selection(data, selected_rows):
    if selected_rows is None:
        selected_rows = []
    df = data[selected_rows[0]]['variant'] if selected_rows else "No selected rows"
    out = str(df)
    return ''.join(out) if df else "No Variant Selected"

#now for the first graph

@app.callback(
    Output('LJ_graph2', 'figure'), #outputs to plotly figure
    Input('table','data'), # gets data from the table
    Input('table','selected_rows'), # gets the seleccted active row
    Input('store1', 'children')) #Gets the stored empty-div dataframe

def update_graph(data, selected_rows, tab0):
    if selected_rows is None:
        selected_rows = []
    var = data[selected_rows[0]]['variant'] if selected_rows else "chr1_115256530_G_T_snp_1" #default selection
    data_long = pd.read_json(tab0).copy() #read the json
    is_var = data_long['variant'] == var #get only the active variant
    filt_dat = data_long[is_var]
    mn = filt_dat['afreq_normcount'].tolist() #calculate statistics for plotting
    sdpos1 = st.median(mn) + np.std(mn)
    sdneg1 = st.median(mn) -(np.std(mn))
    sdpos2 = st.median(mn) + 2*(np.std(mn))
    sdneg2 = st.median(mn) - 2*(np.std(mn))
    if len(filt_dat.index) >50: # if more than 50 samples, get the latest 50 (df is sorted by date)
        filt_dat = filt_dat[-50:]
    figure = go.Figure(data = go.Scatter(x = filt_dat['samplename'], y = mn[-50:], mode='lines+markers', name = 'Value'))
    figure.add_hline(y = sdpos1, line_width=2, line_color="green", name = '+1SD')
    figure.add_hline(y = sdneg1, line_width=2, line_color="green", name = '-1SD')
    figure.add_hline(y = sdpos2, line_width=2, line_color="red", name = '+2SD')
    figure.add_hline(y = sdneg2, line_width=2, line_color="red", name = '-2SD')
    figure.update_xaxes(showticklabels=False)
    return figure

#The second sample-wise graphs

@app.callback(
    Output('table-2','data'), #plot table
    Input('drpdown','value'), #get active sample from dropdown
    Input('store1', 'children'), #get data from dummy div
    Input('store2', 'children')) # get table 2 from dummy div

def update_table2(sel_value, tab0, tab1):
    t0 = pd.read_json(tab0)
    t0 = t0.applymap(lambda x: round(x, 5) if isinstance(x, (int, float)) else x) #when reading from json, floats are 15 decimal points long for some reason.... need to round
    t1 = pd.read_json(tab1)
    t1 = t1.applymap(lambda x: round(x, 5) if isinstance(x, (int, float)) else x)
    filt_dat = t0[t0['samplename'] == sel_value].copy()
    filt_dat['upper_bound'] = t1['upper_bound'].tolist()
    filt_dat['lower_bound'] = t1['lower_bound'].tolist()
    return filt_dat.to_dict('records')

#Refresh button
@app.callback(
    Output('output-container-button', 'children'), #refresh button calls serve_layout.. can also hit refresh on browser... 
    Input('apply-button', 'n_clicks'))

def refresh(n_clicks):
    if not n_clicks:
        return dash.no_update
    else:
        app.layout = serve_layout #not 100% sure this works

#dotplots

@app.callback(
        Output('subplot-div','figure'),
        Input('drpdown','value'),
        Input('store1','children'))


def make_dotplots(value, tab0):
    data_dot = pd.read_json(tab0).copy()
    data_dot.loc[data_dot.samplename != value, 'samplename'] = "Group"
    data_dot.loc[data_dot.samplename == value, 'samplename'] = "Current Sample"
    figure = px.box(data_dot, x = "variant", y = "afreq_normcount",color="samplename", facet_col="variant")
    figure.update_yaxes(matches=None, showticklabels=False)
    figure.update_xaxes(matches=None, tickangle=45, title="")
    figure.for_each_annotation(lambda a: a.update(text=""))
    return figure



@app.callback(
    Output('output-container-button2', 'children'),
    Input('drpdown','value')) #this comment manager takes the dropdown value as input. 

def display_notes(value):
    override = pd.read_csv('comments.txt', sep='\t', names=['sample','time','comment']) #read the comments file to obtain any current comments on the sample.
    filt_or = override.loc[override['sample'] == value] #find the selected dropdown value in the frame
    try:
        slist = filt_or['comment'].tolist()
        slist = [str(i) for i in slist]
        print(slist)
        slist = list(filter(('nan').__ne__, slist))
        return '\n'.join(slist) #print comments
    except:
        return 'No comments'

@app.callback(
    Output("dummy1", "children"), #needed for function, output really does to text file
    Input('button','n_clicks'),
    State('drpdown','value'), #see above, comment manager takes dropdown as input
    State('input-box', 'value'))

def update_notes(n_clicks, drpdown ,input_box):
    now = datetime.now().date()
    f = open('comments.txt', 'a') # file to add comments to
    writer = csv.writer(f, delimiter = "\t")
    row = [drpdown, now, input_box] #row has the sample ID, the date and comments. 
    print(row)
    writer.writerow(row)
    f.close()
    return None

@app.callback(
    Output("dummy2", "children"), #needed for function, output really does to text file
    Input('button2','n_clicks'),
    State('drpdown','value')) #see above, comment manager takes dropdown as input

def remove_outlier(n_clicks, drpdown):
    #insert function to delete and refresh here
    conn = mysql.connector.connect(host="localhost", database = 'HD200',user=" ", passwd=" ")
    cur = conn.cursor()
    q = 'DELETE FROM RunInfo WHERE name = %s'
    val = (drpdown,)
    cur.execute(q, val)
    cur.execute('DELETE FROM CallData WHERE sample IS NULL')
    conn.commit()
    conn.close()

#run on this server with IP

if __name__ == '__main__':
    app.run_server(debug=False, host='10.111.243.16')
