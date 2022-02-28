#!/usr/bin/env python3

#works OK for the DataTable portion, could benefit from adding genes and other annotations to fusions (will have to be done in code) 
#Next step will be to add callbacks to generate Levy Jennings on-the-fly

#Dependencies
import dash
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
    mydb = mysql.connector.connect(host="localhost", database = '',user="", passwd="")
    query = "SELECT CallData.pass_filter,  CallData.afreq,  CallData.coverage,  CallData.norm_count,  CallData.sample, VarData.name, RunInfo.IonWF_version, RunInfo.name, RunInfo.filedate, Transcripts.name , HGVS.transcript, HGVS.HGVSc, HGVS.HGVSp, Genes.name FROM VarData LEFT JOIN HGVS ON HGVS.id = VarData.hgvs LEFT JOIN CallData ON VarData.id = CallData.variant LEFT JOIN RunInfo ON CallData.sample = RunInfo.id LEFT JOIN Transcripts ON Transcripts.id = HGVS.transcript LEFT JOIN Genes ON VarData.gene = Genes.id;"
    df = pd.read_sql(query,mydb)
    mydb.close()

    #read-in data and change duplicated column headers

    #df = pd.read_csv('CAP_compare.gui.txt', sep='\t', header=0)
    #df.rename(columns={'name': 'variant', 'name.1': 'samplename', 'name.2': 'trname', 'name.3': 'gene'}, inplace=True)
    df.columns = ['pass_filter','afreq','coverage','norm_count','sample','variant','IonWF_version','samplename','filedate','trname','transcript','HGVSc', 'HGVSp','gene']
    #print(df.head())

    #get only HD200 and seracare samples
    df = df[df['samplename'].str.contains("HD200")]

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
tabledf = get_sql()
tabledf1 = tabledf.groupby('variant', as_index=False).agg({'samplename': pd.Series.nunique, 'coverage': ['mean'], 'afreq_normcount': ['mean'], 'trname': pd.Series.mode, 'HGVSc': pd.Series.mode, 'HGVSp': pd.Series.mode, 'filedate': ['max'], 'gene': pd.Series.mode, 'sd': ['mean'], 'upper_bound': ['mean'], 'lower_bound':['mean']})
tabledf1.columns = tabledf1.columns.droplevel(1)
neworder = ['variant','gene','afreq_normcount','sd', 'upper_bound', 'lower_bound','coverage','trname','HGVSc','HGVSp','samplename','filedate']
tabledf = tabledf[neworder]
tabledf1 = tabledf1[neworder]
tabledf = tabledf.applymap(lambda x: round(x, 5) if isinstance(x, (int, float)) else x)
tabledf1 = tabledf1.applymap(lambda x: round(x, 5) if isinstance(x, (int, float)) else x)


#create the app layout
app = dash.Dash(__name__)

#defining the layout
app.layout = html.Div(children=[
    #First is a title and a refresh button
    #Then a datatable with selectable rows for later graphs with callback
    html.Div([html.Button('Refresh', id='apply-button', n_clicks=0),
              html.Div(id='output-container-button', children='Click the button to update.')]),
    html.H1(children="HD200 / SeraCare controls"),
    dash_table.DataTable(
    id='table',
    columns=[
        {'name': i, 'id': i, 'deletable': False} for i in tabledf1.columns
        # omit the id column
        if i != 'id'
    ],
    #This to_dict function adds a level to column headers, making them tuples if not dropped after using groupby. See above.
    data=tabledf1.to_dict('records'),
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
    dcc.Dropdown(id = 'drpdown', options = [{'label': i, 'value': i} for i in tabledf['samplename'].unique()], value = tabledf['samplename'].tolist()[-1]),
    # add the conditionnal styling to allele frequencies that are out-of-bounds (2SD)
    dash_table.DataTable(id = 'table-2',
        columns=[{'name':i, 'id':i, 'deletable': False} for i in tabledf.columns if i not in ['id', 'afreqsd','normcountsd']],
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
    html.Div(id='dummy1')
])

#Callbacks and functions

@app.callback(
    Output(component_id='LJ_graph', component_property='children'),
    Input('table', 'data'),
    Input('table', 'selected_rows'))

def print_selection(data, selected_rows):
    if selected_rows is None:
        selected_rows = []
    df = data[selected_rows[0]]['variant'] if selected_rows else "No selected rows"
    out = [str(df), str(selected_rows)]
    return ' '.join(out) if df else "No Variant Selected"

#now for the first graph

@app.callback(
    Output('LJ_graph2', 'figure'),
    Input('table','data'),
    Input('table','selected_rows'))

def update_graph(data, selected_rows):
    if selected_rows is None:
        selected_rows = []
    var = data[selected_rows[0]]['variant'] if selected_rows else "chr1_115256530_G_T_snp_1"
    data_long = tabledf.copy()
    is_var = data_long['variant'] == var
    filt_dat = data_long[is_var]
    if len(filt_dat.index) > 50:
        filt_dat = filt_dat[-50:]
    else:
        filt_dat = filt_dat
    mn = filt_dat['afreq_normcount'].tolist()
    sdpos1 = [(st.median(mn[0:i]) + np.std(mn[0:i])) for i in range(1,len(mn))]
    sdpos1.insert(0,st.median(mn) + np.std(mn))
    sdpos1[1] = st.median(mn) + np.std(mn)
    sdneg1 = [(st.median(mn[0:i]) - np.std(mn[0:i])) for i in range(1,len(mn))]
    sdneg1.insert(0, st.median(mn) -(np.std(mn)))
    sdneg1[1] = st.median(mn) - np.std(mn)
    sdpos2 = [(st.median(mn[0:i]) + 2*np.std(mn[0:i])) for i in range(1,len(mn))]
    sdpos2.insert(0, st.median(mn) + 2*(np.std(mn)))
    sdpos2[1] = st.median(mn) + 2*(np.std(mn))
    sdneg2 = [(st.median(mn[0:i]) - 2*(np.std(mn[0:i]))) for i in range(1,len(mn))]
    sdneg2.insert(0,st.median(mn)-2*(np.std(mn)))
    sdneg2[1] = st.median(mn) - 2*(np.std(mn))
    figure = go.Figure(data = go.Scatter(x = filt_dat['samplename'], y = mn, mode='lines+markers', name = 'Value'))
    figure.add_trace(go.Scatter(x = filt_dat['samplename'], y = sdpos1, mode='lines', name = '+1SD'))
    figure.add_trace(go.Scatter(x = filt_dat['samplename'],y = sdneg1, mode='lines', name = '-1SD'))
    figure.add_trace(go.Scatter(x = filt_dat['samplename'],y = sdpos2, mode='lines', name = '+2SD'))
    figure.add_trace(go.Scatter(x = filt_dat['samplename'],y = sdneg2, mode='lines', name = '-2SD'))
    figure.update_xaxes(showticklabels=False)
    return figure

#The second sample-wise graphs

@app.callback(
    Output('table-2','data'),
    Input('drpdown','value')
)

def update_table2(sel_value):
    filt_dat = tabledf[tabledf['samplename'] == sel_value].copy()
    filt_dat['upper_bound'] = tabledf1['upper_bound'].tolist()
    filt_dat['lower_bound'] = tabledf1['lower_bound'].tolist()
    return filt_dat.to_dict('records')

#Refresh button
@app.callback(
    dash.dependencies.Output('output-container-button', 'children'),
    [dash.dependencies.Input('apply-button', 'n_clicks')])

def refresh(n_clicks):
    if not n_clicks:
        return dash.no_update
    else:
        tabledf = get_sql() # maybe also put a call to the scp function here...

#dotplots

@app.callback(
        Output('subplot-div','figure'),
        Input('drpdown','value')
)

def make_dotplots(value):
    data_dot = tabledf.copy()
    data_dot.loc[data_dot.samplename != value, 'samplename'] = "Group"
    data_dot.loc[data_dot.samplename == value, 'samplename'] = "Current Sample"
    figure = px.box(data_dot, x = "variant", y = "afreq_normcount",color="samplename", facet_col="variant")
    #figure = go.Figure(data = go.box(facet_col = data_dot['variant'], y = data_dot['afreq_normcount'], color = data_dot['samplename'], points = "All", facet_col_wrap=6))
    figure.update_yaxes(matches=None, showticklabels=False)
    figure.update_xaxes(matches=None, tickangle=45, title="")
    figure.for_each_annotation(lambda a: a.update(text=""))
    return figure



@app.callback(
    Output('output-container-button2', 'children'),
    Input('drpdown','value'))

def display_notes(value):
    override = pd.read_csv('comments.txt', sep='\t', names=['sample','time','comment'])
    filt_or = override.loc[override['sample'] == value]
    try:
        #return filt_or[['time','comment']]
        return '\n'.join(filt_or['comment'].tolist())
    except:
        return 'No comments'

@app.callback(
    Output("dummy1", "children"),
    Input('button','n_clicks'),
    State('drpdown','value'),
    State('input-box', 'value'))

def update_notes(n_clicks, drpdown ,input_box):
    now = datetime.now().date()
    f = open('comments.txt', 'w')
    writer = csv.writer(f, delimiter = "\t")
    row = [drpdown, now, input_box]
    print(row)
    writer.writerow(row)
    f.close()
    return None


#run on this server with IP

if __name__ == '__main__':
    app.run_server(debug=False, host='')
