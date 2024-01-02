from bokeh.models.widgets import DataTable, TableColumn
from bokeh.plotting import show, output_file
from bokeh.models import ColumnDataSource
from bokeh.layouts import column

import pandas as pd

def create_html(excel_file):

    sheets = ['Pool-Rep1-N7_S4','Pool-Rep2-N7_S5','Pool-Rep3-N7_S6','GT-Rep1-N7_S1','GT-Rep2-N7_S2','GT-Rep3-N7_S3']
    df = pd.read_excel(excel_file, sheet_name='Integration Percent')

    # Keeping only columns that contain the word 'Complete'
    complete_beacon_df = df[[col for col in df.columns if 'Target' in col or 'Complete Beacon' in col]].round(2)
    complete_beacon_df.columns = [i.split('_')[0] for i in complete_beacon_df.columns]
    PGI_df = df[[col for col in df.columns if 'Target' in col or 'Complete P' in col]].round(2)
    PGI_df.columns = [i.split('_')[0] for i in PGI_df.columns]

    print(complete_beacon_df)
    print(PGI_df)

    # Convert DataFrame to Bokeh DataSource
    source = ColumnDataSource(complete_beacon_df)

    # Create Columns for DataTable
    columns = [TableColumn(field=Ci, title=Ci) for Ci in complete_beacon_df.columns]

    # Create DataTable with specified width and fit columns
    data_table = DataTable(columns=columns, source=source, autosize_mode = 'force_fit',aspect_ratio='auto')

    # Wrap the table in a layout with responsive width
    layout = column(data_table, sizing_mode="stretch_width")

    # Specify the output HTML file
    output_file("first_figure.html")

    # Display the layout
    show(layout)

    

if __name__ == "__main__":
    input_excel_file = '../../TMB103_analysis/work/3e/7e60962aa905734e4ab76ad07243eb/TMB103_results.xlsx'
    create_html(input_excel_file)