import os
import sys
from jinja2 import Environment, FileSystemLoader, PackageLoader, Template
import shutil
import re
from PIL import Image
import pandas as pd


def make_report(exp_param_d):
    """
    Creates a html report for a Crispector run.

    Parameters: 

    Returns:
    """
    file_loader = FileSystemLoader('templates')
    env = Environment(loader=file_loader)

    template = env.get_template('report.html')

    df_result_table = pd.DataFrame.from_dict(exp_param_d["RESULT_TABLE"]["TAB_DATA"])
    result_table_html = df_result_table.to_html(classes="table table-hover", justify="center", border="0")

    output = template.render(exp_param_d=exp_param_d, result_table_html=result_table_html)
    print(output)