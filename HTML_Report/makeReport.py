import os
import sys
from jinja2 import Environment, FileSystemLoader, PackageLoader, Template
import shutil
import re
import pandas as pd
import pickle
import constants_and_types

def create_final_html_report(exp_param_d, site_param_d):
    """
    Creates a html report for a Crispector run.

    Parameters: exp_param_d, site_param_d

    """
    file_loader = FileSystemLoader('crispector_output')
    env = Environment(loader=file_loader)

    crispector_output = env.get_template('file_to_render.html')

    edit_section_result_table_html = create_edit_section_result_table(exp_param_d)
    translocations_trans_res_tab_tab_data_html = create_translocations_trans_res_tab_tab_data(exp_param_d)
    
    for site in exp_param_d["site_names"]:
        print(site)
        print(site_param_d[site])
        create_site_page(site, site_param_d[site])

    final_file = open("crispector_output/final_report.html", "w+")
    output = crispector_output.render(exp_param_d=exp_param_d, 
                                    edit_section_result_table_html=edit_section_result_table_html,
                                    translocations_trans_res_tab_tab_data_html=translocations_trans_res_tab_tab_data_html)
    
    final_file.write(output)
    final_file.close()

def create_edit_section_result_table(exp_param_d):
    df_edit_section_result_table = pd.DataFrame.from_dict(exp_param_d["edit_section"]["result_table"]["tabular_data"])
    df_edit_section_result_table = df_edit_section_result_table.transpose()
    df_edit_section_result_table.columns = df_edit_section_result_table.iloc[0]
    df_edit_section_result_table = df_edit_section_result_table.drop([0])
    edit_section_result_table_html = df_edit_section_result_table.to_html(classes="table table-hover", justify="center", border="0")
    return edit_section_result_table_html

def create_translocations_trans_res_tab_tab_data(exp_param_d):
    df_translocations_trans_res_tab_tab_data = pd.DataFrame.from_dict(exp_param_d["translocations"]["translocations_results_tab"]["tabular_data"])
    df_translocations_trans_res_tab_tab_data = df_translocations_trans_res_tab_tab_data.transpose()
    df_translocations_trans_res_tab_tab_data.columns = df_translocations_trans_res_tab_tab_data.iloc[0]
    df_translocations_trans_res_tab_tab_data = df_translocations_trans_res_tab_tab_data.drop([0]) 
    translocations_trans_res_tab_tab_data_html = df_translocations_trans_res_tab_tab_data.to_html(classes="table table-hover", justify="center", border="0")
    return translocations_trans_res_tab_tab_data_html

def create_site_page(site, site_param_d):
    file_loader = FileSystemLoader('crispector_output')
    env = Environment(loader=file_loader)

    report_name = site + "_report.html"
    f = open("crispector_output/" + report_name, "w+")
    crispector_output = env.get_template("site_file_to_render.html")

    output = crispector_output.render(site_param_d=site_param_d, site=site)
    
    f.write(output)
    f.close()