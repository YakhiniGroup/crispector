from typing import Dict
from crispector.utils.constants_and_types import Path, EDIT_SECTION, RESULT_TABLE, TAB_DATA, TRANSLOCATIONS, \
    TRANS_RES_TAB, HTML_SITES, REPORT_PATH, HTML_SITES_NAME_LIST
import click
import os
import sys
from jinja2 import Environment, FileSystemLoader
import pandas as pd
import pickle


def create_final_html_report(html_param_d: Dict, report_output: Path):
    """
    Creates a html report for a Crispector run.

    Parameters: exp_param_d, site_param_d

    """
    html_templates_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'html_templates')
    file_loader = FileSystemLoader(html_templates_path)
    env = Environment(loader=file_loader)

    crispector_output = env.get_template('file_to_render.html')

    edit_section_result_table_html = create_edit_section_result_table(html_param_d)

    if html_param_d['translocations']['translocations_results_tab'] != '':
        translocations_trans_res_tab_tab_data_html = create_translocations_trans_res_tab_tab_data(html_param_d)

    for site in html_param_d[HTML_SITES][HTML_SITES_NAME_LIST]:
        create_site_page(site, html_param_d, html_param_d[HTML_SITES][site], report_output)

    if html_param_d['translocations']['translocations_results_tab'] != '':
        with open(os.path.join(report_output, "report.html"), "w+") as file:
            output = crispector_output.render(exp_param_d=html_param_d,
                                              edit_section_result_table_html=edit_section_result_table_html,
                                              translocations_trans_res_tab_tab_data_html=translocations_trans_res_tab_tab_data_html)
            file.write(output)
            file.close()

    else:
        with open(os.path.join(report_output, "report.html"), "w+") as file:
            output = crispector_output.render(exp_param_d=html_param_d,
                                              edit_section_result_table_html=edit_section_result_table_html)

            file.write(output)
            file.close()


def create_edit_section_result_table(exp_param_d):
    df_edit_section_result_table = pd.DataFrame.from_dict(exp_param_d[EDIT_SECTION][RESULT_TABLE][TAB_DATA])
    df_edit_section_result_table = df_edit_section_result_table.transpose()
    df_edit_section_result_table.columns = df_edit_section_result_table.iloc[0]
    df_edit_section_result_table = df_edit_section_result_table.drop([0])
    edit_section_result_table_html = df_edit_section_result_table.to_html(classes="table table-hover", justify="center",
                                                                          border="0", index=False)
    return edit_section_result_table_html


def create_translocations_trans_res_tab_tab_data(exp_param_d):
    df_translocations_trans_res_tab_tab_data = pd.DataFrame.from_dict(
        exp_param_d[TRANSLOCATIONS][TRANS_RES_TAB][TAB_DATA])
    df_translocations_trans_res_tab_tab_data = df_translocations_trans_res_tab_tab_data.transpose()
    df_translocations_trans_res_tab_tab_data.columns = df_translocations_trans_res_tab_tab_data.iloc[0]
    df_translocations_trans_res_tab_tab_data = df_translocations_trans_res_tab_tab_data.drop([0])
    translocations_trans_res_tab_tab_data_html = df_translocations_trans_res_tab_tab_data.to_html(
        classes="table table-hover", justify="center", border="0", index=False)
    return translocations_trans_res_tab_tab_data_html


def create_site_page(site: str, exp_param_d: Dict, site_param_d: Dict, report_output: Path):
    html_templates_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'html_templates')
    file_loader = FileSystemLoader(html_templates_path)
    env = Environment(loader=file_loader)

    with open(os.path.join(report_output, site_param_d[REPORT_PATH]), "w+") as file:
        crispector_output = env.get_template("site_file_to_render.html")

        output = crispector_output.render(exp_param_d=exp_param_d, site_param_d=site_param_d, site=site)

        file.write(output)
        file.close()


@click.command()
@click.option('--output_path', type=click.Path(), required=True, help="")
def main(output_path):
    with open(os.path.join(output_path, "crispector_output/html_param_d.pkl"), "rb") as file:
        html_param_d = pickle.load(file)

    create_final_html_report(html_param_d, output_path)


if __name__ == '__main__':
    sys.exit(main())  # pragma: no cover
