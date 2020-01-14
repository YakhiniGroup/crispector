from makeReport import make_report
from PIL import Image


exp_param_d = {}
exp_param_d["PAGE_TITLE"] = "I am an experiment title."
exp_param_d["EDITING_ACTIVITY"] = {}
exp_param_d["EDITING_ACTIVITY"]["TITLE"] = "I am an editing activity title"
exp_param_d["EDITING_ACTIVITY"]["PLOT_PATH"] = "C:/Users/Chi Hong Choi/Desktop/Crispector/templates/img/RAG1_sg_HEK293/editing_activity.png"
exp_param_d["EDITING_ACTIVITY"]["PDF_PATH"] = "C:/Users/Chi Hong Choi/Desktop/Crispector/templates/img/RAG1_sg_HEK293/editing_activity.png"
exp_param_d["READING_STATS"] = {}
exp_param_d["READING_STATS"]["TITLE"] = "Reading Stats Title"
exp_param_d["READING_STATS"]["MAPPING_STATS"] = {}
exp_param_d["READING_STATS"]["MAPPING_STATS"]["TITLE"] = "tab 1 title"
exp_param_d["READING_STATS"]["MAPPING_STATS"]["PLOT_PATH"] = "C:/Users/Chi Hong Choi/Desktop/Crispector/templates/img/RAG1_sg_HEK293/reads_statistics.png"
exp_param_d["READING_STATS"]["MAPPING_STATS"]["PDF_PATH"] = "C:/Users/Chi Hong Choi/Desktop/Crispector/templates/img/RAG1_sg_HEK293/reads_statistics.png"
exp_param_d["READING_STATS"]["MAPPING_PER_SITE"] = {}
exp_param_d["READING_STATS"]["MAPPING_PER_SITE"]["TITLE"] = "tab 2 title"
exp_param_d["READING_STATS"]["MAPPING_PER_SITE"]["PLOT_PATH"] = "C:/Users/Chi Hong Choi/Desktop/Crispector/templates/img/RAG1_sg_HEK293/reads_statistics.png"
exp_param_d["READING_STATS"]["MAPPING_PER_SITE"]["PDF_PATH"] = "C:/Users/Chi Hong Choi/Desktop/Crispector/templates/img/RAG1_sg_HEK293/reads_statistics.png"
exp_param_d["READING_STATS"]["DISCARDED_SITES"] = "We are discareded sites."
exp_param_d["READING_STATS"]["FASTP"] = {}
exp_param_d["READING_STATS"]["FASTP_TX_TEXT"] = "FASTP report for treatment"
exp_param_d["READING_STATS"]["FASTP_TX_PATH"] = "C:/Users/Chi Hong Choi/Desktop/Crispector/templates/img/RAG1_sg_HEK293/treatment_fastp/fastp.html"
exp_param_d["READING_STATS"]["FASTP_MOCK_TEXT"] = "FASTP report for mock"
exp_param_d["READING_STATS"]["FASTP_MOCK_PATH"] = "img/RAG1_sg_HEK293/treatment_fastp/fastp.html"
exp_param_d["RESULT_TABLE"] = {}
exp_param_d["RESULT_TABLE"]["TITLE"] = "I am a result table title"
exp_param_d["RESULT_TABLE"]["TAB_DATA"] = {}
exp_param_d["RESULT_TABLE"]["TAB_DATA"]["COL1"] = [1, 2, 3]
exp_param_d["RESULT_TABLE"]["TAB_DATA"]["COL2"] = [4, 5, 6]
exp_param_d["HTML_SITE_NAMES"] = ["site1", "site2", "site3"]
exp_param_d["LOG_PATH"] = "img/RAG1_sg_HEK293/treatment_fastp/fastp.html"

make_report(exp_param_d)