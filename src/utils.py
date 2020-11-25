import os
import nbformat
from nbconvert import HTMLExporter

def convert_notebook(outdir, report_in_dir, report_out_dir):
#     curdir = os.path.abspath(os.getcwd())
#     indir, _ = os.path.split(report_in_dir)
#     outdir, _ = os.path.split(report_out_dir)
    #out = 
    os.mkdir(outdir)
    
    config = {
        "ExecutePreprocessor": {"enabled": True},
        "TemplateExporter": {"exclude_output_prompt": True, "exclude_input": True, "exclude_input_prompt": True}
    }
    
    nb = nbformat.read(open(report_in_dir), as_version=4)
    html_exporter = HTMLExporter(config=config)
    
    #change into directory with notebook, to execute notebook
    os.chdir(report_in_dir)
    body, resources = (
        html_exporter.from_notebook_node(nb)
    )
    
    current = os.chdir(os.path.abspath(os.getcwd()))
    
    
    with open(report_out_path, 'w') as fh:
        fh.write(body)