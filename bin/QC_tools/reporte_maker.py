import pandas as pd
import re
import dominate
from dominate.tags import *
from dominate.util import raw
from datetime import date, datetime
from QC_tools.css_style import body_style, head_style, div_summary_style

import plotly.express as px
df2 = px.data.iris()
fig = px.scatter(df2, x="sepal_width", y="sepal_length", color='petal_length')

def list_to_dict(params):
    return {"--"+k: v for k, v in (s.split(":", 1) for s in params) if v != "null"}


def string_to_dict(params):
    params_list = params[0].split(', ')
    params_list = {k: v for k, v in (s.split(":", 1) for s in params_list) if v != "null"}
    return params_list
    
def parse_versions(versions):
    try:
        versions_list = {k: v for d, k, v in (re.split(r'(Badread|SPARSim)', s) for s in versions)}
    except:
        versions_list = {k:k for k in versions}
    return versions_list


def read_logo(image_path):
    import base64
    with open(image_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
    return encoded_string

def reporter(stats_table, pie_fig, Q_over_time, atgc, gc, conf_params, work_params, project, versions):
    work_params = string_to_dict(work_params)
    conf_params = list_to_dict(conf_params)

    versions = parse_versions(versions)

    conf_params.update(work_params)

    conf_params.update(versions)

    df = pd.DataFrame(list(conf_params.items()), columns=['Params', 'Value'])


    Basic_inputs_df = df[df['Params'].isin(['--matrix', '--bc_counts', '--transcriptome', 
                                            '--features', '--gtf', '--sim_celltypes', '--cell_types_annotation', 
                                            '--error_model', '--fastq_model', '--ref_genome', '--trained_model', 
                                            '--qscore_model'])]

    simulation_parames_df = df[df['Params'].isin(['--full_length', '--umi_duplication', '--pcr_cycles', 
                                                '--pcr_error_rate', '--pcr_dup_rate', '--pcr_total_reads',
                                                '--badread_identity', '--dT_LENGTH', '--ADPTER_SEQ', '--TSO_SEQ'])]

    run_params_df = df[df['Params'].isin(['container', 'runOptions', 'dsl2', 'startTime', 'threads', 'commandLine'])]

    versions_df = df[df['Params'].isin(['nextflow', 'SPARSim', 'Badread'])]

    doc = dominate.document(title='Styled Table with Tabs')
    current_date = raw(str(date.today().strftime('%b-%d-%Y'))+'<br>'+str(datetime.now().strftime('%H:%M:%S')))

    with doc.head:
        style(head_style())
    
    with doc.body:
        style(body_style())
    
    with doc:
        with header():
            with div():
                h1('AsaruSim v0.1.0', cls='title')
                p(current_date, cls='date')
            img_logo = read_logo('asarusim_v2.png')
            img(src=f'data:image/png;base64,{img_logo}',height='60', alt='Logo', style="float: right;")
            
        with div(cls='container'):
            with div():
                titles = ['Simulation parameters', 
                          'Simulation statistics',
                          'Qscore over sequence position', 
                          'Base content over position', 
                          'Per base GC content']
                raw(div_summary_style(titles))

            with div(cls='content'):
                h1('Simulation parameters', cls='title2')
                with div(cls='tabs'):
                    button('Basic inputs', cls='tab-button active', onclick="openTab(event, 'tab1')")
                    button('Simulation config', cls='tab-button', onclick="openTab(event, 'tab2')")
                    button('Run parameters', cls='tab-button', onclick="openTab(event, 'tab3')")
                    button('Package versions', cls='tab-button', onclick="openTab(event, 'tab4')")

                with div(style='margin-bottom: 20px; min-height: 300px;'):
                    with div(id='tab1', cls='tab-content active'):
                        with table():
                            with thead():
                                with tr():
                                    for col in Basic_inputs_df.columns:
                                        th(col)
                            with tbody():
                                for _, row in Basic_inputs_df.iterrows():
                                    with tr():
                                        for cell in row:
                                            td(cell)
                    
                    with div(id='tab2', cls='tab-content'):
                        with table():
                            with thead():
                                with tr():
                                    for col in simulation_parames_df.columns:
                                        th(col)
                            with tbody():
                                for _, row in simulation_parames_df.iterrows():
                                    with tr():
                                        for cell in row:
                                            td(cell)
                                            
                    with div(id='tab3', cls='tab-content'):
                        with table():
                            with thead():
                                with tr():
                                    for col in run_params_df.columns:
                                        th(col)
                            with tbody():
                                for _, row in run_params_df.iterrows():
                                    with tr():
                                        for cell in row:
                                            td(cell)
                    with div(id='tab4', cls='tab-content'):
                        with table():
                            with thead():
                                with tr():
                                    for col in versions_df.columns:
                                        th(col)
                            with tbody():
                                for _, row in versions_df.iterrows():
                                    with tr():
                                        for cell in row:
                                            td(cell)

                h1('Simulation statistics', cls='title2')
                with div(style='display: flex; justify-content: space-between'):
                    with div(style='max-width: 600px; flex: 2; margin-top: 80px'):
                        raw(stats_table.to_html())
                    with div(style='max-width: 800px; flex: 2'):
                        raw(pie_fig.to_html())

                h1('Read QC', cls='title2')  
                with div(style='margin-bottom: 20px'):
                    
                    raw(Q_over_time.to_html())
                    raw(atgc.to_html())
                    raw(gc.to_html())
                
        with footer():
            with div(style='font-size: 12px; padding: 10px 20px;;'):
                app_url="https://alihamraoui.github.io/AsaruSim/introduction/"
                p(raw(f'Generated by <a href="{app_url}">AsaruSim</a> (version 0.1.0)'))


    with doc.body:
        script(raw("""
        function openTab(evt, tabName) {
            var i, tabcontent, tablinks;
            tabcontent = document.getElementsByClassName("tab-content");
            for (i = 0; i < tabcontent.length; i++) {
                tabcontent[i].style.display = "none";
            }
            tablinks = document.getElementsByClassName("tab-button");
            for (i = 0; i < tablinks.length; i++) {
                tablinks[i].className = tablinks[i].className.replace(" active", "");
            }
            document.getElementById(tabName).style.display = "block";
            evt.currentTarget.className += " active";
        }
        """))
    out_name = project+'.html'
    with open(out_name, 'w') as f:
        f.write(doc.render())

if __name__ == '__main__':
    main()
