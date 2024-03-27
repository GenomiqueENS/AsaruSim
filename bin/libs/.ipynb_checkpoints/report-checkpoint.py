import plotly.io as io

def make_report(Q_over_time, atgc, gc, project, pos):
    html_string = '''
    <html>
        <head>
            <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
            <style>body{ margin:0 100; background:whitesmoke; }</style>
        </head>
        <body>
            <h1>''' + project + '''</h1>
    
            <!-- *** Section 1 *** --->
            <h2>Qscore across ''' + str(pos) + ''' bases</h2>
            ''' + io.to_html(Q_over_time) + '''
            <p></p>
            
            <!-- *** Section 2 *** --->
            <h2>Sequence content across ''' + str(pos) + ''' bases</h2>
            ''' + io.to_html(atgc) + '''
            <p></p>
            <h3>Per base GC content</h3>
            ''' + io.to_html(gc) + '''
        </body>
    </html>'''

    out_name = project+'.html'
    f = open(out_name,'w')
    f.write(html_string)
    f.close()