from dominate.util import raw

def _summary(graphs):
    """
    Compose the summary section of the page
    :param graphs:
    :return: a string with HTML code for the module list
    """
    result = "        <ul class=\"menu-vertical\">\n"
    for i, t in enumerate(graphs):
        result += "          <li class=\"mv-item\"><a href=\"#M" + str(i) + "\">" + t + "</a></li>\n"
    result += "        </ul>\n"
    return result


def div_summary_style(titles):
    style = """
            <style>
            #leftCol {
            position: -webkit-sticky;
            position: sticky;
            top: 10px;
            font-size: 0.9em;
            margin: 0.5em 0 0 0.5em;
            width: calc(250px - 0.01em);
            }
        
            .qc-toc {
                display: none;
            }
        
            .qc-report .qc-toc {
                display: flex;
            }
            .menu-vertical {
            right: 1em;
            padding: 1px;
            margin-right : 1em;
            list-style: none;
            text-align: left;
            background: #F2F2F2;
            }
        
            .mv-item, .mv-item a {
            display: block;
            }
        
            .mv-item a {
            margin: 1px 0;
            padding: 8px 20px;
            color: #666;
            background: #FFF;
            text-decoration: none;
            transition: all .3s;
            }
        
            .mv-item a:hover, .mv-item a:focus {
            background: rgba(var(--bs-dark-rgb), var(--bs-bg-opacity)) !important;
            color: #FFF;
            padding-left: 30px;
            }
        
            .mv-item a:enabled {
            font-weight: bold;
            }
            
            </style>
            """
    menu = """
            <div id="leftCol">
                  <!--h2>Summary</h2-->
            {summary_list}
            </div>
            """.format(summary_list=_summary(titles))
    html = style+menu
    return html
    

def head_style():
    style = raw(f"""
                body {{ font-family: Arial, sans-serif; }}
                header {{ display: flex; 
                          justify-content: space-between; 
                          align-items: center; 
                          background-color: #EBEBEB; 
                          padding: 10px 20px; 
                          border-bottom: 1px solid #ddd; }}
                          
                h1.title {{ font-size: 14px; }}
                
                h1.title2 {{
                    font-size: 16px;
                    max-width: 400px;
                    border-bottom: 2px solid #000000; /* Exemple avec une bordure noire de 2px */
                    padding-bottom: 10px; /* Ajoute un peu d'espace entre le texte et la bordure */
                }}
                p.date {{ font-size: 12px; 
                          color: #666; }}
                
                .container {{ display: flex; 
                              justify-content: space-between; 
                              background-color: #FFFFFF}}
                              
                .content {{ flex: 2;
                            padding: 20px;
                            max-width: 1400px;
                            margin-left: 0;
                            margin-right: auto;
                            }}
                            
                .aside {{ flex: 1; border-right: 1px solid #ddd; padding: 20px; }}
                .tabs {{ display: flex; 
                         margin-top: 20px; 
                         margin-bottom: 20px; }}
                
                .tab-button {{
                    background-color: #f4f4f4;
                    padding: 10px;
                    cursor: pointer;
                    border: none; 
                    border-radius: 4px;  
                }}
                
                .tab-button.active {{
                    background-color: #052F61;
                    color: white;
                    border-radius: 10px; 
                }}
                
                .tab-content {{ display: none; }}
                .tab-content.active {{ display: block; }}
                """)
    return style


def body_style():
    style = """
        body {
            font-family: 'Arial', sans-serif;
            background-color: #f4f4f4;
            margin: 0;
            padding: 0;
        }
        .table-container {
            width: 80%;
            margin: 50px auto;
            box-shadow: 0 2px 8px rgba(0,0,0,0.2);
        }
        table {
            width: 100%;
            font-size: 14px;
            border-collapse: collapse;
            background-color: #fff;
            border-radius: 8px;
            overflow: hidden;
        }
        th, td {
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }
        th {
            background-color: #EBEBEB;
            color: #666;
            
        }
        tr:nth-child(even) {
            background-color: #f2f2f2;
        }
        tr:hover {
            background-color: #e9f1fe;
        }
        """
    return style