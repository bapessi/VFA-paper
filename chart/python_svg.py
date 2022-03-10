from bs4 import BeautifulSoup
import numpy as np
import pandas as pd

file_svg = 'data/metabolic_chart_base.svg'
file_chart= 'data/chart_reactions_positions.xlsx'
 
def change_width(file_svg,positions,width_list):

    with open(file_svg) as f:
        svg = BeautifulSoup(f,'xml')
    a = svg.find(id='OXA-PEP')

    for p in svg.find_all('path'):
        if 'stroke-width' in p['style']:
            if p['id'] in positions:
                print(p['id'])
                style_text = p['style'].partition('stroke-width:')
                width_value = style_text[2].split(';')[0]
                right_side = style_text[2].partition(width_value)[2]
                width_value = width_list[np.where(p['id']==positions)][0]
                width_value = str(width_value)
                p['style'] = style_text[0] + style_text [1] + width_value + right_side 

    fileout = 'data/test.svg'
    with open(fileout,'w') as f:
        f.write(svg.prettify())


df = pd.read_excel(file_chart,index_col=0)
positions = df['Position'].values
width_list = df['width'].values
change_width(file_svg,positions,width_list)
