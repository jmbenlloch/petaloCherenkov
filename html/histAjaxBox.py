''' Present a scatter plot with linked histograms on both axes.
    bokeh serve selection_histogram.py
    http://localhost:5006/selection_histogram
'''
import numpy as np
from bokeh.models import BoxSelectTool, LassoSelectTool, Paragraph, Slider, CustomJS
from bokeh.plotting import figure, hplot, vplot, output_file, show
from bokeh.models.sources import AjaxDataSource

source = AjaxDataSource(data_url='http://localhost:5000/data/10', method='GET',polling_interval=10000)

ph = figure(plot_height=600,plot_width=600,x_range=(-100,100),y_range=(0, 1000), title=None, min_border=10, min_border_left=50)
ph.xgrid.grid_line_color = None

ph.select(BoxSelectTool).select_every_mousemove = False
ph.select(LassoSelectTool).select_every_mousemove = False

#hist = ph.quad(bottom=0, left=0, right=0, top=0, color="white", line_color="#3A5785")
hist = ph.quad(bottom=0, left=[0]*80, right=[0]*80, top=[0]*80, color="white", line_color="#3A5785")
txt = ph.text(70, 900, text=['speed: null'],text_color="firebrick",text_align="center",text_font_size="10pt")
txtMinWL = ph.text(70, 850, text=['minWL: null'],text_color="firebrick",text_align="center",text_font_size="10pt")
txtMaxWL = ph.text(70, 800, text=['maxWL: null'],text_color="firebrick",text_align="center",text_font_size="10pt")

ph.min_border_top = 10
ph.min_border_right = 10

#fitLine = ph.line([0], [0], line_color="#D95B43", line_width=8, alpha=0.7, legend="PDF")

sliderPDE = Slider(start=10, end=100, value=100, step=5, title="pde")
sliderMinWL = Slider(start=150, end=450, value=150, step=50, title="min wl")
sliderASIC = Slider(start=0, end=30, value=0, step=10, title="asic")
sliderSPTR = Slider(start=0, end=80, value=0, step=10, title="sptr")

#plot2
cats = map(str,range(100,1300,100))
p = figure(tools="save", background_fill_color="#EFE8E2", title="", x_range=cats, y_range=(-150,150))
sg1 = p.segment(cats,[0],cats,[0], line_width=2, line_color="black")
sg2 = p.segment(cats,[0],cats,[0], line_width=2, line_color="black")
r1 = p.rect(cats,[0],0.7,[0],fill_color="#E08E79", line_width=2, line_color="black")
r2 = p.rect(cats,[0],0.7,[0],fill_color="#3B8686", line_width=2, line_color="black")
r3 = p.rect(cats,[0],0.2,0.01, line_color="black")
r4 = p.rect(cats,[0],0.2,0.01, line_color="black")
cr = p.circle([0], [0], size=6, color="#F38630", fill_alpha=0.6)

p.xgrid.grid_line_color = None
p.ygrid.grid_line_color = "white"
p.grid.grid_line_width = 2
p.xaxis.major_label_text_font_size="12pt"


callbackHist = CustomJS(args=dict(
            source=hist.data_source,
            sg1Src=sg1.data_source,
            sg2Src=sg2.data_source,
            r1Src=r1.data_source,
            r2Src=r2.data_source,
            r3Src=r3.data_source,
            r4Src=r4.data_source,
            crSrc=cr.data_source,
            sPDE=sliderPDE,
            sMinWL=sliderMinWL,
            sASIC=sliderASIC,
            sSPTR=sliderSPTR,
            txtSrc=txt.data_source,
            txtMinSrc=txtMinWL.data_source,
            txtMaxSrc=txtMaxWL.data_source), code="""
(function(){
  var newscript = document.createElement('script');
     newscript.type = 'text/javascript';
     newscript.async = true;
     newscript.src = 'https://ajax.googleapis.com/ajax/libs/jquery/1.6.1/jquery.min.js';
  (document.getElementsByTagName('head')[0]||document.getElementsByTagName('body')[0]).appendChild(newscript);
})();

        var pde = sPDE.get('value')
        var minwl = sMinWL.get('value')
        var sptr = sSPTR.get('value')
        var asic = sASIC.get('value')

url = 'http://localhost:5000/dataBox/210/' + pde + '/' + minwl + '/1500/' + asic + '/' + sptr
$.getJSON(url, function(data) {
        var dataOld = source.get('data');
        for (i = 0; i < dataOld['right'].length; i++) {
            dataOld['right'][i] = data.right[i]
            dataOld['top'][i] = data.top[i]
            dataOld['left'][i] = data.left[i]
        }

        crData = crSrc.get('data');
        var x = crData['x']
        var y = crData['y']
        x.splice(0,x.length)
        y.splice(0,y.length)
        for (i = 0; i < data['outx'].length; i++) {
            x[i] = data.outx[i]
            y[i] = data.outy[i]
        }

        sg1Data = sg1Src.get('data');
        var sg1x0 = sg1Data['x0']
        var sg1x1 = sg1Data['x1']
        var sg1y0 = sg1Data['y0']
        var sg1y1 = sg1Data['y1']
        sg1x0.splice(0,sg1x0.length)
        sg1x1.splice(0,sg1x1.length)
        sg1y0.splice(0,sg1y0.length)
        sg1y1.splice(0,sg1y1.length)
        for (i = 0; i < data['upper'].length; i++) {
            sg1y0[i] = data.upper[i]
            sg1y1[i] = data.q3[i]
        }

        sg2Data = sg2Src.get('data');
        var sg2x0 = sg2Data['x0']
        var sg2x1 = sg2Data['x1']
        var sg2y0 = sg2Data['y0']
        var sg2y1 = sg2Data['y1']
        sg2x0.splice(0,sg2x0.length)
        sg2x1.splice(0,sg2x1.length)
        sg2y0.splice(0,sg2y0.length)
        sg2y1.splice(0,sg2y1.length)
        for (i = 0; i < data['upper'].length; i++) {
            sg2y0[i] = data.lower[i]
            sg2y1[i] = data.q1[i]
        }

        r1Data = r1Src.get('data');
        var r1x = r1Data['x']
        var r1y = r1Data['y']
        var r1Height = r1Data['height']
        r1x.splice(0,r1x.length)
        r1y.splice(0,r1y.length)
        r1Height.splice(0,r1Height.length)
        for (i = 0; i < data['q3'].length; i++) {
            r1y[i] = (data.q3[i] + data.q2[i])/2
            r1Height[i] = data.q3[i] - data.q2[i]
        }

        r2Data = r2Src.get('data');
        var r2x = r2Data['x']
        var r2y = r2Data['y']
        var r2Height = r2Data['height']
        r2x.splice(0,r2x.length)
        r2y.splice(0,r2y.length)
        r2Height.splice(0,r2Height.length)
        for (i = 0; i < data['q2'].length; i++) {
            r2y[i] = (data.q2[i] + data.q1[i])/2
            r2Height[i] = data.q2[i] - data.q1[i]
        }

        r3Data = r3Src.get('data');
        var r3x = r3Data['x']
        var r3y = r3Data['y']
        r3x.splice(0,r3x.length)
        r3y.splice(0,r3y.length)
        for (i = 0; i < data['lower'].length; i++) {
            r3y[i] = data.lower[i]
        }

        r4Data = r4Src.get('data');
        var r4x = r4Data['x']
        var r4y = r4Data['y']
        r4x.splice(0,r4x.length)
        r4y.splice(0,r4y.length)
        for (i = 0; i < data['upper'].length; i++) {
            r4y[i] = data.upper[i]
        }

        for (i = 0; i < data['cats'].length; i++) {
            r1x[i] = data.cats[i]
            r2x[i] = data.cats[i]
            r3x[i] = data.cats[i]
            r4x[i] = data.cats[i]
            sg1x0[i] = data.cats[i]
            sg1x1[i] = data.cats[i]
            sg2x0[i] = data.cats[i]
            sg2x1[i] = data.cats[i]
        }

        txtData = txtSrc.get('data');
        txtData['text'][0] = 'speed: ' + data.speed
        txtMinData = txtMinSrc.get('data');
        txtMinData['text'][0] = 'minWL: ' + data.minwl
        txtMaxData = txtMaxSrc.get('data');
        txtMaxData['text'][0] = 'maxWL: ' + data.maxwl

});

        source.trigger('change');
        sg1Src.trigger('change');
        sg2Src.trigger('change');
        r1Src.trigger('change');
        r2Src.trigger('change');
        r3Src.trigger('change');
        r4Src.trigger('change');
        crSrc.trigger('change');
        txtSrc.trigger('change');
        txtMinSrc.trigger('change');
        txtMaxSrc.trigger('change');
    """)

sliderPDE.callback = callbackHist
sliderMinWL.callback = callbackHist
sliderASIC.callback = callbackHist
sliderSPTR.callback = callbackHist


sliders = hplot(sliderPDE,sliderMinWL,sliderASIC,sliderSPTR)
layout = vplot(vplot(sliders,hplot(ph,p)), width=800, height=800)

output_file("box.html", title="color_scatter.py example")

show(layout)
