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

sliderPDE = Slider(start=10, end=100, value=100, step=5, title="PDE")
sliderMinWL = Slider(start=150, end=450, value=150, step=50, title="min wl")
sliderASIC = Slider(start=0, end=30, value=0, step=10, title="asic")
sliderSPTR = Slider(start=0, end=80, value=0, step=10, title="sptr")

#plot2
p2 = figure(plot_height=600,plot_width=600, y_range=(-100,100), title=None, min_border=10, min_border_left=50, webgl=True)
p2.xgrid.grid_line_color = None
p2.select(BoxSelectTool).select_every_mousemove = False
p2.select(LassoSelectTool).select_every_mousemove = False
scatterPlot = p2.scatter(x=[], y=[], fill_alpha=0.3, line_color=None)

callbackHist = CustomJS(args=dict(
            source=hist.data_source,
            sPlot=scatterPlot.data_source,
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
        var asic = sASIC.get('value')
        var sptr = sSPTR.get('value')

url = 'http://localhost:5000/data/210/' + pde + '/' + minwl + '/1500/' + asic + '/' + sptr
$.getJSON(url, function(data) {
        var dataOld = source.get('data');
        for (i = 0; i < dataOld['right'].length; i++) {
            dataOld['right'][i] = data.right[i]
            dataOld['top'][i] = data.top[i]
            dataOld['left'][i] = data.left[i]
        }

        var sc2 = sPlot.get('data');
        x = sc2['x']
        y = sc2['y']
        x.splice(0,x.length)
        y.splice(0,y.length)
        for (i = 0; i < data['x'].length; i++) {
            x[i] = data['x'][i]
            y[i] = data['y'][i]
        }

        txtData = txtSrc.get('data');
        txtData['text'][0] = 'speed: ' + data.speed
        txtMinData = txtMinSrc.get('data');
        txtMinData['text'][0] = 'minWL: ' + data.minwl
        txtMaxData = txtMaxSrc.get('data');
        txtMaxData['text'][0] = 'maxWL: ' + data.maxwl
});

        source.trigger('change');
        sPlot.trigger('change');
        txtSrc.trigger('change');
        txtMinSrc.trigger('change');
        txtMaxSrc.trigger('change');

    """)

sliderPDE.callback = callbackHist
sliderMinWL.callback = callbackHist
sliderASIC.callback = callbackHist
sliderSPTR.callback = callbackHist


sliders = hplot(sliderPDE,sliderMinWL,sliderASIC,sliderSPTR)
layout = vplot(vplot(sliders,hplot(ph,p2)), width=800, height=800)

output_file("scatter.html", title="color_scatter.py example")
show(layout)
