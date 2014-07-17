import itertools
import ROOT
import varial.tools
import varial.history
import varial.generators as gen


##################################################### convenience functions ###
legend_tags = ["real", "fakeGamma", "fakeOther", "fake"]


def color(key):
    for tag in legend_tags:
        if tag in key:
            return varial.analysis.get_color(tag)
    if "SihihSB" in key:
        return varial.analysis.get_color("fake")
    return 0


def get_legend_name(key):
    if "Template" in key:
        for tag in legend_tags:
            if tag in key:
                return varial.analysis.get_pretty_name(tag)
    else:
        return varial.analysis.get_pretty_name(key)


def cosmetica1(wrps):
    for w in wrps:
        name = get_legend_name(w.analyzer)
        w.histo.SetTitle(name)
        w.histo.SetYTitle(w.histo.GetYaxis().GetTitle() + " / a.u.")
        w.legend = name
        w.histo.SetLineColor(color(w.analyzer))
        w.histo.SetLineWidth(3)
        w.histo.SetFillStyle(0)
        yield w


def cosmetica2(wrps):
    for w in wrps:
        name = get_legend_name(w.analyzer)
        w.histo.SetTitle(name)
        w.legend = name
        w.histo.SetLineColor(1)
        w.histo.SetLineWidth(1)
        w.histo.SetFillColor(color(w.analyzer))
        w.histo.SetFillStyle(1001)
        yield w


def find_x_range(data_hist):
    x_min = data_hist.GetXaxis().GetXmin()
    x_max = data_hist.GetXaxis().GetXmax()
    for i in xrange(data_hist.GetNbinsX()):
        if data_hist.GetBinContent(i):
            x_min = data_hist.GetXaxis().GetBinLowEdge(i)
            break
    for i in xrange(data_hist.GetNbinsX(), 0, -1):
        if data_hist.GetBinContent(i):
            x_max = data_hist.GetXaxis().GetBinUpEdge(i)
            break
    return x_min - 1e-10, x_max + 1e-10


####################################################### Core Fitter Classes ###
class Fitter(object):
    def __init__(self):
        self.x_min = 0.
        self.x_max = 0.
        self.fit_func = None
        self.fitted = None
        self.mc_tmplts = None

    def build_fit_function(self, fitted, mc_tmplts, x_min, x_max):
        templates = [tw.histo for tw in mc_tmplts]
        size = len(templates)
        self.x_min, self.x_max = x_min, x_max
        self.fitted = fitted
        self.mc_tmplts = mc_tmplts

        def fit_func(x, par):
            value = 0.
            for j, hist in enumerate(templates):
                value += par[j] * hist.GetBinContent(hist.FindBin(x[0]))
            return value

        self.fit_func = ROOT.TF1(
            "MyFitFunc",
            fit_func,
            x_min,
            x_max,
            size
        )
        for i in xrange(0, size):
            self.fit_func.SetParameter(i, 1.)

    def do_the_fit(self):
        self.fitted.histo.Fit(
            self.fit_func, "WL M N", "", self.x_min, self.x_max
        )

    def scale_templates_to_fit(self, templates):
        for i in range(0, len(templates)):
            templates[i].histo.Scale(self.fit_func.GetParameter(i))

    def get_val_err(self, i_par):
        return (
            self.fit_func.GetParameter(i_par),
            self.fit_func.GetParError(i_par)
        )

    def get_ndf(self):
        return self.fit_func.GetNDF()

    def make_fit_result(self, result_wrp, mc_tmplts):
        r = result_wrp
        r.Chi2 = gen.op.chi2(
            (gen.op.sum(mc_tmplts), self.fitted),
            self.x_min, self.x_max
        ).float
        r.NDF = self.get_ndf()
        r.FitProb = ROOT.TMath.Prob(r.Chi2, r.NDF)
        r.dataIntegral = self.fitted.histo.Integral()
        r.legend = []
        r.value = []
        r.error = []
        r.binIntegralMC = []
        r.binIntegralScaled = []
        r.binIntegralScaledError = []
        for i, tmplt in enumerate(mc_tmplts):
            r.legend.append(tmplt.legend)
            val, err = self.get_val_err(i)
            r.value.append(val)
            r.error.append(err)
            r.binIntegralMC.append(tmplt.histo.Integral() / r.value[-1])
            r.binIntegralScaled.append(tmplt.histo.Integral())
            r.binIntegralScaledError.append(
                r.binIntegralScaled[-1] * r.error[-1] / r.value[-1]
            )


############################################################### Loading ... ###
class MassHistoSlicer(varial.tools.Tool):

    def __init__(self, name=None, slices=None):
        super(MassHistoSlicer, self).__init__(name)
        self.slices = slices or [(0, 2), (10, 80)]

    def run(self):
        @varial.history.track_history
        def get_slice_from_th2d(wrp, bin_low, bin_high):
            name = wrp.name + 'from%dto%d' % (bin_low, bin_high)
            histo = wrp.histo.ProjectionY(name, bin_low, bin_high)
            histo = histo.Clone()
            histo.SetTitle(wrp.legend)
            return varial.wrappers.HistoWrapper(histo, **wrp.all_info())

        def slice_generator(wrps):
            for wrp in wrps:
                for low, high in self.slices:
                    yield get_slice_from_th2d(wrp, low, high)

        wrps = filter(
            lambda w: isinstance(w.histo, ROOT.TH2D),
            self.lookup('../FSHistoLoader')
        )
        self.result = list(slice_generator(wrps))


class FitHistosCreator(varial.tools.Tool):
    """Result: [fitted_histo, template1, template2, ...]"""

    def run(self):
        tmplts = itertools.ifilter(
            lambda w: 'Template' in w.name,
            self.lookup("../FSHistoLoader")
        )
        tmplts = (
            varial.wrappers.HistoWrapper(w.histo, **w.all_info())
            for w in tmplts
        )
        tmplts = list(cosmetica2(tmplts))
        if not tmplts[0].is_data:
            data_lumi = varial.analysis.data_lumi_sum_wrp()
            tmplts = list(gen.op.prod((w, data_lumi)) for w in tmplts)

        fitted = [self.lookup('../MassHistoSlicer')[0]]

        self.result = fitted + tmplts


################################################################## Fit Tool ###
class TemplateFitTool(varial.tools.FSPlotter):
    def __init__(self, name=None):
        super(TemplateFitTool, self).__init__(name)
        self.mc_tmplts = None
        self.fitted = None
        self.fitbox_bounds = 0.63, 0.93, 0.60
        self.result = varial.wrappers.Wrapper()
        self.n_templates = 0
        self.fitter = Fitter()
        self.x_min = 0.
        self.x_max = 0.
        self.save_name_lambda = lambda w: w.name.split("_")[1]
        self.input_result_path = 'FitHistosCreator'

        def fix_ratio_histo_name(cnvs):
            for c in cnvs:
                d = c.get_decorator(varial.rendering.BottomPlotRatioSplitErr)
                d.dec_par["y_title"] = "Fit residual"
                yield c

        def set_no_exp(cnvs):
            for c in cnvs:
                c.first_drawn.GetYaxis().SetNoExponent()
                c.bottom_hist.GetYaxis().SetTitleSize(0.14)
                c.canvas.Modified()
                c.canvas.Update()
                yield c

        self.hook_pre_canvas_build = fix_ratio_histo_name
        self.hook_post_canvas_build = set_no_exp
        #TODO self.canvas_decorators.append(com.LumiTitleBox)

    def make_fit_textbox(self):
        res = self.result
        x1, x2, y2 = self.fitbox_bounds
        y1 = y2 - ((len(res.legend) + 1) * 0.04)

        textbox = ROOT.TPaveText(x1, y1, x2, y2, "brNDC")
        textbox.SetBorderSize(1)
        textbox.SetLineColor(0)
        textbox.SetLineStyle(1)
        textbox.SetLineWidth(0)
        textbox.SetFillColor(0)
        textbox.SetFillStyle(1001)
        textbox.SetTextSize(varial.settings.box_text_size)
        textbox.AddText(
            "#chi^{2} / NDF = %d / %i" % (round(res.Chi2, 2), res.NDF)
        )

        text = []
        for i, legend in enumerate(res.legend):
            text.append(
                "N_{%s} = %d #pm %d" % (
                    legend,
                    res.binIntegralScaled[i],
                    res.binIntegralScaledError[i]
                )
            )
        for txt in reversed(text):
            textbox.AddText(txt)
        return textbox

    def set_up_content(self):
        self.result.fitter = self.fitter.__class__.__name__

        self.fitted = self.stream_content.pop(0)
        self.mc_tmplts = self.stream_content
        self.n_templates = len(self.mc_tmplts)

        # do fit procedure
        self.x_min, self.x_max = find_x_range(self.fitted.histo)
        self.fitter.build_fit_function(
            self.fitted, self.mc_tmplts, self.x_min, self.x_max
        )
        self.fitter.do_the_fit()
        self.fitter.scale_templates_to_fit(self.mc_tmplts)
        self.fitter.make_fit_result(self.result, self.mc_tmplts)

        fit_textbox = self.make_fit_textbox()
        self.canvas_decorators.append(
            varial.rendering.TextBoxDecorator(
                None,
                True,
                textbox_dict={self.mc_tmplts[0].name: fit_textbox}
            )
        )

        tmplt_stack = gen.op.stack(self.mc_tmplts)
        self.stream_content = [(tmplt_stack, self.fitted)]

        del self.fitter


############################################################ Fit ToolChains ###
test_fit_chain = varial.tools.ToolChain(
    'TestFitChain', [
        varial.tools.FSHistoLoader(                         # HistoLoader
            filter_keyfunc=lambda w: (w.name, w.analyzer) in [
                ('VertexBeeMassTemplate', 'IvfMergedFilt'),
                ('VertexDeeMassTemplate', 'IvfMergedFilt'),
                ('VertexMassVsDr', 'IvfB2cMergedFilt'),
            ]
        ),
        MassHistoSlicer(                                    # MassHistoSlice
            slices=[0, 4]
        ),
        FitHistosCreator(),                                 # FitHistosCreator
        TemplateFitTool(),                                  # TemplateFitTool
    ]
)


fitter_chain = varial.tools.ToolChain(
    'VtxMassFitter', [
        varial.tools.FSHistoLoader(
            filter_keyfunc=lambda w: w.name in [
                'VertexBeeMassTemplate',
                'VertexDeeMassTemplate',
                'VertexMassVsDr'
            ]
        ),
        MassHistoSlicer(),
        varial.tools.FSPlotter(
            'TemplatePlots',
            filter_keyfunc=lambda w: 'Template' in w.name
        ),
        varial.tools.FSPlotter(
            'MassSlicePlots',
            input_result_path="../MassHistoSlicer",
            hook_loaded_histos=gen.sort
        ),
        test_fit_chain
    ]
)