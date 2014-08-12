import itertools
import ROOT
import varial.analysis
import varial.dbio
import varial.tools
import varial.history
import varial.generators as gen
import varial.operations as op
import varial_result


##################################################### convenience functions ###
legend_tags = ["real", "fakeGamma", "fakeOther", "fake"]
re_bins = list(i / 10. for i in xrange(0, 100, 5))
#re_bins = list(i / 10. for i in xrange(10, 60, 2))
#re_bins = list(i - .5 for i in xrange(0, 21))


def gen_set_legend_and_color(wrps, sig_token, sig="Signal", bg="Background"):
    for w in wrps:
        if sig_token in w.name:
            name = sig
            col = ROOT.kSpring - 4
        else:
            name = bg
            col = ROOT.kRed + 2
        w.histo.SetTitle(name)
        w.legend = name
        w.histo.SetFillColor(col)
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
class HistoSlicer(varial.tools.Tool):
    io = varial.dbio

    def __init__(self, slices, name=None):
        super(HistoSlicer, self).__init__(name)
        self.slices = slices

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
        wrps = slice_generator(wrps)
        #wrps = gen.gen_rebin(wrps, re_bins)
        self.result = list(wrps)


class FitHistosCreator(varial.tools.Tool):
    def __init__(self, template_name, signal_name, name=None):
        super(FitHistosCreator, self).__init__(name)
        self.template_name = template_name
        self.signal_name = signal_name

    def run(self):
        # input
        inp = itertools.ifilter(
            lambda w: self.template_name in w.name,
            self.lookup("../FSHistoLoader")
        )
        inp = itertools.chain(inp, self.lookup('../HistoSlicer'))
        inp = gen.sort(inp)
        inp = gen.group(inp)
        inp = gen.mc_stack_n_data_sum(inp)
        inp = itertools.chain.from_iterable(inp)
        inp = (
            varial.wrappers.HistoWrapper(w.histo, **w.all_info())
            for w in inp
        )
        inp = gen.gen_rebin(inp, re_bins)
        inp = list(inp)
        assert(len(inp) == 3)

        # grab templates
        tmplts = inp[1:3]
        tmplts = gen_set_legend_and_color(tmplts, self.signal_name)
        tmplts = list(tmplts)

        # grab fitted histogram
        fitted = inp[0:1]
        if not fitted[0].is_data:
            fitted[0].legend = 'Pseudo-Data'

        # normalize templates to starting values
        integral = fitted[0].histo.Integral()
        for t in tmplts:
            t.histo.Scale(
                integral / t.histo.Integral() / len(tmplts)
            )

        self.result = fitted + tmplts


class FitHistosCreatorSum(varial.tools.Tool):
    def __init__(self, run_on_mc, name=None):
        super(FitHistosCreatorSum, self).__init__(name)
        self.run_on_mc = run_on_mc

    def run(self):
        inp = self.lookup('../HistoSlicer')
        total = gen.sort(inp)
        total = gen.group(total)
        total = gen.mc_stack_n_data_sum(total)
        total = itertools.chain.from_iterable(total)
        total = (
            varial.wrappers.HistoWrapper(w.histo, **w.all_info())
            for w in total
        )
        total = filter(lambda w: w.is_data != self.run_on_mc, total)
        #fitted, highDR = sorted(total, key=lambda w: 'to80' in w.name)
        fitted = total[0]

        bee_tmplt = next(
            gen.gen_norm_to_data_lumi(
                filter(lambda w: 'to80' not in w.name
                                 and w.sample == 'TTbarTwoMatch', inp)))
        dee_tmplt = next(
            gen.gen_norm_to_data_lumi(
                filter(lambda w: 'to80' not in w.name
                                 and w.sample == 'TTbarBDMatch', inp)))
        # bee_highDR = next(
        #     gen.gen_norm_to_data_lumi(
        #         filter(lambda w: 'to80' in w.name
        #                          and w.sample == 'TTbarTwoMatch', inp)))
        # bee_for_diff = op.prod((
        #     op.norm_to_integral(bee_tmplt), op.integral(bee_highDR)
        # ))
        # fke_tmplt = gen.op.diff((highDR, bee_for_diff))
        fke_tmplt = next(
            gen.gen_norm_to_data_lumi(
                filter(lambda w: 'to80' not in w.name
                                 and w.sample == 'TTbarOneMatch', inp)))

        bee_tmplt.legend = '2 Bee Vertices'
        dee_tmplt.legend = 'Bee + Dee Vertex'
        fke_tmplt.legend = 'Bee + Fake Vertex'
        bee_tmplt.histo.SetTitle('2 Bee Vertices')
        dee_tmplt.histo.SetTitle('Bee + Dee Vertex')
        fke_tmplt.histo.SetTitle('Bee + Fake Vertex')
        bee_tmplt.histo.SetFillColor(ROOT.kSpring - 4)
        dee_tmplt.histo.SetFillColor(ROOT.kRed - 7)
        fke_tmplt.histo.SetFillColor(ROOT.kRed + 2)
        fitted.legend = 'Fit Histo'

        # normalize templates to starting values
        integral = fitted.histo.Integral()
        for t in (bee_tmplt, dee_tmplt, fke_tmplt):
            t.lumi = 1.
            t.histo.Scale(
                integral / (t.histo.Integral() or 1.) / 3.
            )

        self.result = list(gen.gen_rebin(
            (fitted, fke_tmplt, dee_tmplt, bee_tmplt), re_bins))


################################################################## Fit Tool ###
import varial.rendering as rnd

class TemplateFitTool(varial.tools.FSPlotter):
    def __init__(self, input_result_path, name=None):
        super(TemplateFitTool, self).__init__(
            name, input_result_path=input_result_path)
        self.mc_tmplts = None
        self.fitted = None
        self.fitbox_bounds = 0.63, 0.93, 0.60
        self.result = varial.wrappers.Wrapper()
        self.n_templates = 0
        self.fitter = Fitter()
        self.x_min = 0.
        self.x_max = 0.
        self.save_name_lambda = lambda w: w.name.split("_")[1]

        def fix_ratio_histo_name(cnvs):
            for c in cnvs:
                d = c.get_decorator(rnd.BottomPlotRatioSplitErr)
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

    def load_content(self):
        pass

    def set_up_content(self):
        self.result.fitter = self.fitter.__class__.__name__

        wrps = self.lookup(self.input_result_path)
        self.fitted = wrps.pop(0)
        self.mc_tmplts = wrps
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

slices = [(0, 6), (20, 80)]

fitter_plots = varial.tools.ToolChain(
    'VtxMassFitterPlots', [
        varial.tools.FSHistoLoader(
            filter_keyfunc=lambda w: w.name in ['VertexBeeMassTemplate',
                                                'VertexDeeMassTemplate',
                                                'VertexMassVsDr']
        ),
        HistoSlicer(slices),
        varial.tools.FSPlotter(
            'TemplatePlots',
            filter_keyfunc=lambda w: 'Template' in w.name
                                     and 'IvfMergedFilt' in w.analyzer
        ),
        varial.tools.FSPlotter(
            'MassSlicePlots',
            input_result_path="../HistoSlicer",
            hook_loaded_histos=gen.sort
        )
    ]
)


fitter_chain = varial.tools.ToolChain(
    'TestFitChain', [
        varial.tools.FSHistoLoader(
            filter_keyfunc=lambda w: (w.name, w.analyzer) in [
                ('VertexDeeMassTemplate', 'IvfMergedFilt'),
                ('VertexMassVsDr', 'IvfB2cMerged'),
            ] and not w.is_data
        ),
        HistoSlicer(slices),
        FitHistosCreator('Template', 'VertexMassVsDr'),
        varial.tools.FSPlotter(
            'TemplatesAndFittedHisto',
            input_result_path='../FitHistosCreator',
            plot_grouper=None,
            plot_setup=lambda w: [list(varial.tools.overlay_colorizer(w, [1, 2, 3]))]
        ),
        TemplateFitTool(input_result_path='../FitHistosCreator'),
    ]
)


def _mkchnsm(name, slice, coll):
    return varial.tools.ToolChain(
        'FitChainSum'+name, [
            varial.tools.FSHistoLoader(
                filter_keyfunc=lambda w: w.name == 'VertexMassSumVsDr'
                                         and w.analyzer == coll,
            ),
            HistoSlicer([slice]),
            FitHistosCreatorSum(True),
            varial.tools.FSPlotter(
                'TemplatesAndFittedHisto',
                input_result_path='../FitHistosCreatorSum',
                plot_grouper=gen.gen_copy,
                plot_setup=lambda w: [list(varial.tools.overlay_colorizer(w, [1, 2, 3, 4]))]
            ),
            varial.tools.FSPlotter(
                'MassSlicePlots',
                input_result_path="../HistoSlicer",
                hook_loaded_histos=gen.sort
            ),
            TemplateFitTool(input_result_path='../FitHistosCreatorSum'),
            FitHistosCreatorSum(False, name='FitHistosCreatorSumData'),
            TemplateFitTool(name='TemplateFitToolData',
                            input_result_path='../FitHistosCreatorSumData'),
        ]
    )

fitter_chain_sum = varial.tools.ToolChain(
    'FitChainSum', list(
        _mkchnsm("from%02dto%02d_%s" % (s[0], s[1], c), s, c)
        for c in ('IvfB2cMerged', 'IvfB2cMergedCuts')
        for s in ((0, 6), (6, 10), (10, 14), (14, 18))
    ) + [varial_result.summary_chain]
)
