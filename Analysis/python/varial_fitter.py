import itertools
import os
import ROOT

import varial.analysis
import varial.dbio
import varial.tools
import varial.history
import varial.generators as gen
import varial.operations as op
import varial_result

import theta_auto
theta_auto.config.theta_dir = os.environ["CMSSW_BASE"] + "/theta"


##################################################### convenience functions ###
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

    def get_chi2(self):
        return 0.

    def make_fit_result(self, result_wrp, mc_tmplts):
        r = result_wrp
        r.Chi2 = self.get_chi2() or gen.op.chi2(
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


class ThetaFitter(Fitter):
    def __init__(self):
        super(ThetaFitter, self).__init__()
        self.model = None
        self.fit_res = None
        self.ndf = 0

    def _store_histos_for_theta(self, wrp):
        filename = os.path.join(varial.analysis.cwd, wrp.name + ".root")
        f = ROOT.TFile.Open(filename, "RECREATE")
        f.cd()
        for key, value in wrp.__dict__.iteritems():
            if isinstance(value, ROOT.TH1):
                value.SetName(key)
                value.Write()
        f.Close()

    def build_fit_function(self, fitted, mc_tmplts, x_min, x_max):
        self.x_min, self.x_max = x_min, x_max
        self.fitted = fitted
        self.mc_tmplts = mc_tmplts

        theta_root_wrp = varial.wrappers.Wrapper(
            name="ThetaHistos",
            histo__DATA=fitted.histo,
        )
        self.template_names = []
        for i, tmplt in enumerate(mc_tmplts):
            name = 'template%02d' % (i+1)
            self.template_names.append(name)
            setattr(theta_root_wrp, 'histo__' + name, tmplt.histo)
        self._store_histos_for_theta(theta_root_wrp)
        theta_auto.config.workdir = varial.analysis.cwd
        self.model = theta_auto.build_model_from_rootfile(
            os.path.join(varial.analysis.cwd, "ThetaHistos.root"),
            include_mc_uncertainties=True
        )
        self.model.set_signal_processes([self.template_names[-1]])
        for tmplt_name in self.template_names[:-1]:
            self.model.add_lognormal_uncertainty(
                "bg_" + tmplt_name, 1., tmplt_name)
            self.model.distribution.set_distribution_parameters(
                "bg_" + tmplt_name,
                width=50.  # theta_auto.inf
            )
        self.ndf_sum = sum(
            1
            for i in xrange(fitted.histo.GetNbinsX())
            if fitted.histo.GetBinContent(i+1) > .05
        ) - len(self.template_names)

    def do_the_fit(self):
        options = theta_auto.Options()
        options.set('minimizer', 'strategy', 'robust')
        self.fit_res = theta_auto.mle(
            self.model,
            "data",
            1,
            chi2=True,
            options=options,
        )[self.template_names[-1]]
        print self.fit_res

    def scale_templates_to_fit(self, templates):
        par_values = {
            "beta_signal": self.fit_res["beta_signal"][0][0],
        }
        for tmplt_name in self.template_names[:-1]:
            bg_name = "bg_" + tmplt_name
            par_values[bg_name] = self.fit_res[bg_name][0][0]

        self.val_err = []
        for tmplt_name in self.template_names[:-1]:
            val = self.model.get_coeff("histo", tmplt_name).get_value(par_values)
            err = val * abs(self.fit_res["bg_" + tmplt_name][0][1]
                            / (self.fit_res["bg_" + tmplt_name][0][0] or 1e10))
            self.val_err.append((val, err))
        self.val_err.append((
            self.fit_res["beta_signal"][0][0],
            abs(self.fit_res["beta_signal"][0][1])
        ))
        for tmplt, val_err in zip(templates, self.val_err):
                tmplt.histo.Scale(val_err[0])

    def get_val_err(self, i_par):
        return self.val_err[i_par]

    def get_ndf(self):
        return self.ndf_sum

    def get_chi2(self):
        return self.fit_res['__chi2'][0]

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
            histo = wrp.histo.ProjectionY('', bin_low, bin_high)
            histo = histo.Clone()
            histo.SetName(name)
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
    def __init__(self, template_name, signal_name, re_bins=None, name=None):
        super(FitHistosCreator, self).__init__(name)
        self.template_name = template_name
        self.signal_name = signal_name
        self.re_bins = re_bins

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
        if self.re_bins:
            inp = gen.gen_rebin(inp, self.re_bins)
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
    def __init__(self, run_on_mc, re_bins=None, name=None):
        super(FitHistosCreatorSum, self).__init__(name)
        self.run_on_mc = run_on_mc
        self.re_bins = re_bins

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
        # fke_tmplt = next(gen.gen_norm_to_data_lumi(
        #     [op.merge(
        #         itertools.ifilter(
        #             lambda w: 'to80' not in w.name
        #                       and not w.is_data
        #                       and w.sample not in ['TTbarBDMatch',
        #                                            'TTbarTwoMatch'],
        #             inp
        #         )
        #     )]
        # ))
        fke_tmplt = next(
            gen.gen_norm_to_data_lumi(
                filter(lambda w: 'to80' not in w.name
                                 and w.sample == 'TTbarOneMatch', inp)))

        bee_tmplt.legend = 'B + B Vertex'
        dee_tmplt.legend = 'B + D Vertex'
        fke_tmplt.legend = 'B + Fake Vertex'
        bee_tmplt.histo.SetTitle('B + B Vertex')
        dee_tmplt.histo.SetTitle('B + D Vertex')
        fke_tmplt.histo.SetTitle('B + Fake Vertex')
        bee_tmplt.histo.SetFillColor(ROOT.kSpring - 4)
        dee_tmplt.histo.SetFillColor(ROOT.kRed - 7)
        fke_tmplt.histo.SetFillColor(ROOT.kRed + 2)
        if not fitted.is_data:
            fitted.legend = 'Fit Histo'

        # normalize templates to starting values
        if False:  # normalization turned off
            integral = fitted.histo.Integral()
            for t in (bee_tmplt, dee_tmplt, fke_tmplt):
                t.lumi = 1.
                t.histo.Scale(
                    integral / (t.histo.Integral() or 1.) / 3.
                )

        res = [fitted, fke_tmplt, dee_tmplt, bee_tmplt]
        if self.re_bins:
            list(gen.gen_rebin(res, self.re_bins))
        self.result = res


class FitHistosCreatorSimultaneous(varial.tools.Tool):
    def __init__(self, run_on_mc, re_bins=None, name=None):
        super(FitHistosCreatorSimultaneous, self).__init__(name)
        self.run_on_mc = run_on_mc
        self.re_bins = re_bins

    def run(self):
        @varial.history.track_history
        def place_next_to(wrps):
            w1, w2 = wrps
            h1, h2 = w1.histo, w2.histo
            h1_n_bins, h2_n_bins = h1.GetNbinsX(), h2.GetNbinsX()
            histo = ROOT.TH1D(
                '%s_%s' % (w1.name, w2.name),
                ';%s;%s' % (h1.GetXaxis().GetTitle(),
                            h1.GetYaxis().GetTitle()),
                h1_n_bins + h2_n_bins,
                h1.GetXaxis().GetXmin(),
                h1.GetXaxis().GetXmax()
                    - h2.GetXaxis().GetXmin()
                    + h2.GetXaxis().GetXmax(),
            )
            #sqrt2 = 2**.5
            for j in range(h1_n_bins):
                i = j + 1
                histo.SetBinContent(i, h1.GetBinContent(i))  # / 2.)
                histo.SetBinError(i, h1.GetBinError(i))  # / sqrt2)
            for j in range(h2_n_bins):
                i = j + h1_n_bins + 2
                histo.SetBinContent(i, h2.GetBinContent(j+1))  # / 2.)
                histo.SetBinError(i, h2.GetBinError(j+1))  # / sqrt2)
            return varial.wrappers.HistoWrapper(
                histo,
                **w1.all_info()
            )

        inp = self.lookup('../HistoSlicer')
        inp_lead = itertools.ifilter(lambda w: 'PtLead' in w.name, inp)
        inp_subl = itertools.ifilter(lambda w: 'PtSubLead' in w.name, inp)
        inp_lead = gen.gen_rebin(inp_lead, self.re_bins)
        inp_subl = gen.gen_rebin(inp_subl, self.re_bins)
        inp_lead = gen.sort(inp_lead)
        inp_subl = gen.sort(inp_subl)

        inp = list(place_next_to(wrps)
                   for wrps in itertools.izip(inp_lead, inp_subl))

        total = itertools.ifilter(lambda w: w.is_data != self.run_on_mc, inp)
        total = gen.group(total)
        total = gen.mc_stack_n_data_sum(total)
        total = itertools.chain.from_iterable(total)
        total = (
            varial.wrappers.HistoWrapper(w.histo, **w.all_info())
            for w in total
        )
        total = list(total)
        fitted = total[0]

        bee_tmplt = next(
            gen.gen_norm_to_data_lumi(
                filter(lambda w: 'to80' not in w.name
                                 and w.sample == 'TTbarTwoMatch', inp)))
        dee_tmplt = next(
            gen.gen_norm_to_data_lumi(
                filter(lambda w: 'to80' not in w.name
                                 and w.sample == 'TTbarBDMatch', inp)))
        fke_tmplt = next(
            gen.gen_norm_to_data_lumi(
                filter(lambda w: 'to80' not in w.name
                                 and w.sample == 'TTbarOneMatch', inp)))
        # fke_tmplt = next(gen.gen_norm_to_data_lumi(
        #     [op.merge(
        #         itertools.ifilter(
        #             lambda w: 'to80' not in w.name
        #                       and not w.is_data
        #                       and w.sample not in ['TTbarBDMatch',
        #                                            'TTbarTwoMatch'],
        #             inp
        #         )
        #     )]
        # ))

        bee_tmplt.legend = 'B + B Vertex'
        dee_tmplt.legend = 'B + D Vertex'
        fke_tmplt.legend = 'B + Fake Vertex'
        bee_tmplt.histo.SetTitle('B + B Vertex')
        dee_tmplt.histo.SetTitle('B + D Vertex')
        fke_tmplt.histo.SetTitle('B + Fake Vertex')
        bee_tmplt.histo.SetFillColor(ROOT.kSpring - 4)
        dee_tmplt.histo.SetFillColor(ROOT.kRed - 7)
        fke_tmplt.histo.SetFillColor(ROOT.kRed + 2)
        fitted.legend = 'Fit Histo'

        # normalize templates to starting values
        if False:  # normalization turned off
            integral = fitted.histo.Integral()
            for t in (bee_tmplt, dee_tmplt, fke_tmplt):
                t.lumi = 1.
                t.histo.Scale(
                    integral / (t.histo.Integral() or 1.) / 3.
                )

        if dee_tmplt.histo.Integral():
            self.result = [fitted, fke_tmplt, dee_tmplt, bee_tmplt]
        else:
            self.result = [fitted, fke_tmplt, bee_tmplt]


################################################################## Fit Tool ###
import varial.rendering as rnd

class TemplateFitTool(varial.tools.FSPlotter):
    def __init__(self, input_result_path, name=None):
        super(TemplateFitTool, self).__init__(
            name, input_result_path=input_result_path)
        self.mc_tmplts = None
        self.fitted = None
        self.fitbox_bounds = 0.3, 0.6, 0.87
        self.result = varial.wrappers.Wrapper()
        self.n_templates = 0
        self.fitter = ThetaFitter()
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
            plot_setup=lambda w: [list(varial.tools.overlay_colorizer(
                w, [1, ROOT.kRed + 2, ROOT.kRed - 7, ROOT.kSpring - 4]))]
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
            FitHistosCreatorSum(False, name='FitHistosCreatorSumData'),
            varial.tools.FSPlotter(
                'TemplatesAndFittedHisto',
                input_result_path='../FitHistosCreatorSumData',
                plot_grouper=gen.gen_copy,
                plot_setup=lambda ws: [list(varial.tools.overlay_colorizer(
                    filter(lambda w: not w.is_data, ws),
                    [ROOT.kRed + 2, ROOT.kRed - 7, ROOT.kSpring - 4]))],
                canvas_decorators=[varial.rendering.Legend],
            ),
            varial.tools.FSPlotter(
                'MassSlicePlots',
                input_result_path="../HistoSlicer",
                hook_loaded_histos=gen.sort
            ),
            TemplateFitTool(name='TemplateFitToolData',
                            input_result_path='../FitHistosCreatorSumData'),
        ]
    )


def _mkchnsmltn(slice, coll, re_bins):
    name = "from%02dto%02d_%s" % (slice[0], slice[1], coll)
    return varial.tools.ToolChain(
        'FitChainSum'+name, [
            varial.tools.FSHistoLoader(
                filter_keyfunc=lambda w: w.name in ['VertexMassPtLeadVsDr',
                                                    'VertexMassPtSubLeadVsDr']
                                         and w.analyzer == coll,
            ),
            HistoSlicer([slice]),
            FitHistosCreatorSimultaneous(
                False, re_bins=re_bins, name='FitHistosCreatorSimultaneous'),
            varial.tools.FSPlotter(
                'TemplatesAndFittedHisto',
                input_result_path='../FitHistosCreatorSimultaneous',
                plot_grouper=gen.gen_copy,
                plot_setup=lambda ws: [list(
                    varial.tools.overlay_colorizer(
                        gen.apply_histo_linewidth(
                            gen.gen_norm_to_integral(
                                filter(lambda w: not w.is_data, ws)
                            )
                        ),
                        [ROOT.kRed + 2, ROOT.kRed - 7, ROOT.kSpring - 4]
                    )
                )],
                canvas_decorators=[varial.rendering.Legend],
            ),
            varial.tools.FSPlotter(
                'MassSlicePlots',
                input_result_path="../HistoSlicer",
                hook_loaded_histos=gen.sort
            ),
            TemplateFitTool(name='TemplateFitToolData',
                            input_result_path='../FitHistosCreatorSimultaneous'),
        ]
    )


re_bins_1 = list(i / 10. for i in xrange(0, 60, 1))
re_bins_2 = list(i / 10. for i in xrange(0, 60, 2))
re_bins_3 = list(i / 10. for i in xrange(0, 60, 3))
re_bins_4 = list(i / 10. for i in xrange(0, 60, 4))
dr_bins = ((0, 10), (10, 14), (14, 17), (17, 20), (20, 22))
c1 = 'IvfB2cMerged'
c2 = 'IvfB2cMergedCuts'
fitter_chain_sum = varial.tools.ToolChain(
    'FitChainSum', [
        _mkchnsmltn((0, 10), c1, re_bins_2),
        _mkchnsmltn((10, 14), c1, re_bins_2),
        _mkchnsmltn((14, 17), c1, re_bins_2),
        _mkchnsmltn((17, 20), c1, re_bins_2),
        _mkchnsmltn((20, 22), c1, re_bins_2),
        _mkchnsmltn((0, 10), c2, re_bins_2),
        _mkchnsmltn((10, 14), c2, re_bins_2),
        _mkchnsmltn((14, 17), c2, re_bins_3),
        _mkchnsmltn((17, 20), c2, re_bins_2),
        _mkchnsmltn((20, 22), c2, re_bins_2),
        varial_result.summary_chain
    ]
)
#TODO self.canvas_decorators.append(com.LumiTitleBox)