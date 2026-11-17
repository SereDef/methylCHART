from pathlib import Path

from shiny import App, Inputs, reactive, render, req, ui
from faicons import icon_svg
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.layout_engine import PlaceHolderLayoutEngine

# ==== Set-up ==================================================================
INPUT_DIR = Path(__file__).parent          # methylCHART
RESULTS_DIR = INPUT_DIR / "results"
DATA_DIR = INPUT_DIR / "data"

cpg_targets = json.loads((INPUT_DIR / "cpg_targets.json").read_text(encoding="utf-8"))
cpgs = list(cpg_targets)
cpg_labs = {c: f"{c} - {v['label']}" for c, v in cpg_targets.items()}
cpg_pheno = {c: [v["pheno"]] if isinstance(v["pheno"], str) else list(v["pheno"] or [])
             for c, v in cpg_targets.items()}

models = sorted((d.name for d in RESULTS_DIR.iterdir() if d.is_dir()), reverse=True)
ABOUT = (INPUT_DIR / "about.md").read_text(encoding="utf-8")

# Where "good" sits for each metric -- drawn as a dashed guide line.
METRIC_REF = {"EXPV": 0.0, "RMSE": 0.0, "SMSE": 1.0, "MSLL": 0.0, "MACE": 0.0,
              "R2": 0.0, "Rho": 0.0, "Rho_p": 0.05, "MLL": None,
              "ShapiroW": 1.0, "MAPE": 0.0}
DEFAULT_METRICS = ["EXPV", "RMSE", "SMSE", "MSLL", "MACE"]
MAX_METRICS = 5

TRAIN_COL, TEST_COL, HILITE = "#9aa0a6", "#1f77b4", "#d9534f"
STRIPE = "#f2f2f2"   # every other CpG row, so rows are easy to follow across panels
ROW_PX = 16        # vertical pixels per CpG row
LABEL_FS = 8       # CpG label font size
GAP = 0.012        # gap between panels, as a fraction of figure width
X_LAB = {"beta": "beta", "M": "M-value", "Z": "Z"}


def fmt(x):
    return "–" if x is None else f"{x:.2f}"


def stripe_rows(ax, n):
    """Alternating light bands, one per CpG row (replaces the horizontal grid)."""
    for i in range(1, n, 2):
        ax.axhspan(i - 0.5, i + 0.5, color=STRIPE, lw=0, zorder=0)


# ==== UI ======================================================================
app_ui = ui.page_fluid(
    ui.div(
        ui.h5("MethylCHART demo", class_="mb-0 me-auto pb-2"),
        ui.input_select("model", "Model", choices=models, width="560px"),
        ui.input_select("cpg", "CpG", choices=cpg_labs, width="420px"),
        ui.a(icon_svg("github", height="1.6em"),
             href="https://github.com/SereDef/methylCHART",
             target="_blank", title="GitHub repository", class_="ms-3 pb-2"),
        class_="d-flex align-items-end gap-3 flex-wrap px-3 pt-2 border-bottom mb-2",
    ),

    ui.navset_card_underline(

        # ---- model level: every CpG at once ---------------------------------
        ui.nav_panel(
            "Model overview",
            ui.output_ui("overview_boxes"),   # equal slots, 1st empty
            ui.output_ui("overview_ui"),      # figure split the same way
            ui.div(                           # controls live below the plots
                ui.input_checkbox_group("metrics", f"Metrics (max {MAX_METRICS})",
                                        choices=list(METRIC_REF),
                                        selected=DEFAULT_METRICS, inline=True),
                ui.input_radio_buttons("order_by", "Row order",
                                       choices={"metric": "by 1st metric", "name": "by CpG"},
                                       selected="metric", inline=True),
                ui.input_radio_buttons("dist_source", "Distribution",
                                       choices={"val": "values (as modelled)",
                                                "Z": "Z-scores (test)"},
                                       selected="val", inline=True),
                class_="d-flex flex-wrap gap-4 align-items-start small mt-3 pt-2 border-top",
            ),
        ),

        # ---- per-CpG panels -------------------------------------------------
        ui.nav_panel(
            "Trajectories",
            ui.div(ui.span("Colour by", class_="me-2"),
                   ui.input_select("traj_by", None, choices=[]),
                   class_="d-flex align-items-baseline"),
            ui.layout_columns(ui.output_image("traj_pre", height="auto"),
                              ui.output_image("traj_pos", height="auto")),
        ),

        ui.nav_panel(
            "Train/test performance",
            ui.output_ui("perf_boxes"),
            ui.layout_columns(ui.output_image("cent_train", height="auto"),
                              ui.output_image("cent_test", height="auto")),
            ui.layout_columns(ui.output_image("qq_train", height="auto"),
                              ui.output_image("qq_test", height="auto")),
        ),

        ui.nav_panel(
            "Batch effects",
            ui.layout_columns(ui.output_image("dens_pre", height="auto"),
                              ui.output_image("dens_pos", height="auto")),
            ui.output_image("dens_array_pre", height="auto"),
            ui.output_image("dens_array_pos", height="auto"),
            ui.layout_columns(ui.output_image("dens_cohort_pre", height="auto"),
                              ui.output_image("dens_cohort_pos", height="auto")),
        ),

        ui.nav_panel(
            "Phenotype prediction",
            ui.layout_columns(ui.output_text_verbatim("model_summary_raw"),
                              ui.output_text_verbatim("model_summary_norm")),
        ),

        ui.nav_panel("Info & strategy", ui.markdown(ABOUT)),

        id="tabs", full_screen=True,
    ),
)


# ==== Server ==================================================================

def server(input: Inputs):

    # ---- plot images --------------------------------------------------------
    def _png(model, filename):
        path = RESULTS_DIR / model / "plots" / filename
        return {"src": str(path), "width": "100%"} if path.exists() else None

    def img(filename):
        """Plot from the selected model."""
        req(input.model())
        return _png(input.model(), filename)

    def shared_img(filename):
        """Raw-value plots are identical across models: take the first one found,
        so they still show for models run with saveplots=False."""
        req(input.model())
        return next((p for m in [input.model(), *models] if (p := _png(m, filename))), None)

    # ---- results ------------------------------------------------------------
    @reactive.calc
    def modelstats():
        """{'train': df, 'test': df} of metrics (rows) by CpG (columns)."""
        req(input.model())
        path = RESULTS_DIR / input.model() / "results"
        return {s: (pd.read_csv(f, index_col=0)
                    if (f := path / f"statistics_{s}.csv").exists() else None)
                for s in ("train", "test")}

    def cpg_cols(df):
        return [c for c in df.columns if c in cpgs] or list(df.columns)

    def perf(split, metric):
        df, cpg = modelstats()[split], input.cpg()
        if df is None or metric not in df.index or cpg not in df.columns:
            return None
        return float(df.loc[metric, cpg])

    def metric_avg(split, metric):
        """Metric averaged over every CpG of this model."""
        df = modelstats()[split]
        if df is None or metric not in df.index:
            return None
        v = pd.to_numeric(df.loc[metric, cpg_cols(df)], errors="coerce").to_numpy(float)
        return float(np.nanmean(v)) if np.isfinite(v).any() else None

    # ---- model overview -----------------------------------------------------
    _prev = reactive.value(tuple(DEFAULT_METRICS))

    @reactive.effect
    @reactive.event(input.metrics)
    def _cap_metrics():
        """Keep at most MAX_METRICS, dropping the oldest rather than the new click."""
        cur, prev = list(input.metrics()), list(_prev())
        if len(cur) > MAX_METRICS:
            keep = set(([m for m in cur if m not in prev]
                        + [m for m in prev if m in cur])[:MAX_METRICS])
            cur = [m for m in METRIC_REF if m in keep]
            ui.update_checkbox_group("metrics", selected=cur)
        _prev.set(tuple(cur))

    @reactive.calc
    def row_order():
        """CpGs in plotting order (shared by the figure and its height)."""
        s = modelstats()
        base = s["test"] if s["test"] is not None else s["train"]
        if base is None:
            return []
        present, metrics = cpg_cols(base), list(input.metrics())
        if input.order_by() == "metric" and metrics and metrics[0] in base.index:
            return (pd.to_numeric(base.loc[metrics[0], present], errors="coerce")
                      .sort_values().index.tolist())
        return sorted(present, reverse=True)

    @reactive.calc
    def dist_data():
        """(dataframe of CpG columns, kind) for the distribution strip."""
        req(input.model())
        if input.dist_source() == "Z":
            f = RESULTS_DIR / input.model() / "results" / "Z_test.csv"
            if not f.exists():
                return None, "Z"
            df = pd.read_csv(f)
            return df[[c for c in df.columns if c in cpgs]], "Z"

        f = DATA_DIR / "dnam_data.csv"
        if not f.exists():
            return None, "beta"
        df = pd.read_csv(f, usecols=[c for c in pd.read_csv(f, nrows=0).columns if c in cpgs])
        if "Mval" in input.model():        # mirror the pipeline's beta -> M step
            b = np.clip(df.to_numpy(float), 1e-3, 1 - 1e-3)
            return pd.DataFrame(np.log2(b / (1 - b)), columns=df.columns), "M"
        return df, "beta"

    @render.ui
    def overview_boxes():
        """One box per metric on a row of equal slots; the first slot is empty so
        every box sits over its own panel in the figure below."""
        def box(m):
            return ui.div(
                ui.div(f"mean {m}", class_="small text-muted lh-1"),
                ui.div(fmt(metric_avg("test", m)), class_="fs-5 fw-semibold lh-1 my-1"),
                ui.div("train: ", ui.tags.b(fmt(metric_avg("train", m))),
                       class_="small text-muted lh-1"),
                class_="border rounded px-1 py-1 text-center", style="flex:1 1 0;",
            )

        return ui.div(ui.div(style="flex:1 1 0;"),            # labels + distribution
                      *[box(m) for m in input.metrics()],
                      class_="d-flex gap-2 align-items-stretch")

    @render.ui
    def overview_ui():
        n = max(len(row_order()), 1)
        return ui.output_plot("overview", height=f"{min(2600, max(260, ROW_PX * n + 55))}px")

    def _ridge(ax, df, kind, order, sel):
        """Small density per row, spanning only the range where that CpG has data."""
        lo, hi = np.nanpercentile(df.to_numpy(float), [0.2, 99.8])
        if kind == "Z":
            lo, hi = min(lo, -4.0), max(hi, 4.0)
        elif kind == "beta":
            lo, hi = 0.0, 1.0

        for i, c in enumerate(order):
            if c not in df.columns:
                continue
            v = df[c].to_numpy(float)
            v = v[np.isfinite(v)]
            if v.size < 5 or v.max() <= v.min():
                continue
            dens, edges = np.histogram(v, bins=40, range=(v.min(), v.max()), density=True)
            dens = np.convolve(dens, np.ones(3) / 3, mode="same")     # light smoothing
            dens = dens / dens.max() * 0.8 if dens.max() > 0 else dens
            ax.fill_between(np.r_[v.min(), 0.5 * (edges[:-1] + edges[1:]), v.max()],
                            i - 0.4, i - 0.4 + np.r_[0.0, dens, 0.0], lw=0, zorder=2,
                            color=HILITE if c == sel else TEST_COL,
                            alpha=0.90 if c == sel else 0.55)
        if kind == "Z":
            for xv, ls in ((0.0, "--"), (-1.96, ":"), (1.96, ":")):
                ax.axvline(xv, color="grey", lw=0.7, ls=ls, zorder=1)
        ax.set_xlim(lo, hi)
        ax.set_xlabel(X_LAB[kind], fontsize=7)

    @render.plot
    def overview():
        """Figure split into len(metrics)+1 equal slots, matching the box row:
        slot 0 = CpG labels + distribution, then one slot per metric."""
        order, metrics = row_order(), list(input.metrics())
        if not order or not metrics:
            return None
        stats, sel = modelstats(), input.cpg()
        y = np.arange(len(order))
        slot = 1 / (len(metrics) + 1)

        fig = plt.figure(figsize=(2.0 * (len(metrics) + 1),
                                  max(2.0, 0.145 * len(order)) + 0.45))
        # Shiny needs every axes to be a subplot (its coordmap reads the subplotspec),
        # and it forces a tight layout on any figure without a layout engine -- which
        # would undo the positions below. A no-op engine keeps our geometry exact.
        fig.set_layout_engine(PlaceHolderLayoutEngine(adjust_compatible=True,
                                                      colorbar_gridspec=False))
        gs = fig.add_gridspec(1, len(metrics) + 1)
        dax = fig.add_subplot(gs[0, 0])
        maxes = [fig.add_subplot(gs[0, j + 1], sharey=dax) for j in range(len(metrics))]

        lab = max(len(c) for c in order) * LABEL_FS * 0.58 / 72 / fig.get_figwidth() + 0.012
        bottom, height = 0.075, 0.915
        dax.set_position([lab, bottom, max(0.03, slot - lab - GAP), height])
        for j, ax in enumerate(maxes):
            ax.set_position([slot * (j + 1) + GAP / 2, bottom, slot - GAP, height])

        df, kind = dist_data()
        if df is not None and not df.empty:
            _ridge(dax, df, kind, order, sel)
        else:
            dax.set_xticks([])

        k = order.index(sel) if sel in order else None
        for ax in [dax, *maxes]:
            stripe_rows(ax, len(order))

        for ax, m in zip(maxes, metrics):
            if k is not None:
                ax.axhspan(k - 0.5, k + 0.5, color="#ffe08a", alpha=0.45, zorder=0.4)
            if METRIC_REF.get(m) is not None:
                ax.axvline(METRIC_REF[m], color="grey", lw=0.8, ls="--", zorder=1)
            for d_, split, col in ((stats["train"], "train", TRAIN_COL),
                                   (stats["test"], "test", TEST_COL)):
                if d_ is None or m not in d_.index:
                    continue
                v = [float(d_.loc[m, c]) if c in d_.columns else np.nan for c in order]
                ax.scatter(v, y, s=16, label=split, zorder=3, linewidths=0.8,
                           edgecolors=col, facecolors="none" if split == "train" else col)
            ax.grid(axis="x", which="major", lw=0.5, alpha=0.5, color="#b5b5b5")
            ax.set_axisbelow(True)          # grid under the dots, over the stripes
            ax.tick_params(labelsize=6, labelleft=False)

        if k is not None:                   # no vertical grid on the distribution
            dax.axhspan(k - 0.5, k + 0.5, color="#ffe08a", alpha=0.45, zorder=0.4)
        dax.set_yticks(y)
        dax.set_yticklabels(order, fontsize=LABEL_FS)
        dax.set_ylim(-0.8, len(order) - 0.2)
        dax.tick_params(axis="x", labelsize=6)
        if k is not None:
            dax.get_yticklabels()[k].set(fontweight="bold", color=HILITE)
        maxes[-1].legend(loc="lower right", fontsize=7, frameon=False)
        return fig

    # ---- trajectories -------------------------------------------------------
    @reactive.effect
    def _traj_choices():
        choices = ["Array", "Cohort"] + cpg_pheno[input.cpg()]
        ui.update_select("traj_by", choices=choices, selected=choices[0])

    @render.image
    def traj_pre():
        req(input.traj_by())
        return shared_img(f"traj_{input.cpg()}_by_{input.traj_by()}.png")

    @render.image
    def traj_pos():
        req(input.traj_by())
        return img(f"traj_{input.cpg()}Z_by_{input.traj_by()}.png")

    # ---- train / test performance -------------------------------------------
    @render.image
    def cent_train():
        return img(f"centiles_{input.cpg()}_train_harmonized.png")

    @render.image
    def cent_test():
        return img(f"centiles_{input.cpg()}_test_harmonized.png")

    @render.image
    def qq_train():
        return img(f"qq_{input.cpg()}_train.png")

    @render.image
    def qq_test():
        return img(f"qq_{input.cpg()}_test.png")

    @render.ui
    def perf_boxes():
        def box(m):
            return ui.div(
                ui.div(m, class_="small text-muted"),
                ui.div(fmt(perf("test", m)), class_="fs-4 fw-semibold lh-1 my-1"),
                ui.div("train: ", ui.tags.b(fmt(perf("train", m))), class_="small text-muted"),
                class_="border rounded p-2 text-center",
            )

        return ui.layout_columns(*[box(m) for m in DEFAULT_METRICS], fill=False, gap="0.5rem")

    # ---- batch effects ------------------------------------------------------
    @render.image
    def dens_pre():
        return shared_img(f"dens_{input.cpg()}_by_period.png")

    @render.image
    def dens_pos():
        return img(f"dens_{input.cpg()}Z_by_period.png")

    @render.image
    def dens_array_pre():
        return shared_img(f"dens_{input.cpg()}_by_array.png")

    @render.image
    def dens_array_pos():
        return img(f"dens_{input.cpg()}Z_by_array.png")

    @render.image
    def dens_cohort_pre():
        return shared_img(f"dens_{input.cpg()}_by_cohort.png")

    @render.image
    def dens_cohort_pos():
        return img(f"dens_{input.cpg()}Z_by_cohort.png")

    # ---- phenotype prediction -----------------------------------------------
    def summary_text(kind):
        summaries = json.loads((INPUT_DIR / "model_summary.json")
                               .read_text(encoding="utf-8")).get(input.cpg())
        if not summaries:
            return "No model was run"
        blocks = []
        for pheno in cpg_pheno[input.cpg()]:
            txt = summaries.get(pheno, {}).get(kind)
            if txt:
                blocks.append(f"=== {pheno} ===\n"
                              + ("\n".join(txt) if isinstance(txt, list) else txt))
        return "\n\n".join(blocks) or "No model was run"

    @render.code
    def model_summary_raw():
        return summary_text("raw")

    @render.code
    def model_summary_norm():
        return summary_text("norm")


app = App(app_ui, server)
