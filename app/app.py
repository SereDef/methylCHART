from pathlib import Path

from shiny import App, Inputs, reactive, render, req, ui
import json
import pandas as pd

# ==== Set-up ==================================================================
INPUT_DIR = Path(__file__).parent.parent # methylCHART

with open(INPUT_DIR / "cpg_targets.json", "r") as file:
    cpg_targets = json.load(file)

cpgs = list(cpg_targets.keys())
cpg_labs = {cpg: f"{cpg} - {v['label']}" for cpg, v in cpg_targets.items()}

cpg_pheno = {cpg: [v['pheno']] if isinstance(v["pheno"], str) else list(v["pheno"] or []) 
             for cpg, v in cpg_targets.items()}

RESULTS_DIR = INPUT_DIR / "results"

models = sorted(d.name for d in RESULTS_DIR.iterdir() if d.is_dir())

# ==== UI ======================================================================
app_ui = ui.page_navbar(
    ui.nav_spacer(),
    
    ui.nav_panel("Data view",
        ui.navset_card_underline(
            ui.nav_spacer(),
            
            ui.nav_panel("Speghetti", 
              ui.div(ui.span("Colour by", class_="me-2"), 
                    ui.input_select("traj_by", None, choices=[]), class_="d-flex align-items-baseline"),
              ui.row(ui.layout_columns(
                ui.output_image("traj_pre"), ui.output_image("traj_pos")))),
                
            ui.nav_panel("Train/test performance",
              ui.output_ui("perf_boxes"),
              ui.layout_columns(
                ui.output_image("cent_train", height='auto'), ui.output_image("cent_test", height='auto')),
              ui.layout_columns(
                ui.output_image("qq_train", height='auto'), ui.output_image("qq_test", height='auto'))
              ),
            title="Trajectories", full_screen=True, #  height="500px"
        ), 
            
        ui.card(
            ui.card_header("Batch effects"),
            ui.layout_columns(
              ui.output_image("dens_pre", height='auto'), ui.output_image("dens_pos", height='auto')),
            ui.output_image("dens_array_pre", height='auto'), 
            ui.output_image("dens_array_pos", height='auto'),
            ui.layout_columns(
              ui.output_image("dens_cohort_pre", height='auto'), ui.output_image("dens_cohort_pos", height='auto'))
          ),
      ),

      ui.nav_panel("Phenotype prediction performance",
        ui.layout_columns(
          ui.output_text_verbatim("model_summary_raw"),  
          ui.output_text_verbatim("model_summary_norm")) #ui.card(ui.output_data_frame("data")),
      ),
    
      sidebar = ui.sidebar(
          ui.input_select("model", "Model", choices=models),
          ui.input_select("cpg", "CpG", choices=cpg_labs)),
      id="tabs",
      title="MethylCHART demo",
      fillable=False,
)


# ==== Server logic ============================================================

def server(input: Inputs):

    def img(filename: str):
        """Image from the selected model's plot folder (nothing if it doesn't exist)."""
        req(input.model())
        path = RESULTS_DIR / input.model() / "plots" / filename
        
        if not path.exists():
          print("MISSING:", path)
          return None
        
        return {"src": str(path), "width": "100%"}
    
    @reactive.calc
    def modelstats():
        """Performance metrics for the selected model, read once per model."""
        req(input.model())
        path = RESULTS_DIR / input.model() / "results"
        out = {}
        for split in ("train", "test"):
            f = path / f"statistics_{split}.csv"
            out[split] = pd.read_csv(f, index_col=0) if f.exists() else None
        return out
    
    def perf(split, metric):
      df, cpg = modelstats()[split], input.cpg()
      if df is None or metric not in df.index or cpg not in df.columns:
          return None
      return float(df.loc[metric, cpg])
    
    # === Trajectories =========================================================
    # Colour options: phenotypes for this CpG (if exists), then Cohort / Array.
    @reactive.effect
    def _():
        choices = cpg_pheno[input.cpg()] + ["Cohort", "Array"]
        ui.update_select("traj_by", choices=choices, selected=choices[0])

    @render.image
    def traj_pre():
        req(input.traj_by())
        return img(f"traj_{input.cpg()}_by_{input.traj_by()}.png")
      
    @render.image
    def traj_pos():
        req(input.traj_by())
        return img(f"traj_{input.cpg()}Z_by_{input.traj_by()}.png")
    
    # === Train / test performance =============================================
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
        fmt = lambda x: "–" if x is None else f"{x:.2f}"
    
        def box(metric):
            return ui.div(
                ui.div(metric, class_="small text-muted"),
                ui.div(fmt(perf("test", metric)), class_="fs-4 fw-semibold lh-1 my-1"),
                ui.div("train: ", ui.tags.b(fmt(perf("train", metric))), class_="small text-muted"),
                class_="border rounded p-2 text-center",
            )
    
        return ui.layout_columns(*[box(m) for m in ["EXPV", "RMSE", "SMSE", "MSLL", "MACE"]],
                                 fill=False, gap="0.5rem")
        
    # === Batch effects ========================================================
    @render.image
    def dens_pre():
        return img(f"dens_{input.cpg()}_by_period.png")
 
    @render.image
    def dens_pos():
        return img(f"dens_{input.cpg()}Z_by_period.png")
 
    @render.image
    def dens_array_pre():
        return img(f"dens_{input.cpg()}_by_array.png")
 
    @render.image
    def dens_array_pos():
        return img(f"dens_{input.cpg()}Z_by_array.png")
      
    @render.image
    def dens_cohort_pre():
        return img(f"dens_{input.cpg()}_by_cohort.png")
 
    @render.image
    def dens_cohort_pos():
        return img(f"dens_{input.cpg()}Z_by_cohort.png")
    
    # === Pehnotype prediction =================================================
    
    def summary_text(kind: str) -> str:
        with open(INPUT_DIR / "model_summary.json") as file:
            summaries = json.load(file).get(input.cpg())
        if not summaries:
            return "No model was run"
        blocks = []
        for pheno in cpg_pheno[input.cpg()]:
            txt = summaries.get(pheno, {}).get(kind)
            if txt:
                blocks.append(f"=== {pheno} ===\n" + ("\n".join(txt) if isinstance(txt, list) else txt))
        return "\n\n".join(blocks) or "No model was run"
      
    @render.code
    def model_summary_raw():
        return summary_text("raw")
 
    @render.code
    def model_summary_norm():
        return summary_text("norm")
        
    # @render.data_frame
    # def data():
    #     return dat()


app = App(app_ui, server)
