from pathlib import Path

from shiny import App, Inputs, reactive, render, ui
import json


INPUT_DIR = Path(__file__).parent.parent

with open(INPUT_DIR / "cpg_targets.json", "r") as file:
    cpg_targets = json.load(file)

cpgs = list(cpg_targets.keys())

cpg_labs = {cpg: f"{cpg} - {v['label']}" for cpg, v in cpg_targets.items()}

app_ui = ui.page_navbar(
    ui.nav_spacer(),
    
    ui.nav_panel(
        "Data view",
        ui.navset_card_underline(
            ui.nav_spacer(),
            ui.nav_panel("Speghetti", 
              ui.layout_columns(
                ui.output_image("traj_pre"), ui.output_image("traj_pos"))),
            ui.nav_panel("Train/test centiles", 
              ui.layout_columns(
                ui.output_image("cent_train"), ui.output_image("cent_test"))),
            title="Trajectories", full_screen=True, #  height="500px"
            ), 
        ui.card(
            ui.card_header("Density"),
            ui.row(ui.layout_columns(
              ui.output_image("dens_pre"), ui.output_image("dens_pos"))),
            ui.row(ui.layout_columns(
              ui.output_image("dens_array_pre"), ui.output_image("dens_array_pos"))),
          ),
      ),

      ui.nav_panel("Regression performance",
        ui.layout_columns(
          ui.value_box(title="Row count", value=ui.output_text("row_count")),
          ui.value_box(title="Mean training score", value=ui.output_text("mean_score")), fill=False,
        ),  ui.layout_columns(
          ui.output_text_verbatim("model_summary_raw"),  
          ui.output_text_verbatim("model_summary_norm")) #ui.card(ui.output_data_frame("data")),
      ),
    
      sidebar=ui.sidebar(
          ui.input_select("cpg", "CpG", choices=cpg_labs)),
      id="tabs",
      title="MethylCHART demo",
      fillable=True,
)


def server(input: Inputs):
    # @reactive.calc()
    # def dat() -> pl.DataFrame:
    #     return scores.filter(pl.col("account") == input.account())
    
    @render.image
    def traj_pre():
        img: ImgData = {"src": str(INPUT_DIR / 'results' / 'plots' / f'{input.cpg()}_trajectories.png'), 
          "width": "100%"}
        return img
      
    @render.image
    def cent_train():
        img: ImgData = {"src": str(INPUT_DIR / 'results' / 'warpedBLR' / 'plots' / f'centiles_{input.cpg()}_cpgdata_train_harmonized.png'), 
          "width": "100%"}
        return img
      
    @render.image
    def cent_test():
        img: ImgData = {"src": str(INPUT_DIR / 'results' / 'warpedBLR' / 'plots' / f'centiles_{input.cpg()}_cpgdata_test_harmonized.png'),
          "width": "100%"}
        return img
      
    @render.image
    def traj_pos():
        img: ImgData = {"src": str(INPUT_DIR / 'results' / 'plots' / f'{input.cpg()} Z_trajectories.png'), 
          "width": "100%"}
        return img
      
    @render.image
    def dens_pre():
        img: ImgData = {"src": str(INPUT_DIR / 'results' / 'plots' / f'{input.cpg()}_dens_by_period.png'), 
          "width": "100%"}
        return img
      
    @render.image
    def dens_pos():
        img: ImgData = {"src": str(INPUT_DIR / 'results' / 'plots' / f'{input.cpg()} Z_dens_by_period.png'), 
          "width": "100%"}
        return img

    @render.image
    def dens_array_pre():
        img: ImgData = {"src": str(INPUT_DIR / 'results' / 'plots' / f'{input.cpg()}_dens_by_array.png'), 
          "width": "100%"}
        return img
      
    @render.image
    def dens_array_pos():
        img: ImgData = {"src": str(INPUT_DIR / 'results' / 'plots' / f'{input.cpg()} Z_dens_by_array.png'), 
          "width": "100%"}
        return img
      
    @render.text
    def row_count():
        return 10

    @render.text
    def mean_score():
        return 10
      
    @render.text
    def model_summary_raw():
        with open(INPUT_DIR / "model_summary.json", "r") as file:
          model_summaries = json.load(file)
        
        if input.cpg() in list(model_summaries.keys()):
          
          phenos = cpg_targets[input.cpg()]['pheno']
          
          if isinstance(phenos, str):
            phenos = [phenos]
    
          for pheno in phenos:
            summ_text = model_summaries[input.cpg()][pheno]['raw']
          
          return "\n".join(summ_text)
        else:
          return 'No model was run'
        
    @render.text
    def model_summary_norm():
        with open(INPUT_DIR / "model_summary.json", "r") as file:
          model_summaries = json.load(file)
        
        if input.cpg() in list(model_summaries.keys()):
          
          phenos = cpg_targets[input.cpg()]['pheno']
          
          if isinstance(phenos, str):
            phenos = [phenos]
          
          for pheno in phenos:
            summ_text = model_summaries[input.cpg()][pheno]['norm']
          
          return "\n".join(summ_text)
        else:
          return 'No model was run'
        
    # @render.data_frame
    # def data():
    #     return dat()


app = App(app_ui, server)
