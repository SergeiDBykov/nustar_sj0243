
Nu_path = "/Users/sdbykov/work/xray_pulsars/nustar_sj0243/" #path to the repository
analysis_path = Nu_path+"nustar_products/results/" #path to the analysis results
plot_path = Nu_path+'data_analysis/final_results/figures_tables/'  #path to the plots
obj_name = "Swift J0243.6+6124" #name of the pulsar 

#list of observations
ObsList = [
    "90302319002",
    "90302319004",
    "90302319006",
    "90302319008",
    "90401308002",
    "90401334002",
    "90501310002",
]  

#list of observations with high count rate
ObsList_bright = [
    "90302319002",
    "90302319004",
    "90302319006",
    "90302319008",
    "90401334002",
]



ObsAlias = {
    "90302319002": 'I (90302319002)',
    "90302319004": 'II (90302319004)',
    "90302319006": 'III (90302319006)',
    "90302319008": 'IV (90302319008)',
    '90401334002': 'V (90401334002)'

}

#dictionary of pulasr  periods
periods_val = {
    "90302319002": 9.85425,
    "90302319004": 9.84435,
    "90302319006": 9.8234,
    "90302319008": 9.801050,
    "90401334002": 9.7918,
}

#dictionary of dates  of observations
mjd_val = {
    "90302319002": 5.803166871740724e04,
    "90302319004": 5.805731338465247e04,
    "90302319006": 5.806712464344907e04,
    "90302319008": 5.809362565317134e04,
    "90401308002": 5.818752280016200e04,
    "90401334002": 5.836886549981482e04,
    "90501310002": 5.855725539275459e04,
    # '90401308001': 5.818751398361111E+04,  # 200 sec expo, ignore!
}
MJD_REF = 58000 #MJD reference
