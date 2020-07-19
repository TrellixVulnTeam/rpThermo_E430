from equilibrator_api import ComponentContribution, Q_

class rpThermo:
    def __init__(self, ph=7.0, ionic_strength=200, pMg=10.0, temp_k=298.15):
        self.cc = ComponentContribution()
        self.cc.p_h = Q_(ph)
        self.cc.ionic_strength = Q_(str(ionic_strength)+' mM')
        self.cc.p_mg = Q_(pMg)
        self.cc.temperature = Q_(str(temp_k)+' K')

    def reaction_thermo(self, reac):
        for rea in reac.getListOfReactants():
        for pro in reac.getListOfProducts():

        #"kegg:C00002 + kegg:C00001 = kegg:C00008 + kegg:C00009"
        reac_string = ''
        rxn = self.cc.parse_reaction_formula(reac_string)

    def pathway_thermo(self, pathway_id='rp_pathway'):
