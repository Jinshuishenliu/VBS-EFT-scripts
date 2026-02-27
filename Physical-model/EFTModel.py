from HiggsAnalysis.CombinedLimit.PhysicsModel import PhysicsModel

class eftModel(PhysicsModel):
    def setPhysicsOptions(self, physOptions):
        # Quadratic polynomial coefficients
        #self.coeffs = [2.6241595039, 0.0011285359, 0.2013591441] #2016
        #self.coeffs = [3.1500031949, 0.0001023059, 0.1632468855] #2016APV
        #self.coeffs = [6.2664532614, 0.0000676736, 0.4441337301] #2017
        self.coeffs = [10.6526566934, 0.0032300313, 0.8506797757] #run2
        
    def doParametersOfInterest(self):
        # Define POI
        self.modelBuilder.doVar("FT8[0,-10,10]")
        self.modelBuilder.doSet("POI", "FT8")

        # Define constant coefficients
        for i, c in enumerate(self.coeffs):
            self.modelBuilder.doVar(f"c{i}[{c}]")
            self.modelBuilder.out.var(f"c{i}").setConstant(True)

        # Define yield expression (quadratic polynomial)
        #norm = 14.612
        self.modelBuilder.factory_(
            "expr::FT8_scaling('(@0 + @1*@3 + @2*pow(@3,2))/@0', c0, c1, c2, FT8)"
        )

    def getYieldScale(self, bin, process):
        # Apply the polynomial scaling only to the EFT signal
        if process == "eft_ZZ2l2nu":
            return "FT8_scaling"
        else: return 1

eftModel = eftModel()
