"""
This script is used to generate k files from the asset model.

Run this model for 1 step to create d3plot
d3plot is used as reference to test Post-process module

wye
"""
import ansys.heart.preprocessor.models as models
import ansys.heart.writer.dynawriter as writers

model = models.HeartModel.load_model("heart_model.pickle")
writer = writers.MechanicsDynaWriter(model, "ConstantPreloadWindkesselAfterload")
writer.update()
writer.export("post")
