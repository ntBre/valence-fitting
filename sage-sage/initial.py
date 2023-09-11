from openff.toolkit import ForceField

sage = ForceField("openff-2.1.0.offxml")
sage.to_file("initial.offxml")
