#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ./runSensitivity3dSimulations.py > ./output/Sensitivity3d.txt

import os
import subprocess
import time

materials = ["LYSO-Proteus", "BGO", "BGO_EJ232_15", "BGO_BaF2_15", "LYSO-Proteus_EJ232_25", "LYSO-Proteus_BaF2_25"]

for material in materials:
    print("Crystal material is", material)
    for z in [0, 15, 45, 75]:
        print("axial_pos:", z)
        for x in [0, 90, 180, 270]:
            print("x_pos:", x)
            output_dir = f"output/Sensitivity3d/{material}/{x};{z}"
            os.makedirs(output_dir, exist_ok=True)
            subprocess.Popen(["Gate", "-a", f"[crystal,{material}][x_pos,{x}][axial_pos,{z}]", "mainMacro-Sensitivity3d.mac"])
            time.sleep(3000)  # sleep for 30 minutes (30m = 30 * 60 seconds)
