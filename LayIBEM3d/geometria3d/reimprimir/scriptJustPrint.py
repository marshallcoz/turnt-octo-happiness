import bpy
import bmesh
import sys,getopt
import os
from math import pi
bpy.ops.object.mode_set(mode='EDIT')
print("hello")


bpy.context.tool_settings.mesh_select_mode = [False, True, False]
scene = bpy.context.scene
scene.layers = [True] * 20 # Show all layers

contador = 0
maxiter = 5
G = "8.4"
for obj in scene.objects:
    if obj.type == 'MESH':
        bm = bmesh.from_edit_mesh(obj.data)
        print("---------------------------------")
        print(" ")
        print(" ")
        # Presentar los resultados
        for f in bm.faces:
            cen = f.calc_center_bounds()
            nor = f.normal
            are = f.calc_area()
            rad = (are / pi ) ** (0.5)
            print(format(cen.x,G)+" "+format(cen.y,G)+" "+format(cen.z,G)+" "+format(nor.x,G)+" "+format(nor.y,G)+" "+format(nor.z,G)+" "+format(rad,G))
            #print("center: " + str(f.calc_center_bounds()))
            #print("normal: " + str(f.normal))
            #print("area: " + str(f.calc_area()))
            contador = contador + 1
        print(" ")
        print("There are "+ str(contador) + " colocation points")

