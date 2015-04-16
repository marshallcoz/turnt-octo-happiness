import bpy
import bmesh
import sys,getopt
import os
from math import pi
bpy.ops.object.mode_set(mode='EDIT')
print("hello, imprimir vetices del poligono en Z=0")

bpy.context.tool_settings.mesh_select_mode = [True, False, False]
scene = bpy.context.scene
scene.layers = [True] * 20 # Show all layers

G = "8.4"
for obj in scene.objects:
    if obj.type == 'MESH':
        bm = bmesh.from_edit_mesh(obj.data)
        print("------PoligonoEnlaSuperficie--(vértices)-----")
        # Presentar los resultados
        for f in bm.verts:
            cen = f.co
            if cen.z == 0.0:
            	print(format(cen.x,G)+" "+format(cen.y,G))

