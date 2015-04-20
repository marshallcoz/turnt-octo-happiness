import bpy
import bmesh
import sys,getopt
import os
from math import pi
bpy.ops.object.mode_set(mode='OBJECT')
lista = sys.argv
print("Argument List:"+ str(lista))
print(str(lista[len(sys.argv)-1]))
profZ = float(lista[len(sys.argv)-1])
print("hello, imprimir vetices del poligono en " + str(profZ))

#bpy.context.tool_settings.mesh_select_mode = [True, False, False]
scene = bpy.context.scene
scene.layers = [True] * 20 # Show all layers

G = "8.4"
contador = 0
currfi = bpy.data.filepath +"z"+str(profZ) + ".txt"
outfile = open(currfi,'w')
outfile.write("---- "+bpy.data.filepath+" vetices del poligono en " + str(profZ)+"---\n")
#obj = scene.objects['Icosphere']
for obj in scene.objects:
    if obj.type == 'MESH':
        bpy.ops.object.mode_set(mode='OBJECT')
        # Dividier a la profundidad profZ
        bpy.ops.mesh.primitive_cube_add(location=(0.0,0.0,1000-profZ))
        ob = bpy.context.object
        ob.name = 'Borra'
        ob.scale=((1000,1000,1000))
        boo = obj.modifiers.new('Booh', 'BOOLEAN')
        boo.object = ob
        boo.operation = 'DIFFERENCE'
        bpy.ops.object.select_all(action="TOGGLE")
        obj.select = True
        bpy.context.scene.objects.active = obj
        bpy.ops.object.modifier_apply(apply_as='DATA', modifier="Booh")
        bpy.ops.object.select_all(action="TOGGLE")
        ob.select = True
        bpy.context.scene.objects.active = ob
        bpy.ops.object.delete()
        #scene.objects.unlink(ob)
        obj.select = True
        bpy.context.scene.objects.active = obj
        bpy.ops.object.mode_set(mode='EDIT')
        bm = bmesh.from_edit_mesh(obj.data)
        # Presentar los resultados
        outfile.write(str(len(bm.verts))+"\n")
        for f in bm.verts:
            cen = f.co
            if cen.z == profZ:
            	if (abs(cen.x) != 1000) and (abs(cen.y) != 1000):
            		print(format(cen.x,G)+" "+format(cen.y,G))
            		outfile.write(format(cen.x,G)+" "+format(cen.y,G)+"\n")
            		contador = contador + 1
        outfile.close()
print("There are "+ str(contador) + " points")

#import bpy
#import bmesh
#scene = bpy.context.scene
#obj = scene.objects['Icosphere']
#profZ = 0
#bpy.ops.mesh.primitive_cube_add(location=(0.0,0.0,1000-profZ))
#ob = bpy.context.object
#ob.name = 'Borra'
#ob.scale=((1000,1000,1000))
#boo = obj.modifiers.new('Booh', 'BOOLEAN')
#boo.object = ob
#boo.operation = 'DIFFERENCE'
#bpy.ops.object.select_all(action="TOGGLE")
#obj.select = True
#bpy.context.scene.objects.active = obj
#bpy.ops.object.modifier_apply(apply_as='DATA', modifier="Booh")
#bpy.ops.object.select_all(action="TOGGLE")
#ob.select = True
#bpy.context.scene.objects.active = ob
#bpy.ops.object.delete()
