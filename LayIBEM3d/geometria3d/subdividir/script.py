import bpy
import bmesh
import sys,getopt
import os
from math import pi
bpy.ops.object.mode_set(mode='EDIT')
print("hello")

lista = sys.argv
print("Argument List:"+ str(lista))
print(str(lista[len(sys.argv)-1]))

bpy.context.tool_settings.mesh_select_mode = [False, True, False]
scene = bpy.context.scene
scene.layers = [True] * 20 # Show all layers

radioMax = float(lista[len(sys.argv)-1])
print("Radio maximo = "+str(radioMax))

diamMax = radioMax * 2.
factor = 1
maxfactor = 5
contador = 0
keepgoing = True
iter = 0
maxiter = 5
G = "8.4"
minicont = 0
for obj in scene.objects:
    if obj.type == 'MESH':
        bm = bmesh.from_edit_mesh(obj.data)
        while (keepgoing == True):
            maxlen = 0.0
            bpy.ops.mesh.select_all(action = 'DESELECT')
            for f in bm.faces:
                for ed in f.edges:
                    length = ed.calc_length()
                    factor = 1
                    if (length > diamMax):
                        ed.select = True
                        factor = max(factor,int(length/diamMax)+1)
                        factor = min(factor,maxfactor)
            bpy.ops.mesh.subdivide(factor)
            for f in bm.faces:
                minicont = 0
                for ed in f.edges:
                    maxlen = max(maxlen,ed.calc_length())
                    minicont = minicont + 1
                if (minicont >= 5):
                    bpy.ops.mesh.select_all(action = 'DESELECT')
                    f.select = True
                    bpy.ops.mesh.quads_convert_to_tris()
            if (maxlen < radioMax):
                keepgoing = False
            iter = iter + 1
            if (iter >= maxiter):
                keepgoing = False 
        print("---------------------------------")
        currfi = bpy.data.filepath + str(radioMax) + ".txt"
        outfile = open(currfi,'w')
        outfile.write("----"+bpy.data.filepath+str(radioMax)+"----\n")
        outfile.write(str(len(bm.faces))+"\n")
        # Presentar los resultados
        for f in bm.faces:
            cen = f.calc_center_bounds()
            nor = f.normal
            are = f.calc_area()
            rad = (are / pi ) ** (0.5)
            outfile.write(format(cen.x,G)+" "+format(cen.y,G)+" "+format(cen.z,G)+" "+
            	format(nor.x,G)+" "+format(nor.y,G)+" "+format(nor.z,G)+" "+format(rad,G)+"\n")
            outfile.write(str(len(f.verts))+"\n")
            for vf in f.verts:
            	cen = vf.co
            	outfile.write(format(cen.x,G)+" "+format(cen.y,G)+" "+format(cen.z,G)+"\n")
            #print("center: " + str(f.calc_center_bounds()))
            #print("normal: " + str(f.normal))
            #print("area: " + str(f.calc_area()))
            contador = contador + 1
        outfile.close()
        print("There are "+ str(contador) + " colocation points")
#bpy.ops.wm.save_mainfile()
filepath = bpy.data.filepath
directory = os.path.dirname(filepath)
bpy.ops.wm.save_as_mainfile(filepath=directory+"/radio"+str(radioMax)+"m.blend")

