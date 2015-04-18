# hacer una esfera con los cubitos
import bpy

scene = bpy.context.scene
radio = 6.0
radioin = 5.0

for obj in scene.objects:
    if obj.type == 'MESH':
        x = obj.location[0]
        y = obj.location[1]
        z = obj.location[2]
        r = sqrt(x**2+y**2+z**2)
        if (r > radio):
            obj.select = True
            bpy.ops.object.delete(use_global=False)
        if (r < radioin):
            obj.select = True
            bpy.ops.object.delete(use_global=False)
    bpy.ops.object.select_all(action='DESELECT')

