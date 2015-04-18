# hacer una esfera con los cubitos
import bpy

scene = bpy.context.scene

obj1 =  scene.objects['Cube']
obj1.select = True
n = 15
radioin = 2.7
radio = 3
delta = 0.2

for i in range(n):
    for j in range(n):
        for k in range(n):
            x = delta * i
            y = delta * j
            z = delta * k
            r = sqrt(x**2+y**2+z**2)
            obj1.select = True
            if ((r > radioin) and (r <= radio)):
                bpy.ops.object.duplicate_move(\
                    OBJECT_OT_duplicate={"linked":False, "mode":'TRANSLATION'},\
                    TRANSFORM_OT_translate={\
                        "value":(x,y,z), \
                        "constraint_axis":(False, False, False), \
                        "constraint_orientation":'GLOBAL', \
                        "mirror":False, \
                        "proportional":'DISABLED', \
                        "proportional_edit_falloff":'SMOOTH', \
                        "proportional_size":1, \
                        "snap":False, \
                        "snap_target":'CLOSEST', \
                        "snap_point":(0, 0, 0), \
                        "snap_align":False, \
                        "snap_normal":(0, 0, 0), \
                        "texture_space":False, \
                        "remove_on_cancel":False, \
                        "release_confirm":False}\
                )
            bpy.ops.object.select_all(action='DESELECT')