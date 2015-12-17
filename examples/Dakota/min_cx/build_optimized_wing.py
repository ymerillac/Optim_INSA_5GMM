import mesh_utils

file = open('optim_res.dat','r')
tmp = file.readlines()
tmp = tmp[-1]
tmp = tmp.split()
# attention : a changer quand il y aura plus de variables de design. Dernier argument : la corde
root = [float(tmp[2]),float(tmp[4]),float(tmp[6]),1.]
tip = [float(tmp[3]),float(tmp[5]),float(tmp[7]),1.]
mesh_utils.create_mesh_linear_interp('optimized_wing',root,tip,10.,31,50)
file.close()
