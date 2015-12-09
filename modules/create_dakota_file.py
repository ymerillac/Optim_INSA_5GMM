import sys

if len(sys.argv)!=5:
    print 'Erreur sur le nombre de parametres.\nUsage : python script_create_dakota_file.py [NOM DU FICHIER] [METHODE] [NOMBRE DE VARIABLES] [NOM DE LA FONCTION]'
    print 'Exemple d\'utilisation : python create_dakota_file.py minimisation optpp_q_newton 6 myfunc.py'
    raise TypeError('Nb arguments')

fid = open(sys.argv[1]+'.in','w')
fid.write('environment,\n    graphics\n    tabular_graphics_data\n    tabular_graphics_file = \'optim_res.dat\'\n\n')
fid.write('method,\n    '+sys.argv[2]+'\n')
fid.write('    max_iterations = 50\n    convergence_tolerance = 1e-4\n\n')
fid.write('variables,\n    continuous_design = '+sys.argv[3]+'\n')
fid.write('    initial_point\n    lower_bound\n    upper_bound\n    descriptors \n\n')
fid.write('interface,\n    system\n    analysis_driver = \'python '+sys.argv[4]+'\n\n')
fid.write('responses,\n    objective_functions = 1\n    descriptors = \'my_function\'\n    numerical_gradients\n    #no_gradients\n    no_hessians\n    #numerical_hessians')
fid.close()
