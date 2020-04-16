#!/home/lucia3/anaconda/bin/python
# 14.8.14
# simple version, only names and folder change, .out files and Nphot constant
# reads in all model names in a folder
# creates MCMax commands to start the RT in one file, .run
# creates MCMax commands to start the observations in another file, .ray

# path and names of .run .ray files
path = '/home/klarmann/Desktop/MCMax_examples/contineous_evaporation_dust_composition/' #always space as first character!
name = 'evap_models'


#Read in neccesary strings
#image name can be the name of any observation file, not necessarily an IMAGE
common_path = ' /home/klarmann/Desktop/MCMax_examples/contineous_evaporation_dust_composition/'
image_name = 'Image_asy_noscat'
image_name2 = 'Image_asy_scat'
image_name3 = 'Interferometry4'
image_name4 = 'InterferometryN'
photon_number = ' 5e5'
photon_number_raytrace = ' 0'
output_folder = 'Output'
program_name = 'MCMax17'
flag = ' -o'

model_name = [ 'a1e-1_la1e2_r0.0001_L40', 'f_1e0_r0.01_L8']


number_of_models = len(model_name)
#print (number_of_models)

# make sure there are nor trailing lines
line_ending = ['\n'] * (number_of_models - 1)
line_ending.append('')
#print (line_ending)


# open and write file for runs
f_run = open(path + name+'.run', 'w')

# print list of commands for run
for i in range(number_of_models):
    f_run.write(program_name + common_path + model_name[i] +
                '/'+model_name[i]+'.in' + photon_number + flag +
                common_path + model_name[i] + '/'+output_folder + line_ending[i])

f_run.close()

#open and write file for raytrace
f_ray = open(path + name+'.ray', 'w')
print (path + name)

# print list of commands for raytrace
for i in range(number_of_models):
    f_ray.write(program_name + common_path + model_name[i] + '/'+model_name[i]+'.in' +
                photon_number_raytrace + flag + common_path + model_name[i] + '/'+output_folder +
                common_path + model_name[i] + '/'+image_name+'.out' +  # comment line if only spectrum
                common_path + model_name[i] + '/'+image_name2+'.out' +  # comment line if only spectrum
                common_path + model_name[i] + '/'+image_name3+'.out' +  # comment line if only spectrum
                common_path + model_name[i] + '/'+image_name4+'.out' +  # comment line if only spectrum
                line_ending[i])

f_ray.close()

print ('There will be %i models calculated and raytraced' % len(model_name))
print (model_name)
