## Script to launch ensemble of simulations for one parameter


# one variable, more parameters

fldr = "/home/marius/Models/Sicopolis/runs/headers"

file_in = "sico_specs_Run_CALV_UW.h"


params = ["CALV_UW_COEFF"]

values = ["5.0d-01","5.0d-05","5.0d-04","5.0d-03","5.0d-02","5.0d-06"]


for i in range(len(values)):

    #print(i)

    file_out = "sico_specs_Run_CALV_UW_" + str(i+1) + ".h"

    with open(file_in, 'r') as fin, open(file_out, 'w') as fout:

        for line in fin:

            s = line.split()

            if line.startswith("#define") and s[1] in params:

                print(i)

                val = values[i]    

                line = "#define " + params[0] + " " + val + "\n"

            fout.write(line)
