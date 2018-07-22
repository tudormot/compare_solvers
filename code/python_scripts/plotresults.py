import fileinput
import sys
import numpy as np
import matplotlib.pyplot as plt

##arg  1 of command line should be name of graph
##arg  2 of command line should be "--iterations" or "--time" ,which specifies
##          whether the x-axis of our plots should display number of iterations
##          or time taken
##arg  3 of command line should be pardiso time to display for comparison purposes,
##          this parameter only matter if script called with --time option
##next args should be the files

if __name__ == "__main__":

    #quick check to see if args are correct
    if sys.argv[2] != "--iterations" and sys.argv[2] != "--time":
        print "Unrecognised command line arguments"
        raise Exception(' ')

    #we need some arrays to store the data to be ploted:
    no_files=len(sys.argv)-4;
    #print "No of files is "+ str(no_files)
    iter_number=np.zeros((no_files,10005),dtype=int)
    true_tolerance= np.zeros((no_files,10005))
    prec_residual= np.zeros((no_files,10005))
    true_residual=np.zeros((no_files,10005))
    no_iterations = np.zeros(no_files,dtype=int)
    petsc_ksp_time = np.zeros(no_files)
    petsc_setup_time = np.zeros(no_files)
    f_i = fileinput.input(sys.argv[4:])
    filenames = []
    last_token_found = [False for x in xrange(no_files)]
    file_number =-1
    for line in f_i:
        if f_i.filelineno() == 1:
            #print "firstline"
            filenames.append(f_i.filename())
            file_number = file_number + 1

            #print "file number is "+ str(file_number)
        l = line.split()
        #print l
        if l and l[0][0].isdigit() and l[1] == 'KSP' and no_iterations[file_number] <10000:
            #print l
            iter_number[file_number][no_iterations[file_number]] = int(l[0])
            prec_residual[file_number][no_iterations[file_number]] = float(l[5])
            true_residual[file_number][no_iterations[file_number]] = float(l[9])
            true_tolerance[file_number][no_iterations[file_number]] = float (l[11])
            no_iterations[file_number] = no_iterations[file_number] + 1
            #print no_iterations[file_number]
        if l and l[0] == 'Timing':
            if l[4] == 'routine':
                petsc_ksp_time[file_number] = petsc_ksp_time[file_number] + float(l[5])
                last_token_found[file_number] = True
            else:
                petsc_setup_time[file_number] = petsc_setup_time[file_number] + float(l[4])

    #determine file which had the least number of iterations completed:
    min_iterations = np.min(no_iterations)
    print "min_iterations is " + str(min_iterations) + "in file" + str(sys.argv[np.argmin(no_iterations)+4])



    #generate some better legend labels from command line input
    legend_label = [x.replace('_', '').replace('petsc', '').replace('job','').replace('.out','').replace('results', '') for x in filenames]

    if sys.argv[2] == "--iterations":
        plt.figure(1)
        for i in xrange(no_files):
            plt.plot(iter_number[i][0:no_iterations[i]:5],true_tolerance[i][0:no_iterations[i]:5],label=legend_label[i])
        plt.xlabel('Iteration Number')
        plt.ylabel('||r(i)||/||b||')
        plt.title(str(sys.argv[1])+"true tolerance")
        plt.grid(True)
        plt.xscale('log')
        plt.yscale('log')
        plt.legend(loc='best' )
    elif sys.argv[2] == '--time':
        plt.figure(1)
        for i in xrange(no_files):
            if last_token_found[i] is False:
                print "Error: Timing of KSP routine not found in input file: " +str(filenames[i])
                print "Can't plot timing data for this file"
                print "Did the program terminate early (due to cancel?)"
            else:
                time_per_iteration = petsc_ksp_time[i]/no_iterations[i]
                plt.plot([x * time_per_iteration + petsc_setup_time[i] for x in iter_number[i][0:no_iterations[i]:5]],true_tolerance[i][0:no_iterations[i]:5],label=legend_label[i])
        plt.axvline(x=float(sys.argv[3]), color='k', label='pardiso timing') #this introduces a vertical line to show pardiso time
        plt.xlabel('Time (s)')
        plt.ylabel('||r(i)||/||b||')
        plt.title(str(sys.argv[1])+"true tolerance")
        plt.grid(True)
        plt.xscale('log')
        plt.yscale('log')
        plt.legend(loc='best' )

    #plt.legend(legend_label,) #loc='upper left'


    #plt.figure(2)
    #for i in xrange(no_files):
    #    plt.plot(iter_number[i][0:min_iterations:5],true_residual[i][0:min_iterations:5]);
    #for i in xrange(no_files):
    #    pass
    #    plt.plot(iter_number[i][0:min_iterations:5],prec_residual[i][0:min_iterations:5]);

    #plt.xscale('log')
    #plt.xlabel('Iteration Numbe')
    #plt.ylabel('norm of residual')
    #plt.title(str(sys.argv[1])+"comparison between true and preconditioned residual")
    #plt.grid(True)
    #plt.legend([s + "_true_residual" for s in sys.argv[2:]] +[s + "_prec_residual" for s in sys.argv[2:]])

    plt.show()
