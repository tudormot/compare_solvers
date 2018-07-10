import fileinput
import sys
import numpy as np
import matplotlib.pyplot as plt

##arg  one of command line should be name of graph
##next args should be the files

if __name__ == "__main__":

    #we need some arrays to store the data to be ploted:
    no_files=len(sys.argv)-2;
    print "No of files is "+ str(no_files)
    iter_number=np.zeros((no_files,10000),dtype=int)
    true_tolerance= np.zeros((no_files,10000))
    prec_residual= np.zeros((no_files,10000))
    true_residual=np.zeros((no_files,10000))
    no_iterations = np.zeros(no_files,dtype=int)
    f_i = fileinput.input(sys.argv[2:])

    file_number =-1
    for line in f_i:
        if f_i.filelineno() == 1:
            #print "firstline"
            file_number = file_number + 1
            #print "file number is "+ str(file_number)
        l = line.split()
        if l[0][0].isdigit():
            #print l
            iter_number[file_number][f_i.filelineno()-1] = int(l[0])
            prec_residual[file_number][f_i.filelineno()-1] = float(l[5])
            true_residual[file_number][f_i.filelineno()-1] = float(l[9])
            true_tolerance[file_number][f_i.filelineno()-1] = float (l[11])
            no_iterations[file_number] = no_iterations[file_number] + 1
        if no_iterations[file_number] > 10000:
            f_i.nextfile()

    #determine file which had the least number of iterations completed:
    min_iterations = np.min(no_iterations)
    print "min_iterations is " + str(min_iterations)


    plt.figure(1)
    for i in xrange(no_files):
        plt.plot(iter_number[i][0:min_iterations:50],true_tolerance[i][0:min_iterations:50]);

    print iter_number[0:min_iterations:50]

    plt.xlabel('Iteration Number')
    plt.ylabel('||r(i)||/||b||')
    plt.title(str(sys.argv[1])+"true tolerance")
    plt.grid(True)
    plt.legend(sys.argv[2:] ) #loc='upper left'
    plt.figure(2)


    plt.figure(2)
    for i in xrange(no_files):
        plt.plot(iter_number[i][0:min_iterations:50],true_residual[i][0:min_iterations:50]);
    for i in xrange(no_files):
        plt.plot(iter_number[i][0:min_iterations:50],prec_residual[i][0:min_iterations:50]);

    plt.xlabel('Iteration Number')
    plt.ylabel('norm of residual')
    plt.title(str(sys.argv[1])+"comparison between true and preconditioned residual")
    plt.grid(True)
    plt.legend([s + "true_residual" for s in sys.argv[2:]] + [s + "prec_residual" for s in sys.argv[2:]])

    plt.show()
    #plt.plot(t, s)
    #
    #
    #
    #plt.savefig("test.png")
    #

    #plt.plot(x, x)
    #plt.plot(x, 2 * x)
    #plt.plot(x, 3 * x)
    #plt.plot(x, 4 * x)

    #plt.legend(['y = x', 'y = 2x', 'y = 3x', 'y = 4x'], loc='upper left')

    #plt.show()

    #for i in xrange(int(no_iterations[0])):
        #print str(iter_number[0][i]) +"###" + str(true_tolerance[0][i])
            #print l

        #if not mp3filename or mp3filename.startswith('#'):
        #    continue
        #item = SubElement(rss, 'item')
        #itle = SubElement(item, 'title')
        #title.text = mp3filename
        #encl = SubElement(item, 'enclosure', {'type':'audio/mpeg', 'url':mp3filename})
