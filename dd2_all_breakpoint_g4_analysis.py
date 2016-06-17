# Import useful Library tools
import csv
import sys
import random
import numpy
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter,ScalarFormatter,FuncFormatter
from pylab import rcParams

rcParams['figure.figsize']=10,12

# Generate a starting random number from the clock so all trials are different
random.seed()

num_fudge=0

# Create some useful data structures
g4s = []
recombs = []

in_trans_vars = []
in_trans_non_vars = []
in_trans_telomeric_vars = []
in_dels = []
t_in_trans_vars = []
t_in_trans_telomeric_vars = []
t_in_trans_non_vars = []
t_in_dels = []
data_all=[]
t_data_all=[]

# 3d7
#chromosome_lengths = [640851, 947102, 1067971, 1200490, 1343557, 1418242, 1445207, 1472805, 1541735, 1687656, 2038340, 2271494, 2925236, 3291936]
# dd2
chromosome_lengths = [586228,858032,1033953,1113948,1437425,1352663,1406825,1365744,1475315,1661861,1993040,2167057,2932102,3257645]

telomeric_only = False

# Calculate total length of the genome from lengths above
total_length = reduce(lambda x, y: x+y, chromosome_lengths)

# Lists to store delta-bp values for statistical analysis
t_randoms = []
t_recombs = []
t_indels = []

# Function to find a chromosome on the whole genome (for dist. based on length, above)
def which_chromosome(p):
	global num_fudge
	sum=0
	chromosome=1
	for l in chromosome_lengths:
		sum = sum + l
		if sum > p: return chromosome
		chromosome = chromosome + 1

	num_fudge = num_fudge+1
	return chromosome-1

# Function that displays a self-overwriting progress indicator (so we see feedback on long waits)
def progress(percent):
	sys.stdout.write("\r           \r"+str(round(percent,2))+"%")
	sys.stdout.flush()

# The main function that samples distances from G4s over a number of trials.
# It returns a cumulative moving average of the distance, and stores all of the
# distances in the t_randoms list (which is used later to calculate the t-test)
def find_average_distance(number_of_trials, null_type):
	global t_randoms
	t_randoms = []
	n=0
	# Cumulative moving average
	cma=0.0
	onetick=number_of_trials/1000

	for trial in range(0, number_of_trials):
		# Choose a chromosome, either from uniform distribution or else
		# from a distribution reflecting the distribution of recombination
		# sites across the whole genome.
		#if chromosome_distribution_uniform_across_chromosomes:
		if null_type==0:
			print "SHOULD NOT BE USING THIS DISTRIBUTION"
			# Chrom. chosen based on 1/14 chance
			chromosome = random.randint(1, 14)
		
		elif null_type==1:
			# Chrom. chosen based on distribution of recomb sites across chroms.
			c = random.choice(recombs)
			chromosome = c[0]
		
		elif null_type==2:
			# Chromosome based on size of genome
			random_point = random.randint(0, total_length )
			chromosome = which_chromosome( random_point )

		# Choose a random position on the chosen chromosome
		clen = chromosome_lengths[chromosome-1]
		position = random.randint(0, chromosome_lengths[chromosome-1])
		
		# if telomeric_only:
		# 	if bool(random.getrandbits(1)):
		# 		position = random.randint(0,65000)
		# 	else:
		# 		position = random.randint(clen-65000, clen)

		# Find nearest DSS region on the chromosome
		nearest = find_nearest_g4(chromosome, position)

		# Save the value in a list to do a t-test later
		t_randoms.append(nearest[1]);

		# Update moving average with new distance to DSS
		# cma = (nearest[1]+(n*cma))/(n+1)
		# n = n + 1

		# Progress indicator as many trials takes time
		if number_of_trials>1000 and trial%onetick==0:
			progress(100.0*float(trial)/number_of_trials)

# Returns a filtered list of G4s, only the ones on the specified chromosome
def g4s_on_chromosome(chromosome):
	return filter(lambda x: x[0]==chromosome, g4s)

# Searches for the nearest G4 to the specified chromosome/position.
def find_nearest_g4(chromosome, position):
	closest = sys.maxint
	closest_site = ["no g4 on this chromosome"]

	possible_sites = g4s_on_chromosome(chromosome)
	if len(possible_sites) < 1: exit()
	for site in possible_sites:
		# Check the start of the G4 site
		if abs(position-site[1]) < closest:
			closest = abs(position-site[1])
			closest_site = site
		# Check the end of the G4 site
		if abs(position-site[2]) < closest:
			closest = abs(position-site[2])
			closest_site = site

	return [closest_site, closest]

# Loads up the data from the two .csv files
def load_data():
	global g4s, recombs, in_trans_vars, in_trans_telomeric_vars, in_trans_non_vars, in_dels;
	global t_in_trans_vars, t_in_trans_telomeric_vars, t_in_trans_non_vars, t_in_dels;
	global data_all, t_data_all

	with open('data/dd2/dd2-g4s.csv', 'r') as g4csv:
	#with open('../../october-2015/g4_locations.csv','r') as g4csv:
		r = csv.reader(g4csv, delimiter=',')
		for row in r:
			chm = int(row[0])
			sta = int(row[1])
			end = int(row[2])
			g4s = g4s + [[chm,sta,end]]

	g4s=sorted(g4s,key=lambda x: x[0])

	with open('data/dd2/dd2-vars.csv', 'r') as recombscsv:
		r = csv.reader(recombscsv,delimiter=',')
		for row in r:
			chm=int(row[0]) # chromosome
			pos=int(row[1]) # start
			fin=int(row[2]) # end
			# bln=int(row[3]) # length
			# typ=int(row[4]) # type, 1=trans, 2=indel
			# var=row[5] # is var? Y/N
			# tel=row[6] # telomeric Y/N

			nearest=find_nearest_g4(chm,pos)
			# in_trans_vars=in_trans_vars+[[chm,pos,nearest[1],nearest[0]]]
			# t_in_trans_vars = t_in_trans_vars+[nearest[1]]
			print str(chm)+","+str(pos)+","+str(fin)+","+str(nearest)
			data_all=data_all+[[chm,pos,nearest[1],nearest[0]]]
			t_data_all=t_data_all+[nearest[1]]

	# 		if typ==1:

	# 			if tel=="Y":
	# 				in_trans_telomeric_vars=in_trans_telomeric_vars+[[chm,pos,nearest[1],nearest[0]]]
	# 				t_in_trans_telomeric_vars=t_in_trans_telomeric_vars+[nearest[1]]

	# 			if var=="N":
	# 				in_trans_non_vars=in_trans_non_vars+[[chm,pos,nearest[1],nearest[0]]]
	# 				t_in_trans_non_vars = t_in_trans_non_vars+[nearest[1]]
	# 			else:
	# 				in_trans_vars=in_trans_vars+[[chm,pos,nearest[1],nearest[0]]]
	# 				t_in_trans_vars = t_in_trans_vars+[nearest[1]]

	# 		if typ==2:
	# 			in_dels=in_dels+[[chm,pos,nearest[1],nearest[0]]]
	# 			t_in_dels = t_in_dels+[nearest[1]]

	# in_trans_non_vars = sorted(in_trans_non_vars, key=lambda x: x[0])
	# in_trans_vars = sorted(in_trans_vars,key = lambda x: x[0])
	# in_trans_telomeric_vars = sorted(in_trans_telomeric_vars,key=lambda x:x[0])
	# in_dels=sorted(in_dels,key=lambda x: x[0])
	data_all=sorted(data_all,key=lambda x:x[0])

# Prints out some data
def display():
	print "\nRecombination Sites\n-------------------\n\n"
	print "Chrom.\tPosition\t(delta bp)\tNearest G4 Location"
	for r in recombs:
		print str(r[0])+"\t"+str(r[1])+"\t\t"+str(r[2])+"\t\t"+str(r[3])

def do_experiment( breakpoint_type, recombinations, t_recombinations, trials, null_type ):
	global t_randoms, recombs, t_recombs

	recombs = recombinations
	t_recombs = t_recombinations
#	for r in recombs:
#		print ", ".join(map(lambda x: str(x), r))

	find_average_distance( trials, null_type )
	
	# Calculate statistics
	u_recombs = numpy.mean(t_recombs)
	u_randoms = numpy.mean(t_randoms)

	x_recombs = numpy.median(t_recombs)
	x_randoms = numpy.median(t_randoms)

	t_test = stats.ttest_ind(t_recombs, t_randoms, False)
	t_tst1 = stats.ttest_1samp(t_recombs, u_randoms)

	if null_type==0: null_comparison = "uniform chromosome distribution"
	if null_type==1: null_comparison = "proportional to # breaks on chromosome"
	if null_type==2: null_comparison = "every basepair equal chance"	

	print "\r"+breakpoint_type+", "+str(len(t_recombs))+", "+str(u_recombs)+", "+\
		null_comparison+", "+str(u_randoms)+", "+str(u_recombs-u_randoms)+", "+str(t_test[1])+", "+str(t_tst1[1])+", "+str(x_recombs)+", "+str(x_randoms)

	return [[t_recombs, t_randoms]]

def run_experiments():
	e = []
	# data_all = in_trans_vars+in_trans_non_vars+in_dels
	# t_data_all = t_in_trans_vars+t_in_trans_non_vars+t_in_dels
	# data_in_trans = in_trans_vars+in_trans_non_vars
	# t_data_in_trans = t_in_trans_vars+t_in_trans_non_vars

	num_trials = 1<<20

	print "\r"+"breakpoint_type,num_breakpoints,bp_mean,null_comparison,rand_mean,mean_diff,p-value(twosamp),p-value(onesamp),bp_median,rand_median"

	e=e+do_experiment( "all", data_all, t_data_all, num_trials, 1)
	e=e+do_experiment( "all", data_all, t_data_all, num_trials, 2)
	# e=e+do_experiment( "in-trans non-var", in_trans_non_vars, t_in_trans_non_vars, num_trials, 1)
	# e=e+do_experiment( "in-trans non-var", in_trans_non_vars, t_in_trans_non_vars, num_trials, 2)
	# e=e+do_experiment( "in-trans", data_in_trans, t_data_in_trans, num_trials, 1)
	# e=e+do_experiment( "in-trans", data_in_trans, t_data_in_trans, num_trials, 2)
	# e=e+do_experiment( "in-dels", in_dels, t_in_dels, num_trials, 1)
	# e=e+do_experiment( "in-dels", in_dels, t_in_dels, num_trials, 2)
	# e=e+do_experiment( "var only", in_trans_vars, t_in_trans_vars, num_trials, 1)
	# e=e+do_experiment( "var only", in_trans_vars, t_in_trans_vars, num_trials, 2)
	# e=e+do_experiment( "telomeric var only", in_trans_telomeric_vars, t_in_trans_telomeric_vars, num_trials, 1)
	# e=e+do_experiment( "telomeric var only", in_trans_telomeric_vars, t_in_trans_telomeric_vars, num_trials, 2)
	return e

# Program starts here
load_data()

e = run_experiments()

#############################################################################
###
### Figures
###
#############################################################################
def megabases(x, pos):
    'The two args are the value and tick position'
    return '%1.1f' % (x*1e-6)
    #return '%1.1fMbp' % (x*1e-6)

labels=['Recombination Sites','Random Sites']
xlabel = "Site Type"
ylabel = "Distance from G4"

f, ((ax1,ax2), (ax3,ax4), (ax5,ax6)) = plt.subplots(3, 2)
f.subplots_adjust(hspace=0.35,wspace=0.4,top=0.9,bottom=0.1,right=0.9,left=0.18)
f.title="Breakpoint-frequency dependent sampling"
	
axes=[ax1, ax2, ax3, ax4, ax5, ax6]

for ax in axes:
	ax.set_ylabel("Distance from nearest G4 (Mbp)")
	ax.set_yscale('linear')
	ax.set_ylim(0,1600000)
	ax.yaxis.set_major_formatter(FuncFormatter(megabases))

ax1.boxplot(e[0],labels=labels,showmeans=True)
ax1.set_title("All breakpoints (n="+str(len(e[1][0]))+")")

# ax2.boxplot(e[2],labels=labels,showmeans=True)
# ax2.set_title("In-trans non-var (n="+str(len(e[3][0]))+")")

# ax3.boxplot(e[4],labels=labels,showmeans=True)
# ax3.set_title("In-trans (n="+str(len(e[5][0]))+")")

# ax4.boxplot(e[6],labels=labels,showmeans=True)
# ax4.set_title("In-dels (n="+str(len(e[7][0]))+")")

# ax5.boxplot(e[8],labels=labels,showmeans=True)
# ax5.set_title("Var only (n="+str(len(e[9][0]))+")")

# ax6.boxplot(e[10],labels=labels,showmeans=True)
# ax6.set_title("Telomeric var only (n="+str(len(e[11][0]))+")")

plt.suptitle("Measured breakpoint distance from G4 vs random sample distance from G4\n(probability of chromosome selection proportional to # breakpoints on chromosome)")
plt.savefig("output/dd2-nocrosses-1.png")
plt.figure(2)
f, ((ax1,ax2), (ax3,ax4), (ax5,ax6)) = plt.subplots(3, 2)
f.subplots_adjust(hspace=0.35,wspace=0.4,top=0.9,bottom=0.1,right=0.9,left=0.18)
f.title="Chromosome-length dependent sampling"
axes=[ax1, ax2, ax3, ax4, ax5, ax6]

for ax in axes:
	ax.set_ylabel("Distance from nearest G4 (Mbp)")
	ax.set_yscale('linear')
	ax.set_ylim(0,1600000)
	ax.yaxis.set_major_formatter(FuncFormatter(megabases))

ax1.boxplot(e[1],labels=labels,showmeans=True)
ax1.set_title("All breakpoints (n="+str(len(e[0][0]))+")")

# ax2.boxplot(e[3],labels=labels,showmeans=True)
# ax2.set_title("In-trans non-var (n="+str(len(e[2][0]))+")")

# ax3.boxplot(e[5],labels=labels,showmeans=True)
# ax3.set_title("In-trans (n="+str(len(e[4][0]))+")")

# ax4.boxplot(e[7],labels=labels,showmeans=True)
# ax4.set_title("In-dels (n="+str(len(e[6][0]))+")")

# ax5.boxplot(e[9],labels=labels,showmeans=True)
# ax5.set_title("Var only (n="+str(len(e[8][0]))+")")

# ax6.boxplot(e[11],labels=labels,showmeans=True)
# ax6.set_title("Telomeric var only (n="+str(len(e[10][0]))+")")
plt.suptitle("Measured breakpoint distance from G4 vs random sample distance from G4\n(probability of chromosome selection proportional to size)")
plt.savefig("output/dd2-nocrosses-2.png")

plt.show()
