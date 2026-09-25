# Modified 2026-09-24 by Valis Sowilo (github.com/ValisSowilo): updated for NEURON 9, LFPy 2.3, numpy 2 and pandas 2+. See CHANGES_NEURON9.md.
#===========================================================================
# Import, Set up MPI Variables, Load Necessary Files
#===========================================================================
from mpi4py import MPI
import time
tic_0 = time.perf_counter()
import os
from os.path import join
import sys
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
from scipy import stats as st
import neuron
from neuron import h, gui
import LFPy
from LFPy import NetworkCell, Network, Synapse, RecExtElectrode, StimIntElectrode, CurrentDipoleMoment

agegroup = 'y'

#MPI variables:
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()
GLOBALSEED = int(sys.argv[1])
STIMSEED = int(sys.argv[2])

# Create new RandomState for each RANK
SEED = GLOBALSEED*10000
np.random.seed(SEED + RANK)
local_state = np.random.RandomState(SEED + RANK)
halfnorm_rv = st.halfnorm
halfnorm_rv.random_state = local_state
uniform_rv = st.uniform
uniform_rv.random_state = local_state

# Separate seed for background and stimulus inputs
SEED_FOR_STIM = STIMSEED*10000
local_state_stim = np.random.RandomState(SEED_FOR_STIM + RANK)
uniform_rv_stim = st.uniform
uniform_rv_stim.random_state = local_state_stim

from net_params import *

#Mechanisms and files
print('Mechanisms found: ', os.path.isfile('mod/x86_64/special')) if RANK==0 else None
neuron.h('forall delete_section()')
neuron.load_mechanisms('mod/')
h.load_file('net_functions.hoc')
h.load_file('models/biophys_HL5PN1' + agegroup + '.hoc')
h.load_file('models/biophys_HL5MN1.hoc')
h.load_file('models/biophys_HL5BN1.hoc')
h.load_file('models/biophys_HL5VN1.hoc')


print('Importing, setting up MPI variables and loading necessary files took ', str((time.perf_counter() - tic_0)/60)[:5], 'minutes') if RANK==0 else None

#===========================================================================
# Simulation, Analysis, and Plotting Controls
#===========================================================================
#Sim
TESTING = False # i.e.g generate 1 cell/pop, with 0.1 s runtime
no_connectivity = False
saveConMat = False # Generates many files when run in parallel and matrix needs to be rebuilt post-hoc

stimulate = True # Add a stimulus
MDD = False #decrease PN GtonicApic and MN2PN weight by 40%

rec_LFP = True #record LFP from center of layer
rec_DIPOLES = True #record population - wide dipoles

run_circuit_functions = True


#===========================================================================
# Params
#===========================================================================
dt = 0.025 #for both cell and network
tstart = 0.
tstop = 7000.
celsius = 34.
v_init = -80. #for both cell and network

N_cells = 1000.
N_HL5PN = int(0.70*N_cells)
N_HL5MN = int(0.15*N_cells)
N_HL5BN = int(0.10*N_cells)
N_HL5VN = int(0.05*N_cells)

if TESTING:
	OUTPUTPATH = 'Circuit_output_testing'
	N_HL5PN = 1
	N_HL5MN = 1
	N_HL5BN = 1
	N_HL5VN = 1
	tstop = 100
	print('Running test...') if RANK ==0 else None

else:
	OUTPUTPATH = 'Circuit_output'
	print('Running full simulation...') if RANK==0 else None

cellnums = [N_HL5PN, N_HL5MN, N_HL5BN, N_HL5VN] # set after TESTING so test runs use the reduced population sizes

COMM.Barrier()

networkParams = {
	'dt' : dt,
	'tstart': tstart,
	'tstop' : tstop,
	'v_init' : v_init,
	'celsius' : celsius,
	'OUTPUTPATH' : OUTPUTPATH,
	'verbose': False}

def L5depth_shift(depth):
	depth_shift = L5upper + L5PNmaxApic # shift depth to avoid PN crossings at pia
	new_depth = depth - depth_shift

	return new_depth

L5PNmaxApic = 1900. # ~max length of PN apic
L5upper = -1600. # shallow border
L5lower = -2300. # deep border

L5rangedepth = abs(L5lower-L5upper)/2 # Sets +/- range of soma positions
L5meandepth = L5upper - L5rangedepth
L5minSynLoc = L5depth_shift(L5meandepth)-L5rangedepth*2 # Sets minimum depth & range (loc to loc+scale) for synapses

both = {'section' : ['apic', 'dend'],
		'fun' : [uniform_rv, halfnorm_rv],
		'funargs' : [{'loc':L5minSynLoc, 'scale':abs(L5minSynLoc)},{'loc':L5minSynLoc, 'scale':abs(L5minSynLoc)}],
		'funweights' : [1, 1.]}
apic = {'section' : ['apic'],
		'fun' : [uniform_rv],
		'funargs' : [{'loc':L5minSynLoc, 'scale':abs(L5minSynLoc)}],
		'funweights' : [1.]}
dend = {'section' : ['dend'],
		'fun' : [uniform_rv],
		'funargs' : [{'loc':L5minSynLoc, 'scale':abs(L5minSynLoc)}],
		'funweights' : [1.]}
PNdend = {'section' : ['dend'],
		'fun' : [halfnorm_rv],
		'funargs' : [{'loc':L5minSynLoc, 'scale':abs(L5minSynLoc)}],
		'funweights' : [1.]}

pos_args = {
	'none' : dend,
	'HL5PN1HL5PN1' : both,
	'HL5PN1HL5MN1' : dend,
	'HL5PN1HL5BN1' : dend,
	'HL5PN1HL5VN1' : dend,

	'HL5MN1HL5PN1' : apic,
	'HL5MN1HL5MN1' : dend,
	'HL5MN1HL5BN1' : dend,
	'HL5MN1HL5VN1' : dend,

	'HL5BN1HL5PN1' : PNdend,
	'HL5BN1HL5MN1' : dend,
	'HL5BN1HL5BN1' : dend,
	'HL5BN1HL5VN1' : dend,

	'HL5VN1HL5PN1' : PNdend,
	'HL5VN1HL5MN1' : dend,
	'HL5VN1HL5BN1' : dend,
	'HL5VN1HL5VN1' : dend}

L5_pop_args = {'radius':250,
				'loc':L5depth_shift(L5meandepth),
				'scale':L5rangedepth*4,
				'cap': L5rangedepth}

rotations = {'HL5PN1':{'x':1.57,'y':3.72},
			 'HL5MN1':{'x':1.77,'y':2.77},
			 'HL5BN1':{'x':1.26,'y':2.57},
			 'HL5VN1':{'x':-1.57,'y':3.57}}

# class RecExtElectrode parameters:
L5_size = abs(L5upper - L5lower)
e1 = L5depth_shift(L5upper-L5_size*.5)

LFPelectrodeParameters = dict(
	x=np.zeros(1),
	y=np.zeros(1),
	z=[e1],
	N=np.array([[0., 1., 0.] for _ in range(1)]),
	r=5.,
	n=50,
	sigma=0.3,
	method="root_as_point") # see SomaAsPointElectrode below


#method Network.simulate() parameters

simargs = {'rec_imem': False,
		   'rec_vmem': False,
		   'rec_ipas': False,
		   'rec_icap': False,
		   'rec_isyn': False,
		   'rec_vmemsyn': False,
		   'rec_istim': False,
		   'rec_pop_contributions': rec_DIPOLES, # per-population contributions, used for the per-population dipole moments
		   'rec_variables': [],
		   'to_memory': True,
		   'to_file': False,
		   'file_name':'OUTPUT.h5'}

#===========================================================================
# Functions
#===========================================================================
def generateSubPop(popsize,mname,popargs,Gou,Gtonic,GtonicApic):
	print('Initiating ' + mname + ' population...') if RANK==0 else None
	morphpath = 'morphologies/' + mname + '.swc'
	templatepath = 'models/NeuronTemplate.hoc'
	templatename = 'NeuronTemplate'
	
	pt3d = True
	
	cellParams = {
		'morphology': morphpath,
		'templatefile': templatepath,
		'templatename': templatename,
		'templateargs': morphpath,
		'v_init': v_init, #initial membrane potential, d=-65
		'passive': False,#initialize passive mechs, d=T, should be overwritten by biophys
		'dt': dt,
		'tstart': 0.,
		'tstop': tstop,#defaults to 100
		'nsegs_method': None,
		'pt3d': pt3d,#use pt3d-info of the cell geometries switch, d=F
		'delete_sections': False,
		'verbose': False}#verbose output switch, for some reason doens't work, figure out why}
	
	if mname in rotations.keys():
		rotation = rotations.get(mname)
	
	popParams = {
		'CWD': None,
		'CELLPATH': None,
		'Cell' : LFPy.NetworkCell, #play around with this, maybe put popargs into here
		'POP_SIZE': popsize,
		'name': mname,
		'cell_args' : {**cellParams},
		'pop_args' : popargs,
		'rotation_args' : rotation}
	
	network.create_population(**popParams)
	
	# Add biophys, OU processes, & tonic inhibition to cells
	for cellind in range(0,len(network.populations[mname].cells)):
		rseed = int(local_state_stim.uniform()*SEED_FOR_STIM)
		biophys = 'h.biophys_' + mname + '(network.populations[\'' + mname + '\'].cells[' + str(cellind) + '].template)'
		exec(biophys)
		h.createArtificialSyn(rseed,network.populations[mname].cells[cellind].template,Gou)
		h.addTonicInhibition(network.populations[mname].cells[cellind].template,Gtonic,GtonicApic)

class SomaAsPointElectrode(RecExtElectrode):
	# Reproduces LFPy 2.0's method='soma_as_point', where every soma segment is a point source.
	# LFPy >= 2.1 renamed it 'root_as_point' but only treats segment 0 as the point source, which
	# in a Network (all cells on a RANK merged into one cell) would be the first cell's soma only.
	def get_transformation_matrix(self):
		self.kwargs['rootinds'] = self.cell.get_idx('soma')
		return super().get_transformation_matrix()

def gaussian(x, mu, sig):
	return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def addStimulus(stim_list, cell_nums):
	cell_names = ['HL5PN1','HL5MN1','HL5BN1','HL5VN1']
	for stim in stim_list:
		for gid, cell in zip(network.populations[stim['mname']].gids, network.populations[stim['mname']].cells):
			ferror = 0.3 # +/- fractional error (%/100) of Gaussian peak (i.e. max synapse numbers) to shift stimulus input magnitude by
			x_values = np.linspace(0, cell_nums[cell_names.index(stim['mname'])]-1, cell_nums[cell_names.index(stim['mname'])])
			mu_ind = int(cell_nums[cell_names.index(stim['mname'])]*stim['Orientation']*0.5/180) # find index with peak input depending on orientation
			sigma = int(cell_nums[cell_names.index('HL5PN1')]*42*0.5/180) # Set sigma to 42 degrees (Bensmaia et al 2008) and scale by number of PN neurons
			gfunc = gaussian(x_values, mu_ind, sigma)+uniform_rv_stim.rvs(loc=-ferror,scale=ferror*2) # Gaussian vector (y=[0,1]) of length = cellnum
			gfunc[gfunc<0] = 0 # remove negative indices from random shift
			gfunc[x_values>cell_nums[cell_names.index('HL5PN1')]*0.5] = 0 # Ensures only 50% of PNs receive orientation-tuning
			subcount = sum(cell_nums[:cell_names.index(stim['mname'])])
			numsyn = stim['loc_num']*gfunc[gid-subcount]
			if (numsyn>=1):
				numsyn = int(np.round(numsyn)) # because int() rounds down always
				idx = cell.get_rand_idx_area_norm(section=stim['loc'], nidx=numsyn)
				for i in idx:
					time_d=0
					syn = Synapse(cell=cell, idx=i, syntype=stim['stim_type'], weight=1,**stim['syn_params'])
					while time_d <= 0:
						time_d = uniform_rv_stim.rvs(loc = stim['delay'], scale = stim['delay_std'])
					syn.set_spike_times_w_netstim(noise=0.3, start=(stim['start_time']+time_d), number=stim['num_stim'], interval=stim['interval'], seed=SEED_FOR_STIM)

#===========================================================================
# Sim
#===========================================================================
network = Network(**networkParams)

if MDD:
	syn_params['HL5MN1HL5PN1']['gmax'] = syn_params['HL5MN1HL5PN1']['gmax']*0.6
	syn_params['HL5MN1HL5MN1']['gmax'] = syn_params['HL5MN1HL5MN1']['gmax']*0.6
	syn_params['HL5MN1HL5BN1']['gmax'] = syn_params['HL5MN1HL5BN1']['gmax']*0.6
	syn_params['HL5MN1HL5VN1']['gmax'] = syn_params['HL5MN1HL5VN1']['gmax']*0.6
	Gtonic_MN -= (Ncont_MN2MN*con_MN2MN*connection_prob['HL5MN1HL5MN1'])/(Ncont_BN2MN*con_BN2MN*connection_prob['HL5BN1HL5MN1']+Ncont_VN2MN*con_VN2MN*connection_prob['HL5VN1HL5MN1']+Ncont_MN2MN*con_MN2MN*connection_prob['HL5MN1HL5MN1'])*Gtonic_MN*0.4
	Gtonic_BN -= (Ncont_MN2BN*con_MN2BN*connection_prob['HL5MN1HL5BN1'])/(Ncont_BN2BN*con_BN2BN*connection_prob['HL5BN1HL5BN1']+Ncont_VN2BN*con_VN2BN*connection_prob['HL5VN1HL5BN1']+Ncont_MN2BN*con_MN2BN*connection_prob['HL5MN1HL5BN1'])*Gtonic_BN*0.4
	Gtonic_VN -= (Ncont_MN2VN*con_MN2VN*connection_prob['HL5MN1HL5VN1'])/(Ncont_BN2VN*con_BN2VN*connection_prob['HL5BN1HL5VN1']+Ncont_VN2VN*con_VN2VN*connection_prob['HL5VN1HL5VN1']+Ncont_MN2VN*con_MN2VN*connection_prob['HL5MN1HL5VN1'])*Gtonic_VN*0.4
	GtonicApic_PN = Gtonic_PN*0.6
	print('MN tonic reduced by ',(Ncont_MN2MN*con_MN2MN*connection_prob['HL5MN1HL5MN1'])/(Ncont_BN2MN*con_BN2MN*connection_prob['HL5BN1HL5MN1']+Ncont_VN2MN*con_VN2MN*connection_prob['HL5VN1HL5MN1']+Ncont_MN2MN*con_MN2MN*connection_prob['HL5MN1HL5MN1'])*0.4*100, '%') if RANK==0 else None
	print('BN tonic reduced by ',(Ncont_MN2BN*con_MN2BN*connection_prob['HL5MN1HL5BN1'])/(Ncont_BN2BN*con_BN2BN*connection_prob['HL5BN1HL5BN1']+Ncont_VN2BN*con_VN2BN*connection_prob['HL5VN1HL5BN1']+Ncont_MN2BN*con_MN2BN*connection_prob['HL5MN1HL5BN1'])*0.4*100, '%') if RANK==0 else None
	print('VN tonic reduced by ',(Ncont_MN2VN*con_MN2VN*connection_prob['HL5MN1HL5VN1'])/(Ncont_BN2VN*con_BN2VN*connection_prob['HL5BN1HL5VN1']+Ncont_VN2VN*con_VN2VN*connection_prob['HL5VN1HL5VN1']+Ncont_MN2VN*con_MN2VN*connection_prob['HL5MN1HL5VN1'])*0.4*100, '%') if RANK==0 else None

# Generate Populations
tic = time.perf_counter()
generateSubPop(N_HL5PN,'HL5PN1',L5_pop_args,Gou_PN,Gtonic_PN,GtonicApic_PN)
generateSubPop(N_HL5MN,'HL5MN1',L5_pop_args,Gou_MN,Gtonic_MN,Gtonic_MN)
generateSubPop(N_HL5BN,'HL5BN1',L5_pop_args,Gou_BN,Gtonic_BN,Gtonic_BN)
generateSubPop(N_HL5VN,'HL5VN1',L5_pop_args,Gou_VN,Gtonic_VN,Gtonic_VN)
COMM.Barrier()

print('Instantiating all populations took ', str((time.perf_counter() - tic_0)/60)[:5], 'minutes') if RANK==0 else None

tic = time.perf_counter()

# Synaptic Connection Parameters
E_syn = neuron.h.ProbAMPANMDA
I_syn = neuron.h.ProbUDFsyn

weightFunction = local_state.normal

WP = {'loc':1, 'scale':0}
MW=1


weightArguments = [[WP,WP,WP,WP],[WP,WP,WP,WP],[WP,WP,WP,WP],[WP,WP,WP,WP]]
minweightArguments = [[MW,MW,MW,MW],[MW,MW,MW,MW],[MW,MW,MW,MW],[MW,MW,MW,MW]]

delayFunction = local_state.normal
delayParams = {'loc':0.5, 'scale':0}
mindelay=0.5
delayArguments = np.full([4, 4], delayParams)

connectionProbability = [[connection_prob['HL5PN1HL5PN1'],connection_prob['HL5PN1HL5MN1'],connection_prob['HL5PN1HL5BN1'],connection_prob['HL5PN1HL5VN1']],
						 [connection_prob['HL5MN1HL5PN1'],connection_prob['HL5MN1HL5MN1'],connection_prob['HL5MN1HL5BN1'],connection_prob['HL5MN1HL5VN1']],
						 [connection_prob['HL5BN1HL5PN1'],connection_prob['HL5BN1HL5MN1'],connection_prob['HL5BN1HL5BN1'],connection_prob['HL5BN1HL5VN1']],
						 [connection_prob['HL5VN1HL5PN1'],connection_prob['HL5VN1HL5MN1'],connection_prob['HL5VN1HL5BN1'],connection_prob['HL5VN1HL5VN1']]]
if no_connectivity:
	connectionProbability = np.zeros_like(connectionProbability)

multapseFunction = local_state.normal
multapseArguments = [[mult_syns['HL5PN1HL5PN1'],mult_syns['HL5PN1HL5MN1'],mult_syns['HL5PN1HL5BN1'],mult_syns['HL5PN1HL5VN1']],
					 [mult_syns['HL5MN1HL5PN1'],mult_syns['HL5MN1HL5MN1'],mult_syns['HL5MN1HL5BN1'],mult_syns['HL5MN1HL5VN1']],
					 [mult_syns['HL5BN1HL5PN1'],mult_syns['HL5BN1HL5MN1'],mult_syns['HL5BN1HL5BN1'],mult_syns['HL5BN1HL5VN1']],
					 [mult_syns['HL5VN1HL5PN1'],mult_syns['HL5VN1HL5MN1'],mult_syns['HL5VN1HL5BN1'],mult_syns['HL5VN1HL5VN1']]]

synapseParameters = [[syn_params['HL5PN1HL5PN1'],syn_params['HL5PN1HL5MN1'],syn_params['HL5PN1HL5BN1'],syn_params['HL5PN1HL5VN1']],
					 [syn_params['HL5MN1HL5PN1'],syn_params['HL5MN1HL5MN1'],syn_params['HL5MN1HL5BN1'],syn_params['HL5MN1HL5VN1']],
					 [syn_params['HL5BN1HL5PN1'],syn_params['HL5BN1HL5MN1'],syn_params['HL5BN1HL5BN1'],syn_params['HL5BN1HL5VN1']],
					 [syn_params['HL5VN1HL5PN1'],syn_params['HL5VN1HL5MN1'],syn_params['HL5VN1HL5BN1'],syn_params['HL5VN1HL5VN1']]]

synapsePositionArguments = [[pos_args['HL5PN1HL5PN1'],pos_args['HL5PN1HL5MN1'],pos_args['HL5PN1HL5BN1'],pos_args['HL5PN1HL5VN1']],
							[pos_args['HL5MN1HL5PN1'],pos_args['HL5MN1HL5MN1'],pos_args['HL5MN1HL5BN1'],pos_args['HL5MN1HL5VN1']],
							[pos_args['HL5BN1HL5PN1'],pos_args['HL5BN1HL5MN1'],pos_args['HL5BN1HL5BN1'],pos_args['HL5BN1HL5VN1']],
							[pos_args['HL5VN1HL5PN1'],pos_args['HL5VN1HL5MN1'],pos_args['HL5VN1HL5BN1'],pos_args['HL5VN1HL5VN1']]]

for i, pre in enumerate(network.population_names):
	for j, post in enumerate(network.population_names):
		
		connectivity = network.get_connectivity_rand(
			pre=pre,
			post=post,
			connprob=connectionProbability[i][j])
		
		if saveConMat:
			gids = np.array(network.populations[post].gids).astype(int)
			temp_cmat = zip(np.transpose(connectivity),gids)
			np.save(os.path.join(OUTPUTPATH,'conmat_'+pre+post+'_'+str(GLOBALSEED)+'StimSeed'+str(STIMSEED)+'_Rank'+str(RANK)+'.npy'), temp_cmat)
		
		(conncount, syncount) = network.connect(
			pre=pre, post=post,
			connectivity=connectivity,
			syntype=E_syn if pre=='HL5PN1' else I_syn,
			synparams=synapseParameters[i][j],
			weightfun=weightFunction,
			weightargs=weightArguments[i][j],
			minweight=minweightArguments[i][j],
			delayfun=delayFunction,
			delayargs=delayArguments[i][j],
			mindelay=mindelay,
			multapsefun=multapseFunction,
			multapseargs=multapseArguments[i][j],
			syn_pos_args=synapsePositionArguments[i][j])

print('Connecting populations took ', str((time.perf_counter() - tic_0)/60)[:5], 'minutes') if RANK==0 else None

# NEURON 9 aborts in nrn_calc_fast_imem ("Assertion `vec_sav_d' failed") on an MPI rank that owns
# no sections, which happens when there are more ranks than cells (e.g. TESTING with > 4 ranks).
# An unconnected placeholder section on such ranks avoids this; it belongs to no population, so it
# is not part of the recorded spikes, LFP or dipoles. Ranks that own cells are unaffected.
if sum(len(pop.cells) for pop in network.populations.values()) == 0:
	empty_rank_section = h.Section(name='empty_rank_section')

# Setup Extracellular Recording Device
COMM.Barrier()
if stimulate:
	addStimulus(cells_to_stim, cellnums)
COMM.Barrier()


# Run Simulation

tic = time.perf_counter()
# LFPy >= 2.1: the LFP electrode and the current dipole moment are both "probes" whose
# data is filled in by network.simulate(), which now only returns the spikes
probes = []
if rec_LFP:
	LFPelectrode = SomaAsPointElectrode(cell=None, **LFPelectrodeParameters)
	probes.append(LFPelectrode)
if rec_DIPOLES:
	DIPOLEprobe = CurrentDipoleMoment(cell=None)
	probes.append(DIPOLEprobe)

if rec_LFP and not rec_DIPOLES:
	print('Simulating, recording SPIKES and LFP ... ') if RANK==0 else None
elif rec_LFP and rec_DIPOLES:
	print('Simulating, recording SPIKES, LFP, and DIPOLEMOMENTS ... ') if RANK==0 else None
elif rec_DIPOLES and not rec_LFP:
	print('Simulating, recording SPIKES DIPOLEMOMENTS ... ') if RANK==0 else None
elif not rec_LFP and not rec_DIPOLES:
	print('Simulating, recording SPIKES ... ') if RANK==0 else None
SPIKES = network.simulate(probes=probes if probes else None, **simargs)

# Repackage the probe data (summed across RANKs onto RANK 0) in the LFPy 2.0 output
# format used by circuit_functions.py and the L5Circuit Analyses scripts:
#   OUTPUT[0]['imem'] -> (n_contacts, n_timesteps) LFP in mV
#   DIPOLEMOMENT[pop] -> (n_timesteps, 3) dipole moment per population in nA*um
OUTPUT = None
DIPOLEMOMENT = None
if RANK==0:
	if rec_LFP:
		OUTPUT = np.zeros((1,) + LFPelectrode.data.shape, dtype=[('imem', float)])
		OUTPUT[0]['imem'] = LFPelectrode.data['imem']
	if rec_DIPOLES:
		DIPOLEMOMENT = np.zeros(DIPOLEprobe.data.shape[::-1], dtype=[(name, float) for name in network.population_names])
		for name in network.population_names:
			DIPOLEMOMENT[name] = DIPOLEprobe.data[name].T

print('Simulation took ', str((time.perf_counter() - tic_0)/60)[:5], 'minutes') if RANK==0 else None


COMM.Barrier()
if RANK==0:
	tic = time.perf_counter()
	print('Saving simulation output...')
	np.save(os.path.join(OUTPUTPATH,'SPIKES_CircuitSeed'+str(GLOBALSEED)+'StimSeed'+str(STIMSEED)+'.npy'), SPIKES)
	np.save(os.path.join(OUTPUTPATH,'OUTPUT_CircuitSeed'+str(GLOBALSEED)+'StimSeed'+str(STIMSEED)+'.npy'), OUTPUT) if rec_LFP else None
	np.save(os.path.join(OUTPUTPATH,'DIPOLEMOMENT_CircuitSeed'+str(GLOBALSEED)+'StimSeed'+str(STIMSEED)+'.npy'), DIPOLEMOMENT) if rec_DIPOLES else None
	print('Saving simulation took', str((time.perf_counter() - tic_0)/60)[:5], 'minutes')


#===========================================================================
# Plotting
#===========================================================================
if run_circuit_functions and not TESTING: # the plots assume the full-length simulation (2 s transient, stimulus at 6.5 s)

	tstart_plot = 2000
	tstop_plot = tstop

	print('Creating/saving plots...') if RANK==0 else None
	exec(open("circuit_functions.py").read())

#===============
#Final printouts
#===============
if TESTING:
	print('Test complete, switch TESTING to False for full simulation') if RANK==0 else None
elif not TESTING:
	print('Simulation complete') if RANK==0 else None

print('Script completed in ', str((time.perf_counter() - tic_0)/60)[:5], 'minutes') if RANK==0 else None
