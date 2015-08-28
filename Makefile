# Cluster options
Cluster_Name = master
Cluster_GNode = gnode70
Cluster_Node = node10
Cluster_Nb_Nodes = 01
Cluster_ppn = 1
Cluster_Id_Computation = 40
Cluster_walltime = 300:00:00
Cluster_pvmem = 31gb
Cluster_pvmem_bigmem = 63gb
Cluster_List_Gnodes = 70 71 72 73 74 75 76 77 78
# Cluster_Reseau = IB (Infiniband), TCP (Ethernet)
Cluster_Reseau = TCP
Script_Name = launch_make
# Node_Dir
Node_Dir = /usrtmp/$(USER)
# Cluster_Dir = /data1/$(USER), /data2/$(USER), /nutmp/$(Cluster_Node)/$(USER)
Cluster_Dir  = /data1/$(USER)
Cluster_List_Data = 1 2 3 4 5 6
# LMT_Dir, METIL_Dir, DIC_Dir
ERROR_ESTIMATION_Dir   = https://github.com/fpled/ERROR_ESTIMATION.git
LMT_Dir   = git://gitosis.lmt.ens-cachan.fr/LMTpp-test.git
METIL_Dir = git://gitosis.lmt.ens-cachan.fr/METIL.git
DIC_Dir = git://gitosis.lmt.ens-cachan.fr/dic.git
qsub_Dir = /usr/local/torque-current/bin/
scons_file = SConstruct
scons_file_pgd = SConstruct_pgd
scons_file_homog = SConstruct_homog
# Options for the compilation
debug = 0 
opt = 0
timdavis = 0
nb_pro = 8
machine_arch = ''
# Name of machine
List_Machine_Names = pommard chablis
Machine_Name = pommard
Machine_Name_Results = pommard
# Name of station
Station_Name = stationmeca22
# Name of problem
Pb_Name = error_estimation
Pb_Name_Space = space
Pb_Name_Parameter = parameter

# Default ---------------------------
default: compile
	time ./main
	make move_results

# PGD ---------------------------
all_pgd: compile_pgd
	time ./main_pgd
	make move_results

# HOMOG ---------------------------
all_homog: compile_homog
	time ./main_homog
	make move_results

# Codegen ---------------------------
codegen:
	cd LMT/include/codegen; scons -j 1

codegen_cluster:
	ssh $(Cluster_Name) "cd $(Cluster_Dir)/ERROR_ESTIMATION_$(Cluster_Id_Computation)/; make codegen"

# Compile ---------------------------
compile: codegen
#	export METILPATH=../METIL/MET; ../METIL-install/bin/metil formulation.met
	scons --sconstruct=$(scons_file) -j $(nb_pro) arch=$(machine_arch) debug=$(debug) opt=$(opt) timdavis=$(timdavis)

compile_pgd: codegen
#	export METILPATH=../METIL/MET; ../METIL-install/bin/metil formulation_space.met
#	export METILPATH=../METIL/MET; ../METIL-install/bin/metil formulation_parameter.met
	scons --sconstruct=$(scons_file_pgd) -j $(nb_pro) arch=$(machine_arch) debug=$(debug) opt=$(opt) timdavis=$(timdavis)

compile_homog: codegen
#	export METILPATH=../METIL/MET; ../METIL-install/bin/metil formulation_homog.met
	scons --sconstruct=$(scons_file_homog) -j $(nb_pro) arch=$(machine_arch) debug=$(debug) opt=$(opt) timdavis=$(timdavis)

compile_test: codegen
	cd TESTS; scons --sconstruct=$(scons_file) -j $(nb_pro) arch=$(machine_arch) debug=$(debug) opt=$(opt) timdavis=$(timdavis)

# Test ---------------------------
test_metil:
	export METILPATH=../METIL/MET; ../METIL-install/bin/metil TESTS/test.met

test: compile_test
	cd TESTS; ./test

# Move ---------------------------
move_results:
	-find . -maxdepth 1 -name "*.vtu" -exec mv {} RESULTS \;
	-find . -maxdepth 1 -name "*.pvd" -exec mv {} RESULTS \;
	-find . -maxdepth 1 -name "*.log" -exec mv {} RESULTS \;
	-find . -maxdepth 1 -name "*.out" -exec mv {} RESULTS \;
	-find . -maxdepth 1 -name "*.fig" -exec mv {} RESULTS \;
	-find . -maxdepth 1 -name "*.eps" -exec mv {} RESULTS \;
	-find . -maxdepth 1 -name "*.pdf" -exec mv {} RESULTS \;
	-find . -maxdepth 1 -name "*.png" -exec mv {} RESULTS \;
	-find . -maxdepth 1 -name "*.svg" -exec mv {} RESULTS \;
	-find . -maxdepth 1 -name "*.dvi" -exec mv {} RESULTS \;
	-find . -maxdepth 1 -name "*.tex" -exec mv {} RESULTS \;
	-find . -maxdepth 1 -name "*.aux" -exec mv {} RESULTS \;

# Clean ---------------------------
clean:
	scons -c
	-rm -r build

clean_codegen:
	cd LMT/include/codegen; scons -c

clean_results:
	-rm -fr RESULTS/*

clean_all: clean_codegen clean clean_results
	-find . -maxdepth 1 -name "*.o" -exec rm -fr {} \;

clean_metil:
	-rm -fr ../METIL
	-rm -fr ../METIL-build
	-rm -fr ../METIL-install

clean_lmtpp:
	-rm -fr ../LMTpp

# Clean Cluster---------------------------
clean_cluster:
	ssh $(Cluster_Name) "cd $(Cluster_Dir)/ERROR_ESTIMATION_$(Cluster_Id_Computation)/; make clean"

clean_codegen_cluster:
	ssh $(Cluster_Name) "cd $(Cluster_Dir)/ERROR_ESTIMATION_$(Cluster_Id_Computation)/; make clean_codegen"

clean_all_cluster:
	ssh $(Cluster_Name) "cd $(Cluster_Dir)/ERROR_ESTIMATION_$(Cluster_Id_Computation)/; make clean_all"

# Clean GNode---------------------------
clean_gnode:
	ssh $(Cluster_GNode) "cd $(Node_Dir)/ERROR_ESTIMATION/; make clean"

clean_codegen_gnode:
	ssh $(Cluster_Gnode) "cd $(Node_Dir)/ERROR_ESTIMATION/; make clean_codegen"

clean_all_gnode:
	ssh $(Cluster_GNode) "cd $(Node_Dir)/ERROR_ESTIMATION/; make clean_all"

# Gdb ---------------------------
debug_gdb:
	gdb --args make

# Clone ---------------------------
clone_lmtpp: clean_lmtpp
	-git clone $(LMT_Dir) ../LMTpp

clone_metil: clean_metil
	-git clone $(METIL_Dir) ../METIL
	mkdir ../METIL-build; cd ../METIL-build; cmake -i ../METIL

install_metil:
	cd ../METIL-build; make -j4 install

# Excludes for rsync ---------------------------

exclude = --exclude '*.o' --exclude 'core' --exclude '*.log' --exclude '*.eps' --exclude '*.dvi' --exclude '*.pdf' --exclude '*.tex' --exclude '*.aux' --exclude '*.patch' --exclude '*.kdevelop' --exclude '*.kdev4' --exclude '*.filelist' --exclude '*.pcs' --exclude '*.kdevses' --exclude '*.m' --exclude '*.pyc' --exclude '*.cc' --exclude '*.bz2' --exclude '*.cache' --exclude '*.dblite' --exclude '*.git' --exclude '.sconsign' --exclude 'LMT/.git' --exclude 'Doxyfile' --exclude 'build' --exclude 'TESTS' --exclude 'RESULTS' --exclude 'i686.tok'
exclude_main_formulation = --exclude 'main.cpp' --exclude 'formulation.met'
exclude_Makefile = --exclude 'Makefile'
exclude_LMT = --exclude '*.o' --exclude '*.patch' --exclude '*.git' --exclude 'i686.tok' --exclude '*.cache'
exclude_METIL = --exclude '*.patch' --exclude '*.git' --exclude 'i686.tok' --exclude '*.cache' --exclude '*.bz2'
exclude_generated = --exclude 'i686.tok'

# Machine ---------------------------

lmtpp_to_all_machines:
	$(foreach i, $(List_Machine_Names), \
	rsync -auv /utmp/$(Machine_Name)/$(USER)/LMTpp        $(i):$(Node_Dir) $(exclude_LMT);)

metil_to_all_machines:
	$(foreach i, $(List_Machine_Names), \
	rsync -auv /utmp/$(Machine_Name)/$(USER)/METIL      $(i):$(Node_Dir) $(exclude_METIL);)

error_estimation_to_all_machines:
	$(foreach i, $(List_Machine_Names), \
	rsync -auv /utmp/$(Machine_Name)/$(USER)/ERROR_ESTIMATION        $(i):$(Node_Dir) $(exclude); \
	ssh $(i) "cd $(Node_Dir)/ERROR_ESTIMATION/; rm -fr LMT"; \
	ssh $(i) "cd $(Node_Dir)/ERROR_ESTIMATION/; ln -s $(Node_Dir)/LMTpp LMT";)

main_to_all_machines:
	$(foreach i, $(List_Machine_Names), \
	rsync -auv /utmp/$(Machine_Name)/$(USER)/ERROR_ESTIMATION/main.cpp     $(i):$(Node_Dir)/ERROR_ESTIMATION;)

formulation_to_all_machines:
	$(foreach i, $(List_Machine_Names), \
	rsync -auv /utmp/$(Machine_Name)/$(USER)/ERROR_ESTIMATION/formulation.met     $(i):$(Node_Dir)/ERROR_ESTIMATION;)

install_metil_all_machines:
	$(foreach i, $(List_Machine_Names), \
	ssh $(i) "cd $(Node_Dir)/METIL-build/; make -j4 install";)

# Cluster ---------------------------

all_to_cluster:
	rsync -auv /utmp/$(Machine_Name)/$(USER)/LMTpp        $(Cluster_Name):$(Cluster_Dir) $(exclude_LMT)
	rsync -auv /utmp/$(Machine_Name)/$(USER)/METIL      $(Cluster_Name):$(Cluster_Dir) $(exclude_METIL)
	rsync -auv /utmp/$(Machine_Name)/$(USER)/ERROR_ESTIMATION/*     $(Cluster_Name):$(Cluster_Dir)/ERROR_ESTIMATION_$(Cluster_Id_Computation) $(exclude)
	ssh $(Cluster_Name) "cd $(Cluster_Dir)/ERROR_ESTIMATION_$(Cluster_Id_Computation)/; rm -fr LMT"
	ssh $(Cluster_Name) "cd $(Cluster_Dir)/ERROR_ESTIMATION_$(Cluster_Id_Computation)/; ln -s $(Cluster_Dir)/LMTpp LMT"

cluster:
	ssh $(Cluster_Name) "hostname; $(qsub_Dir)/qsub -k eo -l nodes=$(Cluster_Nb_Nodes):ppn=$(Cluster_ppn):$(Cluster_Reseau),walltime=$(Cluster_walltime),pvmem=$(Cluster_pvmem) -v REP='$(Cluster_Dir)/ERROR_ESTIMATION_$(Cluster_Id_Computation)/' $(Script_Name)"

cluster_node:
	ssh $(Cluster_Name) "hostname; $(qsub_Dir)/qsub -k eo -l host=$(Cluster_Node):ppn=$(Cluster_ppn):$(Cluster_Reseau),walltime=$(Cluster_walltime),pvmem=$(Cluster_pvmem) -v REP='$(Cluster_Dir)/ERROR_ESTIMATION_$(Cluster_Id_Computation)/' $(Script_Name)"

cluster_bigmem:
	ssh $(Cluster_Name) "hostname; $(qsub_Dir)/qsub -k eo -q bigmem -l walltime=$(Cluster_walltime),pvmem=$(Cluster_pvmem_bigmem) -v REP='$(Cluster_Dir)/ERROR_ESTIMATION_$(Cluster_Id_Computation)/' $(Script_Name)"

all_cluster: all_to_cluster cluster

install_metil_cluster:
	ssh $(Cluster_GNode) "cd $(Cluster_Dir)/METIL-build; make -j4 install"

install_metil_all_cluster:
	$(foreach i, $(Cluster_List_Data), \
	ssh $(Cluster_GNode) "cd /data$(i)/$(USER)/METIL-build; make -j4 install";)

results_from_cluster:
	rsync -auv $(Cluster_Name):$(Cluster_Dir)/ERROR_ESTIMATION_$(Cluster_Id_Computation)/*.vtu  /utmp/$(Machine_Name_Results)/$(USER)/RESULTS
	rsync -auv $(Cluster_Name):/u/$(USER)/*.out  /utmp/$(Machine_Name_Results)/$(USER)/RESULTS

# GNode ---------------------------

all_to_gnode:
	rsync -auv /utmp/$(Machine_Name)/$(USER)/LMTpp        $(Cluster_GNode):$(Node_Dir) $(exclude_LMT)
	rsync -auv /utmp/$(Machine_Name)/$(USER)/METIL      $(Cluster_GNode):$(Node_Dir) $(exclude_METIL)
	rsync -auv /utmp/$(Machine_Name)/$(USER)/ERROR_ESTIMATION     $(Cluster_GNode):$(Node_Dir) $(exclude)
	ssh $(Cluster_GNode) "cd $(Node_Dir)/ERROR_ESTIMATION/; rm -fr LMT"
	ssh $(Cluster_GNode) "cd $(Node_Dir)/ERROR_ESTIMATION/; ln -s $(Node_Dir)/LMTpp LMT"

gnode:
	ssh $(Cluster_GNode) "hostname; cd $(Node_Dir)/ERROR_ESTIMATION; make"

all_gnode: all_to_gnode gnode

install_metil_gnode:
	ssh $(Cluster_GNode) "cd $(Node_Dir)/METIL-build; make -j4 install"

install_metil_all_gnodes:
	$(foreach i, $(Cluster_List_Gnodes), \
	ssh gnode$(i) "cd $(Node_Dir)/METIL-build; make -j4 install";)

results_from_gnode:
	rsync -auv $(Cluster_GNode):$(Node_Dir)/ERROR_ESTIMATION/*.vtu  /utmp/$(Machine_Name_Results)/$(USER)/RESULTS/
	rsync -auv $(Cluster_GNode):$(Node_Dir)/ERROR_ESTIMATION/*.pvd  /utmp/$(Machine_Name_Results)/$(USER)/RESULTS/
	rsync -auv $(Cluster_GNode):$(Node_Dir)/ERROR_ESTIMATION/*.log  /utmp/$(Machine_Name_Results)/$(USER)/RESULTS/
