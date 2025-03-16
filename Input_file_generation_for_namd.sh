#!/bin/bash

####
file_name=$1
pdb_name=$(basename "$file_name" .pdb)
mkdir $pdb_name
cp $file_name ./$pdb_name
cd $pdb_name

dir=/home/vvithurshan/0014/vvarenthirarajah/Documents/Parameter
#dir=/home/vvithurshan/Documents/topology
para="parameters $dir/par_all36_carb.prm
parameters $dir/par_all36_cgenff.prm
parameters $dir/par_all36_lipid.prm
parameters $dir/par_all36m_prot.prm
parameters $dir/par_all36_na.prm
parameters $dir/par_interface.prm
parameters $dir/toppar_water_ions_namd.str"

####


cat <<EOF > psf_script.py

##################
### PSF MAKER ####

import MDAnalysis as mda
u = mda.Universe("$file_name")
u = u.select_atoms('protein')
chain_ids = [i.segid for i in u.segments]
#chain_ids = ['A']


# Open the file in write mode
with open("psf_maker.tcl", "w") as f:

    # Write the lines to the file
    f.write(f"mol load pdb ${file_name}\n")

    for chain in chain_ids:
    	f.write(f"set {chain} [atomselect top \"chain {chain} and protein\"]\n")
    	f.write(f"\${chain} writepdb {chain}.pdb\n")

    f.write("resetpsf\n")
    f.write("package require psfgen\n")
    f.write("topology $dir/top_all36_prot.rtf\n")
    f.write("pdbalias residue HIS HSD\n")

    for chain in chain_ids:
	    f.write(f"segment {chain} {{pdb {chain}.pdb}}")
	    f.write("\n")
	    f.write(f"coordpdb {chain}.pdb {chain}\n")

    f.write("guesscoord\n")
    output = "${file_name}".split('.')[0]
    f.write(f"writepsf dry_{output}.psf\n")
    f.write(f"writepdb dry_{output}.pdb\n")

    ## moving processed_files

    #f.write("mkdir processed_files\n")

    #for chain in chain_ids:
    	#f.write(f"mv {chain}.pdb ./processed_files\n")

# Close the file
f.close()


##################

EOF

python3 psf_script.py

vmd -dispdev text -eofexit < psf_maker.tcl > psf_maker.log


### Solvation

cat <<EOF > solvate.tcl
package require solvate
solvate dry_${pdb_name}.psf dry_${pdb_name}.pdb -o solvate -s WT -x 15 -y 15 -z 15 +x 15 +y 15 +z 15 -b 2.4
EOF

vmd -dispdev text -eofexit < solvate.tcl > solvate_c.log


### Ionize

cat <<EOF > ionize.tcl
package require autoionize
autoionize -psf solvate.psf -pdb solvate.pdb -sc 0.15 -o ionized
EOF

#vmd -dispdev text -eofexit < ionize.tcl > ionize.log

# cp ionized.psf ./PSF_${pdb_name}.psf
# cp ionized.pdb ./PDB_${pdb_name}.pdb
cp dry_${pdb_name}.psf ./PSF_${pdb_name}.psf
cp dry_${pdb_name}.pdb ./PDB_${pdb_name}.pdb
# mv solvate.* ./processed_files
# mv ionized.* ./processed_files
# mv pdb* ./processed_files
# mv psf* ./processed_files

### Harmonic Restraint

cat <<EOF > restraint.tcl
mol load psf PSF_${pdb_name}.psf pdb PDB_${pdb_name}.pdb
mkdir restraint
set all [atomselect top all]
\$all set beta 0
set protein_1 [atomselect top "protein and not water and not ions"]
\$protein_1 set beta 1
\$all writepdb protein_1.pdb
set all [atomselect top all]
\$all set beta 0
set protein_noh_2 [atomselect top "protein and noh and not water and not ions"]
\$protein_noh_2 set beta 1
\$all writepdb protein_noh_2.pdb
set all [atomselect top all]
\$all set beta 0
set backbone_3 [atomselect top "backbone"]
\$backbone_3 set beta 1
\$all writepdb backbone_3.pdb
set all [atomselect top all]
\$all set beta 0
set sam [atomselect top "(chain A and resid 115 to 179) or (chain B and resid 115 to 179)"]
\$sam set beta 1
\$all writepdb sam.pdb

mv protein_1.pdb ./restraint
mv protein_noh_2.pdb ./restraint
mv backbone_3.pdb ./restraint
mv sam.pdb ./restraint 


EOF

vmd -dispdev text -eofexit < restraint.tcl > restraint.log


### BOX Size

cat <<EOF > boxsize.tcl
# mol load psf PSF_${pdb_name}.psf pdb PDB_${pdb_name}.pdb
mol load psf solvate.psf pdb solvate.pdb
proc get_cell {{molid top}} {
    set outfile [open "boxsize.txt" w]
    set all [atomselect \$molid all]
    set minmax [measure minmax \$all]
    set vec [vecsub [lindex \$minmax 1] [lindex \$minmax 0]]
    puts \$outfile "cellBasisVector1 [lindex \$vec 0] 0 0"
    puts \$outfile "cellBasisVector2 0 [lindex \$vec 1] 0"
    puts \$outfile "cellBasisVector3 0 0 [lindex \$vec 2]"
    set center [measure center \$all]
    puts \$outfile "cellOrigin \$center"
    \$all delete
    close \$outfile
}
get_cell 

EOF

vmd -dispdev text -eofexit < boxsize.tcl > boxsize.log

boxsize_contents=$(cat boxsize.txt)

########### Simulation #################

### minimization
cat <<EOF > Minimization.conf


############################################################################
##START HERE###
##Simulation Template##
# Simulation conditions
coordinates PDB_${pdb_name}.pdb
structure PSF_${pdb_name}.psf

# Simulation conditions
temperature 0


CUDASOAintegrate off
# Harmonic constraints

constraints on
consref ./restraint/protein_1.pdb
conskfile ./restraint/protein_1.pdb
constraintScaling 200
consexp 2
conskcol B


# Output Parameters
set output Minimization
binaryoutput no
outputname \$output
outputenergies 100
outputtiming 100
outputpressure 100
binaryrestart yes
dcdfile \$output.dcd
dcdfreq 1000
XSTFreq 1000
restartfreq 1000
restartname \$output.restart


# Thermostat Parameters
langevin on
langevintemp 0
langevinHydrogen    off
langevindamping 1

# Barostat Parameters

usegrouppressure yes
useflexiblecell no
useConstantArea no
langevinpiston off
langevinpistontarget 1
langevinpistonperiod 200
langevinpistondecay 100
langevinpistontemp 300

# Integrator Parameters

timestep 2
firstTimestep 0
fullElectFrequency 2
nonbondedfreq 1
stepspercycle 10

# Force Field Parameters

paratypecharmm on

$para


exclude scaled1-4
1-4scaling 1.0
rigidbonds all
cutoff 12.0
pairlistdist 14.0
switching on
switchdist 10.0
PME on
PMEGridspacing 1
wrapAll on
wrapWater on



#boundary
$boxsize_contents


# Script

minimize 5000


###########
EOF
cat <<EOF > Minimization2.conf


############################################################################
##START HERE###
##Simulation Template##
# Simulation conditions
coordinates PDB_${pdb_name}.pdb
structure PSF_${pdb_name}.psf

set input Minimization
binCoordinates \$input.restart.coor
#binVelocities \$input.restart.vel
extendedSystem \$input.restart.xsc

# Simulation conditions
temperature 0

CUDASOAintegrate off

# Harmonic constraints

constraints on
consref ./restraint/sam.pdb
conskfile ./restraint/sam.pdb
constraintScaling 100
consexp 2
conskcol B


# Output Parameters
set output Minimization2
binaryoutput no
outputname \$output
outputenergies 100
outputtiming 100
outputpressure 100
binaryrestart yes
dcdfile \$output.dcd
dcdfreq 1000
XSTFreq 1000
restartfreq 1000
restartname \$output.restart


# Thermostat Parameters
langevin on
langevintemp 0
langevinHydrogen    off
langevindamping 1

# Barostat Parameters

usegrouppressure yes
useflexiblecell no
useConstantArea no
langevinpiston off
langevinpistontarget 1
langevinpistonperiod 200
langevinpistondecay 100
langevinpistontemp 300

# Integrator Parameters

timestep 2
firstTimestep 0
fullElectFrequency 2
nonbondedfreq 1
stepspercycle 10

# Force Field Parameters

paratypecharmm on

set dir dirr

$para
exclude scaled1-4
1-4scaling 1.0
rigidbonds all
cutoff 12.0
pairlistdist 14.0
switching on
switchdist 10.0
PME on
PMEGridspacing 1
wrapAll on
wrapWater on


#boundary

#cellBasisVector1 53.70100021362305 0 0
#cellBasisVector2 0 46.17199897766113 0
#cellBasisVector3 0 0 41.00600051879883
#cellOrigin -1.5399737358093262 -0.1962127536535263 1.0410041809082031


# Script

minimize 5000


###############

EOF

cat <<EOF > Minimization3.conf

############################################################################
##START HERE###
##Simulation Template##
# Simulation conditions
coordinates PDB_${pdb_name}.pdb
structure PSF_${pdb_name}.psf

set input Minimization2
binCoordinates \$input.restart.coor
#binVelocities \$input.restart.vel
extendedSystem \$input.restart.xsc

# Simulation conditions
temperature 0

CUDASOAintegrate off

# Harmonic constraints

constraints on
consref ./restraint/sam.pdb
conskfile ./restraint/sam.pdb
constraintScaling 100
consexp 2
conskcol B


# Output Parameters
set output Minimization3
binaryoutput no
outputname \$output
outputenergies 100
outputtiming 100
outputpressure 100
binaryrestart yes
dcdfile \$output.dcd
dcdfreq 1000
XSTFreq 1000
restartfreq 1000
restartname \$output.restart


# Thermostat Parameters
langevin on
langevintemp 0
langevinHydrogen    off
langevindamping 1

# Barostat Parameters

usegrouppressure yes
useflexiblecell no
useConstantArea no
langevinpiston off
langevinpistontarget 1
langevinpistonperiod 200
langevinpistondecay 100
langevinpistontemp 300

# Integrator Parameters

timestep 2
firstTimestep 0
fullElectFrequency 2
nonbondedfreq 1
stepspercycle 10

# Force Field Parameters

paratypecharmm on
#set dir /home/vithurshan/2021/toppar_c36_jul20


$para

exclude scaled1-4
1-4scaling 1.0
rigidbonds all
cutoff 12.0
pairlistdist 14.0
switching on
switchdist 10.0
PME on
PMEGridspacing 1
wrapAll on
wrapWater on

# Script

minimize 10000


#######################

EOF

cat <<EOF > Minimization4.conf

############################################################################
##START HERE###
##Simulation Template##
# Simulation conditions
coordinates $pdb_name.pdb
structure $pdb_name.psf

set input Minimization3
binCoordinates \$input.restart.coor
#binVelocities \$input.restart.vel
extendedSystem \$input.restart.xsc

# Simulation conditions
temperature 0

CUDASOAintegrate off

# Harmonic constraints

constraints on
consref ./restraint/sam.pdb
conskfile ./restraint/sam.pdb
constraintScaling 100
consexp 2
conskcol B


# Output Parameters
set output Minimization4
binaryoutput no
outputname \$output
outputenergies 100
outputtiming 100
outputpressure 100
binaryrestart yes
dcdfile \$output.dcd
dcdfreq 1000
XSTFreq 1000
restartfreq 1000
restartname \$output.restart


# Thermostat Parameters
langevin on
langevintemp 0
langevinHydrogen    off
langevindamping 1

# Barostat Parameters

usegrouppressure yes
useflexiblecell no
useConstantArea no
langevinpiston off
langevinpistontarget 1
langevinpistonperiod 200
langevinpistondecay 100
langevinpistontemp 300

# Integrator Parameters

timestep 2
firstTimestep 0
fullElectFrequency 2
nonbondedfreq 1
stepspercycle 10

# Force Field Parameters

paratypecharmm on
#set dir /home/vithurshan/2021/toppar_c36_jul20


$para

exclude scaled1-4
1-4scaling 1.0
rigidbonds all
cutoff 12.0
pairlistdist 14.0
switching on
switchdist 10.0
PME on
PMEGridspacing 1
wrapAll on
wrapWater on

#boundary

# Script

minimize 20000
reinitvels 0 
# do it after minimization

#################

EOF

cat <<EOF > Annealing.conf

############################################################################
##START HERE###
##Simulation Template##
# Simulation conditions
coordinates PDB_${pdb_name}.pdb
structure PSF_${pdb_name}.psf

set input Minimization4
binCoordinates \$input.restart.coor
#binVelocities \$inputname.restart.coor
extendedSystem \$input.restart.xsc

# Simulation conditions
temperature 0

CUDASOAintegrate on

# Harmonic constraints

constraints on
consref ./restraint/sam.pdb
conskfile ./restraint/sam.pdb
constraintScaling 100
consexp 2
conskcol B

# Output Parameters
set output Annealing

binaryoutput no
outputname \$output
outputenergies 500
outputtiming 500
outputpressure 500
binaryrestart yes
dcdfile \$output.dcd
dcdfreq 10000
XSTFreq 1000
restartfreq 1000
restartname \$output.restart


# Thermostat Parameters
langevin on
langevintemp 0
langevinHydrogen    off
langevindamping 1

# Barostat Parameters


langevinpiston off
usegrouppressure yes
useflexiblecell no
useConstantArea no
langevinpistontarget 1
langevinpistonperiod 200
langevinpistondecay 100
langevinpistontemp 300

# Integrator Parameters

timestep 2
firstTimestep 0
fullElectFrequency 2
nonbondedfreq 1
stepspercycle 10

# Force Field Parameters

paratypecharmm on
#set dir /home/vithurshan/2021/toppar_c36_jul20


$para

exclude scaled1-4
1-4scaling 1.0
rigidbonds all
cutoff 12.0
pairlistdist 14.0
switching on
switchdist 10.0
PME on
PMEGridspacing 1
wrapAll on
wrapWater on

# Script
set Temp 1000
set barostat 0
set nSteps  500
for {set t 0} {\$t <= \$Temp} {incr t} {run \$nSteps;langevintemp \$t;if {\$barostat} {langevinpistontemp \$t}}


####################

EOF

## Equilibration
cat <<EOF > Equilibration.conf

############################################################################
##START HERE###
##Simulation Template##
# Simulation conditions
coordinates PDB_${pdb_name}.pdb
structure PSF_${pdb_name}.psf

set input Annealing

binCoordinates \$input.restart.coor
binVelocities \$input.restart.vel
extendedSystem \$input.restart.xsc

# Simulation conditions
#temperature 60

CUDASOAintegrate on

# Harmonic constraints

constraints on
consref ./restraint/sam.pdb
conskfile ./restraint/sam.pdb
constraintScaling 100
consexp 2
conskcol B

# Output Parameters
set output MD_0
binaryoutput no
outputname \$output
outputenergies 100
outputtiming 100
outputpressure 100
binaryrestart yes
dcdfile \$output.dcd
dcdfreq 1000
XSTFreq 1000
restartfreq 1000
restartname \$output.restart


# Thermostat Parameters
langevin on
langevintemp 300
langevinHydrogen    off
langevindamping 1

# Barostat Parameters


langevinpiston on
usegrouppressure yes
useflexiblecell no
useConstantArea no
langevinpistontarget 1
langevinpistonperiod 200
langevinpistondecay 100
langevinpistontemp 300

# Integrator Parameters

timestep 2
firstTimestep 0
fullElectFrequency 2
nonbondedfreq 1
stepspercycle 10

# Force Field Parameters

paratypecharmm on
#set dir /home/vithurshan/2021/toppar_c36_jul20


$para

exclude scaled1-4
1-4scaling 1.0
rigidbonds all
cutoff 12.0
pairlistdist 14.0
switching on
switchdist 10.0
PME on
PMEGridspacing 1
wrapAll on
wrapWater on


run 1000000 
## for 2 ns


###################
EOF

for((MD=1;MD<=10;MD++))
do
## Production
cat <<EOF > MD$MD.conf

############################################################################
##START HERE###
##Simulation Template##
# Simulation conditions
coordinates PDB_${pdb_name}.pdb
structure PSF_${pdb_name}.psf

set input MD_$((MD-1))

binCoordinates \$input.restart.coor
binVelocities \$input.restart.vel
extendedSystem \$input.restart.xsc

# Simulation conditions
#temperature 0

CUDASOAintegrate on

# Harmonic constraints

constraints on
consref ./restraint/sam.pdb
conskfile ./restraint/sam.pdb
constraintScaling 100
consexp 2
conskcol B

# Output Parameters

set output MD_$((MD))

binaryoutput no
outputname \$output
outputenergies 5000
outputtiming 5000
outputpressure 5000
binaryrestart yes
dcdfile \$output.dcd
dcdfreq 5000
XSTFreq 1000
restartfreq 5000
restartname \$output.restart


# Thermostat Parameters
langevin on
langevintemp 300
langevinHydrogen    off
langevindamping 1

# Barostat Parameters

usegrouppressure yes
useflexiblecell no
useConstantArea no
langevinpiston off
langevinpistontarget 1
langevinpistonperiod 200
langevinpistondecay 100
langevinpistontemp 300

# Integrator Parameters

timestep 2
firstTimestep 0
fullElectFrequency 2
nonbondedfreq 1
stepspercycle 10

# Force Field Parameters

paratypecharmm on
$para

exclude scaled1-4
1-4scaling 1.0
rigidbonds all
cutoff 12.0
pairlistdist 14.0
switching on
switchdist 10.0
PME on
PMEGridspacing 1
wrapAll on
wrapWater on

run 5000000


##################3

EOF
done

cat <<EOF >energy.sh
#!/bin/bash

Work=Plots_Auto
if [ -d "$Work" ]; then rm -Rf $Work; fi
mkdir Plots_Auto
cp *.log ./Plots_Auto/plot.log
cd Plots_Auto
file_name=*.log
grep "ENERGY: " $file_name > Energy_temp.dat
grep -v -e "ABSOLUTE" -e "RELATIVE" Energy_temp.dat > Energy.dat
rm Energy_temp.dat
cat Energy.dat |awk '{print $2 " "$11}' > Kinetic.dat
cp Kinetic.dat ./Kinetic.agr
cat Energy.dat |awk '{print $2 " "$12}' > Total.dat
cp Total.dat ./Total.agr
cat Energy.dat |awk '{print $2 " "$13}' > Temp.dat
cp Temp.dat ./Temp.agr
cat Energy.dat |awk '{print $2 " "$14}' > Pot.dat
cp Pot.dat ./Pot.agr
cat Energy.dat |awk '{print $2 " "$17}' > Pressure.dat
cp Pressure.dat ./Pressure.agr
cat Energy.dat |awk '{print $2 " "$19}' > Volume.dat
cp Volume.dat ./Volume.agr
echo Done;

EOF



cat <<EOF > bash_for_gpu.sh

ppn=4
namd3=/home/vvithurshan/0014/vvarenthirarajah/Documents/NAMD_3.0alpha12_Linux-x86_64-multicore-CUDA/namd3
MD_Times=10
pe=240-241
gpu_number=0,1


echo "Minimization";
\$namd3 +p\${ppn} +idlepoll +setcpuaffinity   +devices 0  +pemap \$pe   Minimization.conf > Minimization1.log

echo "Minimization2";
\$namd3 +p\${ppn} +idlepoll +setcpuaffinity   +devices 0  +pemap \$pe   Minimization2.conf > Minimization2.log

echo "Minimization3";
\$namd3 +p\${ppn} +idlepoll +setcpuaffinity   +devices 0  +pemap \$pe   Minimization3.conf > Minimization3.log

echo "Minimization4";
\$namd3 +p\${ppn} +idlepoll +setcpuaffinity   +devices 0  +pemap \$pe  Minimization4.conf > Minimization4.log

for (( min=1; min<=4; min++))
do
	mkdir min_\$min
	mv Minimization\$min.log ./min_\$min
	cd min_\$min
	source energy.sh
	cd ../../
done

echo "Annealing";
\$namd3 +p\${ppn} +idlepoll +setcpuaffinity   +devices \${gpu_number}  +pemap \$pe   Annealing.conf > Annealing.log

echo "Equilibration";
\$namd3 +p\${ppn} +idlepoll +setcpuaffinity   +devices \${gpu_number}  +pemap \$pe  Equilibration.conf > Equilibration.log

for (( X=1; X<=\${MD_Times}; X++))
do
	echo "MD\${X} is running";
	\$namd3 +p\${ppn} +idlepoll +setcpuaffinity   \${gpu_number}  +pemap \$pe MD\${X}.conf > MD\${X}.log
done

echo "Done";
EOF
