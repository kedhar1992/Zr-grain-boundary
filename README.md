# Zr-grain-boundary

#Kedharnath & Bharat 
#Zirconium grain boundary

#Enter 3 orthogonal vectors
#Grain 1
variable f1 equal 0
variable f2 equal 0 
variable f3 equal 1

variable g1 equal 0
variable g2 equal -3 
variable g3 equal 0

variable h1 equal 3
variable h2 equal 0 
variable h3 equal 0

#Grain 2
variable f4 equal 0
variable f5 equal 1
variable f6 equal 0

variable g4 equal 0
variable g5 equal 0 
variable g6 equal 3

variable h4 equal ${h1}
variable h5 equal ${h2}
variable h6 equal ${h3}

#Grain 3
variable f7 equal ${f1}
variable f8 equal ${f2}
variable f9 equal ${f3}

variable g7 equal ${g1}
variable g8 equal ${g2}
variable g9 equal ${g3}

variable h7 equal ${h1}
variable h8 equal ${h2}
variable h9 equal ${h3}

#Constants
variable tempo equal 1000000000
variable latticeconst equal 3.234          # c = 5.168
variable caratio equal 1.598021027
variable minimumenergy equal -6.635
variable pairstyle index "eam/fs"
variable paircoeff index "* * Zr_3.eam.fs Zr Zr Zr"

#Desired filename for LAMMPS to create and dump atoms to
variable filename index "minimum_energy"

# Calculate lattice parameters  
variable xlattice equal ${latticeconst}
variable ylattice equal sqrt(3)
variable zlattice equal ${caratio}

variable fmag equal sqrt((${f1})^2+(${f2})^2+(${f3})^2)
variable gmag equal sqrt((${g1})^2+(${g2})^2+(${g3})^2)
variable hmag equal sqrt((${h1})^2+(${h2})^2+(${h3})^2)

#For Grain 1
variable fx equal ${f1}/${fmag}
variable fy equal ${g1}/${gmag}
variable fz equal ${h1}/${hmag}
variable gx equal ${f2}*${ylattice}/${fmag}
variable gy equal ${g2}*${ylattice}/${gmag}
variable gz equal ${h2}*${ylattice}/${hmag}
variable hx equal ${f3}*${zlattice}/${fmag}
variable hy equal ${g3}*${zlattice}/${gmag}
variable hz equal ${h3}*${zlattice}/${hmag}

#For Grain 2
variable fx2 equal ${f4}/${fmag}
variable fy2 equal ${g4}/${gmag}
variable fz2 equal ${h4}/${hmag}
variable gx2 equal ${f5}*${ylattice}/${fmag}
variable gy2 equal ${g5}*${ylattice}/${gmag}
variable gz2 equal ${h5}*${ylattice}/${hmag}
variable hx2 equal ${f6}*${zlattice}/${fmag}
variable hy2 equal ${g6}*${zlattice}/${gmag}
variable hz2 equal ${h6}*${zlattice}/${hmag}

#For Grain 3
variable fx3 equal ${f7}/${fmag}
variable fy3 equal ${g7}/${gmag}
variable fz3 equal ${h7}/${hmag}
variable gx3 equal ${f8}*${ylattice}/${fmag}
variable gy3 equal ${g8}*${ylattice}/${gmag}
variable gz3 equal ${h8}*${ylattice}/${hmag}
variable hx3 equal ${f9}*${zlattice}/${fmag}
variable hy3 equal ${g9}*${zlattice}/${gmag}
variable hz3 equal ${h9}*${zlattice}/${hmag}

#Simulation parameters
variable etol equal 1.0e-25 
variable ftol equal 1.0e-25 
variable maxiter equal 10000
variable maxeval equal 100000
variable overlapboth equal 2 
variable counter equal 0 

#Calculate simulation box size

variable xlen1 equal ${xlattice}*sqrt((${f1})^2+(${f2})^2+(${f3})^2)
variable ylen1 equal ${xlattice}*sqrt((${g1})^2+(${g2})^2+(${g3})^2)
variable zlen1 equal ${xlattice}*sqrt((${h1})^2+(${h2})^2+(${h3})^2)

variable xsize equal sqrt((${f1})^2+(${f2})^2+(${f3})^2)
variable ysize equal sqrt((${g1})^2+(${g2})^2+(${g3})^2)  
variable zsize equal sqrt((${h1})^2+(${h2})^2+(${h3})^2)

variable xlen equal "v_xlen1*12"
variable ylen equal "v_ylen1*12"
variable zlen equal "v_zlen1*10"

variable grn1 equal "v_ylen-(v_ylen/4)"
variable grn2 equal "v_ylen/4"

print "xlen = ${xlen}, ylen = ${ylen}, zlen = ${zlen}"

# Implement overlap criterion
variable overlapinc equal 545

####################Begin Loop Structure###################################
#This will allow multiple configurations as defined by the
#Xtranslations,Ztranslations, and Overlap distance

variable inc equal "v_latticeconst / 6"
variable xinc equal "2*floor(v_xlen1 / v_inc)"
variable zinc equal "2*floor(v_zlen1 / v_inc)"

label loopa 
variable a loop ${xinc}
print "value of a = ${a}"
variable tx equal "(v_a - 1)/v_xinc*v_xsize" 
print "trans in x = ${tx}"
    label loopb 
    variable b loop ${zinc}
    print "value of b = ${b}"
    variable tz equal "(v_b - 1)/v_zinc*v_zsize" 
    print "trans in z = ${tz}"
        label loopd 
        variable d loop ${overlapboth} 
            label loopc 
            variable c loop ${overlapinc} 
            variable overlapdist equal ".500 + 0.005 * (v_c-1)" 

#Calculate counter for each loop
            variable ctemp equal ${counter}+1 
            variable counter equal ${ctemp} 
            variable ctemp delete 
            print "Counter: ${counter}" 
#Create directory to dump atoms
shell mkdir ${filename}
            
#####################Initiallize Simulation################################

            clear
            units metal
            dimension	3
            boundary	p	p	p
            atom_style	atomic
            atom_modify map array

#######################Create Lattice Structure############################

            region box block 0 ${xlen} 0 ${ylen} 0 ${zlen} units box
            create_box	3 box
	
            lattice custom ${xlattice} a1 ${fx3} ${fy3} ${fz3} a2 ${gx3} ${gy3} ${gz3} a3 ${hx3} ${hy3} ${hz3} basis 0.0 0.0 0.0 basis 0.5 0.5 0 basis 0.5 0.83333333 0.5 basis 0 0.33333333 0.5 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
	    region grain1 block INF INF INF ${grn2} INF INF units box
            create_atoms 1 region grain1 basis 1 1 basis 2 1 basis 3 1 basis 4 1  


            lattice custom ${xlattice} a1 ${fx2} ${fy2} ${fz2} a2 ${gx2} ${gy2} ${gz2} a3 ${hx2} ${hy2} ${hz2} basis 0.0 0.0 0.0 basis 0.5 0.5 0 basis 0.5 0.83333333 0.5 basis 0 0.33333333 0.5 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1 
	    region grain2 block INF INF ${grn2} ${grn1} INF INF units box
            create_atoms 2 region grain2 basis 1 2 basis 2 2 basis 3 2 basis 4 2  


            lattice custom ${xlattice} a1 ${fx} ${fy} ${fz} a2 ${gx} ${gy} ${gz} a3 ${hx} ${hy} ${hz} basis 0.0 0.0 0.0 basis 0.5 0.5 0 basis 0.5 0.83333333 0.5 basis 0 0.33333333 0.5 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
	    region grain3 block INF INF ${grn1} INF INF INF units box
            create_atoms 3 region grain3 basis 1 3 basis 2 3 basis 3 3 basis 4 3 

            group lower type 1
            group middle type 2
	    group upper type 3

#########################Induce Pair Potentials############################

            pair_style	${pairstyle}
            pair_coeff	${paircoeff}
            neighbor 2.0 bin
            neigh_modify delay 10 check yes

###############Translate atoms and delete overlapping atoms################

displace_atoms middle move ${tx} 0 ${tz} units lattice
print "trans in x = ${tx}, trans in z = ${tz} "

if "$c == 1" then "variable atomprev equal 1"
if "$d == 1" then "delete_atoms overlap ${overlapdist} middle upper" #middle atom is deleted
if "$d == 1" then "delete_atoms overlap ${overlapdist} middle lower" 
if "$d == 2" then "delete_atoms overlap ${overlapdist} upper middle"
if "$d == 2" then "delete_atoms overlap ${overlapdist} lower middle"

#Ensure that equivalent configuration has not been tested
#If the same number of atoms is present as before, keep looping
#until overlap distance increases enough to delete more atoms

            variable natoms equal "count(all)" 
            print "Previous: ${atomprev}, Present: ${natoms}" 
            if "${atomprev} == ${natoms}" then "jump in.zr_ebsd loopend" 
 
#############################Computes######################################

            compute csym all centro/atom 12 
            compute eng all pe/atom 
            compute eatoms all reduce sum c_eng 

########################Minimization###############################

            reset_timestep 0 
            thermo 100
            thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms
            min_style cg
            minimize 1e-15 1e-15 ${maxiter} ${maxeval}

# ---------- Run Minimization 2---------------------
# Now allow the box to expand/contract perpendicular to the grain boundary
reset_timestep 0
thermo 10
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms
fix 1 all box/relax y 0.0 vmax 0.001
#Relax only in y directions
min_style cg
minimize ${etol} ${ftol} ${maxiter} ${maxeval}
# ---------- Run Minimization 3---------------------
reset_timestep 0
thermo 10
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms
min_style cg
minimize ${etol} ${ftol} ${maxiter} ${maxeval}
            
#Calculate GB Energy after minimization 

            variable gbarea equal "lx * lz * 2" 
            variable gbe equal "(c_eatoms - (v_minimumenergy * count(all)))/v_gbarea" 
            variable gbemJm2 equal ${gbe}*16021.7733 
            variable gbernd equal round(${gbemJm2}) 
            print "After second minimization:" 
            print "GB energy is ${gbemJm2} mJ/m^2" 

# Store number of atoms for overlap criterion for use in comparing next configurations 
            variable atomprev equal ${natoms}

#####################Dump data into Data file##############################

	  shell cd ${filename}  
          reset_timestep 0 
          timestep 0.001 
velocity all create 20 95812384
fix values all print 1 "${tx} ${tz} ${d} ${c} ${gbemJm2}" append Values.txt title "" 
            
#Store only if it is minimum energy than previous GB
if "${gbemJm2} <= ${tempo}" then "dump  1 all  custom 1000 dump.${tx}_${tz}_${d}_${c}_${gbemJm2} id type x y z c_csym c_eng"
if "${gbemJm2} <= ${tempo}" then "variable tempo equal ${gbemJm2}" 
if "${gbemJm2} <= ${tempo}" then "variable transx equal ${tx}" 
if "${gbemJm2} <= ${tempo}" then "variable transz equal ${tz}"
print "trans in x = ${tx}, trans in z = ${tz} "
            run 1             
            shell cd ..

#####################Finish Loop Structure#################################

            label loopend 
            next c 
            jump in.zr_ebsd loopc 
        variable c delete 
        next d 
        jump in.zr_ebsd loopd 
    variable d delete 
    next b 
    jump in.zr_ebsd loopb 
variable b delete 
next a 
jump in.zr_ebsd loopa 

print "All done" 
