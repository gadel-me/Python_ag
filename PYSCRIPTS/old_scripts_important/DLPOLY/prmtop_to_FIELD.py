#===============import data==============

file_name_1 = raw_input("Select .prmtop file: ")

lines = (open('%s' % file_name_1).readlines())

file_name_2 = raw_input("Select charge input or press 'Enter' to get charges from .prmtop: ")

if file_name_2:
	charges = (open('%s' % file_name_2).readlines())
	case = 0
if not file_name_2:
	case = 1

#================create lists===============

atom_name = []; amber_charges = []; masses = []
atom_type_index = []; nonbonded_parm_index = []; bond_force_constant = []
bond_equil_value = []; angle_force_constant = []; angle_equil_value = []
dihedral_force_constant = []; dihedral_periodicity = []; dihedral_phase = []
scee_scale_factor = []; scnb_scale_factor = []; LJ_acoef = []; LJ_bcoef = []
bonds = []; angles = []; dihedrals = []; amber_atom_type = []

for i in xrange(0,len(lines)):
	if len(lines[i].split()) == 2:
		
		if lines[i].split()[1] == "POINTERS":
			atom_count = int(lines[i+2].split()[0])
			atom_LJ_types = int(lines[i+2].split()[1])
			bond_count = int(lines[i+2].split()[2]) + int(lines[i+2].split()[3])
			angle_count = int(lines[i+2].split()[4]) + int(lines[i+2].split()[5])
			dihedral_count = int(lines[i+2].split()[6]) + int(lines[i+2].split()[7])
		
		if lines[i].split()[1] == "ATOM_NAME":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "CHARGE":
					break
				atom_name.extend(lines[j].split())

		if lines[i].split()[1] == "CHARGE":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "ATOMIC_NUMBER":
					break
				amber_charges.extend(lines[j].split())

		if lines[i].split()[1] == "MASS":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "ATOM_TYPE_INDEX":
					break
				masses.extend(lines[j].split())

		if lines[i].split()[1] == "ATOM_TYPE_INDEX":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "NUMBER_EXCLUDED_ATOMS":
					break
				atom_type_index.extend(lines[j].split())

		if lines[i].split()[1] == "NONBONDED_PARM_INDEX":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "RESIDUE_LABEL":
					break
				nonbonded_parm_index.extend(lines[j].split())

		if lines[i].split()[1] == "BOND_FORCE_CONSTANT":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "BOND_EQUIL_VALUE":
					break
				bond_force_constant.extend(lines[j].split())

		if lines[i].split()[1] == "BOND_EQUIL_VALUE":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "ANGLE_FORCE_CONSTANT":
					break
				bond_equil_value.extend(lines[j].split())

		if lines[i].split()[1] == "ANGLE_FORCE_CONSTANT":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "ANGLE_EQUIL_VALUE":
					break
				angle_force_constant.extend(lines[j].split())

		if lines[i].split()[1] == "ANGLE_EQUIL_VALUE":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "DIHEDRAL_FORCE_CONSTANT":
					break
				angle_equil_value.extend(lines[j].split())

		if lines[i].split()[1] == "DIHEDRAL_FORCE_CONSTANT":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "DIHEDRAL_PERIODICITY":
					break
				dihedral_force_constant.extend(lines[j].split())

		if lines[i].split()[1] == "DIHEDRAL_PERIODICITY":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "DIHEDRAL_PHASE":
					break
				dihedral_periodicity.extend(lines[j].split())

		if lines[i].split()[1] == "DIHEDRAL_PHASE":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "SCEE_SCALE_FACTOR":
					break
				dihedral_phase.extend(lines[j].split())

		if lines[i].split()[1] == "SCEE_SCALE_FACTOR":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "SCNB_SCALE_FACTOR":
					break
				scee_scale_factor.extend(lines[j].split())

		if lines[i].split()[1] == "SCNB_SCALE_FACTOR":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "SOLTY":
					break
				scnb_scale_factor.extend(lines[j].split())

		if lines[i].split()[1] == "LENNARD_JONES_ACOEF":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "LENNARD_JONES_BCOEF":
					break
				LJ_acoef.extend(lines[j].split())

		if lines[i].split()[1] == "LENNARD_JONES_BCOEF":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "BONDS_INC_HYDROGEN":
					break
				LJ_bcoef.extend(lines[j].split())

		if lines[i].split()[1] == "BONDS_INC_HYDROGEN":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "BONDS_WITHOUT_HYDROGEN":
					break
				bonds.extend(lines[j].split())

		if lines[i].split()[1] == "BONDS_WITHOUT_HYDROGEN":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "ANGLES_INC_HYDROGEN":
					break
				bonds.extend(lines[j].split())

		if lines[i].split()[1] == "ANGLES_INC_HYDROGEN":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "ANGLES_WITHOUT_HYDROGEN":
					break
				angles.extend(lines[j].split())

		if lines[i].split()[1] == "ANGLES_WITHOUT_HYDROGEN":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "DIHEDRALS_INC_HYDROGEN":
					break
				angles.extend(lines[j].split())

		if lines[i].split()[1] == "DIHEDRALS_INC_HYDROGEN":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "DIHEDRALS_WITHOUT_HYDROGEN":
					break
				dihedrals.extend(lines[j].split())

		if lines[i].split()[1] == "DIHEDRALS_WITHOUT_HYDROGEN":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "EXCLUDED_ATOMS_LIST":
					break
				dihedrals.extend(lines[j].split())

		if lines[i].split()[1] == "AMBER_ATOM_TYPE":
			for j in xrange(i+2,len(lines)):
				if len(lines[j].split()) == 2 and lines[j].split()[1] == "TREE_CHAIN_CLASSIFICATION":
					break
				amber_atom_type.extend(lines[j].split())

bonds = map(int,bonds); bonds = [bonds[3*i:3*(i+1)] for i in range(len(bonds)/3)]
bonds = sorted(bonds,key=lambda item: item[1]); bonds = sorted(bonds,key=lambda item: item[0])

angles = map(int,angles); angles = [angles[4*i:4*(i+1)] for i in range(len(angles)/4)]
angles = sorted(angles,key=lambda item: item[2]); angles = sorted(angles,key=lambda item: item[1]); angles = sorted(angles,key=lambda item: item[0])

dihedrals = map(int,dihedrals); dihedrals = [dihedrals[5*i:5*(i+1)] for i in range(len(dihedrals)/5)]
dihedrals = sorted(dihedrals,key=lambda item: abs(item[3])); dihedrals = sorted(dihedrals,key=lambda item: abs(item[2]))
dihedrals = sorted(dihedrals,key=lambda item: item[1]); dihedrals = sorted(dihedrals,key=lambda item: item[0])

if case == 0:

	gaussian_charges = []; var = 0
	
	for i in xrange(0,len(charges)):
		if charges[i].startswith(' Charges from ESP fit with hydrogens summed into heavy atoms:'):
			break
		if charges[i].startswith(' Fitting point charges to electrostatic potential'):
			var = i
		if var > 0 and i - var >= 4 and len(charges[i].split()) == 3:
			gaussian_charges.append(charges[i].split()[2])
	
	if len(gaussian_charges) != atom_count:
		print "Something went wrong with the GAUSSIAN charges! Check the format of the GAUSSIAN charge output!"

#===============generate FIELD file================

system_name = raw_input("Set system name: ")
molecule_name = raw_input("Set molecule name: ")

open('FIELD', 'w').write('%s\n' % system_name)
open('FIELD', 'a+').write('UNITS kJ / mol\n')
open('FIELD', 'a+').write('MOLECULES 1\n')
open('FIELD', 'a+').write('%s\n' % molecule_name)
open('FIELD', 'a+').write('NUMMOLS 1\n')
open('FIELD', 'a+').write('ATOMS %d\n' % atom_count)

if case == 1:
	for i in xrange(0,atom_count):
		open('FIELD', 'a+').write('%s\t%8.4f\t%8.4f\t1\t0\n' % (amber_atom_type[i].upper(),float(masses[i]),float(amber_charges[i])/18.222615))

if case == 0:
	for i in xrange(0,atom_count):
		open('FIELD', 'a+').write('%s\t%8.4f\t%8.4f\t1\t0\n' % (amber_atom_type[i].upper(),float(masses[i]),float(gaussian_charges[i])))

open('FIELD', 'a+').write('BONDS %d\n' % bond_count)

for i in xrange(0,len(bonds)):
	open('FIELD', 'a+').write('harm\t%2.d\t%2.d\t' % (1+bonds[i][0]/3,1+bonds[i][1]/3))
	open('FIELD', 'a+').write('%8.5f\t%8.5f\n' % (float(bond_force_constant[bonds[i][2]-1])*4.184*2,float(bond_equil_value[bonds[i][2]-1])))

open('FIELD', 'a+').write('ANGLES %d\n' % angle_count)

for i in xrange(0,len(angles)):
	open('FIELD', 'a+').write('harm\t%2.d\t%2.d\t%2.d\t' % (1+angles[i][0]/3,1+angles[i][1]/3,1+angles[i][2]/3))
	open('FIELD', 'a+').write('%8.4f\t%8.4f\n' % (float(angle_force_constant[angles[i][3]-1])*4.184*2,float(angle_equil_value[angles[i][3]-1])*57.29577951))

open('FIELD', 'a+').write('DIHEDRALS %d\n' % dihedral_count)

for i in xrange(0,len(dihedrals)):

	if dihedrals[i][3] >= 0:
		open('FIELD', 'a+').write('cos\t%2.d\t%2.d\t%2.d\t%2.d\t' % (1+dihedrals[i][0]/3,1+dihedrals[i][1]/3,1+abs(dihedrals[i][2])/3,1+abs(dihedrals[i][3])/3))
		open('FIELD', 'a+').write('%8.4f\t%6.3f\t' % (float(dihedral_force_constant[dihedrals[i][4]-1])*4.184,float(dihedral_phase[dihedrals[i][4]-1])*57.29577951))
		open('FIELD', 'a+').write('%6.3f\t' % (float(dihedral_periodicity[dihedrals[i][4]-1])))
		open('FIELD', 'a+').write('%6.4f\t%6.4f\n' % (1/1.2,1/2.0))

bonds_indices = []

for i in xrange(0,len(bonds)):
	bonds_indices.append([bonds[i][0],bonds[i][1]])

for i in xrange(0,len(dihedrals)):

	if dihedrals[i][3] < 0:

		if dihedrals[i][0] != 0 and dihedrals[i][1] != 0:
			open('FIELD', 'a+').write('cos\t%2.d\t%2.d\t%2.d\t%2.d\t' % (1+abs(dihedrals[i][2])/3,1+dihedrals[i][0]/3,1+dihedrals[i][1]/3,1+abs(dihedrals[i][3])/3))
			open('FIELD', 'a+').write('%8.4f\t%6.3f\t' % (float(dihedral_force_constant[dihedrals[i][4]-1])*4.184,float(180)))
			open('FIELD', 'a+').write('%6.3f\t' % (float(1)))
			open('FIELD', 'a+').write('%6.4f\t%6.4f\n' % (0.0,0.0))

		if dihedrals[i][0] == 0 or dihedrals[i][1] == 0:

			if [abs(dihedrals[i][2]),dihedrals[i][0]] in bonds_indices or [dihedrals[i][0],abs(dihedrals[i][2])] in bonds_indices:
				if [abs(dihedrals[i][2]),dihedrals[i][1]] in bonds_indices or [dihedrals[i][1],abs(dihedrals[i][2])] in bonds_indices:
					if [abs(dihedrals[i][2]),abs(dihedrals[i][3])] in bonds_indices or [abs(dihedrals[i][3]),abs(dihedrals[i][2])] in bonds_indices:
						open('FIELD', 'a+').write('cos\t%2.d\t%2.d\t%2.d\t%2.d\t' % (1+abs(dihedrals[i][2])/3,1+dihedrals[i][0]/3,1+dihedrals[i][1]/3,1+abs(dihedrals[i][3])/3))
						open('FIELD', 'a+').write('%8.4f\t%6.3f\t' % (float(dihedral_force_constant[dihedrals[i][4]-1])*4.184,float(180)))
						open('FIELD', 'a+').write('%6.3f\t' % (float(1)))
						open('FIELD', 'a+').write('%6.4f\t%6.4f\n' % (0.0,0.0))
						continue

			open('FIELD', 'a+').write('cos\t%2.d\t%2.d\t%2.d\t%2.d\t' % (1+dihedrals[i][1]/3,1+abs(dihedrals[i][3])/3,1+abs(dihedrals[i][2])/3,1+dihedrals[i][0]/3))
			open('FIELD', 'a+').write('%8.4f\t%6.3f\t' % (float(dihedral_force_constant[dihedrals[i][4]-1])*4.184,float(180)))
			open('FIELD', 'a+').write('%6.3f\t' % (float(1)))
			open('FIELD', 'a+').write('%6.4f\t%6.4f\n' % (0.0,0.0))
					
unique_atom_type = sorted(list(set(amber_atom_type)))
unique_atom_type_and_index = []

for i in xrange(0,len(unique_atom_type)):
	for j in xrange(0,len(amber_atom_type)):
		if unique_atom_type[i] == amber_atom_type[j]:
			unique_atom_type_and_index.append([unique_atom_type[i],atom_type_index[j]])
			break

open('FIELD', 'a+').write('FINISH\n')
open('FIELD', 'a+').write('VDW %d\n' % int((len(unique_atom_type)**2+len(unique_atom_type))/2))

i_1 = 0; i_2 = 0; a_t_i_1 = 0; a_t_i_2 = 0; n_p_i = 0

for i in xrange(0,len(unique_atom_type_and_index)):
	for j in xrange(i,len(unique_atom_type_and_index)):

		# get atom type indices from atom_type_index
		a_t_i_1 = int(unique_atom_type_and_index[i][1])
		a_t_i_2 = int(unique_atom_type_and_index[j][1])

		# get nonbonden parm indices from nonbonded_parm_index
		n_p_i = int(nonbonded_parm_index[int(atom_LJ_types*(a_t_i_1-1)+a_t_i_2)-1])

		open('FIELD', 'a+').write('%s\t%s\tlj\t' % (unique_atom_type[i].upper(),unique_atom_type[j].upper()))
		
		if float(LJ_acoef[n_p_i-1]) == 0:

			open('FIELD', 'a+').write('%8.5f\t' % 0.0)

		if float(LJ_acoef[n_p_i-1]) != 0:

			open('FIELD', 'a+').write('%8.5f\t' % float( 4.184*((float(LJ_bcoef[n_p_i-1]))**2)/(4*float(LJ_acoef[n_p_i-1]))))

		if float(LJ_bcoef[n_p_i-1]) == 0:

			open('FIELD', 'a+').write('%8.5f\n' % 0.0)

		if float(LJ_bcoef[n_p_i-1]) != 0:

			open('FIELD', 'a+').write('%8.5f\n' % float((float(LJ_acoef[n_p_i-1])/float(LJ_bcoef[n_p_i-1]))**(1/6.0)))

open('FIELD', 'a+').write('CLOSE')
			
