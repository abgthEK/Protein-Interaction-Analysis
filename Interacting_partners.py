# Using inteface residues as input (interface chain known) , finding the interacting partners, node & edge types

f = open ( '1k8w_201706091115_binding_residues1.tbl' , 'r' )
g = open ( '1k8w.pdb' , 'r' )
h = open ( 'chothia.radii' , 'r' )
r = open ('1k8w_out.dat' ,'w' )

# To obtain the Van der waal radii data from text file into the dictionary

Van_radi = {}

for line in h:
	datasplit = line.split()
	Van_radi[datasplit[0]] = ",".join(datasplit[1:])

# To obtain the interacting residue numbers

lines = f.readlines()

glines = g.readlines()

thresh = 0.5

a = 0

# Interface chain known to be 'A'

# Outer loop checks the input file line by line and selects the 1st residue that needs to be checked for interaction

for i in range(0, len(lines)):
	line = lines[i]
	data = line.split()
	if line[0] != '_' and data[0] == 'A':      
		x = data[1]
		for j in range(0, len(lines)):			# This loop is to select the residue to which the initial one needs to be checked 
			a = 0
			xline = lines[j]
			data_in = xline.split()
			if line[0] != '_' and data[0] == 'A':
				if data[0] == 'A' and data_in[1] != x:
					y = data_in[1]
					for k in range(0, len(glines)):  # This loop opens the .pdb file to obtain the coordinates of the individual atoms in each residue for distance calculation
						z = 0
						yline = glines[k]
						gdata = yline.split()
						if gdata[0] == 'ATOM' and gdata[5] == x:
							m = float(gdata[6])	# These are the coordinate values
							n = float(gdata[7])
							o = float(gdata[8])
							res_name1 = gdata[3]
							atom_name1 = gdata[2]
							res_name_dict_list1 = [gdata[3], '_', gdata[2]]	# Joins the residue and atom name to suit our van der waal radii dictionary
							res_name_dict1 = ''.join(res_name_dict_list1)
							for l in range(0, len(glines)):		# This loop obtain the coordinates for the atoms of the interacting partner residue
								zline = glines[l]
								ldata = zline.split()
								if ldata[0] == 'ATOM' and ldata[5] == y:
									t = float(ldata[6])
									u = float(ldata[7])
									v = float(ldata[8])
									res_name2 = ldata[3]
									atom_name2 = ldata[2]	
									res_name_dict_list2 = [ldata[3], '_', ldata[2]]
									res_name_dict2 = ''.join(res_name_dict_list2)
									euc_distance = (((m)-(t))**2+((n)-(u))**2+((o)-(v))**2)**(0.5) # Calculating euclidean distance
									thresh_distance = float(Van_radi[res_name_dict1]) + float(Van_radi[res_name_dict2]) + thresh	# Calculating threshold distance based on our definition
									if euc_distance < thresh_distance:	# To check our interaction condition
										a = a + 1
					for k in range(0, len(glines)):
						z = 0
						yline = glines[k]
						gdata = yline.split()
						if gdata[0] == 'ATOM' and gdata[5] == x:
							res_name1 = gdata[3]
							atom_name1 = gdata[2]
							for l in range(0, len(glines)):
								zline = glines[l]
								ldata = zline.split()
								if ldata[0] == 'ATOM' and ldata[5] == y:
									res_name2 = ldata[3]
									atom_name2 = ldata[2]	
								if a >=1:		# As per our condition for interaction
									r.write(x)
									r.write(' ')
									r.write(y)																						
									r.write(' ')
									if int(x)-int(y) == 1 or int(x)-int(y)== -1:	# To find out the edge type of each interaction pairs
										r.write('1')
										r.write('\n')
									else:
										r.write('2')
										r.write('\n')
									break
							break


					