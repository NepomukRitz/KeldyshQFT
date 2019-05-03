import re
import numpy as np
import cbin_to_py 


def get_python_class_object(metadata):
	mmmcd = re.compile('matrix.*?matrix.*?matrix.*?complex<double>')
	mmmd = re.compile('matrix.*?matrix.*?matrix.*?double')
	mmscd = re.compile('matrix.*?.*?matrix.*?.*?syma.*?.*?complex<double>')
	mmcd = re.compile('matrix.*?.*?matrix.*?.*?complex<double>')
	mmd =  re.compile('matrix.*?.*?matrix.*?.*?double')
	mscd = re.compile('matrix.*?.*?syma.*?.*?complex<double>')
	mcd = re.compile('matrix.*?.*?complex<double>')
	mi = re.compile('matrix.*?.*?int')
	md = re.compile('matrix.*?.*?double')
	scd = re.compile('syma.*?.*?complex<double>')
	si = re.compile('syma.*?.*?int')
	sd = re.compile('syma.*?.*?double')
	
	liste = [mmmcd, mmmd, mmscd, mmcd, mmd, mscd, mcd, mi, md, scd, si, sd]

	i=0
	for x in liste:
		if x.match(metadata):
			break
		i=i+1
	
	return i
	
	

def load_as_np(foldername, filename):
	filename_extended_meta = foldername + "/" + filename + ".meta"
	with open(filename_extended_meta,'r') as data:
		metadata = data.read() 
	
	dimensions = [int(s) for s in metadata.split() if s.isdigit()]
	
	Ap = np.zeros(dimensions, dtype=complex)

	typenumber = get_python_class_object(metadata)	

	#print(typenumber)
	
	if typenumber==1:
		A = cbin_to_py.MMMd_py()
		A.load_components(foldername,filename)
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				for i2 in range(dimensions[2]):
					for i3 in range(dimensions[3]):
						for i4 in range(dimensions[4]):
							for i5 in range(dimensions[5]):
							 	Ap[i0,i1,i2,i3,i4,i5] = A.component(i0,i1,i2,i3,i4,i5)

	elif typenumber==0:
		A = cbin_to_py.MMMcd_py()
		A.load_components(foldername,filename)
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				for i2 in range(dimensions[2]):
					for i3 in range(dimensions[3]):
						for i4 in range(dimensions[4]):
							for i5 in range(dimensions[5]):
							 	Ap[i0,i1,i2,i3,i4,i5] = A.component(i0,i1,i2,i3,i4,i5)
	elif typenumber==2:
		A = cbin_to_py.MMScd_py()
		A.load_components(foldername,filename)
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				for i2 in range(dimensions[2]):
					for i3 in range(dimensions[3]):
						for i4 in range(dimensions[4]):
							for i5 in range(dimensions[5]):
							 	Ap[i0,i1,i2,i3,i4,i5] = A.component(i0,i1,i2,i3,i4,i5)
	elif typenumber==4:
		A = cbin_to_py.MMd_py()
		A.load_components(foldername,filename)
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				for i2 in range(dimensions[2]):
					for i3 in range(dimensions[3]):
					 	Ap[i0,i1,i2,i3] = A.component(i0,i1,i2,i3)
	elif typenumber==3:
		A = cbin_to_py.MMcd_py()
		A.load_components(foldername,filename)
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				for i2 in range(dimensions[2]):
					for i3 in range(dimensions[3]):
					 	Ap[i0,i1,i2,i3] = A.component(i0,i1,i2,i3)
	elif typenumber==5:
		A = cbin_to_py.MScd_py()
		A.load_components(foldername,filename)
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				for i2 in range(dimensions[2]):
					for i3 in range(dimensions[3]):
					 	Ap[i0,i1,i2,i3] = A.component(i0,i1,i2,i3)
	elif typenumber==6:
		A = cbin_to_py.Mcd_py()
		A.load_components(foldername,filename)
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				Ap[i0,i1] = A.component(i0,i1)
	elif typenumber==7:
		A = cbin_to_py.Mi_py()
		A.load_components(foldername,filename)
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				Ap[i0,i1] = A.component(i0,i1)
	elif typenumber==8:
		A = cbin_to_py.Md_py()
		A.load_components(foldername,filename)
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				Ap[i0,i1] = A.component(i0,i1)
	elif typenumber==9:
		A = cbin_to_py.Scd_py()
		A.load_components(foldername,filename)
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				Ap[i0,i1] = A.component(i0,i1)
	elif typenumber==10:
		A = cbin_to_py.Si_py()
		A.load_components(foldername,filename)
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				Ap[i0,i1] = A.component(i0,i1)
	elif typenumber==11:
		A = cbin_to_py.Sd_py()
		A.load_components(foldername,filename)
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				Ap[i0,i1] = A.component(i0,i1)
	elif typenumber==12:
		print("Error: unknown structure:")
		print(metadata)
	print(typenumber)
	return Ap

def save_as_cbin(A, matrix_type, foldername, variablename):
	dimensions = np.shape(A)
	number_of_dimensions = np.size(dimensions)
	typenumber = 99;
	if matrix_type == "mmmd":
		typenumber = 0;
	if matrix_type == "mmmcd":
		typenumber = 1;
	elif matrix_type =="mmscd":
		typenumber = 2;
	elif matrix_type =="mmd":
		typenumber = 3;
	elif matrix_type =="mmcd":
		typenumber = 4;
	elif matrix_type =="mscd":
		typenumber = 5;
	elif matrix_type =="mcd":
		typenumber = 6;
	elif matrix_type =="mi":
		typenumber = 7;
	elif matrix_type =="md":
		typenumber = 8;
	elif matrix_type =="scd":
		typenumber = 9;
	elif matrix_type =="si":
		typenumber = 10;
	elif matrix_type =="sd":
		typenumber = 11;

	
	if typenumber==0:
		A_c = cbin_to_py.MMMd_py()
		A_c.resize(dimensions[0], dimensions[1], dimensions[2], dimensions[3], dimensions[4], dimensions[5])
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				for i2 in range(dimensions[2]):
					for i3 in range(dimensions[3]):
						for i4 in range(dimensions[4]):
							for i5 in range(dimensions[5]):
								A_c.set_component(i0,i1,i2,i3,i4,i5,A[i0,i1,i2,i3,i4,i5])
	if typenumber==1:
		A_c = cbin_to_py.MMMcd_py()
		A_c.resize(dimensions[0], dimensions[1], dimensions[2], dimensions[3], dimensions[4], dimensions[5])
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				for i2 in range(dimensions[2]):
					for i3 in range(dimensions[3]):
						for i4 in range(dimensions[4]):
							for i5 in range(dimensions[5]):
								A_c.set_component(i0,i1,i2,i3,i4,i5,A[i0,i1,i2,i3,i4,i5])
	if typenumber==2:
		A_c = cbin_to_py.MMScd_py()
		A_c.resize(dimensions[0], dimensions[1], dimensions[2], dimensions[3], dimensions[4])
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				for i2 in range(dimensions[2]):
					for i3 in range(dimensions[3]):
						for i4 in range(dimensions[4]):
							for i5 in range(i4+1):
								A_c.set_component(i0,i1,i2,i3,i4,i5,A[i0,i1,i2,i3,i4,i5])
	if typenumber==3:
		A_c = cbin_to_py.MMd_py()
		A_c.resize(dimensions[0], dimensions[1], dimensions[2], dimensions[3])
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				for i2 in range(dimensions[2]):
					for i3 in range(dimensions[3]):
						A_c.set_component(i0,i1,i2,i3,A[i0,i1,i2,i3])
	if typenumber==4:
		A_c = cbin_to_py.MMcd_py()
		A_c.resize(dimensions[0], dimensions[1], dimensions[2], dimensions[3])
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				for i2 in range(dimensions[2]):
					for i3 in range(dimensions[3]):
						A_c.set_component(i0,i1,i2,i3,A[i0,i1,i2,i3])
	if typenumber==5:
		A_c = cbin_to_py.MScd_py()
		A_c.resize(dimensions[0], dimensions[1],dimensions[2])
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				for i2 in range(dimensions[2]):
					for i3 in range(i2+1):
						A_c.set_component(i0,i1,i2,i3,A[i0,i1,i2,i3])
	if typenumber==6:
		A_c = cbin_to_py.Mcd_py()
		A_c.resize(dimensions[0], dimensions[1])
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				A_c.set_component(i0,i1,A[i0,i1])
	if typenumber==7:
		A_c = cbin_to_py.Mi_py()
		A_c.resize(dimensions[0], dimensions[1])
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				A_c.set_component(i0,i1,A[i0,i1])
	if typenumber==8:
		A_c = cbin_to_py.Md_py()
		A_c.resize(dimensions[0], dimensions[1])
		for i0 in range(dimensions[0]):
			for i1 in range(dimensions[1]):
				A_c.set_component(i0,i1,A[i0,i1])
	if typenumber==9:
		A_c = cbin_to_py.Scd_py()
		A_c.resize(dimensions[0])
		for i0 in range(dimensions[0]):
			for i1 in range(i0+1):
				A_c.set_component(i0,i1,A[i0,i1])
	if typenumber==10:
		A_c = cbin_to_py.Si_py()
		A_c.resize(dimensions[0])
		for i0 in range(dimensions[0]):
			for i1 in range(i0+1):
				A_c.set_component(i0,i1,A[i0,i1])
	if typenumber==11:
		A_c = cbin_to_py.Sd_py()
		A_c.resize(dimensions[0])
		for i0 in range(dimensions[0]):
			for i1 in range(i0+1):
				A_c.set_component(i0,i1,A[i0,i1])

	A_c.save_components(foldername,variablename)

		
