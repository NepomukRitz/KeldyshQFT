
def analyse_channel(l,k,Nf,L,L_structure,pos_feedback,a_dyn,a_stat):
	a_r = np.empty([Nf]) 
	a_i = np.empty([Nf]) 
	for i in range(Nf):
		Li = L_structure[i]
		a_r[i] = 0 
		a_i[i] = 0 
		if(abs(l)<=Li and abs(k)<=Li):
			a_r[i] = np.amax(np.absolute(np.real(a_dyn[i][l+Li,k+Li]))) 
			a_i[i] = np.amax(np.absolute(np.imag(a_dyn[i][l+Li,k+Li]))) 
	a_r[pos_feedback] = np.amax(np.absolute(np.real(a_stat[l+L,k+L]))) 
	a_i[pos_feedback] = 0.0 
	return [a_r,a_i]

def analyse_self(Nff,self):
	s_r = np.empty([Nff])
	s_i = np.empty([Nff])
	for i in range(Nff):
		#s_r[i] = np.amax(np.absolute(np.real(self[0,i,:,:])))
		#s_i[i] = np.amax(np.absolute(np.imag(self[0,i,:,:])))
		s_r[i] = np.real(self[0,i,30,30])
		s_i[i] = np.imag(self[0,i,30,30])
	return [s_r, s_i]

#Analyse P structure:
[apud00_r0,apud00_i0] = analyse_channel(0,0,NfbP0,L0,Lp_structure0,pos_NfbP_2mu0,aPud_dyn0,aPud_stat0)
[apud00_r1,apud00_i1] = analyse_channel(0,0,NfbP1,L1,Lp_structure1,pos_NfbP_2mu1,aPud_dyn1,aPud_stat1)
[apud00_r2,apud00_i2] = analyse_channel(0,0,NfbP2,L2,Lp_structure2,pos_NfbP_2mu2,aPud_dyn2,aPud_stat2)
#[apud10_r1,apud10_i1] = analyse_channel(1,0,NfbP1,L1,Lp_structure1,pos_NfbP_2mu1,aPud_dyn1,aPud_stat1)

[axud00_r0,axud00_i0] = analyse_channel(0,0,NfbX0,L0,Lx_structure0,pos_NfbX_00,aXud_dyn0,aXud_stat0)
[axud00_r1,axud00_i1] = analyse_channel(0,0,NfbX1,L1,Lx_structure1,pos_NfbX_01,aXud_dyn1,aXud_stat1)
[axud00_r2,axud00_i2] = analyse_channel(0,0,NfbX2,L2,Lx_structure2,pos_NfbX_02,aXud_dyn2,aXud_stat2)
#[axud10_r1,axud10_i1] = analyse_channel(1,0,NfbX1,L1,Lx_structure1,pos_NfbX_01,aXud_dyn1,aXud_stat1)

[s_r0,s_i0] = analyse_self(Nff0,ERetu0)	
[s_r1,s_i1] = analyse_self(Nff1,ERetu1)	
[s_r2,s_i2] = analyse_self(Nff2,ERetu2)	
#[s_r2,s_i2] = analyse_self(Nff2,ERetu2)	
	

fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
#axis.plot(wbP,apud00_r,color='blue',marker='o')
#axis.plot(wbP,apud00_i,color='red',marker='o')
#axis.plot(wbP,apud10_r*9,color='green',marker='o')
#axis.plot(wbP,apud10_i*9,color='orange',marker='o')

#axis.plot(wbX,axud00_r,color='blue',marker='o')
#axis.plot(wbX,axud00_i,color='red',marker='o')
#axis.plot(wbX,axud10_r*9,color='green',marker='o')
#axis.plot(wbX,axud10_i*9,color='orange',marker='o')

#axis.plot(wf,s_r,color='blue',marker='o')
#axis.plot(wf0,s_i0,color='green',marker='o')
#axis.plot(wf1,s_i1,color='blue',marker='o')
#axis.plot(wf2,s_i2,color='red',marker='o')

#cs = axis.contourf(np.imag(ERetu1[0,pos_Nff_mu1,:,:]- ERetu2[0,pos_Nff_mu2,:,:]))
#cbar = fig.colorbar(cs)

#axis.plot(wbP0,apud00_r0,color='blue',marker='o')
#axis.plot(wbP0,apud00_i0,color='red',marker='o')
#axis.plot(wbP1,apud00_r1,color='green',marker='o')
#axis.plot(wbP1,apud00_i1,color='orange',marker='o')
axis.plot(wbX2,axud00_r2,color='black',marker='o')
axis.plot(wbX2,axud00_i2,color='magenta',marker='o')

#axis.set_xlim([-4,-2])
axis.set_xlim([-5,+5])

plt.show()

