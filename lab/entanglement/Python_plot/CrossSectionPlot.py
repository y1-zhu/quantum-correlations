import numpy as np
import matplotlib.font_manager
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

matplotlib.matplotlib_fname()


#import matplotlib.pyplot as plt
#import numpy as np
# define subplots with subplot size
#fig, ax = plt.subplots(2, 2, figsize=(8,6))
# define data
#x = [0, 1, 2, 3, 4, 5]
#y = [i**2 for i in x]
# create subplots
#ax[0, 0].plot(x, y, color='black')
#ax[0, 1].plot(x, y, color='green')
#ax[1, 0].plot(x, y, color='red')
#ax[1, 1].plot(x, y, color='blue')
# display the plot
#plt.show()

# Some example data to display
#x = np.linspace(0, 2 * np.pi, 400)
#y = np.sin(x ** 2)


#fig.suptitle('Horizontally stacked subplots')
#axs[0].plot(x, y)
#axs[1].plot(x, -y)

fig, axs = plt.subplots(2,2,figsize=(25,25))

#axs[0].figure(figsize=(10,10), layout='constrained')
#tb= np.loadtxt("ThreeQubit-GHZ+W.txt",delimiter="\t")

#load data 
tb= np.loadtxt("3QCS_W+GHZ(250V).txt", delimiter="\t")
add_tb = np.loadtxt("W-GHZ-allStates.txt", delimiter="\t")
angle_list = np.loadtxt("angle-list-W-GHZ.txt", delimiter="\t")
W_angle = np.loadtxt("angle-list-allStates-W-GHZ.txt", delimiter="\t")
tb_as = np.concatenate((tb, add_tb))
#add_W = np.loadtxt("W-state.txt", delimiter="\t")
#coeff = np.loadtxt("W-state-coefficients.txt", delimiter="\t")
##PPT cutoff
for i in range(len(tb[:,0])): 
	if tb[i,2] > tb[i,0]:
        	tb[i,2] = tb[i,0]



FSEP_COLOR = "skyblue"

tt= np.linspace(0,2*np.pi,250);
## add angle list for all states bound in vicinity of W state
tt_as = np.concatenate((tt, angle_list))
xx=np.cos(tt);
yy=np.sin(tt);
xx_as = np.cos(tt_as);
yy_as = np.sin(tt_as);

arr=[]
a = 1
FONT = 60

## 0: all states  
## 1 : PPT mixture 
## 2 : PPT 
## 5 : FbiSEP (inner) 
## 6 : FbiSEP (outer)
## 13 : FSEP(inner) 
## 14 : FSEP(outer)
## 15 : BiSEP(inner) 
## 16 : BiSEP(outer) 


axs[0,0].axis('off')
#axs[0,0].scatter(xx[0]*tb[0,0],yy[0]*tb[0,0],s=200,color="blue")
#axs[0,0].scatter(np.multiply(coeff[0],add_W[0]),np.multiply(coeff[1],add_W[0]),s=50,color="blue")

#test_array_1 = np.concatenate((np.multiply(xx,tb[:,0]),np.multiply(coeff[0],add_W)))
#test_array_2 = np.concatenate((np.multiply(yy,tb[:,0]),np.multiply(coeff[1],add_W)))
#x = np.random.rand(1)
#y = np.random.rand(1)
#axs[0,0].scatter(x,y,s = 2,color='green')
axs[0,0].fill(np.multiply(xx_as,tb_as[:,0]),np.multiply(yy_as,tb_as[:,0]),label='all states',linewidth=a)
#axs[0,0].fill(test_array_1,test_array_2,label='all states',linewidth=3)
axs[0,0].fill(np.multiply(xx,tb[:,15]),np.multiply(yy,tb[:,15]),label='BSEP',linewidth=a)
#white background for Bi separability
axs[0,0].fill(np.multiply(xx,tb[:,5]),np.multiply(yy,tb[:,5]),color='white')
#axs[0,0].plot(np.multiply(xx,tb[:,5]),np.multiply(yy,tb[:,5]),linewidth=a,color='green')
#axs[0,0].fill(np.multiply(xx,tb[:,5]),np.multiply(yy,tb[:,5]),linewidth=a,color='green',alpha=0.2)
#axs[0,0].fill(np.multiply(xx,tb[:,5]),np.multiply(yy,tb[:,5]),linewidth=a,color='green',alpha=0.2)
axs[0,0].fill(np.multiply(xx,tb[:,5]),np.multiply(yy,tb[:,5]),label='SEP A$\\vert$BC SEP AB$\\vert$C',linewidth=a,color='green',alpha=0.2)
#axs[0,0].plot(np.multiply(xx,tb[:,5]),np.multiply(yy,tb[:,5]),linewidth=a,color='green')
axs[0,0].fill(np.multiply(xx,tb[:,5]),np.multiply(yy,tb[:,5]),linewidth=a,color='green',alpha=0.2)
axs[0,0].fill(np.multiply(xx,tb[:,5]),np.multiply(yy,tb[:,5]),linewidth=a,color='green',alpha=0.2)
#axs[0,0].fill(arr,arr,label='SEP A|BC',linewidth=a,alpha=0.2,color='green')
#axs[0,0].fill(arr,arr,label='SEP B$\\vert$AC',linewidth=a,color='green',alpha=0.2)
#axs[0,0].fill(arr,arr,label='SEP AB$\\vert$C',linewidth=a,color='green',alpha=0.2)
axs[0,0].fill(np.multiply(xx,tb[:,13]),np.multiply(yy,tb[:,13]),label='FSEP',linewidth=a,color=FSEP_COLOR)

# plotting the W state 

#axs[0,0].fill(np.multiply(coeff[0],add_W[0]),np.multiply(coeff[1],add_W[0]),label='all states',linewidth=3)
#axs[0,0].fill(np.multiply(coeff[0],add_W[15]),np.multiply(coeff[1],add_W[15]),label='biSEP',linewidth=3)
#axs[0,0].fill(np.multiply(coeff[0],add_W[5]),np.multiply(coeff[1],add_W[5]),label='SEP A$\\vert$BC',linewidth=3,color='green',alpha=0.2)
#axs[0,0].plot(np.multiply(coeff[0],add_W[5]),np.multiply(coeff[1],add_W[5]),linewidth=a,color='green')
#axs[0,0].fill(np.multiply(coeff[0],add_W[5]),np.multiply(coeff[1],add_W[5]),linewidth=3,color='green',alpha=0.2)
#axs[0,0].fill(np.multiply(coeff[0],add_W[5]),np.multiply(coeff[1],add_W[5]),linewidth=3,color='green',alpha=0.2)
#axs[0,0].fill(arr,arr,label='SEP B$\\vert$AC',linewidth=a,color='green',alpha=0.2)
#axs[0,0].fill(arr,arr,label='SEP AB$\\vert$C',linewidth=a,color='green',alpha=0.2)
#axs[0,0].fill(np.multiply(coeff[0],add_W[13]),np.multiply(coeff[1],add_W[13]),label='FSEP',linewidth=3,color='lightseagreen')
#axs[0,0].plot(arr,arr,label='SEP A|BC',linewidth=a,alpha=0.2,color='green')
#axs[0,0].plot(arr,arr,label='SEP B|AC',linewidth=a,color='green',alpha=0.2)
#axs[0,0].plot(arr,arr,label='SEP AB|C',linewidth=a,color='green',alpha=0.2)
#axs[0].plot(np.multiply(xx,tb[:,1]),np.multiply(yy,tb[:,1]))
#axs[0].plot(np.multiply(xx,tb[:,2]),np.multiply(yy,tb[:,2]))


#axs[0].plot(np.multiply(xx,tb[:,5]),np.multiply(yy,tb[:,5]))
#axs[0].plot(np.multiply(xx,tb[:,6]),np.multiply(yy,tb[:,6]))
#axs[0].plot(np.multiply(xx,tb[:,7]),np.multiply(yy,tb[:,7]))


axs[0,0].xaxis.set_ticklabels([])
axs[0,0].yaxis.set_ticklabels([])

axs[0,0].text(xx[0]*tb[0,0]-0.07,yy[0]*tb[0,0]+0.05,'GHZ',fontsize=FONT)
axs[0,0].text(0.9,-0.2,'(a)',fontsize=FONT)
axs[0,0].text(W_angle[0],W_angle[1],'W',fontsize=FONT)
axs[0,0].scatter(W_angle[0],W_angle[1],s=100,color="blue")
axs[0,0].scatter(np.multiply(xx_as[0],tb_as[0,0]),np.multiply(yy_as[0],tb_as[0,0]),s=100,color="blue")

axs[0,0].legend(loc='upper center', bbox_to_anchor=(1.02, 1.27),
          ncol=3,fontsize=FONT,frameon=False)

#axs[0,0].legend(loc='upper center', bbox_to_anchor=(0.8, 1.0),
#          ncol=2,fontsize=36,frameon=False)

#axs[0,0].legend(loc='upper center', bbox_to_anchor=(0.6, 1.0),
         # ncol=2,fontsize=28,frameon=False)

#tb= np.loadtxt("ThreeQubit-Random+GHZ.txt",delimiter="\t")
tb= np.loadtxt("3QCS_Random+GHZ(250V).txt",delimiter="\t")

tt= np.linspace(0,2*np.pi,250);
xx=np.cos(tt);
yy=np.sin(tt);

## 0: all states  
## 1 : PPT mixture 
## 2 : PPT A|BC
## 3 : PPT B|AC
## 4 : PPT AB|C
## 5 : biSEP A|BC (inner) 
## 6 : biSEP A|BC (outer)
## 7 : biSEP B|AC (inner)
## 8 : biSEP B|AC (outer)
## 9 : biSEP AB|C (inner)
## 10 : biSEP AB|C (outer)
## 11 : FbiSEP (inner)
## 12 : FbiSEP (outer)
## 13 : FSEP(inner) 
## 14 : FSEP(outer)
## 15 : BiSEP(inner) 
## 16 : BiSEP(outer) 


## PPT cutoff 
state_num = len(tb[:,0])
ppt = np.zeros(state_num)
for i in range(state_num): 
	if tb[i,2] > tb[i,0]:
        	tb[i,2] = tb[i,0]
	if tb[i,3] > tb[i,0]:
		tb[i,3] = tb[i,0]
	if tb[i,4] > tb[i,0]:
		tb[i,4] = tb[i,0]
	ppt[i] = np.min([tb[i,2],tb[i,3],tb[i,4]])


axs[0,1].axis('off')
axs[0,1].scatter(xx[0]*tb[0,0],yy[0]*tb[0,0],s=200,color="blue")
#axs[0,1].scatter(xx[67]*tb[67,0],yy[62]*tb[67,0],s=200,color="blue")
axs[0,1].fill(np.multiply(xx,tb[:,0]),np.multiply(yy,tb[:,0]),label='all states',linewidth=a)
#axs[1].plot(np.multiply(xx,tb[:,1]),np.multiply(yy,tb[:,1]),label='biSEP outer',linewidth=a)
#axs[1].plot(np.multiply(xx,tb[:,13]),np.multiply(yy,tb[:,13]),label='PPT mixture',linewidth=a)
axs[0,1].fill(np.multiply(xx,tb[:,15]),np.multiply(yy,tb[:,15]),label='biSEP',linewidth=a)
#axs[1].plot(np.multiply(xx,tb[:,7]),np.multiply(yy,tb[:,7]),label='FbiSEP outer',linewidth=a)
#axs[1].plot(np.multiply(xx,ppt),np.multiply(yy,ppt),label='PPT',linewidth=a)
#axs[1].plot(np.multiply(xx,tb[:,12]),np.multiply(yy,tb[:,12]),label='FbiSEP',linewidth=a)
axs[0,1].fill(np.multiply(xx,tb[:,5]),np.multiply(yy,tb[:,5]),color='white')
axs[0,1].fill(np.multiply(xx,tb[:,7]),np.multiply(yy,tb[:,7]),color='white')
axs[0,1].fill(np.multiply(xx,tb[:,9]),np.multiply(yy,tb[:,9]),color='white')

axs[0,1].fill(np.multiply(xx,tb[:,5]),np.multiply(yy,tb[:,5]),label='SEP A$\\vert$BC',linewidth=a,alpha=0.2,color='green')
axs[0,1].fill(np.multiply(xx,tb[:,7]),np.multiply(yy,tb[:,7]),label='SEP B$\\vert$AC',linewidth=a,color='green',alpha=0.2)
axs[0,1].fill(np.multiply(xx,tb[:,9]),np.multiply(yy,tb[:,9]),label='SEP AB$\\vert$C',linewidth=a,color='green',alpha=0.2)
axs[0,1].plot(np.multiply(xx,tb[:,5]),np.multiply(yy,tb[:,5]),linewidth=a,color='green')
axs[0,1].plot(np.multiply(xx,tb[:,7]),np.multiply(yy,tb[:,7]),linewidth=a,color='green')
axs[0,1].plot(np.multiply(xx,tb[:,9]),np.multiply(yy,tb[:,9]),linewidth=a,color='green')
#axs[1].plot(np.multiply(xx,tb[:,14]),np.multiply(yy,tb[:,14]),label='FSEP outer',linewidth=a)
axs[0,1].fill(np.multiply(xx,tb[:,13]),np.multiply(yy,tb[:,13]),label='FSEP',linewidth=a,color=FSEP_COLOR)

axs[0,1].xaxis.set_ticklabels([])
axs[0,1].yaxis.set_ticklabels([])

axs[0,1].text(0.85,-0.05,'GHZ',fontsize=FONT)
axs[0,1].text(0.9,-0.2,'(b)',fontsize=FONT)

#tb= np.loadtxt("3QCS_W+Random(250V).txt",delimiter="\t")
tb= np.loadtxt("3QCS_W+Random(250V)_new.txt", delimiter="\t")


tt= np.linspace(0,2*np.pi,250);
xx=np.cos(tt);
yy=np.sin(tt);

## 0: all states  
## 1 : PPT mixture 
## 2 : PPT A|BC
## 3 : PPT B|AC
## 4 : PPT AB|C
## 5 : biSEP A|BC (inner) 
## 6 : biSEP A|BC (outer)
## 7 : biSEP B|AC (inner)
## 8 : biSEP B|AC (outer)
## 9 : biSEP AB|C (inner)
## 10 : biSEP AB|C (outer)
## 11 : FbiSEP (inner)
## 12 : FbiSEP (outer)
## 13 : FSEP(inner) 
## 14 : FSEP(outer)
## 15 : BiSEP(inner) 
## 16 : BiSEP(outer) 

## PPT cutoff 
state_num = len(tb[:,0])
ppt = np.zeros(state_num)
for i in range(state_num): 
	if tb[i,2] > tb[i,0]:
        	tb[i,2] = tb[i,0]
	if tb[i,3] > tb[i,0]:
		tb[i,3] = tb[i,0]
	if tb[i,4] > tb[i,0]:
		tb[i,4] = tb[i,0]
	ppt[i] = np.min([tb[i,2],tb[i,3],tb[i,4]])



axs[1,0].axis('off')
axs[1,0].scatter(xx[0]*tb[0,0],yy[0]*tb[0,0],s=200,color="blue")
#axs[1,0].scatter(xx[62]*tb[62,0],yy[62]*tb[62,0],s=200,color="blue")
axs[1,0].fill(np.multiply(xx,tb[:,0]),np.multiply(yy,tb[:,0]),label='all states',linewidth=a)
#axs[1].plot(np.multiply(xx,tb[:,1]),np.multiply(yy,tb[:,1]),label='biSEP outer',linewidth=a)
#axs[1].plot(np.multiply(xx,tb[:,13]),np.multiply(yy,tb[:,13]),label='PPT mixture',linewidth=a)
axs[1,0].fill(np.multiply(xx,tb[:,15]),np.multiply(yy,tb[:,15]),label='biSEP',linewidth=a)
#axs[1].plot(np.multiply(xx,tb[:,7]),np.multiply(yy,tb[:,7]),label='FbiSEP outer',linewidth=a)
#axs[1].plot(np.multiply(xx,ppt),np.multiply(yy,ppt),label='PPT',linewidth=a)
#axs[1].plot(np.multiply(xx,tb[:,12]),np.multiply(yy,tb[:,12]),label='FbiSEP',linewidth=a)
axs[1,0].fill(np.multiply(xx,tb[:,5]),np.multiply(yy,tb[:,5]),color='white')
axs[1,0].fill(np.multiply(xx,tb[:,7]),np.multiply(yy,tb[:,7]),color='white')
axs[1,0].fill(np.multiply(xx,tb[:,9]),np.multiply(yy,tb[:,9]),color='white')

axs[1,0].fill(np.multiply(xx,tb[:,5]),np.multiply(yy,tb[:,5]),label='SEP A$\\vert$BC',linewidth=a,alpha=0.2,color='green')
axs[1,0].fill(np.multiply(xx,tb[:,7]),np.multiply(yy,tb[:,7]),label='SEP B$\\vert$AC',linewidth=a,color='green',alpha=0.2)
axs[1,0].fill(np.multiply(xx,tb[:,9]),np.multiply(yy,tb[:,9]),label='SEP AB$\\vert$C',linewidth=a,color='green',alpha=0.2)
axs[1,0].plot(np.multiply(xx,tb[:,5]),np.multiply(yy,tb[:,5]),linewidth=a,color='green')
axs[1,0].plot(np.multiply(xx,tb[:,7]),np.multiply(yy,tb[:,7]),linewidth=a,color='green')
axs[1,0].plot(np.multiply(xx,tb[:,9]),np.multiply(yy,tb[:,9]),linewidth=a,color='green')
#axs[1].plot(np.multiply(xx,tb[:,14]),np.multiply(yy,tb[:,14]),label='FSEP outer',linewidth=a)
axs[1,0].fill(np.multiply(xx,tb[:,13]),np.multiply(yy,tb[:,13]),label='FSEP',linewidth=a,color=FSEP_COLOR)

axs[1,0].xaxis.set_ticklabels([])
axs[1,0].yaxis.set_ticklabels([])

#axs[1,0].legend(loc='upper center', bbox_to_anchor=(0.7, 1.0),
 #         ncol=2,fontsize=32,frameon=False)

#axs[1,0].text(np.cos(1.7115926535897932)*tb[62,0],np.sin(1.7115926535897932)*tb[62,0],'W',fontsize=40)
axs[1,0].text(0.9,-0.2,'(c)',fontsize=FONT)
axs[1,0].text(xx[0]*tb[0,0],yy[0]*tb[0,0], 'W', fontsize=FONT)

tb= np.loadtxt("3QCS_Random+Random(250V).txt",delimiter="\t")



tt= np.linspace(0,2*np.pi,250);
xx=np.cos(tt);
yy=np.sin(tt);

##1 FSEP inner
##2 biSEP outer 
##3 SEP1 outer 
##4 SEP2 inner 
##5 SEP3 outer 
##6 all states 
##7 SEP1 inner 
##8 Min SEP outer
##9 PPT1
##10 PPT3
##11 SEP3 inner 
##12 SEP2 outer 
##13 Min SEP inner 
##14 PPT mixture 
##15 FSEP outer 
##16 biSEP inner 
##17 PPT2

## PPT cutoff 
state_num = len(tb[:,0])
ppt = np.zeros(state_num)
for i in range(state_num): 
	if tb[i,2] > tb[i,0]:
        	tb[i,2] = tb[i,0]
	if tb[i,3] > tb[i,0]:
		tb[i,3] = tb[i,0]
	if tb[i,4] > tb[i,0]:
		tb[i,4] = tb[i,0]
	ppt[i] = np.min([tb[i,2],tb[i,3],tb[i,4]])



axs[1,1].axis('off')
#axs[1,1].scatter(xx[0]*tb[0,0],yy[0]*tb[0,0],s=200,color="blue")
#axs[1,1].scatter(xx[62]*tb[62,0],yy[62]*tb[62,0],s=200,color="blue")
axs[1,1].fill(np.multiply(xx,tb[:,0]),np.multiply(yy,tb[:,0]),label='all states',linewidth=a)
#axs[1].plot(np.multiply(xx,tb[:,1]),np.multiply(yy,tb[:,1]),label='biSEP outer',linewidth=a)
#axs[1].plot(np.multiply(xx,tb[:,13]),np.multiply(yy,tb[:,13]),label='PPT mixture',linewidth=a)
axs[1,1].fill(np.multiply(xx,tb[:,15]),np.multiply(yy,tb[:,15]),label='biSEP',linewidth=a)
#axs[1].plot(np.multiply(xx,tb[:,7]),np.multiply(yy,tb[:,7]),label='FbiSEP outer',linewidth=a)
#axs[1].plot(np.multiply(xx,ppt),np.multiply(yy,ppt),label='PPT',linewidth=a)
#axs[1].plot(np.multiply(xx,tb[:,12]),np.multiply(yy,tb[:,12]),label='FbiSEP',linewidth=a)
axs[1,1].fill(np.multiply(xx,tb[:,5]),np.multiply(yy,tb[:,5]),color='white')
axs[1,1].fill(np.multiply(xx,tb[:,7]),np.multiply(yy,tb[:,7]),color='white')
axs[1,1].fill(np.multiply(xx,tb[:,9]),np.multiply(yy,tb[:,9]),color='white')

axs[1,1].fill(np.multiply(xx,tb[:,5]),np.multiply(yy,tb[:,5]),label='SEP A$\\vert$BC',linewidth=a,alpha=0.2,color='green')
axs[1,1].fill(np.multiply(xx,tb[:,7]),np.multiply(yy,tb[:,7]),label='SEP B$\\vert$AC',linewidth=a,color='green',alpha=0.2)
axs[1,1].fill(np.multiply(xx,tb[:,9]),np.multiply(yy,tb[:,9]),label='SEP AB$\\vert$C',linewidth=a,color='green',alpha=0.2)
axs[1,1].plot(np.multiply(xx,tb[:,5]),np.multiply(yy,tb[:,5]),linewidth=a,color='green')
axs[1,1].plot(np.multiply(xx,tb[:,7]),np.multiply(yy,tb[:,7]),linewidth=a,color='green')
axs[1,1].plot(np.multiply(xx,tb[:,9]),np.multiply(yy,tb[:,9]),linewidth=a,color='green')
#axs[1].plot(np.multiply(xx,tb[:,14]),np.multiply(yy,tb[:,14]),label='FSEP outer',linewidth=a)
axs[1,1].fill(np.multiply(xx,tb[:,13]),np.multiply(yy,tb[:,13]),label='FSEP',linewidth=a,color=FSEP_COLOR)

axs[1,1].set_xlim(-0.17,0.4)
axs[1,1].set_ylim(-0.21,0.35)
axs[1,1].xaxis.set_ticklabels([])
axs[1,1].yaxis.set_ticklabels([])

axs[1,1].text(0.35,-0.18,'(d)',fontsize=FONT)

#fig.tight_layout()
#plt.subplot_tool()
#plt.show()

#axs[1,0].legend(loc='upper center', bbox_to_anchor=(0.7, 1.0),
 #         ncol=2,fontsize=32,frameon=False)

#axs[1,0].text(-0.16,0.85,'W',fontsize=40)

plt.subplots_adjust(wspace=0.05, hspace=0.0)


plt.savefig("test3.pdf",transparent=True)

