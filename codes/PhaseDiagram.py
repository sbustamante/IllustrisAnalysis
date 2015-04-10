#PhaseDiagram.py
#
#This code makes a video of the evolution of a slide of the simulation
#Usage: PhaseDiagram.py <TYPE OF DIAGRAM, format = XYC (X property, Y property, Coloured property)>
#			<Plot (0) or video (1)> <number of snaps>
#	R - density
#	P - pressure
#	T - temperature
#	U - internal energy
#	M - mass
#Example: PhaseDiagram.py RPT 1
#
#by: Sebastian Bustamante

execfile('_head.py')

#==================================================================================================
#			PARAMETERS
#==================================================================================================
#Simulation
simulation = "VPH_064"
#Box lenght [kpc h^-1]
Box = 20000
#Snapbase
snapbase = "snap"
#Number of snapshots
snaps = int(sys.argv[3])
#Number of files per snap
snapfiles = 1
#Sampling rate
sampling = 1
#Type of print <0-all or 1-gas or 2-DM> (This should be always set to 1)
type = 1

#Gas properties
prop_dict = {"U":"Energy", "R":"Density", "P":"Pressure", "T":"Temperature", "M":"Mass", "D":"Rdomain"}
units = { \
"Energy":		"(km/sec)$^2$",\
"Density":		"$10^{10}$ h$^{-1}$M$_{\odot}$/(h$^{-1}$ kpc)$^3$",\
"Pressure":		"Pa",\
"Temperature":		"K",\
"Mass":			"$10^{10}$ h$^{-1}$M$_{\odot}$",\
"Rdomain":		"R$_{vir}$"}
notat = { \
"Energy":		"u",\
"Density":		"$\\rho$",\
"Pressure":		"P",\
"Temperature":		"T",\
"Mass":			"M",\
"Rdomain":		"R$_{dom}$"}	
ranges = { \
"Energy":		[ 1, 6 ],\
"Density":		[ -12, 0 ],\
"Pressure":		[ -22, -11 ],\
"Temperature":		[ 2, 8 ],\
"Mass":			[ -2, 2 ],\
"Rdomain":		[ -2.5, 2.5 ]}
indexes = { "Energy":8, "Density":9, "Pressure":10, "Temperature":11, "Mass":7, "Rdomain":14 }
properties = [prop_dict[prop] for prop in sys.argv[1]]	#Specific order of properties

#Window parameters
prop1 = ranges[properties[0]]	#X
prop2 = ranges[properties[1]]	#Y
prop3 = ranges[properties[2]]	#Color

prop1_hist = 5
prop2_hist = 5
#Bins of 1D histograms
bins = 28*2.
#Bins of 2D histogram
bins2D = 300.

#==================================================================================================
#			MAKING VIDEO OF EVOLUTION
#==================================================================================================
if sys.argv[2] == "0":		#Plot
    snaprange = [snaps]
else:				#Video
    snaprange = xrange( snaps )
    
for snap in snaprange:
    #Formating plot................................................................................
    #no labels
    nullfmt = NullFormatter()

    #definitions for the axes
    left, width = 0.1, 0.65
    bottom_v, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02

    rect_hist2D = [left, bottom_v, width, height]
    rect_histx = [left, bottom_h, 1.24*width, 0.2]
    rect_histy = [left_h, bottom_v, 0.2, height]

    #start with a rectangular Figure
    plt.figure(1, figsize=(8,8))

    axHist2D = plt.axes(rect_hist2D)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    #no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    
    #Loading data of current snap..................................................................
    data = ascii( "%s/%s/%s"%(datafold,simulation,snapbase), snap, snapfiles, sampling, type )
    
    #Calculating domains...........................................................................
    for prop in sys.argv[1]:
	if prop == "D":
	    data_domains = gasdomains( data[:,1], data[:,2], data[:,3], simulation, snapbase, snap, Box )
	    data2 = np.zeros( (len(data), len(data[0])+1 ) )
	    data2[:,:-1] = data
	    data2[:,-1] = data_domains
	    data = data2
	    del data2
	    break
    
    #Ordered properties
    property1 = np.log10(data[:,indexes[ properties[0] ]])
    property2 = np.log10(data[:,indexes[ properties[1] ]])
    property3 = data[:,indexes[ properties[2] ]]

    #Phase diagram
    ##Scatter.......................................................................................
    #scatter2d = axHist2D.scatter( property1, property2, c=property3, s=2, marker='.',\
    #linewidth=0.01, cmap='jet', vmin = prop3[0], vmax = prop3[1] );
    ##Create the colorbar
    #axc, kw = matplotlib.colorbar.make_axes( axHistx,\
    #orientation = "vertical", shrink=1., pad=.04, aspect=10, anchor=(.5,1.5) )
    #cb = matplotlib.colorbar.Colorbar( axc, scatter2d, orientation = "vertical" )
    #cb.set_label( "log$_{10}$(%s  [%s])"%(notat[properties[2]],units[properties[2]]), \
    #labelpad=-80, fontsize=12 )
  
    #Histogram2D...................................................................................
    Hist_phase = np.transpose(np.histogram2d( property1, property2, weights = property3,
    bins = bins2D, range = (prop1, prop2) )[0][::,::-1] )
    Hist_phase_num = np.transpose(np.histogram2d( property1, property2,
    bins = bins2D, range = (prop1, prop2) )[0][::,::-1] )
    
    map2d = axHist2D.imshow( np.log10(Hist_phase[::,::]/Hist_phase_num[::,::]), interpolation='nearest', aspect = 'auto',
    cmap = 'jet', extent = (prop1[0],prop1[1],prop2[0],prop2[1]), vmin = prop3[0], vmax = prop3[1] )
    #Create the colorbar
    axc, kw = matplotlib.colorbar.make_axes( axHistx,\
    orientation = "vertical", shrink=1., pad=.04, aspect=10, anchor=(.5,1.5) )
    cb = matplotlib.colorbar.Colorbar( axc, map2d, orientation = "vertical" )
    cb.set_label( "log$_{10}$(%s  [%s])"%(notat[properties[2]],units[properties[2]]), \
    labelpad=-70, fontsize=12 )
    ticks = np.linspace( prop3[0],prop3[1],9 )
    cb.set_ticks( ticks )
    cb.set_ticklabels( ['{0:3.1f}'.format(t) for t in ticks] )
    
    #Histogram X (Prop1)
    histx = np.histogram( property1, bins=bins, normed=True, range=prop1 )
    axHistx.bar( histx[1][:-1], histx[0], width = (prop1[1]-prop1[0])/bins,\
    linewidth=2.0, color="gray" )
    #Histogram Y (Prop2)
    histy = np.histogram( property2, bins=bins, normed=True, range=prop2 )
    axHisty.barh( histy[1][:-1], histy[0], height = (prop2[1]-prop2[0])/bins,\
    linewidth=2.0, color="gray" )
    
    #Plot format...................................................................................
    #Upper histogram
    axHistx.set_xlim( prop1 )
    axHistx.set_xticks( np.linspace( prop1[0],prop1[1],bins/4.+1 ) )
    yticks = np.linspace( 0, np.max(histx[0]), prop1_hist+1 )
    axHistx.set_yticks( yticks )
    axHistx.set_yticklabels( ['{0:1.2f}'.format(yt) for yt in yticks] )
    axHistx.grid( color='black', linestyle='--', linewidth=1., alpha=0.3 )
    axHistx.set_ylabel( "Normed distribution" )
    
    #Right histogram
    axHisty.set_ylim( prop2 )
    axHisty.set_yticks( np.linspace( prop2[0],prop2[1],bins/4.+1 ) )
    xticks = np.linspace( 0, np.max(histy[0]), prop2_hist+1 )
    axHisty.set_xticks( xticks )
    axHisty.set_xticklabels( ['{0:1.2f}'.format(xt) for xt in xticks], rotation=-90 )
    axHisty.grid( color='black', linestyle='--', linewidth=1., alpha=0.3 )
    axHisty.set_xlabel( "Normed distribution" )
    
    #2D Scatter
    axHist2D.grid( color='black', linestyle='--', linewidth=1., alpha=0.3 )
    axHist2D.set_xlim( prop1 )
    axHist2D.set_ylim( prop2 )
    xticks = np.linspace( prop1[0],prop1[1],bins/4.+1 )
    axHist2D.set_xticks( xticks, ['{0:1.2f}'.format(xt) for xt in xticks] )
    yticks = np.linspace( prop2[0],prop2[1],bins/4.+1 )
    axHist2D.set_yticks( yticks, ['{0:1.2f}'.format(yt) for yt in yticks] )
    axHist2D.set_xlabel( "log$_{10}$(%s  [%s])"%(notat[properties[0]],units[properties[0]]) )
    axHist2D.set_ylabel( "log$_{10}$(%s  [%s])"%(notat[properties[1]],units[properties[1]]) )
    axHist2D.legend( loc='upper right', fancybox=True, shadow=True, ncol = 1, prop={'size':10} )
    axHist2D.text( prop1[1] - 0.95*(prop1[1]-prop1[0]),prop2[1] - 0.08*(prop2[1]-prop2[0]),\
    "%s:snapshot %03d\nz = %1.2f"%(simulation,snap,abs(data[0,12])), fontweight="bold", fontsize = 10 )

    if sys.argv[2] == "1":		#Video
	fname='_tmp-%03d.png'%snap
	plt.savefig(fname)
	plt.close()


#Making the video
if sys.argv[2] == "0":		#Plot
    plt.savefig("Phase_diagram_%s_%s.png"%(sys.argv[1], simulation))
    plt.show()
else:				#Video
    print 'Making movie animation.mpg - this make take a while'
    os.system("ffmpeg -qscale 1 -r 10 -b 9600 -i _tmp-%03d.png  video.mp4")
    #Deleting temporal images
    os.system('rm -rf *.png')