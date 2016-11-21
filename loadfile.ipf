#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//function loadWaves(aaa,bbb)
//	variable aaa
//	variable bbb
//	variable i;
//	
//	for(i = aaa; i <= bbb; i = i+1)
//		AutoRunV3("Flea3","J:\\JB's experiment\\2016Sep26_AnalysisOfData\\data\\Flea3_23Sep2016_0"+num2str(i)+".ibw") 
//	endfor
//end


// "ListFigures" is used for automatically generating figures of data
// Input:  aaa is the starting file number ex. aaa is "18" for Flea3_03Oct2016_0018
//		bbb is the ending file number ex. bbb is "33" for Flea3_03Oct2016_0033
//		ccc is either "insitu" or "tof"
// Output: none
// 2016/10/17, written by JB

function ListFigures(dd, fn1, fn2, x1, x2, y1, y2, sorting)
       string dd // date of data, format should be like "16Nov2016"
       variable fn1 // start filenumber
       variable fn2 // end filenumber
       variable x1 // left position of the region to be shown
       variable x2 // right position of the region to be shown
       variable y1 // top position of the region to be shown
       variable y2 // bottom position of the region to be shown
       variable sorting //sorting or not

	string filepath = "J:\\JB's experiment\\rawdata\\"+"Flea3_"
	string sss // formated filenumber
	string figname // name of the figures
	
	variable i
	variable/D n_file = fn2 - fn1 + 1
	wave p = root:fittings:parameters
	//make/O/N=(1,6) ODs
	
	wave optdepth = root:Flea3:optdepth
	wave pulsetime = root:Flea3:IndexedWaves:pulsetime
	wave holdtime = root:Flea3:IndexedWaves:holdtime
	wave sortidx = root:Flea3:sortidx
	wave para = root:fittings:parameters
	wave/T filenames = root:Flea3:IndexedWaves:FileNames
	
	//for taking image
	//NewDataFolder/O root:fittings
	//NewPath path1, "J:"
	//NewMovie/A/F=2/O/P=path1
	
	NewDataFolder/O root:gfigures
	SetDataFolder root:Flea3
	if(sorting == 1)
		redimension/N=(n_file) pulsetime
		for(i = fn1; i <= fn2; i = i+1)
			sprintf sss,"%04.0f", i
			AutoRunV3("Flea3",filepath+sss+".ibw") 
			duplicate/O optdepth $("raw"+sss)
		endfor
		duplicate/O pulsetime sortidx
		MakeIndex pulsetime sortidx
		make/T/O/N=(n_file) filenumber // sorted filenumber
		for(i=0; i < n_file; i=i+1)
			splitstring /E="_([[:digit:]]+).ibw" filenames[sortidx[i]], sss
			filenumber[i] = sss
		endfor	
	else
		make/T/O/N=(n_file) filenumber
		for(i = fn1; i <= fn2; i = i+1)
			sprintf sss,"%04.0f", i
			filenumber[i-fn1] = sss
		endfor
	endif
	
	for(i = 0; i < n_file; i = i+1)
		AutoRunV3("Flea3",filepath+dd+"_"+filenumber[i]+".ibw") 
		SetDataFolder root:gfigures
			figname = "a" + dd + "_" + filenumber[i] + "_insitu"
			// Modify the output pictures here
			duplicate/O/R=(x1, x2)(y1, y2) root:Flea3:optdepth $figname
			//R=(-105,-45)(90,150) 
		
			//Fittings for atom number and adiabatic radius
			//SetDataFolder root:fittings
			//imagetransform/G=407 getcol root:Flea3:optdepth
			//duplicate/o root:fittings:W_ExtractedCol root:gfigures:$("h_a" + dd + "_" +sss+"_insitu")
			//imagetransform/G=190 getrow root:Flea3:optdepth
			//duplicate/o root:fittings:W_ExtractedRow root:gfigures:$("h_v" + dd + "_" +sss+"_insitu")
			//SetDataFolder root:fittings
			//n_ring(root:Flea3:optdepth)
			
			SetDataFolder root:gfigures
			//Create image
			NewImage  $figname
			SetAxis/R left y1,y2
			ModifyGraph height={Aspect,1}
			ModifyImage $figname ctab= {-0.1,2,Rainbow,1}
			
			//Edit the textbox here
			TextBox/C/N=text0 "\\Z07"+dd+"_"+sss+" mf=+1, in situ\rraman peak at 100us, detune 250khz\rramptime="+num2str(holdtime[i]*10^3)+"ms"
			TextBox/C/N=text0/X=5/Y=80
		endfor
	//DeletePoints 0,1, ODs
	//KillWaves temp_w
	//CloseMovie
//Execute "TileWindows/G=30/O=1/W=(5,25,1427,673)"
end

// "peak_radius" is used for calculating adiabatic radius of ring BEC
// Usage: 	Set the working dirrectory to be the folder of your data, i.e. pictures of ring BEC.
//			Run peak_radius()
// 2016/10/17, written by JB

function peak_radius()
	string list = WaveList("*",";","")
	string filename
	String savedDF= GetDataFolder(1)
	//String foldername= GetDataFolder(0)
	DFREF dfr_o = GetDataFolderDFR()
	variable numFolder = CountObjectsDFR(dfr_o, 1)
	SetDataFolder root:temp
	DFREF dfr_w = GetDataFolderDFR()
	make/O/N=1 p_r
	InsertPoints 0,numFolder-1, p_r
	make/O/N=1 root:temp:W_coef
	wave w = root:temp:W_coef
	wave lineprofile =  root:temp:W_ImageLineProfile
	wave lineprofilex =  root:temp:W_LineProfileX
	wave lineprofiley =  root:temp:W_LineProfileY
//Make/D/N=1/O W_coef
//W_coef[0] = {1}
//FuncFit/NTHR=0 fff W_coef  fit_W_ImageLineProfile /D 
	variable i;
	variable xcenter, ycenter
	for (i=0; i < numFolder; i=i+1)
		filename = StringFromList(i, list)
		// Find center of the cloud
		CurveFit/M=2/W=0/Q Gauss2D, dfr_o:$filename
		xcenter = w[2] 
		ycenter = w[4]
		
		// Use Lineprofile + Gaussfit to find position of peaks along x- and y-direction 
		make/O/N=2 xtrace1 = {w[2]-15,w[2]+15}, ytrace1 = {w[4],w[4]}, xtrace2 = {w[2],w[2]}, ytrace2 = {w[4]-15,w[4]+15}, xmax={0,0}, ymax={0,0}
		ImageLineProfile/SC srcWave=dfr_o:$filename, xWave=xTrace1, yWave=yTrace1
		//CurveFit/NTHR=0/Q gauss  lineprofile[,11] /X=lineprofilex
		FuncFit/NTHR=0 fff w lineprofile[11,15] /X=lineprofilex /D
		xmax[0] = w[2]
		CurveFit/NTHR=0/Q gauss  lineprofile[11,] /X=lineprofilex
		xmax[1] = w[2]
		ImageLineProfile/SC srcWave=dfr_o:$filename, xWave=xTrace2, yWave=yTrace2
		CurveFit/NTHR=0/Q gauss  lineprofile[,11] /X=lineprofiley
		ymax[0] = w[2]
		CurveFit/NTHR=0/Q gauss   lineprofile[11,] /X=lineprofiley
		ymax[1] = w[2]

		// Calculate mean of distance btween center and peaks
		print xcenter, ycenter, xmax, ymax 
		p_r[i] = ((abs(ymax[0] - ycenter)+abs(ymax[1] - ycenter))/2 +  (abs(xmax[0] - xcenter)+abs(xmax[1] - xcenter))/2)/2
	endfor
	SetDataFolder savedDF
end


// "N-Ring" is used for calculating number of atoms of ring BEC, which, due to its profile, is hard to fitting using AutoRunV3
// Input: the path to the RAW image of ring BEC, which should be a XY image, w/ magnification 4, and camera's pixel size to be 5.6um
// Output: wave "parameter" is generated, stored with
//		  row position of maximum, column position of maximum, range of worked out ROI, background OD, sum OD, sum OD subtracted from background OD, atom #
// 2016/10/17, written by JB

function N_ring(ww)
	wave ww // original raw image
	wave stats = root:fittings:M_WaveStats
	variable rMaxLoc // location of maximum in row direction
	variable cMaxLoc // location of maximum in column direction
	variable rMaxLocO // location of maximum in row direction at raw image
	variable cMaxLocO // location of maximum in column direction at raw image
	variable OD = 0 // Summation of OD of image of interest
	variable avg_OD // average of background OD
	variable range = 50 // the initial region of interest
	variable step = 2 // this is the diffrenece between summation area of temp_inner_sum and of  temp_outer_sum
//	variable r_sum =0
//	variable ri_sum = 0
//	variable c_sum = 0
//	variable ci_sum = 0
	variable temp_inner_sum = 0 // for summation of OD
	variable temp_outer_sum = 0 // for summation of OD
	variable range4back = 100 // region for background
	variable back_sum = 0 // for summation of background OD
	variable magnification = 4 // magnification of image
	variable pixelsize = 5.6 // pixel size of the camera 
	variable i
	variable j

	NewDataFolder/O root:fittings
	SetDataFolder root:fittings
	make/O/N=7 parameters
	
	// find center of the cloud
	Duplicate/o/R=(-300,100)(-150,270) ww, w
	matrixop/o sr = sumrows(w)
	matrixop/o sc = sumcols(w)^t
	WaveStats/Q/W sr
	rMaxLoc = stats[11]
	WaveStats/Q/W sc
	cMaxLoc = stats[11]
	duplicate/o/R=(DimOffset(w, 0) + (rMaxLoc-range)*1.4, DimOffset(w, 0) + (rMaxLoc+range)*1.4)(DimOffset(w, 1) + (cMaxLoc-range)*1.4,DimOffset(w, 1) + (cMaxLoc+range)*1.4) w IOR0
//	r_sum = 0
//	ri_sum = 0
//	c_sum = 0
//	ci_sum = 0
	// calculating the center of mass of the selected area to find center of the cloud more accurately
	// first finding the center in vertical direction
//	for (i=-range; i<=range; i=i+1)
//		for (j=-range; j<=range; j=j+1)
//			r_sum = r_sum + w[rMaxLoc+i][cMaxLoc+j]*(cMaxLoc+j)
//			ri_sum = ri_sum + w[rMaxLoc+i][cMaxLoc+j]
//		endfor
//		c_sum = r_sum / ri_sum * (rMaxLoc+i)
//		print r_sum / ri_sum
//		ci_sum = ci_sum + (rMaxLoc+i)
//	endfor
//	rMaxLoc = c_sum / ci_sum

//	r_sum = 0
//	ri_sum = 0
//	c_sum = 0
//	ci_sum = 0
//	for (j=-range; j<=range; j=j+1)
//		for (i=-range; i<=range; i=i+1)
//			c_sum = c_sum + w[rMaxLoc+i][cMaxLoc+j]*(rMaxLoc+i)
//			ci_sum = ci_sum + (rMaxLoc+i)
//		endfor
//		r_sum = c_sum / ci_sum * (cMaxLoc+j)
//		ri_sum = ri_sum + (cMaxLoc+j)
//	endfor
//	cMaxLoc = r_sum / ri_sum
	
	// Find region of interest
	do 
		temp_inner_sum = 0
		temp_outer_sum = 0
		for (i=-range; i<=range; i=i+1)
			for (j=-range; j<=range; j=j+1)
				temp_inner_sum = temp_inner_sum + w[rMaxLoc+i][cMaxLoc+j]
			endfor
		endfor
		for (i=-range-step; i<=range+step; i=i+1)
			for (j=-range-step; j<=range+step; j=j+1)
				temp_outer_sum = temp_outer_sum + w[rMaxLoc+i][cMaxLoc+j]
			endfor
		endfor
		range = range - 1
	while (temp_inner_sum > temp_outer_sum * 0.95)

	// Duplicate image of interest to IOR
	duplicate/o/R=(DimOffset(w, 0) + (rMaxLoc-range)*1.4, DimOffset(w, 0) + (rMaxLoc+range)*1.4)(DimOffset(w, 1) + (cMaxLoc-range)*1.4,DimOffset(w, 1) + (cMaxLoc+range)*1.4) w IOR

	// Store sum of OD in ROI, here we use  temp_outer_sum as sum OD
	OD = temp_outer_sum
	range = range + 1 +step

	// Calculate average of backgroud OD by selecting a region next to the image of interest 
	back_sum = 0
	rMaxLocO = ((rMaxLoc*DimDelta(w,0)+DimOffset(w,0))-DimOffset(ww,0))/DimDelta(ww,0)
	cMaxLocO = ((cMaxLoc*DimDelta(w,1)+DimOffset(w,1))-DimOffset(ww,1))/DimDelta(ww,1)
	//print rMaxLoc, cMaxLoc, rMaxLocO, cMaxLocO
	for (i=0; i<2*range4back; i=i+1)
		for (j=0; j<2*range4back; j=j+1)
			back_sum = back_sum + ww[rMaxLocO+range+i][cMaxLocO+range+j]
		endfor
	endfor
	avg_OD = back_sum / range4back^2
	
	// Calculate OD of interest
	make/O/N=6 parameters={rMaxLoc, cMaxLoc, range, avg_OD, OD, (OD - avg_OD * range^2)*(pixelsize/magnification)^2/0.23}
end

//function N_rings()
//	string list = WaveList("*",";",""), savedDF= GetDataFolder(1)
//	string filename
//	DFREF dfr_o = GetDataFolderDFR()
//	variable numFolder = CountObjectsDFR(dfr_o, 1)
//	
//	SetDataFolder root:nring
//	DFREF dfr_w = GetDataFolderDFR()
//	make/O/N=1 p_r
//	InsertPoints 0,numFolder-1, p_r
//	make/O/N=(1,6) ODs
//	wave p = root:nring:parameters
//	variable i
//
//	for (i=0 ; i < numfolder ; i = i+1)
//		filename = StringFromList(i, list)
//		n_ring(dfr_o:$filename)
//		matrixop/o p = p^t
//		print p
//		concatenate/O/NP=0 {ODs,p}, temp_w
//		duplicate/O temp_w, ODs
//	endfor
//	DeletePoints 0,1, ODs
//	KillWaves temp_w
//	SetDataFolder savedDF
//end
//
//function ff(x)
//	variable x
//	variable y = (1-exp(-0.00000011*x^2))^2
//
//	return y
//end

//function dimd(x,y,lambda,theta,D)
//	variable x,y,lambda,theta,D
//	
//	return x
//
//end