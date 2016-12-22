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

	string filepath = "J:\\JB\'s experiment\\rawdata\\"+"Flea3_"
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
	wave wcoef = root:fittings:W_coef
	wave havehole = root:fittings:'havingHoles_img1125_2001-2030'
	wave adrffreq = root:Flea3:IndexedWaves:adrffreq
	
	//for taking image
	//NewDataFolder/O root:fittings
	//NewPath path1, "J:"
	//NewMovie/A/F=5/O/P=path1
	
	NewDataFolder/O root:gfigures
	SetDataFolder root:Flea3
	
	if(sorting == 1)
		redimension/N=(n_file) pulsetime
		for(i = fn1; i <= fn2; i = i+1)
			sprintf sss,"%04.0f", i
			//print filepath+dd+"_"+sss+".ibw"
			AutoRunV3("Flea3",filepath+dd+"_"+sss+".ibw") 
			//duplicate/O optdepth $("raw"+sss)
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
	//SetDataFolder root:fittings
	//make/O/N=(n_file,4) centers=nan
	//SetDataFolder root:gfigures
	for(i = 0; i < n_file; i = i+1)
		AutoRunV3("Flea3",filepath+dd+"_"+filenumber[i]+".ibw") 
		//SetDataFolder root:gfigures
		//figname = "a" + dd + "_" + filenumber[i]
		// Modify the output pictures here
		//duplicate/O/R=(x1, x2)(y1, y2) root:Flea3:optdepth $figname
		//R=(-105,-45)(90,150) 
		
		//SetDataFolder root:gfigures
		//Create image
		//NewImage  $figname
		//SetAxis/R left y1,y2
		//ModifyGraph width=288, height={Aspect,1}
		//ModifyGraph height={Aspect,1}
		//ModifyImage $figname ctab= {-0.1,2,Rainbow,1}
		
		//Edit the textbox here
		//TextBox/C/N=text0 "\\Z07"+dd+"_"+filenumber[i]+" mf=+1, TOF 24ms\rraman peak at 100us, holdtime 13ms\rramptime 7.5ms, detune="+num2str((adrffreq[i]-adrffreq[0])*10^3)+" kHz"
		//TextBox/C/N=text0 "\\Z07"+dd+"_"+filenumber[i]+" pulsetime="+num2str(pulsetime[sortidx[i]]*10^6)+" us\r number = "+ num2str(n_ring($figname, x1, x2, y1, y2,0))
		//TextBox/C/N=text0 "\\Z07"+dd+"_"+filenumber[i]+" mf=-1, TOF 24ms\rN="+num2str(para[5])
		//sprintf sss,"%05.1f", holdtime[i]*10^3
		//TextBox/C/N=text0 "\Z18Hold time = "+sss+" ms"
		//TextBox/C/N=text0/X=1/Y=90
		//TextBox/C/N=text0 "\\Z07"+dd+"_"+filenumber[i]

		//Fittings for atom number and adiabatic radius				
		//SetDataFolder root:fittings
		//CurveFit/NTHR=0/Q Gauss2D  root:gfigures:$figname
		//K6 = 0;
		//CurveFit/Q/H="0000001"/NTHR=0 Gauss2D   root:gfigures:$figname


		//do
		//Make/D/N=13/O W_coef=0
		//if ( havehole[i] == 5)
			//W_coef[0] = {0.0007,13,-32,4.4,159,4.4,0,-13,-31,7,159.8,6,0}
			//W_coef[0] = {0.1,0.4,-35,5.5,163,5.25,0,-5,-40,4,160,4,0} // with BECwithHole2 for Dec05
			//W_coef[0] = {0.1,4,-35,5.5,163,5.25,0,-5,-38,4,162,4,0} // with BECwithHole2 for Dec05
		//	W_coef[0] = {-0.01,1,-35,5,155,7,0,5,-31,10,159,5,0} // with BECwithHole2 for Dec06
			//W_coef[0] = {0.004,5+enoise(5),-39.8+enoise(0.2),4.8,154+enoise(0.2),5.78,0,-6+enoise(1),-35.8,5+enoise(0.2),156.8,5+enoise(0.2),0}
		//	W_coef={-0.00075845,0.48922,-40.32+enoise(0.2),4.0706,154.48+enoise(0.2),4.328,0.027159,44.071,-36.961+enoise(0.2),8.6943+enoise(0.2),155.57+enoise(0.2),44.991+enoise(0.2),-0.46734}
		 // 	W_coef={0.012218,7.3559+enoise(1),-35.627+enoise(1),4.2407,153.75+enoise(1),3.8687,0,-2.4862+enoise(1),-.436+enoise(1),14.659,154+enoise(1),3.9898,0}	
		  //	FuncFitMD/H="0000001000001"/Q/NTHR=0 BECwithHole2 W_coef  root:gfigures:$figname
			//centers[i][0] = wcoef[2]
			//centers[i][1] = wcoef[4]
		//	centers[i][2] = wcoef[8]
		//	centers[i][3] = wcoef[10]
		//while(centers [i][2] > -41)
		//else
		//	CurveFit/NTHR=0/Q Gauss2D  root:gfigures:$figname
		//	centers[i][0] = wcoef[2]
		//	centers[i][1] = wcoef[4]
		//	centers[i][2] = nan
		//	centers[i][3] = nan
		//endif
		//imagetransform/G=407 getcol root:Flea3:optdepth
		//duplicate/o root:fittings:W_ExtractedCol root:gfigures:$("h_a" + dd + "_" +sss+"_insitu")
		//imagetransform/G=190 getrow root:Flea3:optdepth
		//duplicate/o root:fittings:W_ExtractedRow root:gfigures:$("h_v" + dd + "_" +sss+"_insitu")
		//DoUpdate
		//AddMovieFrame
	endfor
	//DeletePoints 0,1, ODs
	//KillWaves temp_w
	//CloseMovie
//Execute "TileWindows/G=30/O=1/W=(5,25,1427,673)"
//Execute "TileWindows/G=5/O=1/W=(5,25,1427,673)"
//Execute "TileWindows/G=5/A=(4,0)/O=1/W=(5,29,1427,736)"
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

function/D N_ring(optdepth, x1, x2, y1, y2, resize)
	wave optdepth // original raw image
	variable resize
	// XY image
	variable x1 //200//100
	variable x2 //240//330
	variable y1 //410//290
	variable y2 //450//550
	
	wave wcoef = root:fittings:W_coef
	variable rMaxLoc // location of maximum in row direction
	variable cMaxLoc // location of maximum in column direction
	variable range=100  // the initial region of interest
	variable step = 2 // this is the diffrenece between summation area of temp_inner_sum and of  temp_outer_sum
	variable temp_inner_sum = 0 // for summation of OD
	variable temp_outer_sum = 0 // for summation of OD
	variable range4back = 75 // region for background
	variable back_sum = 0 // for summation of background OD
	variable magnification = 2 // magnification of image
	variable pixelsize = 5.6 // pixel size of the camera 
	variable OD = 0 // Summation of OD of image of interest
	variable avg_OD = 0 //
	variable i, j
	variable threshold
	

	
	// background range
//	variable x1_back = 180
//	variable x2_back = 380
//	variable y1_back = 80
//	variable y2_back = 280
	
	//Histogram bin setting
	variable hist_offset = -2
	variable hist_binsize = 0.001
	variable hist_numbin = 7000
			
	NewDataFolder/O root:fittings
	SetDataFolder root:fittings
	make/O/N=7 parameters

	//Histogram algo
	//make/O/N=(DimSize(optdepth,0),DimSize(optdepth,1)) roiwave
	//roiwave[x1,x2][y1,y2] = 0
	if(resize==1)
		duplicate/o/R=(x1,x2)(y1,y2) optdepth ROI
	else
		duplicate/o optdepth ROI
	endif
	Make/N=(hist_numbin)/O ROI_hist
	Histogram/B={hist_offset,hist_binsize,hist_numbin} ROI, ROI_hist
	
	make/O/N=(DimSize(ROI_hist,0)) ROI_hist_multi
	ROI_hist_multi = ROI_Hist[p]*(p*hist_binsize+hist_offset)
	Integrate ROI_hist_multi/D=ROI_hist_multi_INT
	
	Make/N=(hist_numbin)/O ROI_hist_count
	ROI_hist_count = (ROI_hist[x]>0) ? ROI_hist[x] : 0
	Integrate ROI_hist_count/D=ROI_hist_count_INT
	
	K0 = 0;
	CurveFit/H="1000"/Q/NTHR=0/TBOX=768 gauss  ROI_hist_multi_INT[0,-hist_offset/hist_binsize+50]
	variable borderpoint = floor(sqrt(-ln((-0.1-wcoef[0])/wcoef[1])*wcoef[3]^2)+wcoef[2])
	
	//variable yy=2100
	if(WaveMax(ROI_hist_multi_INT) > 0)
		variable pixel_counted = WaveMax(ROI_hist_count_INT) - ROI_hist_count_INT[borderpoint]
		variable background_point_counted = ROI_hist_count_INT[borderpoint] - ROI_hist_count_INT[-hist_offset/hist_binsize-1]
		variable avg_backOD = ROI_hist_multi_INT[borderpoint]/background_point_counted
		OD = WaveMax(ROI_hist_multi_INT) - ROI_hist_multi_INT[borderpoint] - pixel_counted * avg_backOD
	else
		OD = 0
	endif
	//CurveFit/Q/NTHR=0 gauss ROI_hist[-hist_offset/hist_binsize*3/4,-hist_offset/hist_binsize*5/4]
	//avg_OD = wcoef[2]
	//print pixel_counted, background_point_counted, avg_backOD
	
	// Find center of the cloud
//	ImageStats/G={x1,x2,y1,y2} root:Flea3:optdepth // prevent back dead point in CCD's effect
//	if(V_max > 1)
//		rMaxLoc = V_maxRowLoc
//		cMaxLoc = V_maxColLoc
//		threshold = 0.95
//	else
		//CurveFit/NTHR=0/Q Gauss2D  optdepth[x1,x2][y1,y2] 	
		//rMaxLoc = ceil((wcoef[2]- DimOffset(optdepth, 0))/DimDelta(optdepth,0))
		//cMaxLoc = ceil((wcoef[4]- DimOffset(optdepth, 1))/DimDelta(optdepth,1))
		//if (rMaxLoc < x1 || rMaxLoc > x2 || cMaxLoc < y1 || cMaxLoc > y2)
//			rMaxLoc = 220
//			cMaxLoc = 433
		//endif
//		threshold = 0.86
//	endif
//	print rMaxLoc, cMaxLoc, threshold
	// Find region of interest
//	do 
//		temp_inner_sum = 0
//		temp_outer_sum = 0
//		for (i=-range; i<=range; i=i+1)
//			for (j=-range; j<=range; j=j+1)
//				temp_inner_sum = temp_inner_sum + optdepth[rMaxLoc+i][cMaxLoc+j]
//			endfor
//		endfor
//		for (i=-range-step; i<=range+step; i=i+1)
//			for (j=-range-step; j<=range+step; j=j+1)
//				temp_outer_sum = temp_outer_sum + optdepth[rMaxLoc+i][cMaxLoc+j]
//			endfor
//		endfor
//		range = range - 1
//	while (temp_inner_sum > temp_outer_sum * threshold)
//	range = range + 1 +step
	// Duplicate image of interest to ROI
//	duplicate/o/R=(DimOffset(optdepth, 0) + (rMaxLoc-range)*1.4, DimOffset(optdepth, 0) + (rMaxLoc+range)*1.4)(DimOffset(optdepth, 1) + (cMaxLoc-range)*1.4,DimOffset(optdepth, 1) + (cMaxLoc+range)*1.4) optdepth ROI

	// Store sum of OD in ROI, here we use  temp_outer_sum as sum OD
//	OD = temp_outer_sum

	// Calculate average of backgroud OD by selecting a region next to the image of interest 
//	back_sum = 0
//	for (i=x1_back; i<x2_back; i=i+1)
//		for (j=y1_back; j<y1_back; j=j+1)
//			back_sum = back_sum + optdepth[i][j]
//		endfor
//	endfor
//	avg_OD =  back_sum/((x2_back-x1_back)*(y2_back-y1_back))
	//print rMaxLoc, cMaxLoc, threshold, rMaxLoc+range+2*range4back, cMaxLoc-range-2*range4back
	
	// Calculate OD of interest
	// row maximum location, column maximum location, range, average background OD, summation of OD, atom number
	//make/O/N=6 parameters={DimOffset(optdepth, 0)+rMaxLoc*1.4, DimOffset(optdepth, 1)+cMaxLoc*1.4, range*1.4, avg_OD, OD, (OD - avg_OD * range^2)*(pixelsize/magnification)^2/0.23}
	//make/O/N=3 parameters={avg_OD, pixel_counted, OD, (OD - avg_OD*pixel_counted)*(pixelsize/magnification)^2/0.23}
	
	//return atom number
	//return (OD - avg_OD *pixel_counted)*(pixelsize/magnification)^2/0.23
	//return OD*(pixelsize/magnification)^2/0.23
	return OD*(pixelsize/magnification)^2/(3*0.78^2/2/3.14159/(1+1.9^2*4)/(1+0.5/(1+1.9^2*4)))
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
