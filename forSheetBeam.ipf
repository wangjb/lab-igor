function ListFigures(aaa,bbb,ccc)
	variable aaa starting file number
	variable bbb ending file number
	string ccc either insitu or tof
	string sss  formated filenumber
	string figname  name of the figures
	variable i
		
	wave ramptime = rootFlea3IndexedWavesramptime
	wave pos = rootsheetpos
	wave power = rootsheetpower
	wave tilt= rootsheettilt
	wave para = rootsheetparameters
	
	NewDataFolderO rootgfigures
	NewDataFolderO rootfittings
	
	 If want to include atom number, uncomment these lines
	wave p = rootfittingsparameters
	makeON=(1,6) ODs
	
	for(i = aaa-1000; i = bbb-1000; i = i+1)
		print i
		if(i  29) 
			sprintf sss,%04.0f, (i+1000)
		elseif(i == 29)
			continue
		else
			sprintf sss,%04.0f, (i+1000)
			i = i-1
		endif
		print sss
		AutoRunV3(Flea3,JJB's experimentrawdataFlea3_21Oct2016_+sss+.ibw) 
		N_ring(rootFlea3optdepth)
		SetDataFolder rootgfigures
		if(stringmatch(ccc,insitu) && floor((i-1)16) == 0)
			figname = a20161021+sss+_insitu
			 Modify the output pictures here
			duplicateoR=(-100,-30)(70,140) rootFlea3optdepth $figname
			duplicateoR=(-105,-45)(90,150) rootFlea3optdepth $figname
			duplicateoR=(-400,400)(-200,-600) rootFlea3optdepth $figname
			
			Fittings for atom number and adiabatic radius
			SetDataFolder rootfittings
			imagetransformG=400 getcol rootFlea3optdepth
			duplicateo rootfittingsW_ExtractedCol rootgfigures$(h_a20160921+sss+_insitu)
			imagetransformG=197 getrow rootFlea3optdepth
			duplicateo rootfittingsW_ExtractedRow rootgfigures$(v_a20160921+sss+_insitu)
			SetDataFolder rootfittings
			Create image
			NewImage  $figname
			SetAxisR left 70,140
			SetAxisR left 90,150
			SetAxisR left -200,-600
			ModifyGraph height={Aspect,1}
			ModifyGraph height={Aspect,0.5}
			ModifyImage $figname ctab= {-0.1,2,Rainbow,1}
			
			Edit the textbox here
			TextBoxCN=text0 Z0720161014+sss+ mf=-1, in situr ramptime=+num2str(ramptime[i-aaa]10^3)+ms

			TextBoxCN=text0 Z0720161021+sss+ pos=+num2str(pos[floor((i-1)16)])+rpower=+num2str(power[mod(i,4)])+, tilt=+num2str(tilt[mod(floor((i-1)4),4)])

			TextBoxCN=text0 Z0720161021+sss+ pos=+num2str(pos[floor((i-1-1)16)])+rpower=+num2str(power[mod(i-1,4)])+, tilt=+num2str(tilt[mod(floor((i-1-1)4),4)])

			print mod(i-1000-1,4)
			print mod(floor((i-1000-1)4),4)
			TextBoxCN=text0X=10Y=80
			ModifyGraph grid(top)=2,gridRGB(top)=(0,0,0)
		elseif(stringmatch(ccc,tof))
			figname = a20161014+sss+_tof
			
			 Modify the output pictures here
			duplicateoR=(-200,0)(0,200) rootFlea3optdepth $figname
			duplicateoR=(-230,-30)(0,200) rootFlea3optdepth $figname
			
			Create image
			NewImage $figname
			SetAxisR left 0,200
			ModifyGraph height={Aspect,1}
			ModifyImage $figname ctab= {-0.1,2,Rainbow,1}
			Edit the textbox here
			TextBoxCN=text0 Z0720161014+sss+ mf=-1, TOF 16msr ramptime=+num2str(ramptime[i-aaa]10^3)+ms
			TextBoxCN=text0X=15.00Y=80
			TextBoxCN=text0X=2Y=80
		endif
		N_ring(rootFlea3optdepth)
		matrixopo p = p^t
		concatenateONP=0 {ODs,p}, temp_w
		duplicateO temp_w, ODs
		if(i = 29) 
			i = i+1
		endif
	endfor
	DeletePoints 0,1, ODs
	KillWaves temp_w
	Execute TileWindowsG=30O=1W=(14,15,1382,695)
end
