#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later
#include <Waves Average>
#include "PXPUtils"

Function AnalyseMF041043()
	LoadWave/O/A/J/K=0/L={0,1,0,0,0}/W
	// we want to grab
	// ATG9 and TPD54
	// rapalog and norapalog
	// need to cycle through these four combos and extract data for each uid that matches
	// just need to decide what column we take (tt is always column 0, data in column 1
	String condAstr = "ATG9;TPD54;"
	String condBstr = "rapalog;norapalog;"
	WAVE/Z/T id, Channel, Condition
	FindDuplicates/FREE/RT=uid id
	Make/O/N=(numpnts(uid),5)/T theKeyW = ""
	theKeyW[][0] = uid[p]
	Wave/Z frac_mean_mito, tt
	String wName, condA, condB, plotName
	String yList, xList, avName, errName, avList
	
	Variable i,j,k
	
	for(i = 0; i < ItemsInList(condAstr); i += 1)
		condA = StringFromList(i,condAstr)
		for(j = 0; j < ItemsInList(condBstr); j += 1)
			condB = StringFromList(j,condBstr)
			plotName = "p_" + condA + "_" + condB
			KillWindow/Z $plotName
			Display/N=$plotName
			for(k = 0; k < numpnts(uid); k += 1)
				wName = condA + "_" + condB + "_" + num2str(k)
				Duplicate/O frac_mean_mito, $(wName + "_y")
				Duplicate/O tt, $(wName + "_x")
				Wave yW = $(wName + "_y")
				Wave xW = $(wName + "_x")
				yW[] = (cmpstr(Channel[p],condA) == 0 && cmpstr(Condition[p],condB) == 0 && cmpstr(id[p],uid[k]) == 0) ? yW[p] : NaN
				xW[] = (cmpstr(Channel[p],condA) == 0 && cmpstr(Condition[p],condB) == 0 && cmpstr(id[p],uid[k]) == 0) ? xW[p] : NaN
				WaveTransform zapnans xW
				WaveTransform zapnans yW
				if(numpnts(xW) == 0 || k == 67)
					KillWaves/Z xW,yW
				else
					AppendToGraph/W=$plotName yW vs xW
					// now normalise to mean of points 1,2,3
					WaveStats/Q/M=1/RMD=[1,3] yW
					yW[] = yW[p] / V_avg
					theKeyW[k][1 + (2 * i) + (j)] = wName
				endif
			endfor
			yList = Wavelist("*_y",";","WIN:"+ plotName)
			xList = ReplaceString("_y",yList,"_x")
			avName = ReplaceString("p_", plotName, "W_Ave_")
			errName = ReplaceString("Ave", avName, "Err")
			fWaveAverage(yList, xList, 3, 1, avName, errName)
			AppendToGraph/W=$plotName $avName
			Label/W=$plotName left "Mitochondrial signal (Fold-change)"
			Label/W=$plotName bottom "Time (s)"
			ErrorBars/W=$plotName $avName SHADE= {0,0,(0,0,0,0),(0,0,0,0)},wave=($ErrName,$ErrName)
			ModifyGraph/W=$plotName lsize($avName)=2,rgb($avName)=(0,0,0)
			SetAxis/W=$plotName/A/N=1 left
		endfor
	endfor
End

// do batch curve fit and then

Function PlotTheFits()
	Wave/Z matA = batchCurveFitRun0_ResultsCopy
	String wList = WaveList("ATG9_rapa*_y",";","")
	String resultList = PXPUtils#GetDimLabelsFromWave(matA, 0)
	Variable nWaves = ItemsInList(wList)
	Make/O/N=(nWaves)/T Results_Names
	Make/O/N=(nWaves,2) Results_Taus
	Make/O/N=(nWaves,2) Results_Sigmas
	String wName
	Variable theRow
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		Results_Names[i] = wName
		theRow = WhichListItem(wName, resultList)
		if(theRow != -1)
			Results_Taus[i][0] = matA[theRow][%tau]
			Results_Sigmas[i][0] = matA[theRow][%'tau Sigma']
		else
			Results_Taus[i][0] = NaN
			Results_Sigmas[i][1] = NaN
		endif
		theRow = WhichListItem(ReplaceString("ATG9",wName,"TPD54"), resultList)
		if(theRow != -1)
			Results_Taus[i][1] = matA[theRow][%tau]
			Results_Sigmas[i][1] = matA[theRow][%'tau Sigma']
		else
			Results_Taus[i][1] = NaN
			Results_Sigmas[i][1] = NaN
		endif
	endfor
	Results_Taus[][] = (Results_Taus[p][q] > 200) ? NaN : Results_Taus[p][q]
	KillWindow/Z p_theFits
	Display/N=p_theFits Results_Taus[][1] vs Results_Taus[][0]
	ModifyGraph mode=3
	SetAxis left 0,200
	SetAxis bottom 0,200
End	