;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;Editor: Eric Xu;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;Date: 2015.04.01;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
;**************************************************************************
; Define parameters
;**************************************************************************
   if(.not.isvar("fi_data")) then
   	fi_data = "/mnt/CMIP5/cmip5_tmp/1010/ori_data"
   end if
   if(.not.isvar("fi_name")) then
   	fi_name = "sic_OImon_FGOALS-g2_historicalGHG_r1i1p1_192701-199112_ori.nc"
   end if
   if(.not.isvar("fig_data")) then
   	fig_data = "/mnt/CMIP5/cmip5_tmp/1010/fig_data"
   end if
   if(.not.isvar("fig_name")) then
   	fig_name = "sic_OImon_FGOALS-g2_historicalGHG_r1i1p1_192701-199112_EOF_NH_60_90"
   end if
   if(.not.isvar("nc_data")) then
   	nc_data = "/mnt/CMIP5/cmip5_tmp/1010/nc_data"
   end if
   if(.not.isvar("nc_name")) then
   	nc_name = "sic_OImon_FGOALS-g2_historicalGHG_r1i1p1_192701-199112_EOF_NH_60_90.nc"
   end if
   if(.not.isvar("txt_data")) then
   	txt_data = "/mnt/CMIP5/cmip5_tmp/1010/txt_data"
   end if
   if(.not.isvar("txt_name")) then
   	txt_name = "sic_OImon_FGOALS-g2_historicalGHG_r1i1p1_192701-199112_EOFTS_NH_60_90.txt"
   end if
   if(.not.isvar("var_name")) then
   	var_name="sic"
   end if
   if(.not.isvar("fig_type")) then
   	fig_type= "png"
   end if
   if(.not.isvar("latmin")) then
   	latmin=60
   end if
   if(.not.isvar("latmax")) then
   	latmax=90
   end if
;**************************************************************************
; read in netCDF file and computation
;**************************************************************************
   fi  =   addfile(fi_data+"/"+fi_name,"r")
   var =   fi->$var_name$                       ; Y0.01-YN.12
   lat =   fi->lat
   lon =   fi->lon
   ntime   = dimsizes(var&time)
   ny      = dimsizes(var&lat)
   nx      = dimsizes(var&lon)
   
   ;printVarSummary(var)
   ;TIME     = fi->time(12:ntime-1) ; exclude the past year, Y1.01-YN.12
   ;YYYYMM   = cd_calendar(TIME,-1)
   ;YYYY=YYYYMM/100  ; Y1.01-YN.12
   
;**************************************************************************
;  Monthly data pre-processing
;**************************************************************************  
   var0_reordered = var(time|12:ntime-1,lat|:,lon|:) ;Y1.01-YN.12
   x = var0_reordered(:,{latmin:latmax},:)
   x = dtrend_leftdim(x,False)    ; 消除最左边一维的最小方差线性变化趋势的影响
   xClm      = clmMonTLL(x)       ; 计算多年的月平均值,xClm{J, F, M, A, M, J, J, A, S, O, N, D}
   ;printVarSummary(xClm)
;printVarSummary(x)
   
   xAnom           = calcMonAnomTLL (x, xClm)      ; 异常场(消除季节循环), Y1d.01-YNd.12
   xAnom@long_name = "ANOMALIES: "+var@long_name
   ;printVarSummary(xAnom)
   
   delete(var0_reordered)
   delete(xClm)  
; =================================================================
; create weights:  sqrt(cos(lat))   [or sqrt(gw) ]
; =================================================================
  rad    = 4.*atan(1.)/180.    ;pi/180(rad/度)
  clat   = x&lat
  clat   = cos(rad*clat)
  clat   = where(clat.lt.0, 0.0, clat)  ; avoid a potential numerical issue at pole  
  clat   = sqrt( clat )
  printVarSummary(clat)
; =================================================================
; weight all observations 
; =================================================================
  clat0=tofloat(clat)
  xw     = xAnom*conform(xAnom, clat0, 1)  ;对异常场进行加权
  copy_VarMeta(xAnom, xw)
  xw@long_name = "Wgt: "+var@long_name
  printVarSummary(xw)
  ;xw!0="time"
  ;xw!1="lat"
  ;xw!2="lon"
  ;printVarSummary(xw)
  wx1 = xw(lat|:,lon|:,time|:)      ; Reorder (lat,lon,time), Y1w.01-YNw.12
  printVarSummary(wx1)
  delete(xAnom)
  
;**************************************************************************
; Summer data pre-processing
;**************************************************************************
  xw2 = month_to_season(xw, "JJA")  ; ntime = NYEAR
  wx2 = xw2(lat|:,lon|:,time|:)     ; Reorder (lat,lon,time), Y1w.JJA-YNw.JJA
  ;printVarSummary(wx2)
  delete(xw)
  
;**************************************************************************
; Winter data pre-processing
;************************************************************************** 
   var1_reordered=var(time|11:ntime-2,lat|:,lon|:) ;Y0.12-YN.11
   x=var1_reordered(:,{latmin:latmax},:)
   x=dtrend_leftdim(x,False)
   xClm      = clmMonTLL(x)      ; xClm{D, J, F, M, A, M, J, J, A, S, O, N}
   ;printVarSummary(xClm)
   
   xAnom           = calcMonAnomTLL (x, xClm)    ; Y0d.12-YNd.11
   xAnom@long_name = "ANOMALIES: "+var@long_name
   ;printVarSummary(xAnom)
   
   delete(var1_reordered)
   delete(xClm)
; =================================================================
; weight all observations 
; =================================================================
  clat0=tofloat(clat)
  xw     = xAnom*conform(xAnom, clat0, 1)  ; Y0d.12-YNd.11 (time,lat,lon)
  copy_VarMeta(xAnom, xw)
  xw@long_name = "Wgt: "+var@long_name
  ;xw!0="time"
  ;xw!1="lat"
  ;xw!2="lon"
  printVarSummary(xw)
  dim = dimsizes(xw)
; =================================================================
; winter average
; =================================================================  
  xw_1d = ndtooned(xw)
  xw_4d = onedtond(xw_1d,(/dim(0)/12,12,dim(1),dim(2)/))
  xw_4d_winter   = dim_avg_n(xw_4d(:,0:2,:,:),1) ; ntime = NYEAR
  copy_VarMeta(xw2, xw_4d_winter)
  ;xw_4d_winter!0="time"
  ;xw_4d_winter!1="lat"
  ;xw_4d_winter!2="lon"
  wx3 = xw_4d_winter(lat|:,lon|:,time|:)      ; Y1w.DJF-YNw.DJF
  printVarSummary(wx3)
  
  delete(xAnom)
  delete(xw)

;**************************************************************************
;  Monthly EOF fields
;************************************************************************** 
   optEof = False                 ; optEOF@jopt = 0使用协方差矩阵
   optETS = False                 ; optETS@jopt = 0表示使用元数据wx和eof
   neof = 4
   
   eof1    = eofunc_Wrap(wx1, neof, optEof)
  ; printVarSummary(eof1)
   eof1_ts = eofunc_ts_Wrap( wx1, eof1,  optETS )   ; Y1ts.01-YNts.12
   eof1_ts = dim_standardize_n( eof1_ts, 0, 1) ;主成分时间序列标准化
   printVarSummary(eof1_ts)
   ;print(eof1_ts&time)
   
   YYYYMM = cd_calendar(eof1_ts&time, -1)
   YYYY = YYYYMM/100                                ; ntime = NYEAR*12
   year = YYYY(::12)
   ;print(year)
;**************************************************************************
;  Summer EOF fields
;************************************************************************** 
   eof2    = eofunc_Wrap(wx2, neof, optEof)
  ;printVarSummary(eof2)
   eof2_ts = eofunc_ts_Wrap( wx2, eof2,  optETS )   ; Y1ts.JJA-YNts.JJA
   eof2_ts = dim_standardize_n( eof2_ts, 0, 1)
  ; printVarSummary(eof2_ts)   

;**************************************************************************
;  Winter EOF fields
;************************************************************************** 
   eof3    = eofunc_Wrap(wx3, neof, optEof)
  ; printVarSummary(eof3)
   eof3_ts = eofunc_ts_Wrap( wx3, eof3,  optETS )   ; Y1ts.DJF-YNts.DJF
   eof3_ts = dim_standardize_n( eof3_ts, 0, 1)
   printVarSummary(eof3_ts)

;**************************************************************************
; Creating plots
;**************************************************************************
  wks = gsn_open_wks(fig_type, fig_data+"/"+fig_name) 
  gsn_define_colormap(wks,"BlueRed")
  plot = new(neof,graphic)
  
  res                      = True         
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
  res@gsnPolar             = "NH"
  res@gsnAddCyclic         = False
  res@gsnSpreadColors      = True         ; spread out color table

  res@mpFillOn             = False        ; turn off map fill
  res@mpMinLatF            = latmin         ; zoom in on map
  res@mpMaxLatF            = latmax
  res@mpGeophysicalLineThicknessF =2
  res@gsnLeftStringFontHeightF =0.022
  res@gsnRightStringFontHeightF =0.022
  res@gsnPolarLabelFontHeightF = .02

  res@cnFillOn             = True         ; turn on color fill
  res@cnLinesOn            = False        ; True is default
  res@cnLineLabelsOn       = False        ; True is default
  res@lbLabelBarOn         = False        ; turn off individual lb's
  
  symMinMaxPlt(eof1, 16, False, res)       ; set symmetric plot min/max
  ; panel plot only resources
  resP                     = True         ; modify the panel plot
  resP@gsnMaximize         = True         ; large format
  resP@gsnPanelLabelBar    = True         ; add common colorbar
  resP@lbLabelAutoStride   = True         ; auto stride on labels
  resP@txString             = var@long_name+": "+tostring(min(YYYY))+"-"+tostring(max(YYYY))+ " ("+tostring(latmin)+"N - "+tostring(latmax)+"N)"

;=================================================================
;  Monthly EOFs Spatial Patterns
;=================================================================
  do n=0,neof-1
	res@gsnLeftString  = "EOF "+(n+1)
	res@gsnRightString = sprintf("%5.1f", eof1@pcvar(n)) +"%"
	plot(n)=gsn_csm_contour_map_polar(wks,eof1(n,:,:),res)
  end do
  gsn_panel(wks,plot(0:3),(/2,2/),resP)     ; draw all 'neof' as one plot

;=================================================================
;  Summer EOFs Spatial Patterns
;=================================================================
  symMinMaxPlt(eof2, 16, False, res) 
  resP@txString             = var@long_name+": JJA "+tostring(min(YYYY))+"-"+tostring(max(YYYY))+ " ("+tostring(latmin)+"N - "+tostring(latmax)+"N)" 
  do n=0,neof-1
	res@gsnLeftString  = "EOF "+(n+1)
	res@gsnRightString = sprintf("%5.1f", eof2@pcvar(n)) +"%"
	plot(n)=gsn_csm_contour_map_polar(wks,eof2(n,:,:),res)
  end do
  gsn_panel(wks,plot(0:3),(/2,2/),resP)	 

;=================================================================
;  Winter EOFs Spatial Patterns
;=================================================================
  symMinMaxPlt(eof3, 16, False, res) 
  resP@txString             = var@long_name+": DJF "+tostring(min(YYYY))+"-"+tostring(max(YYYY))+ " ("+tostring(latmin)+"N - "+tostring(latmax)+"N)" 
  do n=0,neof-1
	res@gsnLeftString  = "EOF "+(n+1)
	res@gsnRightString = sprintf("%5.1f", eof3@pcvar(n)) +"%"
	plot(n)=gsn_csm_contour_map_polar(wks,eof3(n,:,:),res)
  end do
  gsn_panel(wks,plot(0:3),(/2,2/),resP)
  
;**************************************************************************
; Time Series (principal component) plots
;**************************************************************************
  rts           = True
  rts@gsnDraw   = False       ; don't draw yet
  rts@gsnFrame  = False       ; don't advance frame yet
  rts@gsnScale  = True        ; force text scaling
  rts@gsnMaximize  = True 

  rts@vpHeightF = 0.40        ; Changes the aspect ratio
  rts@vpWidthF  = 0.90
  rts@vpXF      = 0.05        ; change start locations
  rts@vpYF      = 0.75        ; the plot
  
; =================================================================
; these four resources allow the user to stretch the plot size, and
; decide exactly where on the page to draw it.
; =================================================================
  rts@gsnYRefLine           = 0.              ; reference line   
  ;rts@gsnXYBarChart         = True            ; create bar chart 
  rts@gsnAboveYRefLineColor = "red"           ; above ref line fill red
  rts@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue
  rts@gsnXYBarChartOutlineThicknessF  = 0.
  rts@tmXTOn  = False
  rts@tmYROn  = False
  
  rts@trXMinF = min(YYYY)
  rts@trXMaxF = max(YYYY)
  ;rts@tmXBMode  = "Manual"
  ;rts@tmXBTickStartF = min(year)
  ;rts@tmXBTickEndF = max(year)
  ;rts@tmXBTickSpacingF= 5

  rtsP                     = True             ; modify the panel plot
  rtsP@gsnMaximize         = True             ; large format
  
  rtsP@gsnPanelFigureStrings = (/"a)","b)","c)","d)"/)
  rtsP@gsnPanelFigureStringsPerimOn  = False
  rtsP@amJust     = "BottomRight"
  rtsP@gsnPanelFigureStringsFontHeightF  = 0.01
  
  rtsP@txString            = var@long_name+": "+tostring(min(YYYY))+"-"+tostring(max(YYYY))+ " ("+tostring(latmin)+"N - "+tostring(latmax)+"N)"  
;=================================================================
;  Monthly EOFs Time Series
;=================================================================
  eof1_ts@long_name = "Standardized"
  do n=0,neof-1
     rts@gsnLeftString  = "EOF "+(n+1)
     rts@gsnRightString = sprintf("%5.1f",eof1@pcvar(n)) +"%"
     plot(n) = gsn_csm_xy (wks,YYYY,eof1_ts(n,:),rts)
  end do
  gsn_panel(wks,plot(0:3),(/2,2/),rtsP)        ; draw all 'neof' as one plot
  
;=================================================================
;  Summer EOFs Time Series
;=================================================================
  eof2_ts@long_name = "Standardized"
  rts@gsnXYBarChart         = True            ; create bar chart 
  rtsP@txString             = var@long_name+": JJA "+tostring(min(YYYY))+"-"+tostring(max(YYYY))+ " ("+tostring(latmin)+"N - "+tostring(latmax)+"N)" 
  do n=0,neof-1
     rts@gsnLeftString  = "EOF "+(n+1)
     rts@gsnRightString = sprintf("%5.1f",eof2@pcvar(n)) +"%"
     plot(n) = gsn_csm_xy (wks,year,eof2_ts(n,:),rts)
  end do
  gsn_panel(wks,plot(0:3),(/2,2/),rtsP)

;=================================================================
;  Winter EOFs Time Series
;=================================================================
  eof3_ts@long_name = "Standardized"
  rtsP@txString             = var@long_name+": DJF "+tostring(min(YYYY))+"-"+tostring(max(YYYY))+ " ("+tostring(latmin)+"N - "+tostring(latmax)+"N)" 
  do n=0,neof-1
     rts@gsnLeftString  = "EOF "+(n+1)
     rts@gsnRightString = sprintf("%5.1f",eof3@pcvar(n)) +"%"
     plot(n) = gsn_csm_xy (wks,year,eof3_ts(n,:),rts)
  end do
  gsn_panel(wks,plot(0:3),(/2,2/),rtsP)

;**************************************************************************
; export EOF nc data
;**************************************************************************
  system("rm -f " + nc_data+"/"+nc_name)
  fout= addfile(nc_data+"/"+nc_name,"c")
  setfileoption(fout,"DefineMode",True)
  
  fAtt               = True            ; assign file attributes
  fAtt@title         = "EOF of "+var@long_name+": NH : "+tostring(min(YYYY))+"-"+tostring(max(YYYY)) 
  fAtt@source_file   = fi_name
  fAtt@Conventions   = "None"   
  fAtt@creation_date = systemfunc ("date")        
  fileattdef( fout, fAtt )            ; copy file attributes  
  
  nlon=dimsizes(eof1&lon)
  nlat=dimsizes(eof1&lat)
  ;evn=eof1&evn
  ;print(evn)
  ;evn2=eof2&evn
  ;evn3=eof3&evn

  dimNames = (/"evn", "lat", "lon"/)  
  dimSizes = (/neof,  nlat,  nlon/)
  dimUnlim = (/False, False, False/)
  ;print(dimSizes)
  filedimdef(fout,dimNames,dimSizes,dimUnlim)
  ;printVarSummary(eof1)

  ;filevardef(fout, "evn"  ,typeof(evn),getvardims(evn))
  filevardef(fout, "lat"  ,typeof(eof1&lat),getvardims(eof1&lat))
  filevardef(fout, "lon"  ,typeof(eof1&lon),getvardims(eof1&lon))
  filevardef(fout, "eof1"  ,typeof(eof1)  ,getvardims(eof1))
  filevardef(fout, "eof2"  ,typeof(eof2)  ,getvardims(eof2))
  filevardef(fout, "eof3"  ,typeof(eof3)  ,getvardims(eof3))

  ;filevarattdef(fout,"evn"  ,evn)
  filevarattdef(fout,"lat"  ,eof1&lat)
  filevarattdef(fout,"lon"  ,eof1&lon)
  eof1@long_name = "monthly EOFs"
  filevarattdef(fout,"eof1"  ,eof1)
  ;printVarSummary(eof1)
  eof2@long_name = "summer EOFs"
  filevarattdef(fout,"eof2"  ,eof2)
  eof3@long_name = "winter EOFs"
  filevarattdef(fout,"eof3"  ,eof3)
  setfileoption(fout,"DefineMode",False)

  fout->lat    = (/eof1&lat/)
  fout->lon    = (/eof1&lon/)
  fout->eof1      = (/eof1/)
  fout->eof2      = (/eof2/)
  fout->eof3      = (/eof3/)
   
;**************************************************************************
; export time series txt data
;************************************************************************** 
  header1 = (/"-------------------------------------------", \
               "time    EOF1    EOF2    EOF3    EOF4", \
               "-------------------------------------------"/)
  footer1 = (/"--------------------------------------------------", \
            "EOF1:"+sprintf("%5.1f", eof1@pcvar(0)) +"%"+"  EOF2:"+sprintf("%5.1f", eof1@pcvar(1)) +"%"+"  EOF3:"+sprintf("%5.1f", eof1@pcvar(2)) +"%"+"  EOF4:"+sprintf("%5.1f", eof1@pcvar(3)) +"%", \
            "--------------------------------------------------"/)
  alist1 = [/YYYYMM,eof1_ts(0,:),eof1_ts(1,:),eof1_ts(2,:),eof1_ts(3,:)/]
  hlist1 = [/header1/]
  flist1 = [/footer1/]
  write_table(txt_data+"/"+txt_name, "w", hlist1, "%s")
  write_table(txt_data+"/"+txt_name, "a", alist1, "%d %8.3f %8.3f %8.3f %8.3f")
  write_table(txt_data+"/"+txt_name, "a", flist1, "%s")

  header2 = (/"-------------------------------------------", \
               "JJA     EOF1    EOF2    EOF3    EOF4", \
               "-------------------------------------------"/)
  footer2 = (/"--------------------------------------------------", \
            "EOF1:"+sprintf("%5.1f", eof2@pcvar(0)) +"%"+"  EOF2:"+sprintf("%5.1f", eof2@pcvar(1)) +"%"+"  EOF3:"+sprintf("%5.1f", eof2@pcvar(2)) +"%"+"  EOF4:"+sprintf("%5.1f", eof2@pcvar(3)) +"%", \
            "--------------------------------------------------"/)
  alist2 = [/year,eof2_ts(0,:),eof2_ts(1,:),eof2_ts(2,:),eof2_ts(3,:)/]
  hlist2 = [/header2/]
  flist2 = [/footer2/]  
  write_table(txt_data+"/"+txt_name, "a", hlist2, "%s")
  write_table(txt_data+"/"+txt_name, "a", alist2, "%d %8.3f %8.3f %8.3f %8.3f")
  write_table(txt_data+"/"+txt_name, "a", flist2, "%s")

  header3 = (/"-------------------------------------------", \
               "DJF     EOF1    EOF2    EOF3    EOF4", \
               "-------------------------------------------"/)
  footer3 = (/"--------------------------------------------------", \
            "EOF1:"+sprintf("%5.1f", eof3@pcvar(0)) +"%"+"  EOF2:"+sprintf("%5.1f", eof3@pcvar(1)) +"%"+"  EOF3:"+sprintf("%5.1f", eof3@pcvar(2)) +"%"+"  EOF4:"+sprintf("%5.1f", eof3@pcvar(3)) +"%", \
            "--------------------------------------------------"/)
  alist3 = [/year,eof3_ts(0,:),eof3_ts(1,:),eof3_ts(2,:),eof3_ts(3,:)/]
  hlist3 = [/header3/]
  flist3 = [/footer3/]  
  write_table(txt_data+"/"+txt_name, "a", hlist3, "%s")
  write_table(txt_data+"/"+txt_name, "a", alist3, "%d %8.3f %8.3f %8.3f %8.3f")
  write_table(txt_data+"/"+txt_name, "a", flist3, "%s")
  print("END")
end
