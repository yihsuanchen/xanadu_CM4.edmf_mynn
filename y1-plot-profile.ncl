load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/ncrc/home2/Yi-hsuan.Chen//script/ncl//ignorantNCL.ncl"
;load "./ignorantNCL.ncl"
load "/ncrc/home2/Yi-hsuan.Chen//script/ncl//ignorantNCL-colormaps.ncl"
;load "zz-read_parameters.ncl"

begin
  kx = 33
  tx = 4

  th = new( (/tx,kx/), "float"  )

  z_full = (/   37724.035  , 32227.250  , 28718.551  , 26018.232  , 23780.551  , 21804.475  , 19977.420  , 18260.461  , 16643.578  , 15120.780  , 13686.597  , 12335.868  , 11041.139  ,  9768.828  ,  8514.327  ,  7294.708  ,  6141.825  ,  5095.068  ,  4182.057  ,  3407.400  ,  2758.814  ,  2218.686  ,  1770.327  ,  1399.202  ,  1092.666  ,   839.890  ,   632.109  ,   462.201  ,   324.340  ,   214.095  ,   127.775  ,    62.942  ,    17.861  /)
  
  th (0,:) = (/    1154.470  ,   869.131  ,   730.368  ,   639.500  ,   573.057  ,   520.270  ,   475.841  ,   437.562  ,   404.337  ,   375.358  ,   349.968  ,   327.631  ,   318.730  ,   316.250  ,   313.803  ,   311.425  ,   309.177  ,   307.136  ,   305.355  ,   303.845  ,   302.579  ,   301.527  ,   300.652  ,   299.929  ,   299.331  ,   298.838  ,   298.433  ,   298.101  ,   297.832  ,   297.618  ,   297.449  ,   297.322  ,   297.234  /) 
  th (1,:) = (/    1154.470  ,   869.131  ,   730.368  ,   639.500  ,   573.057  ,   520.270  ,   475.841  ,   437.562  ,   404.337  ,   375.358  ,   349.968  ,   327.631  ,   318.730  ,   316.250  ,   313.803  ,   311.425  ,   309.177  ,   307.136  ,   305.355  ,   303.845  ,   302.579  ,   301.527  ,   300.652  ,   299.929  ,   299.331  ,   298.838  ,   298.433  ,   298.101  ,   297.832  ,   297.618  ,   297.449  ,   297.347  ,   303.268  /)
  
  th (2,:) = (/    1154.470  ,   869.131  ,   730.368  ,   639.500  ,   573.057  ,   520.270  ,   475.841  ,   437.562  ,   404.337  ,   375.358  ,   349.968  ,   327.631  ,   318.730  ,   316.250  ,   313.803  ,   311.425  ,   309.177  ,   307.136  ,   305.355  ,   303.845  ,   302.579  ,   301.527  ,   300.652  ,   299.929  ,   299.331  ,   298.838  ,   298.433  ,   298.101  ,   297.832  ,   297.618  ,   297.465  ,   298.128  ,   308.233  /)
  th (3,:) = (/    1154.470  ,   869.131  ,   730.368  ,   639.500  ,   573.057  ,   520.270  ,   475.841  ,   437.562  ,   404.337  ,   375.358  ,   349.968  ,   327.631  ,   318.730  ,   316.250  ,   313.803  ,   311.425  ,   309.177  ,   307.136  ,   305.355  ,   303.845  ,   302.579  ,   301.527  ,   300.652  ,   299.929  ,   299.331  ,   298.838  ,   298.433  ,   298.101  ,   297.832  ,   297.618  ,   297.581  ,   298.783  ,   313.284  /)

  th@long_name = "Theta"
  th@units = "K"

  z_full@long_name = "height"
  z_full@units = "m"

;------------------
; general settings
;------------------

  wks_name = "test111"
  wks_type = "png"  ; or x11,pdf,ncgm
  wks = gsn_open_wks(wks_type, wks_name)
  res = True

;==================
  ks = 28
  ke = kx-1
  yy = z_full(ks:ke)
  xx = th(:,ks:ke)
  labels_xy1 = new ( tx, "string" )

  ss1 = "dry SCM, Hour "
  do tt=0,tx-1
    labels_xy1(tt) = ss1+tt
  end do

  res_xy1 = res
    res_xy1@tiMainString = "Offline MYNN - DCBL_C1"
  opt_xy1 = True
    opt_xy1@gsnString = "x"
    opt_xy1@legend_mean = False

;--------
; plot
;--------
  plot_xy1 = ignorant_xy(wks, xx, yy, labels_xy1, res_xy1, opt_xy1)
  opt1 = True
  trim_fig(wks_type, wks_name, opt1)


end
