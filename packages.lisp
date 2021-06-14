
;; (in-package #:CL-USER)

(defpackage :vector-ops
  (:use #:common-lisp)
  (:nicknames #:vo #:vop #:vops)
  (:export
   #:vec
   #:blit
   #:vlog
   #:vexp
   #:vabs
   #:vsub
   #:vadd
   #:vmul
   #:vdiv
   #:vscale
   #:voffset
   #:vmin
   #:vmax
   #:vsum
   #:vprod
   #:vround
   #:vfloor
   #:vceiling
   #:vsqr
   #:vsqrt
   #:vlog10
   #:vpow10
   #:vexpt
   #:destructure-vector
   #:vlimit<
   #:vlimit>
   #:extrema
   #:vextrema
   #:clip
   #:vclip
   #:wrap
   #:vwrap
   ))

(defpackage :fft
  (:use #:common-lisp)
  (:export
   #:$fftw-forward
   #:$fftw-inverse

   #:d2zfft
   #:z2dfft
   #:z2zfft
   #:unsafe-z2zfft

   #:r2cfft
   #:c2rfft
   #:c2cfft
   #:unsafe-c2cfft
   
   #:d2zfft2d
   #:z2dfft2d
   #:z2zfft2d
   #:unsafe-z2zfft2d
   
   #:r2cfft2d
   #:c2rfft2d
   #:c2cfft2d
   #:unsafe-c2cfft2d

   #:make-fft-buffer
   #:fft-buffer-p
   #:fft-buffer-r
   #:fft-buffer-roff
   #:fft-buffer-i
   #:fft-buffer-ioff
   #:fft-buffer-hr
   #:fft-buffer-nx
   #:get-real
   #:get-imag
   #:set-real
   #:set-imag
   #:r2c
   #:c2r
   #:c2c
   #:fwd
   #:inv
   #:fwd-magnitude
   #:fwd-magnitude-db
   #:fwd-power
   #:fwd-phase
   #:fwd-phase-deg
   #:get-align16-offset
   #:siglab_sbFFT

   #:symmetric-fill
   #:symmetric-replace
   #:center-implant

   #:get-stwids
   #:get-dtwids
   #:get-process-split-tmp
   #:get-stmp
   #:get-dtmp
   #:get-c-address
   ))

#+:MACOSX
(defpackage :fft2d
  (:use #:common-lisp)
  (:export
   #:make-fft-buffer
   #:fft-buffer-p
   #:fft-buffer-r
   #:fft-buffer-roff
   #:fft-buffer-i
   #:fft-buffer-ioff
   #:fft-buffer-hr
   #:fft-buffer-ny
   #:fft-buffer-nx
   #:get-real
   #:get-imag
   #:set-real
   #:set-imag
   #:r2c
   #:c2r
   #:fwd
   #:inv
   #:fwd-magnitude
   #:fwd-magnitude-db
   #:fwd-power
   #:fwd-phase
   #:fwd-phase-deg
   #:siglab_sbFFT))

#+:MACOSX
(defpackage :sfft
  (:use #:common-lisp)
  (:export
   #:make-fft-buffer
   #:fft-buffer-p
   #:fft-buffer-r
   #:fft-buffer-roff
   #:fft-buffer-i
   #:fft-buffer-ioff
   #:fft-buffer-hr
   #:fft-buffer-nx
   #:get-real
   #:get-imag
   #:set-real
   #:set-imag
   #:r2c
   #:c2r
   #:c2c
   #:fwd
   #:inv
   #:fwd-magnitude
   #:fwd-magnitude-db
   #:fwd-power
   #:fwd-phase
   #:fwd-phase-deg))

#+:MACOSX
(defpackage :dfft
  (:use #:common-lisp)
  (:export
   #:make-fft-buffer
   #:fft-buffer-p
   #:fft-buffer-r
   #:fft-buffer-roff
   #:fft-buffer-i
   #:fft-buffer-ioff
   #:fft-buffer-hr
   #:fft-buffer-nx
   #:get-real
   #:get-imag
   #:set-real
   #:set-imag
   #:d2z
   #:z2d
   #:z2z
   #:fwd
   #:inv
   #:fwd-magnitude
   #:fwd-magnitude-db
   #:fwd-power
   #:fwd-phase
   #:fwd-phase-deg))

#+:MACOSX
(defpackage :aquaterm
  (:use #:common-lisp)
  (:nicknames #:aqt)
  (:export
   #:def-proxy-fli-function
   #:terminate
   #:open-plot
   #:select-plot
   #:set-plot-size
   #:set-plot-title
   #:save-context
   #:restore-context
   #:set-cliprect
   #:set-affine-transform
   #:render-plot
   #:clear-plot
   #:close-plot
   #:save-plot
   #:set-color
   #:set-background-color
   #:set-font-name
   #:set-font-size
   #:add-label
   #:set-line-width
   #:set-line-cap-style
   #:move-to
   #:add-line-to
   #:add-poly-line
   #:move-to-vertex
   #:add-edge-to-vertex
   #:add-polygon
   #:add-filled-rect
   #:erase-rect
   #:set-image-transform
   #:reset-image-transform
   #:add-image-with-bitmap
   #:add-transformed-image-with-bitmap

   #:set-accepting-events
   #:wait-next-event
   #:get-last-event

   #:butt-line-cap-style
   #:round-line-cap-style
   #:square-line-cap-style

   #:align-left
   #:align-center
   #:align-right

   #:align-middle
   #:align-baseline
   #:align-bottom
   #:align-top ))

(defpackage :scigraph
  (:use #:common-lisp)
  (:nicknames #:sg)
  (:export
   #:tvscl
   #:plot
   #:oplot
   #:draw-text
   #:plot-axes
   #:plot-image
   #:draw-grid
   #:window
   #:wset
   #:wshow
   #:werase
   #:wkill
   #:wattach
   #:wdetach
   #:wsave
   #:delay-update
   #:update
   #:with-delayed-update
   #:plot-polys
   #:set-cmap
   #:get-mouse
   #:get-coord-values
   #:get-wsize
   #:copy-graphic-to-clipboard
   #:$sym-circle
   #:$sym-square
   #:$sym-box
   #:$sym-dot
   #:$sym-cross
   #:$sym-triangle
   #:$sym-histogram
   #:$penpat-solid
   #:$penpat-dash
   #:$penpat-dot
   #:$penpat-dashdot
   #:$penpat-dashdotdot
   #:$penpat-null
   #:$penpat-inside-frame
   #:rgb
   #:$black
   #:$white
   #:$red
   #:$green
   #:$blue
   #:$darkred
   #:$darkgreen
   #:$darkblue
   #:$skyblue
   #:$purple
   #:$magenta
   #:$cyan
   #:$yellow
   #:$orange
   #:$gray50
   #:$gray25
   #:$gray75

   #:plotfft
   #:tvfft
   #:histo
   #:log-histo

   #:load-cursor-from-file
   #:direct-redraw
   
   #:prompt-for-filenames
   #:move-memory
   #:clear-memory
   #:shift-image-left
   #:set-image-cursor-transform

   #:set-heat-colormap
   #:set-gray-colormap

   #:save-plot

   #:fplot
   #:ofplot
   #:paramplot
   #:oparamplot

   #:stamp-logo

   #:get-next-event
   ))

(defpackage :vectorized-math
  (:use #:common-lisp #:vector-ops)
  (:nicknames #:vmath #:vm)
  (:export
   #:iramp
   #:framp
   #:dramp
   #:bipolar-framp
   #:unoise
   #:hypot
   #:sinc
   #:logabs
   #:vector-of
   #:gasdev
   #:gnoise
   #:median
   #:mad
   #:negmad
   #:total
   #:mean
   #:stdev
   #:variance
   #:wtmean
   #:wtvariance
   #:wtstdev
   #:percentile
   #:percentiles
   #:standard-percentiles
   #:copy
   #:histogram

   #:fzeros
   #:fones
   #:izeros
   #:iones
   
   #:inner-prod
   #:general-inner-prod
   #:outer-prod
   #:general-outer-prod
   #:dist
   #:shift
   #:shifth
   #:slice
   #:reshape
   #:map-dimension
   #:reduce-dimension
   #:vrotl
   #:vrotr
   #:vshifth
   #:split
   #:transpose
   #:acanon
   #:xplane
   #:yplane
   #:bipolar-xplane
   #:bipolar-yplane
   #:shifted-bipolar-framp
   #:shifted-bipolar-xplane
   #:shifted-bipolar-yplane

   #:vectorwise
   #:elementwise
   #:elementwise-reduce
   #:gensyms-for-args
   #:ixsort

   #:require-same-shape
   #:require-same-size
   #:defun-ffp

   #:dgesvd-solve
   #:dgesvd
   #:dgesvd-bksb
   #:dgesvd-predict

   #:simplex
   
   #:without-denormals
   ))
   
(defpackage :surface-plots
  (:nicknames #:surf #:surfplots)
  (:use #:common-lisp)
  (:export
   #:plot-surface
   #:lego-plot
   #:rot
   #:red-color-fn
   #:gray-color-fn
   #:lamps-color-fn
   #:*lamps*
   #:make-lamp
   #:lamp-dir-vector
   #:lamp-rgb-triple
   #:lamp-intensity
   ))

(defpackage :image-processing
  (:use #:common-lisp)
  (:nicknames #:img)
  (:export
   #:<image-array>
   #:<matrix-array>
   #:image-array-arena
   "(SETF IMAGE-ARRAY-ARENA)"
   #:make-image
   #:make-matrix
   #:make-similar
   #:pixelwise
   #:pixelwise-reduce
   #:as-image
   #:as-matrix
   #:<subimage-array>
   #:subimage-centered
   #:toroidal-coords
   #:subimage
   #:place-subimage
   #:fill-subimage
   #:shift
   #:shifth
   #:vshift
   #:col-vector
   #:row-vector
   #:get-col
   #:get-row
   #:get-column-vector
   #:get-row-vector
   #:x-slice
   #:y-slice
   #:flipv
   #:fliph
   #:transpose
   #:matrix-diagonal
   #:tvscl
   #:max-ix
   #:fft
   #:ifft
   #:xplane
   #:yplane
   #:where
   #:total
   #:mean
   #:plot-surface
   #:lego-plot
   #:indices-of
   #:find-peak
   #:show-peak))

(defpackage :photometry
  (:use #:common-lisp)
  (:nicknames #:phot)
  (:export
   #:photom))

(defpackage #:interpolation
  (:use #:common-lisp)
  (:nicknames #:interp)
  (:export
   #:linint
   #:polint
   #:ratint
   #:locate
   #:locate-subtable
   #:spline
   #:splint
   #:monotonic-spline
   #:monotonic-splint
   ))

(defpackage #:integrate
  (:use #:common-lisp)
  (:nicknames #:integ)
  (:export
   #:trapm
   #:mtriple
   #:qromb
   #:qrombm
   #:sub-array
   ))

(defpackage #:roots
  (:use #:common-lisp)
  (:export
   #:find-root
   #:newton
   #:deriv
   #:zbrac
   #:rtbis
   #:rtsec
   #:zbrent
   #:rtsafe
   ))

(defpackage :kaiser
  (:use :cl)
  (:export
   :kaiser-window
   :sinc
   :beta-factor
   :compute-fir-filter-order
   :bessel-i0
   ))

#|
;; generate HTML documentation
(in-package :vmath)
(doctools:gen-docs
 :asdf-system-name :vmath
 :package-name     :vmath
 :directory        (translate-logical-pathname "PROJECTS:LISP;vmath;")
 :subtitle         "a library for vectorized math")
|#
