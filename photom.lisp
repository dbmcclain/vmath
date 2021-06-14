;; photom.lisp -- Interactive stellar photometry
;; DM/MCFA  12/99
;;

(in-package "PHOTOMETRY")

;; set the colormap to grays
(defun init-colormap ()
  (let ((r (vm:iramp 256)))
    (plt:set-cmap r r r)))

(topgui:define-toplevel-app-interface interface-1 ()
  ((image  :accessor interface-image :initform nil))
  (:panes

   (image-display
    capi:output-pane
    :accessor image-display
    :cursor :crosshair
    :input-model
    '(
      ("Control-c"    copy-to-clipboard)
      ("Left"          move-left 1)
      ("Right"         move-right 1)
      ("Up"            move-up 1)
      ("Down"          move-down 1)
      ("Control-left"  move-left 10)
      ("Control-right" move-right 10)
      ("Control-up"    move-up 10)
      ("Control-down"  move-down 10)
      (#\z             set-ref-magnitude-from-kbd)
      (:motion         show-position)
      ((:button-3 :press)     set-ref-magnitude)
      ;((:button-1 :motion)   redraw-crosshairs)
      ;((:button-1 :release)  undraw-crosshairs)
      #|
      (#\control-\c    copy-to-clipboard)
      (#\left          move-left 1)
      (#\right         move-right 1)
      (#\up            move-up 1)
      (#\down          move-down 1)
      (#\control-left  move-left 10)
      (#\control-right move-right 10)
      (#\control-up    move-up 10)
      (#\control-down  move-down 10)
      (#\z             set-ref-magnitude-from-kbd)
      (:motion         show-position)
      ((:button-3 :press)     set-ref-magnitude)
      ;((:button-1 :motion)   redraw-crosshairs)
      ;((:button-1 :release)  undraw-crosshairs)
      |#
      )
    :display-callback 'redisplay-image
    :min-width 640
    :min-height 480)

   (color-bar-display
    capi:output-pane
    :accessor color-bar-display
    :display-callback 'redraw-color-bar
    :min-width 9
    :max-width 9
    :min-height 256
    :max-height 256)
   
   (magn-image-display
    capi:output-pane
    :accessor magn-image-display
    ;; :cursor :arrow
    :min-width (* 11 8)
    :max-width (* 11 8)
    :min-height (* 11 8)
    :max-height (* 11 8)
    :display-callback 'draw-magn-frame-fiducials)

   (x-display-title
    capi:title-pane
    :accessor x-display-title
    :text "X Position")
   (x-display
    capi:title-pane
    :accessor x-display
    :foreground :white
    :background :gray50
    :min-width 50
    :max-width 50)

   (y-display-title
    capi:title-pane
    :accessor y-display-title
    :text "Y Position")
   (y-display
    capi:title-pane
    :accessor y-display
    :foreground :white
    :background :gray50
    :min-width 50
    :max-width 50)

   (z-display-title
    capi:title-pane
    :accessor z-display-title
    :text "Z Value")
   (z-display
    capi:title-pane
    :accessor z-display
    :foreground :white
    :background :gray50
    :min-width 50
    :max-width 50)

   (bg-display-title
    capi:title-pane
    :accessor bg-display-title
    :text "BG Median")
   (bg-display
    capi:title-pane
    :accessor bg-display
    :foreground :white
    :background :gray50
    :min-width 50
    :max-width 50)

   (m-display-title
    capi:title-pane
    :accessor m-display-title
    :text "Magnitude")
   (m-display
    capi:title-pane
    :accessor m-display
    :foreground :white
    :background :gray50
    :min-width 50
    :max-width 50)

   (offset-slider
    capi:slider
    :accessor offset-slider
    :start 0
    :end 500
    :slug-start 0
    :min-width 200
    :max-width 200
    :callback 'slider-refresh-image
    :title "Offset")

   (range-slider
    capi:slider
    :accessor range-slider
    :start 0
    :end 500
    :slug-start 500
    :min-width 200
    :max-width 200
    :callback 'slider-refresh-image
    :title "Range")

   (neg-button
    capi:check-button
    :accessor neg-button
    :text "Negative"
    :selection-callback 'refresh-image
    :retract-callback 'refresh-image
    :min-width 80
    :max-width 80
    :enabled t)

   (logscale-button
    capi:check-button
    :accessor logscale-button
    :selection-callback 'refresh-image
    :retract-callback 'refresh-image
    :text "Log Scale"
    :min-width 80
    :max-width 80
    :enabled t)
   )
  (:layouts
   (row-layout-1
    capi:row-layout
    '(image-display
      color-bar-display
      column-layout-1))
   (column-layout-1
    capi:column-layout
    '(row-layout-3   ;; magn image + readouts
      row-layout-4   ;; check boxes
      offset-slider
      range-slider
      ))
   (row-layout-3
    capi:row-layout
    '(magn-image-display
      column-layout-2
      column-layout-3))
   (column-layout-2
    capi:column-layout
    '(x-display-title
      y-display-title
      z-display-title
      bg-display-title
      m-display-title))
   (column-layout-3
    capi:column-layout
    '(x-display
      y-display
      z-display
      bg-display
      m-display))
   (row-layout-4
    capi:row-layout
    '(neg-button
      logscale-button)))

  (:menu-bar menu-1)
  (:menus
   (menu-1
    "File"
    (("Open Image..." :selection-callback 'get-image)
     ("Exit"          :selection-callback 'kill-interface))))

  (:default-initargs
   :layout 'row-layout-1
   :title "Stellar Photometry"
   ))

(defun kill-interface (data intf)
  (declare (ignore data))
  (capi:destroy intf))

#|
(defvar *bullseye-cursor*
  (lazy:make (plt:load-cursor-from-file "bullseye.cur")))

(defun photom ()
  (labels
      ((startup (intf)
                (init-colormap)
                (let* ((pane (image-display intf))
                       (repr (slot-value pane 'capi-internals:representation)))
                  (setf (slot-value repr 'capi-win32-lib::cursor)
                        #|(lazy:force *bullseye-cursor*)|#
                        :cross))
                ))
    (topgui:run-toplevel-app-interface 'interface-1
                                       :after #'startup)))
|#

(defun photom ()
  (labels
      ((startup (intf)
                (init-colormap)
                ))
    (topgui:run-toplevel-app-interface 'interface-1
                                       :after #'startup)))
(defun copy-to-clipboard (&rest args)
  (declare (ignore args))
  (plt:copy-graphic-to-clipboard)
  (capi:display-message "image captured"))

(defclass <image-array-x> ()
  ((data :accessor image-array-arena :initarg :arena)
   (tint :accessor integration-time   :initarg  :tint)))

(defun make-image-x (arr tint)
  (if (= (ca:carray-rank arr) 2)
      (if (ca:is-float-array arr)
          (make-instance '<image-array-x>
		         :arena arr
                         :tint  (or tint 1.0))
        (let ((xarr (ca:convert-to-float arr :type :float)))
          (ca:discard-carray arr)
          (make-image-x xarr tint)))
    (progn
      (ca:discard-carray arr)
      (error "flat 2-D array required"))))

#|
(defmethod img:make-similar ((a <image-array-x>) arr)
  (make-instance '<image-array-x>
		 :arena arr
                 :tint  (integration-time a)))
|#

(defun set-window-title (intf str)
  (setf (capi:interface-title intf) str))

(defun set-interface-image (intf img)
  (setf (interface-image intf) img))

(defun select-kth (vec ixs k left right)
  (declare 
   #+:lispworks 
   (optimize (debug 1) (speed 3) (safety 0)) ;; what is this trying to achieve?
   (type fixnum k left right)
   (type (array fixnum 1) ixs))
  (ca:with-row-major-access (p vec)
    (labels ((ind-aref (ix)
               ;; Corman Lisp does not yet handle DECLARE inside of (LABELS ...) (v/1.5 12/01)
               #-:cormanlisp 
               (declare (type fixnum ix))
               (fli:dereference p :index (aref ixs ix))))
      (let ((l   left)
            (r   right))
        (declare (type fixnum l r))
        (tagbody
         again-1
         (let ((v (ind-aref k))
               (i l)
               (j r))
           (declare (type fixnum i j))
           (tagbody
            again-2
            (let ((ip (do ((ix i (1+ ix)))
                          ((<= v (ind-aref ix)) ix)
                        (declare (type fixnum ix))
                      ))
                  (jp (do ((ix j (1- ix)))
                          ((>= v (ind-aref ix)) ix)
                        (declare (type fixnum ix))
                        )))
              (declare (type fixnum ip jp))
              (if (<= ip jp)
                  (let ((ipp (1+ ip))
                        (jpp (1- jp)))
                    (declare (type fixnum ipp jpp))
                    (rotatef (aref ixs ip) (aref ixs jp))
                    (setf i ipp
                          j jpp)
                    (if (<= ipp jpp)
                        (go again-2)))
                (setf i ip
                      j jp))
              ))
           (if (< j k) (setf l i))
           (if (< k i) (setf r j))
           (if (< l r) (go again-1))
           )))
      (values (ind-aref k) ixs))))
  
(defun percentiles (arr)
  (let* ((len   (ca:carray-total-size arr))
         (limit (1- len))
         (ixs   (vm:iramp len))
         (n50   (round (* 50 limit) 100))
         (n01   (round limit 100))
         (n05   (round (*  5 limit) 100))
         (n10   (round (* 10 limit) 100))
         (n25   (round (* 25 limit) 100))
         (n75   (round (* 75 limit) 100))
         (n90   (round (* 90 limit) 100))
         (n95   (round (* 95 limit) 100))
         (n99   (round (* 99 limit) 100))
         (pc50  (select-kth arr ixs n50 0 limit))
         (pc25  (select-kth arr ixs n25 0 n50))
         (pc10  (select-kth arr ixs n10 0 n25))
         (pc05  (select-kth arr ixs n05 0 n10))
         (pc01  (select-kth arr ixs n01 0 n05))
         (pc75  (select-kth arr ixs n75 n50 limit))
         (pc90  (select-kth arr ixs n90 n75 limit))
         (pc95  (select-kth arr ixs n95 n90 limit))
         (pc99  (select-kth arr ixs n99 n95 limit)))
    (list
     :pc01 pc01
     :pc05 pc05
     :pc10 pc10
     :pc25 pc25
     :pc50 pc50
     :pc75 pc75
     :pc90 pc90
     :pc95 pc95
     :pc99 pc99)))

(defun set-slider-bounds (intf img)
  (let* ((pcs   (percentiles (image-array-arena img)))
         (pc01  (truncate (getf pcs :pc01)))
         (pc99  (max (truncate (getf pcs :pc99))
                     (+ pc01 500)))
         (range (* 2 (- pc99 pc01)))
         (offset-slider (offset-slider intf))
         (range-slider  (range-slider intf)))
    (setf (capi:range-end range-slider)  range
          (capi:range-end offset-slider) pc99
          (capi:range-slug-start offset-slider) pc01
          (capi:range-slug-start range-slider)  range)))

(defmethod float-value (v &optional (default 1.0))
  (declare (ignore v))
  default)

(defmethod float-value ((v number) &optional default)
  (declare (ignore default))
  (float v))

(defmethod float-value ((v string) &optional default)
  (declare (ignore default))
  (float (read-from-string v)))

#+:LISPWORKS4.2
(defun get-image (data intf)
  (handler-case
      (multiple-value-bind (img attrs fname)
          (scids:getvar :attrs '("exptime") :fast t)
        (let* ((arr  (sa:make-variant-carray img))
               (tint (float-value (first attrs)))
               (ximg (make-image-x arr tint)))
          (set-interface-image intf ximg) ;; (img:flipv ximg))
          (set-window-title intf (namestring fname))
          (set-slider-bounds intf ximg)
          (refresh-image data intf)))
    (scids:USER-CANCEL ())))

#-:LISPWORKS4.2
(defun get-image (data intf)
  (handler-case
      (multiple-value-bind (img attrs fname)
          (scids:getvar :attrs '("exptime"))
        (let* ((arr  (ca:make-carray :float (array-dimensions img)))
               (tint (float-value (first attrs)))
               (ximg (make-image-x arr tint)))
          (ca:copy-lisp-array-to-float-carray img arr)
          (set-interface-image intf ximg) ;; (img:flipv ximg))
          (set-window-title intf (namestring fname))
          (set-slider-bounds intf ximg)
          (refresh-image data intf)))
    (scids:USER-CANCEL ())))

(defvar *color-bar*
  ;; allocate only once an use over and over again...
  (let ((img (ca:make-carray :float '(256 9))))
    (dotimes (iy 256)
      (dotimes (ix 9)
        (setf (ca:caref img iy ix) (float iy))))
    img))

(defun redraw-color-bar (pane &rest args)
  (declare (ignore args))
  (plt:direct-redraw pane))

(defun draw-color-bar (intf)
  (let ((wpane (color-bar-display intf)))
    (plt:wset wpane)
    (let ((neg-img (capi:button-selected (neg-button intf))))
      (plt:tvscl *color-bar* :neg neg-img))
    ))

(defvar *magnification*  2)

(defun draw-image-display (intf img)
  (let ((wpane           (image-display intf))
        (offs-pane       (offset-slider intf))
        (range-pane      (range-slider intf))
        (negate-button   (neg-button intf))
        (logscale-button (logscale-button intf)))
    (plt:window wpane)
    (plt:clear  wpane)
    (let* ((display-min (capi:range-slug-start offs-pane))
           (display-max (+ display-min
                           (capi:range-slug-start range-pane)))
           (neg-img     (capi:button-selected negate-button))
           (use-log     (capi:button-selected logscale-button)))
      (if use-log
          (show-log-stretch img
                            display-min
                            display-max
                            neg-img)
        (show-linear-stretch img
                             display-min
                             display-max
                             neg-img)))
    ))

(defmacro with-fast-unsafe-code (&body body)
  `(locally
     (declare (optimize (speed  3)
                        (safety 0)
                        (debug  0)))
     ,@body))

(defun show-log-stretch (img minv maxv neg)
  (let* ((arr  (image-array-arena img))
         (minv (float minv))
         (maxv (float maxv))
         (fn   #'(lambda (v)
                   (declare (type float v))
                   (log (max 1.0 (- v minv))))
               ))
    (declare (type float minv maxv))
    (ca:with-dynamic-carray (ximg :float
                                  (ca:carray-dimensions arr))
      (ca:with-row-major-access (pdst ximg)
        (ca:with-row-major-access (psrc arr)
          (dotimes (ix (ca:carray-total-size ximg))
            (setf (fli:dereference pdst :index ix)
                  (funcall fn (fli:dereference psrc :index ix))))
          ))
      (plt:tvscl ximg
                :range (list 0 (funcall fn maxv))
                :magn  *magnification*
                :neg   neg))
    ))

(defun show-linear-stretch (img minv maxv neg)
  (plt:tvscl (image-array-arena img)
            :range (list minv maxv)
            :magn  *magnification*
            :neg   neg))

(defun slider-refresh-image (slider val action)
  (when (eq action :move)
    (refresh-image val (slot-value slider 'capi:interface))))

(defun refresh-image (data intf)
  (declare (ignore data))
  (let ((img (interface-image intf)))
    (when img
      (draw-color-bar intf)
      (draw-image-display intf img))))

(defun redisplay-image (pane &rest args)
  (declare (ignore args))
  (plt:direct-redraw pane))

(defvar *sub-image*
  ;; allocate once and use over and over again...
  (ca:make-carray :float '(11 11)))

(um:defun* get-subimage-centered (arr (ctrx ctry))
  (destructuring-bind (dimy dimx)
      (ca:carray-dimensions arr)
    (let* ((minval (ca:caref arr ctry ctrx))
           (maxval minval))
      (declare (type float minval maxval))
      (if (and (< 4 ctry (- dimy 5))
               (< 4 ctrx (- dimx 5)))
          (loop for iy from (- ctry 5) to (+ ctry 5)
                and jy from 0 below 11 do
                (loop for ix from (- ctrx 5) to (+ ctrx 5)
                      and jx from 0 below 11 do
                      (let ((v (ca:caref arr iy ix)))
                        (declare (type float v))
                        (setf minval (min minval v)
                              maxval (max maxval v))
                        (setf (ca:caref *sub-image* jy jx) v))))
        (let ((rpos -1)
              (rvec #.(make-array 121)))
          (loop for iy from (- ctry 5) to (+ ctry 5)
                and jy from 0 below 11 do
                (loop for ix from (- ctrx 5) to (+ ctrx 5)
                      and jx from 0 below 11 do
                      (if (and (< -1 iy dimy)
                               (< -1 ix dimx))
                          (let ((v (ca:caref arr iy ix)))
                            (declare (type float v))
                            (setf minval (min minval v)
                                  maxval (max maxval v))
                            (setf (ca:caref *sub-image* jy jx) v))
                        (let ((pos (+ jx (* 11 jy))))
                          (setf (aref rvec pos) rpos)
                          (setf rpos pos)))
                      ))
          (ca:with-row-major-access (p *sub-image*)
            (do ((rpos rpos))
                ((minusp rpos))
              (let ((rpos2 (aref rvec rpos)))
                (setf (fli:dereference p :index rpos) minval)
                (setf rpos rpos2))
              ))
          ))
      (values *sub-image* minval maxval)
      )))
          
(defun show-magnified-selection (intf img x y #|med|# )
  (let* ((pane  (magn-image-display intf))
         (arr   (image-array-arena img))
	 (negate-button   (neg-button intf))
	 (neg-img     (capi:button-selected negate-button)))
    (multiple-value-bind (ximg minval maxval)
        (get-subimage-centered arr (list y x))
      (plt:wset pane)
      (plt:tvscl ximg
                :magn  8
                :range (list minval
                             (max (ca:caref ximg 5 5) ;; pick out ctr pixel
                                  (+ 100 (/ (+ minval maxval) 2)))) ;; + med 100)))
                :neg   neg-img)
      )))

(defun draw-magn-frame-fiducials (pane x y width height)
  (declare (ignore x y width height))
  (plt:direct-redraw pane)
  (gp:with-graphics-state (pane :operation boole-1
                                :foreground :red)
    (gp:draw-line pane 0 43 38 43)
    (gp:draw-line pane 48 43 87 43)
    (gp:draw-line pane 43 0 43 38)
    (gp:draw-line pane 43 48 43 87)
    (gp:draw-rectangle pane 23 23 40 40)
    (gp:draw-circle pane 43 43 10)))

(defparameter *wmask*
  (let ((mask nil))
    ;; a list of (iy ix) index pairs defining the "moat"
    (loop for iy from -20 to 20 do
          (let ((ysq (* iy iy)))
            (loop for ix from -20 to 20 do
                  (if (> (+ (* ix ix) ysq) #.(* 24 24)) ;; gives 112 pts
                      (push (list iy ix) mask))
                  )))
    mask))

(defparameter *moat*
  (make-array (length *wmask*)
              :element-type 'float))

(defun compute-moat-median (img ctry ctrx)
  (let ((arr (image-array-arena img)))
    (destructuring-bind (dimy dimx) (ca:carray-dimensions arr)
      (let* ((moat (let ((cnt 0))
                     (if (and (< 19 ctry (- dimy 20))
                              (< 19 ctrx (- dimx 20)))
                         (loop for (iy ix) in *wmask* do
                               (let ((ypos (+ iy ctry))
                                     (xpos (+ ix ctrx)))
                                 (setf (aref *moat* cnt)
                                       (ca:caref arr ypos xpos))
                                 (incf cnt)))
                       (loop for (iy ix) in *wmask* do
                             (let ((ypos (+ iy ctry))
                                   (xpos (+ ix ctrx)))
                               (when (and (< -1 ypos dimy)
                                          (< -1 xpos dimx))
                                 (setf (aref *moat* cnt)
                                       (ca:caref arr ypos xpos))
                                 (incf cnt))
                               )))
                     (make-array cnt
                                 :displaced-to *moat*
                                 :element-type 'float)))
             (med (vm:median moat))
             (tot (let ((tot 0))
                    (if (and (< 4 ctry (- dimy 5))
                             (< 4 ctrx (- dimx 5)))
                        (loop for iy from (- ctry 5) to (+ ctry 5) do
                              (loop for ix from (- ctrx 5) to (+ ctrx 5) do
                                    (incf tot (max 0 (- (ca:caref arr iy ix) med)))
                                    ))
                      (loop for iy from (- ctry 5) to (+ ctry 5) do
                            (loop for ix from (- ctrx 5) to (+ ctrx 5) do
                                  (when (and (< -1 iy dimy)
                                             (< -1 ix dimx))
                                    (incf tot (max 0 (- (ca:caref arr iy ix) med))))
                                  )))
                    tot)))
        (values med tot))
      )))

(defvar *ref-magn*  (+ 14.31 -0.84 0.89 -1.83 0.12))

(defun set-bgmed-display (intf med)
  (setf (capi:title-pane-text (bg-display intf))
        (format nil "~D" (round med))))

(defun compute-magnitude (x y img intf)
  (show-magnified-selection intf img x y)
  (multiple-value-bind (med tot)
      (compute-moat-median img y x)
    (set-bgmed-display intf med)
    (if (plusp tot)
        (+ (* -2.5 (- (log tot 10)
                      (log (integration-time img) 10)))
           *ref-magn*)
      99)))

(defparameter *lastx* 0)
(defparameter *lasty* 0)

(defun set-ref-magnitude-from-kbd (pane &rest args)
  (declare (ignore args))
  (set-ref-magnitude pane *lastx* *lasty*))

#|
(defun cvt-display-to-image-coords (img x y)
  (let ((ylimit (ca:carray-dimension (image-array-arena img) 0)))
    (values (truncate x *magnification*)
            (truncate (- (* ylimit *magnification*) 1 y) *magnification*))))
|#

(defun cvt-display-to-image-coords (img x y)
  (declare (ignore img))
  (values (truncate x *magnification*)
          (truncate y *magnification*)))

(defun set-ref-magnitude (pane x y)
  (setf *lastx* x
        *lasty* y)
  (let* ((intf (capi:element-interface pane))
         (img  (interface-image intf)))
    (when img
      (multiple-value-bind (x y) (cvt-display-to-image-coords img x y)
        (let* ((magn      (compute-magnitude x y img intf))
               (user-magn (capi:prompt-for-string
                           "Enter reference magnitude"
                           :initial-value (format nil "~,2F" magn))))
          (when user-magn
            (setf *ref-magn* (+ *ref-magn*
                                (- (float-value user-magn) magn))))
          )))))

(defun set-x-display (intf x)
  (setf (capi:title-pane-text (x-display intf))
        (format nil "~D" x)))

(defun set-y-display (intf y)
  (setf (capi:title-pane-text (y-display intf))
        (format nil "~D" y)))

(defun set-z-display (intf z)
  (setf (capi:title-pane-text (z-display intf))
        (format nil "~D" (round z))))

(defun set-m-display (intf m)
  (setf (capi:title-pane-text (m-display intf))
        (format nil "~,2F" m)))

(defun clear-z-display (intf)
  (setf (capi:title-pane-text (z-display intf)) ""))

(defun clear-m-display (intf)
  (setf (capi:title-pane-text (m-display intf)) ""))

(defun show-position (pane x y)
  (setf *lastx* x
        *lasty* y)
  (let* ((intf (capi:element-interface pane))
         (img  (interface-image intf)))
    (when img
      (let* ((a     (image-array-arena img))
             (xlim  (1- (ca:carray-dimension a 1)))
             (ylim  (1- (ca:carray-dimension a 0))))
        (multiple-value-bind (x y) (cvt-display-to-image-coords img x y)
          (set-x-display intf x)
          (set-y-display intf y)
          (if (and (>= x 0)
                   (<= x xlim)
                   (>= y 0)
                   (<= y ylim))
              (progn
                (set-z-display intf (ca:caref a y x))
                (set-m-display intf (compute-magnitude x y img intf)))
            (progn
              (clear-z-display intf)
              (clear-m-display intf)))
          )))))

(defun adjust-cursor (dx dy)
  (fli:with-dynamic-foreign-objects ()
    (let ((pos (fli:allocate-dynamic-foreign-object
                :type :int
                :nelems 2)))
      (win32:get-cursor-pos pos)
      (win32:set-cursor-pos (+ dx (fli:dereference pos :index 0))
                            (+ dy (fli:dereference pos :index 1)))
      )))

(labels
    ((move-cursor (pane x y dx dy)
                  (let ((dx (* dx *magnification*))
                        (dy (* dy *magnification*)))
                    (show-position pane (+ x dx) (+ y dy))
                    (adjust-cursor dx dy))))

  (defun move-up (pane x y key data)
    (declare (ignore key))
    (move-cursor pane x y 0 (- data)))
  
  (defun move-down (pane x y key data)
    (declare (ignore key))
    (move-cursor pane x y 0 data))
  
  (defun move-left (pane x y key data)
    (declare (ignore key))
    (move-cursor pane x y (- data) 0))
  
  (defun move-right (pane x y key data)
    (declare (ignore key))
    (move-cursor pane x y data 0)))

  
;; ----------------------------------------------------------

;; -- end of photom.lisp -- ;;
