;; fft.lisp -- interface to FFTX 2-D FFT Routines
;;
;; DM/MCFA  08/99
;; --------------------------------------------------

(in-package #:FFT)

;; ------------------------------------------------------------------
;; Show us the FFT Library Version during loading...
;;
(multiple-value-bind (ans str)
    (getFFTVersionString "")
  (declare (ignore ans))
  (print str))

;; ------------------------------------------------------------------

(defstruct (fft-buffer
            (:constructor make-fft-buf))
  nx
  r i
  roff ioff
  hr)

(defun half-dim (n)
  (1+ (truncate n 2)))

(defun nfloats (nb type-size)
  (assert (zerop (logand nb (1- type-size)))) ;; assure that our offsets are multiples of float size
  (truncate nb type-size))

(defun make-fft-buffer (nx type)
  (let* ((type-size (ecase type
                      (single-float 4)
                      (double-float 8)))
         (nxa   (max 8 (um:ceiling-pwr2 nx)))
         (rarr  (make-array (+ nxa 3) :element-type type :allocation :static))
         (roff  (nfloats (get-align16-offset rarr) type-size))
         (iarr  (make-array (+ nxa 7) :element-type type :allocation :static))
         (ioff  (nfloats (get-align16-offset iarr) type-size))
         (hrarr (make-array (half-dim nxa)
                            :element-type type
                            :displaced-to rarr
                            :displaced-index-offset roff)))

    (if (zerop (logand (- (+ (* type-size ioff) (sys:object-address iarr))
                          (+ (* type-size roff) (sys:object-address rarr)))
                       (1- 1024)))
        ;; offset ioff by another 4 to avoid the Pentium quirk when two buffer addresses differ
        ;; by multiple of 1024 bytes.
        (incf ioff 4)) ;; bump by another 16 (or 32) bytes
                          
    (make-fft-buf
     :nx   nxa
     :r    rarr
     :roff roff
     :i    iarr
     :ioff ioff
     :hr   hrarr)
    ))

(defun get-real (fftbuf)
  (values (fft-buffer-r fftbuf) (fft-buffer-roff fftbuf)))

(defun get-imag (fftbuf)
  (values (fft-buffer-i fftbuf) (fft-buffer-ioff fftbuf)))

(defmethod set-real (fftbuf (arr vector))
  (replace (fft-buffer-r fftbuf) arr
           :start1 (fft-buffer-roff fftbuf)))

(defmethod set-real (fftbuf (val real))
  (multiple-value-bind (buf off) (get-real fftbuf)
    (fill buf (coerce val (array-element-type buf)) :start off)))

(defmethod set-imag (fftbuf arr)
  (replace (fft-buffer-i fftbuf) arr
           :start1 (fft-buffer-ioff fftbuf)))

(defmethod set-imag (fftbuf (val real))
  (multiple-value-bind (buf off) (get-imag fftbuf)
    (fill buf (coerce val (array-element-type buf)) :start off)))

(defun copy-fft-buffer-contents (src dst)
  (replace (fft-buffer-r dst) (fft-buffer-r src)
           :start1 (fft-buffer-roff dst)
           :start2 (fft-buffer-roff src))
  (replace (fft-buffer-i dst) (fft-buffer-i src)
           :start1 (fft-buffer-ioff dst)
           :start2 (fft-buffer-ioff src)))

;; ---------------------------------------------------------------

(defmethod check-dimension ((arr vector))
  ;; return array dimensions (ny, nx) of arr rounded up to next
  ;; power of 2.
  (max 8 (um:ceiling-pwr2 (length arr))))

(defmethod check-dimension ((arr fft-buffer))
  (fft-buffer-nx arr))
  
(defmethod effective-array-dimension ((arr vector))
  ;; return a pair (ny,nx) of array dimensions for arr.
  (length arr))

(defmethod effective-array-dimension ((arr fft-buffer))
  (fft-buffer-nx arr))

(defun effective-ctype (precision)
  (ecase precision
    ((:float :single :single-float :altivec) 
     (values :float 'single-float))
    ((:double :double-float :fftw)
     (values :double 'double-float))))

(defun vec (arr &key (size (array-total-size arr)) (offset 0))
  (make-array size
              :displaced-to arr
              :displaced-index-offset offset
              :element-type (array-element-type arr)))

(defun copy-array-to-real-cvect (arr cdst nx &key precision)
  (let ((nxa (effective-array-dimension arr))
        (pad (coerce 0 precision)))
    (dotimes (ix nxa)
      (setf (fli:dereference cdst :index ix)
            (coerce (aref arr ix) precision)))
    ;; pad remainder of row with zeros to get to power-of-two size
    (dotimes (ix (- nx nxa))
      (setf (fli:dereference cdst :index (+ nxa ix)) pad))
    ))

(defun copy-array-to-complex-cvect (arr cdst nx &key precision)
  (let ((nxa (effective-array-dimension arr))
        (pad (coerce 0 precision)))
    (dotimes (ix nxa)
      (let ((v (aref arr ix)))
        (setf (fli:dereference cdst :index (+ ix ix))
              (coerce (realpart v) precision)
              (fli:dereference cdst :index (+ ix ix 1))
              (coerce (imagpart v) precision))))
    ;; pad remainder of row with zeros to get to power-of-two size
    (let ((coff (+ nxa nxa)))
      (dotimes (ix (* 2 (- nx nxa)))
        (setf (fli:dereference cdst :index (+ coff ix)) pad))
      )))

(defun make-result-array (nx precision)
  (make-array nx :element-type precision))

(defun make-complex-result-array (nx precision)
  (declare (ignore precision))
  (make-array nx :element-type 'complex))

(defun convert-real-cvect-to-array (csrc nx
                                         &key
                                         dest
                                         precision)
  (let ((rslt (or dest (make-result-array nx precision))))
    (fli:replace-foreign-array rslt csrc)
    ))

(defun convert-complex-cvect-to-array (csrc nx
                                            &key
                                            dest
                                            precision)
  (let ((rslt  (or dest (make-complex-result-array nx precision))))
    (do ((cix 0 (+ cix 2))
         (rix 0 (1+ rix)))
        ((>= rix nx) rslt)
      (setf (aref rslt rix)
            (complex (fli:dereference csrc :index cix)
                     (fli:dereference csrc :index (1+ cix))))
      )))

(defun db10 (x)
  (if (plusp x)
      (* 10.0 (log x 10.0))
    -140.0))

(defun dtor (x)
  (* #. (/ (atan 1 1) 45) x))

(defun rtod (x)
  (* #. (/ 45 (atan 1 1)) x))

(defun phase-deg (x)
  (rtod (phase x)))

(defun convert-complex-cvect-magnitudes-to-array (csrc nx
                                                       &key
                                                       dest
                                                       precision)
  (let ((rslt (or dest (make-result-array nx precision))))
    (do ((cix 0 (+ cix 2))
         (rix 0 (1+ rix)))
        ((>= rix nx) rslt)
      (setf (aref rslt rix)
            (abs (complex (fli:dereference csrc :index cix)
                          (fli:dereference csrc :index (1+ cix)))))
      )))

(defun convert-complex-cvect-power-to-array (csrc nx
                                                  &key
                                                  dest
                                                  precision)
  (let ((rslt (or dest (make-result-array nx precision))))
    (do ((cix 0 (+ cix 2))
         (rix 0 (1+ rix)))
        ((>= rix nx) rslt)
      (setf (aref rslt rix)
            (let ((c (complex (fli:dereference csrc :index cix)
                              (fli:dereference csrc :index (1+ cix)))))
              (realpart (* c (conjugate c)))))
      )))

(defun convert-complex-cvect-magnitudes-db-to-array (csrc nx
							  &key
							  dest
							  precision)
  (let ((rslt (or dest (make-result-array nx precision))))
    (do ((cix 0 (+ cix 2))
         (rix 0 (1+ rix)))
        ((>= rix nx) rslt)
      (setf (aref rslt rix)
            (let ((c (complex (fli:dereference csrc :index cix)
                              (fli:dereference csrc :index (1+ cix)))))
              (db10 (realpart (* c (conjugate c))))))
      )))

(defun convert-complex-cvect-phases-to-array (csrc nx
                                                   &key
                                                   dest
                                                   precision)
  (let ((rslt (or dest (make-result-array nx precision))))
    (do ((cix 0 (+ cix 2))
         (rix 0 (1+ rix)))
        ((>= rix nx) rslt)
      (setf (aref rslt rix)
            (phase (complex (fli:dereference csrc :index cix)
                            (fli:dereference csrc :index (1+ cix)))))
      )))

(defun convert-complex-cvect-phases-deg-to-array (csrc nx
                                                   &key
                                                   dest
                                                   precision)
  (let ((rslt (or dest (make-result-array nx precision))))
    (do ((cix 0 (+ cix 2))
         (rix 0 (1+ rix)))
        ((>= rix nx) rslt)
      (setf (aref rslt rix)
            (phase-deg
             (complex (fli:dereference csrc :index cix)
                      (fli:dereference csrc :index (1+ cix)))))
      )))

#|
(defun tst (arr)
  (let ((nx (check-dimension arr)))
    (fli:with-dynamic-foreign-objects ()
      (let ((cvec (fli:allocate-dynamic-foreign-object
                   :type :double
                   :nelems nx)))
        (copy-array-to-real-cvect arr cvec nx)
        (convert-real-cvect-to-array cvec nx))
      )))

(defun tstz (arr)
  (let ((nx (check-dimension arr)))
    (fli:with-dynamic-foreign-objects ()
      (let ((cvec (fli:allocate-dynamic-foreign-object
                   :type :double
                   :nelems (* 2 nx))))
        (copy-array-to-complex-cvect arr cvec nx)
        (convert-complex-cvect-to-array cvec nx))
      )))
|#

(defun r2c-body (arr finish-fn &key
                     (precision :single-float)
                     dest)
  (um:bind*
      ((nx (check-dimension arr))
       (:values (ctype ltype) (effective-ctype precision)))
    (fli:with-dynamic-foreign-objects ()
      (let ((cdst (fli:allocate-dynamic-foreign-object
                   :type   ctype
                   :nelems (* 2 nx))))
        (fli:with-dynamic-foreign-objects ()
          (let ((csrc (fli:allocate-dynamic-foreign-object
                       :type   ctype
                       :nelems nx)))
            (copy-array-to-real-cvect arr csrc nx
                                      :precision ltype)
            (if (eq ctype :double)
                (d2zfft nx csrc cdst)
              (r2cfft nx csrc cdst))))
        (funcall finish-fn cdst nx
                 :dest dest
                 :precision ltype))
      )))

(defun r2c (arr &key
                (precision :single-float)
                dest)
  (r2c-body arr #'convert-complex-cvect-to-array
            :precision precision
            :dest dest))

(defun slow-fwd-magnitude (arr precision dest)
  (r2c-body arr #'(lambda (cdst nx &key dest precision)
                    (convert-complex-cvect-magnitudes-to-array 
                     cdst
                     (half-dim nx)
                     :dest dest
                     :precision precision))
            :precision precision
            :dest dest))

(defun slow-fwd-power (arr precision dest)
  (r2c-body arr #'(lambda (cdst nx &key dest precision)
                    (convert-complex-cvect-power-to-array 
                     cdst
                     (half-dim nx)
                     :dest dest
                     :precision precision))
            :precision precision
            :dest dest))

(defun slow-fwd-magnitude-db (arr precision dest)
  (r2c-body arr #'(lambda (cdst nx &key dest precision)
                    (convert-complex-cvect-magnitudes-db-to-array 
                     cdst
                     (half-dim nx)
                     :dest dest
                     :precision precision))
            :precision precision
            :dest dest))

(defun slow-fwd-phase (arr precision dest)
  (r2c-body arr #'(lambda (cdst nx &key dest precision)
                    (convert-complex-cvect-phases-to-array 
                     cdst
                     (half-dim nx)
                     :dest dest
                     :precision precision))
            :precision precision
            :dest dest))

(defun slow-fwd-phase-deg (arr precision dest)
  (r2c-body arr #'(lambda (cdst nx &key dest precision)
                    (convert-complex-cvect-phases-deg-to-array 
                     cdst
                     (half-dim nx)
                     :dest dest
                     :precision precision))
            :precision precision
            :dest dest))

;; -----------------------------------------------------------

(defun c2r (arr &key
                (precision :single-float)
                dest)
  (um:bind*
      ((nx (check-dimension arr))
       (:values (ctype ltype) (effective-ctype precision)))
    (fli:with-dynamic-foreign-objects ()
      (let ((cdst (fli:allocate-dynamic-foreign-object
                   :type   ctype
                   :nelems nx)))
        (fli:with-dynamic-foreign-objects ()
          (let ((csrc (fli:allocate-dynamic-foreign-object
                       :type   ctype
                       :nelems (* 2 nx))))
            (copy-array-to-complex-cvect arr csrc nx
                                         :precision ltype)
            (if (eq ctype :double)
                (z2dfft nx csrc cdst)
              (c2rfft nx csrc cdst))
            ))
        (convert-real-cvect-to-array cdst nx
                                     :dest dest
                                     :precision ltype))
      )))

#|
(defun z2z (arr dir
                &key
                (precision :single-float)
                dest)
  (um:bind*
      ((nx (check-dimension arr))
       (:values (ctype ltype) (effective-ctype precision)))
    (fli:with-dynamic-foreign-objects ()
      (let ((cdst (fli:allocate-dynamic-foreign-object
                   :type   ctype
                   :nelems (* 2 nx))))
        (copy-array-to-complex-cvect arr cdst nx 
                                     :precision ltype)
        (if (eq ctype :double)
            (z2zfft nx cdst cdst dir)
          (c2cfft nx cdst cdst dir))
        (convert-complex-cvect-to-array cdst nx
                                        :dest dest
                                        :precision ltype))
      )))
|#

;; --------------------------------------------------------------
;; fast routines for split-complex FFT requirements

(defvar *split-tmp* nil)

(defun get-split-temp-array (nx type)
  (let ((tmp *split-tmp*))
    (unless (and (fft-buffer-p tmp)
                 (>= (array-total-size (fft-buffer-r tmp)) nx)
                 (eql type (array-element-type (fft-buffer-r tmp))))
      (setf tmp (make-fft-buffer nx type)
            *split-tmp* tmp))
    (unless (= (fft-buffer-nx tmp) nx)
      (setf (fft-buffer-nx tmp) nx
            (fft-buffer-hr tmp) (make-array (half-dim nx)
                                            :element-type type
                                            :displaced-to (fft-buffer-r tmp)
                                            :displaced-index-offset (fft-buffer-roff tmp))))
    (values (fft-buffer-r tmp) (fft-buffer-roff tmp)
            (fft-buffer-i tmp) (fft-buffer-ioff tmp))
    ))

;; --------------------------------------------------------------

(defun d-copy-array-to-split-complex-cvect (arr dst-r roff dst-i ioff nx nxa)
  (declare (optimize (float 0) (safety 0) (speed 3)))
  (declare (type (array double-float (*)) dst-r dst-i))
  (declare (type fixnum nx nxa))
  (dotimes (ix nxa)
    (declare (fixnum ix))
    (let ((v (aref arr ix)))
      (setf (aref dst-r (+ ix roff)) (coerce (realpart v) 'double-float)
            (aref dst-i (+ ix ioff)) (coerce (imagpart v) 'double-float))))
  (when (> nx nxa)
    (fill dst-r 0d0 :start (+ nxa roff))
    (fill dst-i 0d0 :start (+ nxa ioff))
    ))

(defun d-convert-split-complex-cvect-to-array (src-r roff src-i ioff nx dest)
  (declare (optimize (float 0) (safety 0) (speed 3) (debug 0)))
  (declare (type (array double-float (*)) src-r src-i))
  (declare (type fixnum roff ioff nx))
  (let ((ans (or dest (make-complex-result-array nx 'double-float))))
    (dotimes (ix nx)
      (declare (fixnum ix))
      (setf (aref ans ix)
            (complex (aref src-r (+ ix roff))
                     (aref src-i (+ ix ioff)))))
    ans))

(defun unsafe-z2z (arr dir nx nxa dest)
  (declare (optimize (float 0) (safety 0) (speed 3)))
  (declare (fixnum nx nxa))
  (um:bind*
      ((:values (tmp-r roff tmp-i ioff) (get-split-temp-array nx 'double-float))
       (:declare (type (array double-float (*)) tmp-r tmp-i)))
    (d-copy-array-to-split-complex-cvect arr tmp-r roff tmp-i ioff nx nxa)
    (unsafe-z2zfft nx tmp-r roff tmp-i ioff dir)
    (d-convert-split-complex-cvect-to-array tmp-r roff tmp-i ioff nx dest)
    ))

;; --------------------------------------------------------------

(defun s-copy-array-to-split-complex-cvect (arr dst-r roff dst-i ioff nx nxa)
  (declare (optimize (float 0) (safety 0) (speed 3) (debug 0)))
  (declare (type (array single-float (*)) dst-r dst-i))
  (declare (type fixnum roff ioff nx nxa))
  (dotimes (ix nxa)
    (declare (fixnum ix))
    (let ((v (aref arr ix)))
      (setf (aref dst-r (+ ix roff)) (coerce (realpart v) 'single-float)
            (aref dst-i (+ ix ioff)) (coerce (imagpart v) 'single-float))))
  (when (> nx nxa)
    (fill dst-r 0e0 :start (+ nxa roff))
    (fill dst-i 0e0 :start (+ nxa ioff))
    ))

(defun s-convert-split-complex-cvect-to-array (src-r roff src-i ioff nx dest)
  (declare (optimize (float 0) (safety 0) (speed 3)))
  (declare (type (array single-float (*)) src-r src-i))
  (declare (type fixnum roff ioff nx))
  (let ((ans (or dest (make-complex-result-array nx 'single-float))))
    (dotimes (ix nx)
      (declare (fixnum ix))
      (setf (aref ans ix)
            (complex (aref src-r (+ ix roff))
                     (aref src-i (+ ix ioff)))))
    ans))

(defun unsafe-c2c (arr dir nx nxa dest)
  (declare (optimize (float 0) (safety 0) (speed 3)))
  (declare (fixnum nx nxa))
  (um:bind*
      ((:values (tmp-r roff tmp-i ioff) (get-split-temp-array nx 'single-float))
       (:declare (type (array single-float (*)) tmp-r tmp-i)))
    (s-copy-array-to-split-complex-cvect arr tmp-r roff tmp-i ioff nx nxa)
    (unsafe-c2cfft nx tmp-r roff tmp-i ioff dir)
    (s-convert-split-complex-cvect-to-array tmp-r roff tmp-i ioff nx dest)
    ))

;; --------------------------------------------------------------

(defun z2z (arr direction precision dest)
  (um:bind*
      ((nx  (check-dimension arr))
       (:declare (type fixnum nx))
       (nxa (effective-array-dimension arr))
       (:declare (type fixnum nxa))
       (:values (ctype ltype) (effective-ctype precision))
       (:declare (ignore ctype)))
    (labels ((setup-dst (dst)
             (unless (eq dst arr)
               (copy-fft-buffer-contents arr dst))
             dst))
      (cond ((and (fft-buffer-p arr)
                  (eql ltype 'single-float)
                  (eql (array-element-type (fft-buffer-r arr)) 'single-float))
             (let ((dst (setup-dst (or dest arr))))
               (unsafe-c2cfft nxa
                              (fft-buffer-r dst)
                              (fft-buffer-roff dst)
                              (fft-buffer-i dst)
                              (fft-buffer-ioff dst)
                              direction)
               dst))

          ((and (fft-buffer-p arr)
                (eql ltype 'double-float)
                (eql (array-element-type (fft-buffer-r arr)) 'double-float))
           (let ((dst (setup-dst (or dest arr))))
               (unsafe-z2zfft nxa
                              (fft-buffer-r dst)
                              (fft-buffer-roff dst)
                              (fft-buffer-i dst)
                              (fft-buffer-ioff dst)
                              direction)
               dst))

          ((eql ltype 'single-float)
           (unsafe-c2c arr direction nx nxa dest))

          ((eql ltype 'double-float)
           (unsafe-z2z arr direction nx nxa dest))

          (t (error "Invalid precision for FFT"))
          ))))

;; --------------------------------------------------------------

(defun fwd (arr &key
                (precision :single-float)
                dest)
  (z2z arr $fftw-forward precision dest))


(defun inv (arr &key
                (precision :single-float)
                dest)
  (z2z arr $fftw-inverse precision dest))

;; --------------------------------------------------------------

(defun fast-capable-p (arr)
  (fft-buffer-p arr))

(defun fast-fft-oper (arr after-fn dest)
  (um:bind*
      ((:values (rarr roff) (get-real arr))
       (:values (iarr ioff) (get-imag arr))
       (type    (array-element-type rarr))
       (dst     (or dest (fft-buffer-hr arr))))
    (set-imag arr 0)
    (funcall (if (eql type 'single-float)
                 #'unsafe-c2cfft
               #'unsafe-z2zfft)
             (fft-buffer-nx arr)
             rarr roff iarr ioff $fftw-forward)
    (dotimes (ix (length dst))
      (setf (aref dst ix) (funcall after-fn
                                   (aref rarr (+ ix roff))
                                   (aref iarr (+ ix ioff)))))
    dst))
    
;; --------------------------------------------------------------

(defun fwd-magnitude (arr &key
                          (precision :single-float)
                          dest)
  (cond ((fast-capable-p arr)
         (fast-fwd-magnitude arr dest))

        (t (slow-fwd-magnitude arr precision dest))
        ))

(defun fast-fwd-magnitude (arr dest)
  (fast-fft-oper arr (lambda (r i)
                       (sqrt (+ (* r r) (* i i))))
                 dest))


;; --------------------------------------------------------------

(defun fwd-power (arr &key
                      (precision :single-float)
                      dest)
  (cond ((fast-capable-p arr)
         (fast-fwd-power arr dest))

        (t (slow-fwd-power arr precision dest))
        ))

(defun fast-fwd-power (arr dest)
  (fast-fft-oper arr (lambda (r i)
                       (+ (* r r) (* i i)))
                 dest))
    
;; --------------------------------------------------------------

(defun fwd-magnitude-db (arr &key
			     (precision :single-float)
			     dest)
  (cond ((fast-capable-p arr)
         (fast-fwd-magnitude-db arr dest))

        (t (slow-fwd-magnitude-db arr precision dest))
        ))

(defun fast-fwd-magnitude-db (arr dest)
  (fast-fft-oper arr (lambda (r i)
                       (db10 (+ (* r r) (* i i))))
                 dest))

;; --------------------------------------------------------------

(defun fwd-phase (arr &key
                      (precision :single-float)
                      dest)
  (cond ((fast-capable-p arr)
         (fast-fwd-phase arr dest))

        (t (slow-fwd-phase arr precision dest))
        ))

(defun fast-fwd-phase (arr dest)
  (fast-fft-oper arr (lambda (r i)
                       (atan i r))
                 dest))

;; --------------------------------------------------------------

(defun fwd-phase-deg (arr &key
			  (precision :single-float)
			  dest)
  (cond ((fast-capable-p arr)
         (fast-fwd-phase-deg arr dest))

        (t (slow-fwd-phase-deg arr precision dest))
        ))

(defun fast-fwd-phase-deg (arr dest)
  (fast-fft-oper arr (lambda (r i)
                       (rtod (atan i r)))
                 dest))

;; --------------------------------------------------------------

#|
(defun doitf (n)
  (let ((x (vm:gnoise 4096))
        (y (make-array 4096 :element-type 'complex)))
    (time
     (dotimes (ix n)
      (fft:fwd x :dest y :precision :float)))))
(compile 'doitf)

(defun doitd (n)
  (let ((x (vm:gnoise 4096))
        (y (make-array 4096 :element-type 'complex)))
    (time
     (dotimes (ix n)
       (fft:fwd x :dest y :precision :double)))))
(compile 'doitd)

(defun doitfs (n)
  (let ((x (make-fft-buffer 4096 'single-float))
        (y (make-fft-buffer 4096 'single-float)))
    (um:move (vm:gnoise 4096) 0 (fft-buffer-r x) 0 4096)
    (um:move (vm:gnoise 4096) 0 (fft-buffer-i x) 0 4096)
    (time
     (dotimes (ix n)
      (fft:fwd x :dest y :precision :float)))))
(compile 'doitf)

(defun doitds (n)
  (let ((x (make-fft-buffer 4096 'double-float))
        (y (make-fft-buffer 4096 'double-float)))
    (map-into (fft-buffer-r x) (um:rcurry #'coerce 'double-float) (vm:gnoise 4096))
    (map-into (fft-buffer-i x) (um:rcurry #'coerce 'double-float) (vm:gnoise 4096))
    (time
     (dotimes (ix n)
       (fft:fwd x :dest y :precision :double)))))
(compile 'doitd)

|#

;; -- end of fft.lisp -- ;;

