;; fft.lisp -- interface to FFTX 2-D FFT Routines
;;
;; DM/MCFA  08/99
;; --------------------------------------------------

(in-package :FFT2D)

;; ------------------------------------------------------------------

(defstruct (fft-buffer
            (:constructor make-fft-buf))
  ny nx
  r i
  roff ioff
  hr
  pr pi)

(defun half-dim (n)
  (1+ (truncate n 2)))

(defun nfloats (nb type-size)
  (assert (zerop (logand nb (1- type-size)))) ;; assure that our offsets are multiples of 4 bytes
  (truncate nb type-size))

(defun make-fft-buffer (ny nx type)
  (let* ((type-size (ecase type
                      (single-float 4)
                      (double-float 8)))
         (nya   (if (= 1 ny) 1 (max 8 (um:ceiling-pwr2 ny))))
         (nxa   (max 8 (um:ceiling-pwr2 nx)))
         (rarr  (make-array (+ (* nxa nya) 3) :element-type type :allocation :static))
         (roff  (nfloats (fft:get-align16-offset rarr) type-size))
         (iarr  (make-array (+ (* nxa nya) 7) :element-type type :allocation :static))
         (ioff  (nfloats (fft:get-align16-offset iarr) type-size))
         (hrarr (make-array (if (> ny 1)
                                (list (half-dim nya) (half-dim nxa))
                              (half-dim nxa))
                            :element-type type
                            :displaced-to rarr
                            :displaced-index-offset roff)))
    (if (zerop (logand (- (+ (* 4 ioff) (sys:object-address iarr))
                          (+ (* 4 roff) (sys:object-address rarr)))
                       (1- 1024)))
        ;; offset ioff by another 4 to avoid the Pentium quirk when two buffer addresses differ
        ;; by multiple of 1024 bytes.
        (incf ioff 4)) ;; bump by another 16 bytes
    (let ((ptr  (+ (* roff type-size) (fft:get-c-address rarr)))
          (pti  (+ (* ioff type-size) (fft:get-c-address iarr))))
      (assert (zerop (logand 15 ptr)))
      (assert (zerop (logand 15 pti)))
                          
      (make-fft-buf
       :ny   nya
       :nx   nxa
       :r    rarr
       :roff roff
       :i    iarr
       :ioff ioff
       :hr   hrarr
       :pr   ptr
       :pi   pti)
      )))

(defun get-real (fftbuf)
  (values (fft-buffer-r fftbuf) (fft-buffer-roff fftbuf) (fft-buffer-pr fftbuf)))

(defun get-imag (fftbuf)
  (values (fft-buffer-i fftbuf) (fft-buffer-ioff fftbuf) (fft-buffer-pi fftbuf)))

(defmethod set-real (fftbuf (arr vector))
  (replace (fft-buffer-r fftbuf) arr
           :start1 (fft-buffer-roff fftbuf)))

(defmethod set-real (fftbuf (val real))
  (multiple-value-bind (buf off) (get-real fftbuf)
    (fill buf (coerce val (array-element-type buf)) :start off)))

(defmethod set-imag (fftbuf (arr vector))
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

;; --------------------------------------------------------------------

(defmethod check-dimensions ((arr array))
  ;; return array dimensions (ny, nx) of arr rounded up to next
  ;; power of 2.
  (case (array-rank arr)
    (1 (values 1 (max 8 (um:ceiling-pwr2 (array-dimension arr 0)))))
    (2 (values (max 8 (um:ceiling-pwr2 (array-dimension arr 0)))
               (max 8 (um:ceiling-pwr2 (array-dimension arr 1)))))
    (otherwise (error "FFT only handles 1 and 2 dimensional arrays"))))


(defmethod check-dimensions ((arr fft-buffer))
  (values (fft-buffer-ny arr)
          (fft-buffer-nx arr)))
  
(defmethod effective-array-dimensions ((arr array))
  ;; return a pair (ny,nx) of array dimensions for arr.
  (ecase (array-rank arr)
    (1 (values 1 (array-dimension arr 0)))
    (2 (values (array-dimension arr 0)
               (array-dimension arr 1)))
    ))

(defmethod effective-array-dimensions ((arr fft-buffer))
  (values (fft-buffer-ny arr)
          (fft-buffer-nx arr)))

(defun effective-type (type)
  (cond ((and (consp type)
              (eql 'complex (car type)))
         (cadr type))
        (t
         type)))

(defmethod effective-array-element-type ((arr vector))
  (let ((type (array-element-type arr)))
    (cond ((eql type 't)
           (effective-type (type-of (aref arr 0))))
          (t
           (effective-type type))
          )))

(defmethod effective-array-element-type ((arr array))
  (let ((type (array-element-type arr)))
    (cond ((eql type 't)
           (effective-type (type-of (row-major-aref arr 0))))
          (t
           (effective-type (array-element-type arr)))
          )))

(defmethod effective-array-element-type ((arr fft-buffer))
  (effective-type (array-element-type (fft-buffer-r arr))))

(defun effective-ctype (precision)
  (ecase precision
    ((:float :single :single-float :altivec single-float) 
     (values :float 'single-float))
    ((:double :double-float :fftw double-float)
     (values :double 'double-float))))

(defun vec (arr &key (size (array-total-size arr)) (offset 0))
  (make-array size
              :displaced-to arr
              :displaced-index-offset offset
              :element-type (array-element-type arr)))

(defun copy-array-to-real-cvect (arr cdst ny nx &key precision)
  (multiple-value-bind (nya nxa) (effective-array-dimensions arr)
    (let ((pad (coerce 0 precision)))
      (dotimes (iy nya)
        (let ((va (vec arr :size nxa :offset (* iy nxa)))
              (coff (* iy nx)))
          (dotimes (ix nxa)
            (setf (fli:dereference cdst :index (+ coff ix))
                  (coerce (aref va ix) precision)))
          ;; pad remainder of row with zeros to get to power-of-two size
          (incf coff nxa)
          (dotimes (ix (- nx nxa))
            (setf (fli:dereference cdst :index (+ coff ix)) pad))
          ))
      ;; pad remaining rows with zeros to get to power-of-two size
      (let ((coff (* nx nya)))
        (dotimes (ix (* nx (- ny nya)))
          (setf (fli:dereference cdst :index (+ coff ix)) pad)))
      )))

(defun copy-array-to-complex-cvect (arr cdst ny nx &key precision)
  (multiple-value-bind (nya nxa) (effective-array-dimensions arr)
    (let ((pad (coerce 0 precision)))
      (dotimes (iy nya)
        (let ((va (vec arr :size nxa :offset (* iy nxa)))
              (coff (* 2 iy nx)))
          (dotimes (ix nxa)
            (let ((v (aref va ix)))
              (setf (fli:dereference cdst :index (+ coff ix ix))
                    (coerce (realpart v) precision)
                    (fli:dereference cdst :index (+ coff ix ix 1))
                    (coerce (imagpart v) precision))))
          ;; pad remainder of row with zeros to get to power-of-two size
          (incf coff (+ nxa nxa))
          (dotimes (ix (* 2 (- nx nxa)))
            (setf (fli:dereference cdst :index (+ coff ix)) pad))
          ))
      ;; pad remaining rows with zeros to get to power-of-two size
      (let ((coff (* 2 nx nya)))
        (dotimes (ix (* 2 nx (- ny nya)))
          (setf (fli:dereference cdst :index (+ coff ix)) pad)))
      )))

(defun reduced-dimensions (ny nx)
  (if (> ny 1)
      (list ny nx)
    (list nx)))

(defun make-result-array (ny nx precision)
  (make-array (reduced-dimensions ny nx)
              :element-type precision))

(defun make-complex-result-array (ny nx precision)
  (declare (ignore precision))
  (make-array (reduced-dimensions ny nx)
              :element-type 'complex))

(defun convert-real-cvect-to-array (csrc ny nx
                                         &key
                                         dest
                                         precision)
  (let* ((rslt  (or dest (make-result-array ny nx precision)))
         (vrslt (vec rslt))
         (arrsiz (* nx ny)))
    (do ((ix 0 (1+ ix)))
        ((>= ix arrsiz) rslt)
      (setf (aref vrslt ix) (fli:dereference csrc :index ix)))
    ))

(defun convert-complex-cvect-to-array (csrc ny nx
                                            &key
                                            dest
                                            precision)
  (let* ((rslt  (or dest (make-complex-result-array ny nx precision)))
         (vrslt (vec rslt))
         (arrsiz (* nx ny)))
    (do ((cix 0 (+ cix 2))
         (rix 0 (1+ rix)))
        ((>= rix arrsiz) rslt)
      (setf (aref vrslt rix)
            (complex (fli:dereference csrc :index cix)
                     (fli:dereference csrc :index (1+ cix))))
      )))

(defun db20 (x)
  (* 20.0 (log (abs x) 10.0)))

(defun db10 (x)
  (* 10.0 (log (abs x) 10.0)))

(defun dtor (x)
  (* #. (/ (atan 1 1) 45) x))

(defun rtod (x)
  (* #. (/ 45 (atan 1 1)) x))

(defun phase-deg (x)
  (rtod (phase x)))

(defun convert-complex-cvect-magnitudes-to-array (csrc ny nx
                                                       &key
                                                       dest
                                                       precision)
  (let* ((rslt  (or dest (make-result-array ny nx precision)))
         (vrslt (vec rslt))
         (arrsiz (* nx ny)))
    (do ((cix 0 (+ cix 2))
         (rix 0 (1+ rix)))
        ((>= rix arrsiz) rslt)
      (setf (aref vrslt rix)
            (coerce
             (abs (complex (fli:dereference csrc :index cix)
                           (fli:dereference csrc :index (1+ cix))))
             precision))
      )))

(defun sq (x)
  (* x x))

(defun convert-complex-cvect-power-to-array (csrc ny nx
                                                       &key
                                                       dest
                                                       precision)
  (let* ((rslt  (or dest (make-result-array ny nx precision)))
         (vrslt (vec rslt))
         (arrsiz (* nx ny)))
    (do ((cix 0 (+ cix 2))
         (rix 0 (1+ rix)))
        ((>= rix arrsiz) rslt)
      (setf (aref vrslt rix)
            (coerce
             (+ (sq (fli:dereference csrc :index cix))
                (sq (fli:dereference csrc :index (1+ cix))))
             precision))
      )))

(defun convert-complex-cvect-magnitudes-db-to-array (csrc ny nx
							  &key
							  dest
							  precision)
  (let* ((rslt  (or dest (make-result-array ny nx precision)))
         (vrslt (vec rslt))
         (arrsiz (* nx ny)))
    (do ((cix 0 (+ cix 2))
         (rix 0 (1+ rix)))
        ((>= rix arrsiz) rslt)
      (setf (aref vrslt rix)
            (coerce
             (db20 (complex (fli:dereference csrc :index cix)
                            (fli:dereference csrc :index (1+ cix))))
             precision))
      )))

(defun convert-complex-cvect-phases-to-array (csrc ny nx
                                                   &key
                                                   dest
                                                   precision)
  (let* ((rslt  (or dest (make-result-array ny nx precision)))
         (vrslt (vec rslt))
         (arrsiz (* nx ny)))
    (do ((cix 0 (+ cix 2))
         (rix 0 (1+ rix)))
        ((>= rix arrsiz) rslt)
      (setf (aref vrslt rix)
            (coerce
             (phase (complex (fli:dereference csrc :index cix)
                             (fli:dereference csrc :index (1+ cix))))
             precision))
      )))

(defun convert-complex-cvect-phases-deg-to-array (csrc ny nx
                                                   &key
                                                   dest
                                                   precision)
  (let* ((rslt  (or dest (make-result-array ny nx precision)))
         (vrslt (vec rslt))
         (arrsiz (* nx ny)))
    (do ((cix 0 (+ cix 2))
         (rix 0 (1+ rix)))
        ((>= rix arrsiz) rslt)
      (setf (aref vrslt rix)
            (coerce
             (phase-deg
              (complex (fli:dereference csrc :index cix)
                       (fli:dereference csrc :index (1+ cix))))
             precision))
      )))

#|
(defun tst (arr)
  (multiple-value-bind (ny nx) (check-dimensions arr)
    (fli:with-dynamic-foreign-objects ()
      (let ((cvec (fli:allocate-dynamic-foreign-object
                   :type :double
                   :nelems (* nx ny))))
        (copy-array-to-real-cvect arr cvec ny nx)
        (convert-real-cvect-to-array cvec ny nx))
      )))

(defun tstz (arr)
  (multiple-value-bind (ny nx) (check-dimensions arr)
    (fli:with-dynamic-foreign-objects ()
      (let ((cvec (fli:allocate-dynamic-foreign-object
                   :type :double
                   :nelems (* 2 nx ny))))
        (copy-array-to-complex-cvect arr cvec ny nx)
        (convert-complex-cvect-to-array cvec ny nx))
      )))
|#

(defun r2c-body (arr finish-fn &key
                     (precision (effective-array-element-type arr))
                     dest)
  (multiple-value-bind (ny nx) (check-dimensions arr)
    (let ((maxdim (max nx ny))
          (maxbuf (* nx ny)))
      (multiple-value-bind (ctype ltype) (effective-ctype precision)
        (fli:with-dynamic-foreign-objects ()
          (let ((cdst (fli:allocate-dynamic-foreign-object
                       :type   ctype
                       :nelems (* 2 nx ny))))
            (fli:with-dynamic-foreign-objects ()
              (let ((csrc (fli:allocate-dynamic-foreign-object
                           :type   ctype
                           :nelems (* nx ny))))
                (copy-array-to-real-cvect arr csrc ny nx
                                          :precision ltype)
                (if (eq ctype :double)
                    (fft:with-dtwids (twids maxdim)
                      (multiple-value-bind (prtmp pitmp) (fft:get-dtmp maxbuf)
                        (fft:d2zfft2d nx ny csrc cdst prtmp pitmp twids)))
                  (fft:with-stwids (twids maxdim)
                    (multiple-value-bind (prtmp pitmp) (fft:get-stmp maxbuf)
                      (fft:r2cfft2d nx ny csrc cdst prtmp pitmp twids))))
                ))
            (funcall finish-fn cdst ny nx
                     :dest dest
                     :precision ltype))
          ))
      )))

(defun r2c (arr &key
                (precision (effective-array-element-type arr))
                dest)
  (r2c-body arr #'convert-complex-cvect-to-array
            :precision precision
            :dest dest))

(defun slow-fwd-magnitude (arr precision dest)
  (r2c-body arr #'(lambda (cdst ny nx &key dest precision)
                    (convert-complex-cvect-magnitudes-to-array 
                     cdst
                     ny ;; (half-dim ny)
                     nx ;; (half-dim nx)
                     :dest dest
                     :precision precision))
            :precision precision
            :dest dest))

(defun slow-fwd-power (arr precision dest)
  (r2c-body arr #'(lambda (cdst ny nx &key dest precision)
                    (convert-complex-cvect-power-to-array 
                     cdst
                     ny ;; (half-dim ny)
                     nx ;; (half-dim nx)
                     :dest dest
                     :precision precision))
            :precision precision
            :dest dest))

(defun slow-fwd-magnitude-db (arr precision dest)
  (r2c-body arr #'(lambda (cdst ny nx &key dest precision)
                    (convert-complex-cvect-magnitudes-db-to-array 
                     cdst
                     ny ;; (half-dim ny)
                     nx ;; (half-dim nx)
                     :dest dest
                     :precision precision))
            :precision precision
            :dest dest))

(defun slow-fwd-phase (arr precision dest)
  (r2c-body arr #'(lambda (cdst ny nx &key dest precision)
                    (convert-complex-cvect-phases-to-array 
                     cdst
                     ny ;; (half-dim ny)
                     nx ;; (half-dim nx)
                     :dest dest
                     :precision precision))
            :precision precision
            :dest dest))

(defun slow-fwd-phase-deg (arr precision dest)
  (r2c-body arr #'(lambda (cdst ny nx &key dest precision)
                    (convert-complex-cvect-phases-deg-to-array 
                     cdst
                     ny ;; (half-dim ny)
                     nx ;; (half-dim nx)
                     :dest dest
                     :precision precision))
            :precision precision
            :dest dest))

;; -----------------------------------------------------------

(defun c2r (arr &key
                (precision (effective-array-element-type arr))
                dest)
  (multiple-value-bind (ny nx) (check-dimensions arr)
    (let ((maxdim (max nx ny))
          (maxbuf (* nx ny)))
      (multiple-value-bind (ctype ltype) (effective-ctype precision)
        (fli:with-dynamic-foreign-objects ()
          (let ((cdst (fli:allocate-dynamic-foreign-object
                       :type   ctype
                       :nelems (* nx ny))))
            (fli:with-dynamic-foreign-objects ()
              (let ((csrc (fli:allocate-dynamic-foreign-object
                           :type   ctype
                           :nelems (* 2 nx ny))))
                (copy-array-to-complex-cvect arr csrc ny nx
                                             :precision ltype)
                (if (eq ctype :double)
                    (fft:with-dtwids (twids maxdim)
                      (multiple-value-bind (prtmp pitmp) (fft:get-dtmp maxbuf)
                        (fft:z2dfft2d nx ny csrc cdst prtmp pitmp twids)))
                  (fft:with-stwids (twids maxdim)
                    (multiple-value-bind (prtmp pitmp) (fft:get-stmp maxbuf)
                      (fft:c2rfft2d nx ny csrc cdst prtmp pitmp twids))))
                ))
            (convert-real-cvect-to-array cdst ny nx
                                         :dest dest
                                         :precision ltype))
          ))
      )))

#|
(defun z2z (arr dir
                &key
                (precision :single-float)
                dest)
  (multiple-value-bind (ny nx) (check-dimensions arr)
    (multiple-value-bind (ctype ltype) (effective-ctype precision)
      (fli:with-dynamic-foreign-objects ()
        (let ((cdst (fli:allocate-dynamic-foreign-object
                     :type   ctype
                     :nelems (* 2 nx ny))))
          (copy-array-to-complex-cvect arr cdst ny nx 
                                       :precision ltype)
          (if (eq ctype :double)
              (fft:z2zfft2d nx ny cdst cdst dir)
            (fft:c2cfft2d nx ny cdst cdst dir))
          (convert-complex-cvect-to-array cdst ny nx
                                          :dest dest
                                          :precision ltype))
        ))
    ))
|#

;; --------------------------------------------------------------
;; fast routines for split-complex FFT requirements

(defun get-split-temp-array-2D (ny nx type)
  (let ((tmp (fft:get-process-split-tmp)))
    (unless (and (fft-buffer-p tmp)
                 (>= (array-total-size (fft-buffer-r tmp)) (* ny nx))
                 (eql type (array-element-type (fft-buffer-r tmp))))
      (setf tmp (make-fft-buffer ny nx type)
            (fft:get-process-split-tmp) tmp))
    (unless (and (= (fft-buffer-nx tmp) nx)
                 (= (fft-buffer-ny tmp) ny))
      (setf (fft-buffer-ny tmp) ny
            (fft-buffer-nx tmp) nx
            (fft-buffer-hr tmp) (make-array (* (half-dim ny) (half-dim nx))
                                            :element-type type
                                            :displaced-to (fft-buffer-r tmp)
                                            :displaced-index-offset (fft-buffer-roff tmp))))
    (values (fft-buffer-r tmp) (fft-buffer-roff tmp)
            (fft-buffer-i tmp) (fft-buffer-ioff tmp))
    ))

(defun d-copy-array-to-split-complex-cvect (arr dst-r roff dst-i ioff ny nx nya nxa)
  (declare (optimize (float 0) (safety 0) (speed 3)))
  (declare (type (array double-float (*)) dst-r dst-i))
  (declare (type fixnum ny nx nya nxa roff ioff))
  (let ((ny*nx   (the fixnum (* ny nx)))
        (nya*nxa (the fixnum (* nya nxa))))
    (declare (fixnum ny*nx nya*nxa))
    (do ((iy        0    (the fixnum (1+ iy)))
         (src-off   0    (the fixnum (+ src-off nxa)))
         (dstr-off  roff (the fixnum (+ dstr-off nx)))
         (dsti-off  ioff (the fixnum (+ dsti-off nx))))
        ((>= iy nya))
      (declare (fixnum iy src-off dst-off))
      (dotimes (ix nxa)
        (declare (fixnum ix))
        (let ((v (row-major-aref arr (the fixnum (+ src-off ix)))))
          (setf (aref dst-r (the fixnum (+ dstr-off ix)))
                (coerce (realpart v) 'double-float)
                (aref dst-i (the fixnum (+ dsti-off ix)))
                (coerce (imagpart v) 'double-float))))
      (when (> nx nxa)
        (fill dst-r 0d0 :start (the fixnum (+ dstr-off nxa)) :end (the fixnum (+ dstr-off nx)))
        (fill dst-i 0d0 :start (the fixnum (+ dsti-off nxa)) :end (the fixnum (+ dsti-off nx)))
        ))
    (when (> ny nya)
      (fill dst-r 0d0 :start (+ nya*nxa roff) :end (+ ny*nx roff))
      (fill dst-i 0d0 :start (+ nya*nxa ioff) :end (+ ny*nx ioff)))
    ))

(defun d-convert-split-complex-cvect-to-array (src-r roff src-i ioff ny nx dest)
  (declare (optimize (float 0) (safety 0) (speed 3) (debug 0)))
  (declare (type (array double-float (*)) src-r src-i))
  (declare (type fixnum ny nx roff ioff))
  (um:bind*
      ((ans   (or dest (make-complex-result-array ny nx 'double-float)))
       (ny*nx (the fixnum (* ny nx)))
       (:declare (fixnum (ny*nx))))
    (dotimes (ix ny*nx)
      (declare (fixnum ix))
      (setf (row-major-aref ans ix)
            (complex (aref src-r (+ ix roff))
                     (aref src-i (+ ix ioff)))))
    ans))

(defun unsafe-z2z (arr dir ny nx nya nxa dest)
  (declare (optimize (float 0) (safety 0) (speed 3)))
  (declare (fixnum ny nx nya nxa))
  (um:bind*
      ((:values (tmp-r roff ptr tmp-i ioff pti) (get-split-temp-array-2D ny nx 'double-float))
       (:declare (type (array double-float (*)) tmp-r tmp-i)))
    (d-copy-array-to-split-complex-cvect arr tmp-r roff tmp-i ioff ny nx nya nxa)
    (fft:with-dtwids (twids (max nx ny))
      (multiple-value-bind (prtmp pitmp) (fft:get-dtmp (* nx ny))
        (fft:unsafe-z2zfft2d nx ny ptr pti dir prtmp pitmp twids)))
    (d-convert-split-complex-cvect-to-array tmp-r roff tmp-i ioff ny nx dest)
    ))

;; --------------------------------------------------------------

(defun s-copy-array-to-split-complex-cvect (arr dst-r roff dst-i ioff ny nx nya nxa)
  (declare (optimize (float 0) (safety 0) (speed 3) (debug 0)))
  (declare (type (array single-float (*)) dst-r dst-i))
  (declare (type fixnum ny nx nya nxa))
  (let ((ny*nx   (* ny nx))
        (nya*nxa (* nya nxa)))
    (declare (fixnum ny*nx nya*nxa))
    (do ((iy       0    (the fixnum (1+ iy)))
         (src-off  0    (the fixnum (+ src-off nxa)))
         (dstr-off roff (the fixnum (+ dstr-off nx)))
         (dsti-off ioff (the fixnum (+ dsti-off nx))))
        ((>= iy nya))
      (declare (fixnum iy src-off dst-off))
      (dotimes (ix nxa)
        (declare (fixnum ix))
        (let ((v (row-major-aref arr (the fixnum (+ src-off ix)))))
          (setf (aref dst-r (the fixnum (+ dstr-off ix)))
                (coerce (realpart v) 'single-float)
                (aref dst-i (the fixnum (+ dsti-off ix)))
                (coerce (imagpart v) 'single-float))))
      (when (> nx nxa)
        (fill dst-r 0e0 :start (the fixnum (+ dstr-off nxa)) :end (the fixnum (+ dstr-off nx)))
        (fill dst-i 0e0 :start (the fixnum (+ dsti-off nxa)) :end (the fixnum (+ dsti-off nx)))
        ))
    (when (> ny nya)
      (fill dst-r 0e0 :start (+ nya*nxa roff) :end (* ny*nx roff))
      (fill dst-i 0e0 :start (+ nya*nxa ioff) :end (* ny*nx ioff)))
    ))

(defun s-convert-split-complex-cvect-to-array (src-r roff src-i ioff ny nx dest)
  (declare (optimize (float 0) (safety 0) (speed 3)))
  (declare (type (array single-float (*)) src-r src-i))
  (declare (type fixnum ny nx roff ioff))
  (um:bind*
      ((ans   (or dest (make-complex-result-array ny nx 'single-float)))
       (ny*nx (the fixnum (* ny nx)))
       (:declare (fixnum (ny*nx))))
    (dotimes (ix ny*nx)
      (declare (fixnum ix))
      (setf (row-major-aref ans ix)
            (complex (aref src-r (+ ix roff))
                     (aref src-i (+ ix ioff)))))
    ans))

(defun unsafe-c2c (arr dir ny nx nya nxa dest)
  (declare (optimize (float 0) (safety 0) (speed 3)))
  (declare (fixnum ny nx nya nxa))
  (um:bind*
      ((:values (tmp-r roff ptr tmp-i ioff pti) (get-split-temp-array-2D ny nx 'single-float))
       (:declare (type (array single-float (*)) tmp-r tmp-i)))
    (s-copy-array-to-split-complex-cvect arr tmp-r roff tmp-i ioff ny nx nya nxa)
    (fft:with-stwids (twids (max nx ny))
      (multiple-value-bind (prtmp pitmp) (fft:get-stmp (* nx ny))
        (fft:unsafe-c2cfft2d nx ny ptr pti dir prtmp pitmp twids)))
    (s-convert-split-complex-cvect-to-array tmp-r roff tmp-i ioff ny nx dest)
    ))

;; --------------------------------------------------------------

(defun z2z (arr direction precision dest)
  (um:bind*
      ((:values (ny nx)       (check-dimensions arr))
       (:values (nya nxa)     (effective-array-dimensions arr))
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
               (fft:with-stwids (twids (max nx ny))
                 (multiple-value-bind (prtmp pitmp) (fft:get-stmp (* nx ny))
                   (fft:unsafe-c2cfft2d nx ny
                                        (fft-buffer-pr dst)
                                        (fft-buffer-pi dst)
                                        direction
                                        prtmp pitmp twids)))
               dst))

          ((and (fft-buffer-p arr)
                (eql ltype 'double-float)
                (eql (array-element-type (fft-buffer-r arr)) 'double-float))
           (let ((dst (setup-dst (or dest arr))))
             (fft:with-dtwids (twids (max nx ny))
               (multiple-value-bind (prtmp pitmp) (fft:get-dtmp (* nx ny))
                 (fft:unsafe-z2zfft2d nx ny
                                      (fft-buffer-pr dst)
                                      (fft-buffer-pi dst)
                                      direction
                                      prtmp pitmp twids)))
             dst))
          
          ((eql ltype 'single-float)
           (unsafe-c2c arr direction ny nx nya nxa dest))

          ((eql ltype 'double-float)
           (unsafe-z2z arr direction ny nx nya nxa dest))

          (t (error "Invalid precision for FFT"))
          ))))

;; --------------------------------------------------------------

(defun fwd (arr &key
                (precision (effective-array-element-type arr))
                dest)
  (z2z arr fft:$fftw-forward precision dest))


(defun inv (arr &key
                (precision (effective-array-element-type arr))
                dest)
  (z2z arr fft:$fftw-inverse precision dest))

;; --------------------------------------------------------------

(defun fast-capable-p (arr)
  (and (fft-buffer-p arr)
       (= 1 (fft-buffer-ny arr))))

(defun fast-fft-oper (arr after-fn dest)
  (um:bind*
      ((:values (rarr roff ptr) (get-real arr))
       (:values (iarr ioff pti) (get-imag arr))
       (type    (array-element-type rarr))
       (dst     (or dest (fft-buffer-hr arr))))
    (let* ((nx (fft-buffer-nx arr))
           (ny (fft-buffer-ny arr))
           (maxdim (max nx ny))
           (maxbuf (* nx ny)))
      (set-imag arr 0)
      (if (eql type 'single-float)
          (fft:with-stwids (twids maxdim)
            (multiple-value-bind (prtmp pitmp) (fft:get-stmp maxbuf)
              (fft:unsafe-c2cfft2d
                 nx ny
                 ptr pti fft:$fftw-forward
                 prtmp pitmp twids)))
          (fft:with-dtwids (twids maxdim)
            (multiple-value-bind (prtmp pitmp) (fft:get-dtmp maxbuf)
              (fft:unsafe-z2zfft2d
                 nx ny
                 ptr pti fft:$fftw-forward
                 prtmp pitmp twids))))
      (if (= 2 (array-rank dst))
          (um:bind*
              (((ny nx) (array-dimensions dst))
               (nnx (fft-buffer-nx arr)))
            (dotimes (iy ny)
              (dotimes (ix nx)
                (setf (aref dst iy ix) (funcall after-fn
                                                (row-major-aref rarr (+ roff ix (* iy nnx)))
                                                (row-major-aref iarr (+ ioff ix (* iy nnx))))
                      ))))
        (dotimes (ix (length dst))
          (setf (aref dst ix) (funcall after-fn
                                       (aref rarr (+ ix roff))
                                       (aref iarr (+ ix ioff))))))
      dst)))
        
;; --------------------------------------------------------------

(defun fwd-magnitude (arr &key
                          (precision (effective-array-element-type arr))
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
                      (precision (effective-array-element-type arr))
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
			     (precision (effective-array-element-type arr))
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
                      (precision (effective-array-element-type arr))
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
			  (precision (effective-array-element-type arr))
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
      (fwd x :dest y :precision :float)))))
(compile 'doitf)

(defun doitd (n)
  (let ((x (vm:gnoise 4096))
        (y (make-array 4096 :element-type 'complex)))
    (time
     (dotimes (ix n)
       (fwd x :dest y :precision :double)))))
(compile 'doitd)

(defun doitfs (n)
  (let ((x (make-fft-buffer 1 4096 'single-float))
        (y (make-fft-buffer 1 4096 'single-float)))
    (set-real x (vm:gnoise 4096))
    (set-imag x (vm:gnoise 4096))
    (time
     (dotimes (ix n)
      (fwd x :dest y :precision :float)))))
(compile 'doitf)

(defun doitds (n)
  (let ((x (make-fft-buffer 1 4096 'double-float))
        (y (make-fft-buffer 1 4096 'double-float)))
    (set-real x (map 'vector (um:rcurry #'coerce 'double-float) (vm:gnoise 4096)))
    (set-imag x (map 'vector (um:rcurry #'coerce 'double-float) (vm:gnoise 4096)))
    (time
     (dotimes (ix n)
       (fwd x :dest y :precision :double)))))
(compile 'doitd)

|#

;; -- end of fft.lisp -- ;;

