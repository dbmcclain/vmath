;; mac-dfft.lisp -- interface to double-precision 1D FFT routines
;;
;; DM/MCFA  08/99
;; --------------------------------------------------

(in-package :DFFT)

;; ------------------------------------------------------------------

(defstruct (fft-buffer
            (:constructor make-fft-buf))
  nx
  r i
  roff ioff
  hr)

(defun half-dim (n)
  (1+ (truncate n 2)))

(defun nfloats (nb)
  (assert (zerop (logand nb 7))) ;; assure that our offsets are multiples of 8 bytes
  (truncate nb 8))

(defun make-fft-buffer (nx)
  (let* ((nxa   (max 8 (um:ceiling-pwr2 nx)))
         (rarr  (make-array (+ nxa 1) :element-type 'double-float :allocation :static))
         (iarr  (make-array (+ nxa 3) :element-type 'double-float :allocation :static))
         (roff  (nfloats (fft:get-align16-offset rarr)))
         (ioff  (nfloats (fft:get-align16-offset iarr)))
         (hrarr (make-array (half-dim nxa)
                            :element-type 'double-float
                            :displaced-to rarr
                            :displaced-index-offset roff)))
    (if (zerop (logand (- (+ (* 8 ioff) (sys:object-address iarr))
                          (+ (* 8 roff) (sys:object-address rarr)))
                       (1- 1024)))
        ;; offset ioff by another 2 to avoid the Pentium quirk when two buffer addresses differ
        ;; by multiple of 1024 bytes.
        (incf ioff 2)) ;; bump by another 16 bytes
                          
    (make-fft-buf
     :nx   nxa
     :r    rarr
     :roff roff
     :i    iarr
     :ioff ioff
     :hr   hrarr)
    ))

(defun get-real (fftbuf)
  (values (fft-buffer-r fftbuf) (fft-buffer-roff fftbuf) (fft-buffer-nx fftbuf)))

(defun get-imag (fftbuf)
  (values (fft-buffer-i fftbuf) (fft-buffer-ioff fftbuf) (fft-buffer-nx fftbuf)))

(defmethod set-real (fftbuf (arr vector))
  (replace (fft-buffer-r fftbuf) arr
           :start1 (fft-buffer-roff fftbuf)))

(defmethod set-real (fftbuf (val double-float))
  (fill (fft-buffer-r fftbuf) val :start (fft-buffer-roff fftbuf)))

(defmethod set-imag (fftbuf (arr vector))
  (replace (fft-buffer-i fftbuf) arr
           :start1 (fft-buffer-ioff fftbuf)))

(defmethod set-imag (fftbuf (val double-float))
  (fill (fft-buffer-i fftbuf) val :start (fft-buffer-ioff fftbuf)))

(defun copy-fft-buffer-contents (src dst)
  (replace (fft-buffer-r dst) (fft-buffer-r src)
           :start1 (fft-buffer-roff dst)
           :start2 (fft-buffer-roff src))
  (replace (fft-buffer-i dst) (fft-buffer-i src)
           :start1 (fft-buffer-ioff dst)
           :start2 (fft-buffer-ioff src)))

;; ------------------------------------------------------------------

(defun pwr (r i)
  (declare (type double-float r i))
  (+ (* r r) (* i i)))

(defun ampl (r i)
  (declare (type double-float r i))
  (sqrt (pwr r i)))

(defun db10 (r i)
  (declare (type double-float r i))
  (* 10d0 (log (pwr r i) 10d0)))

(defun rtod (x)
  (declare (type double-float x))
  (* #.(/ 180d0 pi) x))

(defun phs-deg (r i)
  (declare (type double-float r i))
  (rtod (phs r i)))

(defun phs (r i)
  (declare (type double-float r i))
  (atan i r))

;; -----------------------------------------------------------

(defun d2z (arr)
  ;; in-place routine
  (fill (fft-buffer-i arr) 0d0)
  (fft:unsafe-z2zfft (fft-buffer-nx arr)
                     (fft-buffer-r  arr)
                     (fft-buffer-roff arr)
                     (fft-buffer-i  arr)
                     (fft-buffer-ioff arr)
                     fft:$fftw-forward))

(defun z2d (arr)
  ;; in-place routine
  (fft:unsafe-z2zfft (fft-buffer-nx arr)
                     (fft-buffer-r  arr)
                     (fft-buffer-roff arr)
                     (fft-buffer-i  arr)
                     (fft-buffer-ioff arr)
                     fft:$fftw-inverse))

(defun z2z (arr dir)
  ;; in-place-routine
  (fft:unsafe-z2zfft (fft-buffer-nx arr)
                     (fft-buffer-r  arr)
                     (fft-buffer-roff arr)
                     (fft-buffer-i  arr)
                     (fft-buffer-ioff arr)
                     dir))

;; --------------------------------------------------------------

(defun fwd (arr &key dest)
  (cond (dest
         (copy-fft-buffer-contents arr dest)
         (fwd dest))

        (t (z2z arr fft:$fftw-forward))
        ))


(defun inv (arr &key dest)
  (cond (dest
         (copy-fft-buffer-contents arr dest)
         (inv dest))

        (t (z2z arr fft:$fftw-inverse))
        ))

;; --------------------------------------------------------------

(defun fast-real-fft-oper (arr after-fn dest)
  (um:bind*
      ((:values (rarr roff) (get-real arr))
       (:declare (type fixnum roff))
       (:declare (type (array double-float (*)) rarr))
       (:values (iarr ioff) (get-imag arr))
       (:declare (type fixnum ioff))
       (:declare (type (array double-float (*)) iarr))
       (dst     (or dest (fft-buffer-hr arr)))
       (:declare (type (array double-float (*)) dst)))
    (set-imag arr 0d0)
    (fft:unsafe-c2cfft (fft-buffer-nx arr)
                       rarr roff iarr ioff fft:$fftw-forward)
    (dotimes (ix (length dst))
      (declare (type fixnum ix))
      (setf (aref dst ix) (funcall after-fn
                                   (aref rarr (+ ix roff))
                                   (aref iarr (+ ix ioff)))))
    dst))
    
;; --------------------------------------------------------------

(defun fwd-magnitude (arr &key dest)
  (fast-real-fft-oper arr #'ampl dest))

;; --------------------------------------------------------------

(defun fwd-power (arr &key dest)
  (fast-real-fft-oper arr #'pwr dest))
    
;; --------------------------------------------------------------

(defun fwd-magnitude-db (arr &key dest)
  (fast-real-fft-oper arr #'db10 dest))

;; --------------------------------------------------------------

(defun fwd-phase (arr &key dest)
  (fast-real-fft-oper arr #'phs dest))
    
;; --------------------------------------------------------------

(defun fwd-phase-deg (arr &key dest)
  (fast-real-fft-oper arr #'phs-deg dest))


#|
(defun doitd (n)
  (let* ((nfft 1024)
         (x (make-fft-buffer nfft))
         (y (make-fft-buffer nfft)))
    (set-real x (map 'vector (um:rcurry #'coerce 'double-float) (vm:gnoise nfft)))
    (set-imag x (map 'vector (um:rcurry #'coerce 'double-float) (vm:gnoise nfft)))
    (time
     (dotimes (ix n)
      (fwd x :dest y)))))
(compile 'doitd)

|#

;; -- end of mac-dfft.lisp -- ;;

