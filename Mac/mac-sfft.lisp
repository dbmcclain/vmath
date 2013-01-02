;; mac-sfft.lisp -- interface to single-precision 1D FFT routines
;;
;; DM/MCFA  08/99
;; --------------------------------------------------

(in-package :SFFT)

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
  (assert (zerop (logand nb 3))) ;; assure that our offsets are multiples of 4 bytes
  (truncate nb 4))

(defun make-fft-buffer (nx)
  (let* ((nxa   (max 8 (um:ceiling-pwr2 nx)))
         (rarr  (make-array (+ nxa 3) :element-type 'single-float :allocation :static))
         (iarr  (make-array (+ nxa 7) :element-type 'single-float :allocation :static))
         (roff  (nfloats (fft:get-align16-offset rarr)))
         (ioff  (nfloats (fft:get-align16-offset iarr)))
         (hrarr (make-array (half-dim nxa)
                            :element-type 'single-float
                            :displaced-to rarr
                            :displaced-index-offset roff)))
    (if (zerop (logand (- (+ (* 4 ioff) (sys:object-address iarr))
                          (+ (* 4 roff) (sys:object-address rarr)))
                       (1- 1024)))
        ;; offset ioff by another 4 to avoid the Pentium quirk when two buffer addresses differ
        ;; by multiple of 1024 bytes.
        (incf ioff 4)) ;; bump by another 16 bytes
                          
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

(defmethod set-real (fftbuf (val single-float))
  (fill (fft-buffer-r fftbuf) val :start (fft-buffer-roff fftbuf)))

(defmethod set-imag (fftbuf (arr vector))
  (replace (fft-buffer-i fftbuf) arr
           :start1 (fft-buffer-ioff fftbuf)))

(defmethod set-imag (fftbuf (val single-float))
  (fill (fft-buffer-i fftbuf) val :start (fft-buffer-ioff fftbuf)))

(defun copy-fft-buffer-contents (src dst)
  (replace (fft-buffer-r dst) (fft-buffer-r src)
           :start1 (fft-buffer-roff dst)
           :start2 (fft-buffer-roff src))
  (replace (fft-buffer-i dst) (fft-buffer-i src)
           :start1 (fft-buffer-ioff dst)
           :start2 (fft-buffer-ioff src)))

(defun pwr (r i)
  (declare (type single-float r i))
  (+ (* r r) (* i i)))

(defun ampl (r i)
  (declare (type single-float r i))
  (sqrt (pwr r i)))

(defun db10 (r i)
  (declare (type single-float r i))
  (* 10e0 (log (pwr r i) 10e0)))

(defun rtod (x)
  (declare (type single-float x))
  (float (* #. (/ 180e0 pi) x) 1e0))

(defun phs-deg (r i)
  (declare (type single-float r i))
  (rtod (phs r i)))

(defun phs (r i)
  (declare (type single-float r i))
  (atan i r))

;; -----------------------------------------------------------

(defun r2c (arr)
  ;; in-place routine
  (fill (fft-buffer-i arr) 0e0)
  (fft:unsafe-c2cfft (fft-buffer-nx arr)
                     (fft-buffer-r  arr)
                     (fft-buffer-roff arr)
                     (fft-buffer-i  arr)
                     (fft-buffer-ioff arr)
                     fft:$fftw-forward))

(defun c2r (arr)
  ;; in-place routine
  (fft:unsafe-c2cfft (fft-buffer-nx arr)
                     (fft-buffer-r  arr)
                     (fft-buffer-roff arr)
                     (fft-buffer-i  arr)
                     (fft-buffer-ioff arr)
                     fft:$fftw-inverse))

(defun c2c (arr dir)
  ;; in-place routine
  (fft:unsafe-c2cfft (fft-buffer-nx arr)
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

        (t (c2c arr fft:$fftw-forward))))


(defun inv (arr &key dest)
  (cond (dest
         (copy-fft-buffer-contents arr dest)
         (inv dest))

        (t (c2c arr fft:$fftw-inverse))))

;; --------------------------------------------------------------

(defun fast-real-fft-oper (arr after-fn dest)
  (um:bind*
      ((:values (rarr roff) (get-real arr))
       (:declare (type fixnum roff))
       (:declare (type (array single-float (*)) rarr))
       (:values (iarr ioff) (get-imag arr))
       (:declare (type fixnum ioff))
       (:declare (type (array single-float (*)) iarr))
       (dst     (or dest (fft-buffer-hr arr)))
       (:declare (type (array single-float (*)) dst)))
    ;; (set-imag arr 0e0)
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
(defun doits (n)
  (let* ((nfft 1024)
         (x (make-fft-buffer nfft))
         (y (make-fft-buffer nfft)))
    (replace (fft-buffer-r x) (vm:gnoise nfft))
    (replace (fft-buffer-i x) (vm:gnoise nfft))
    (time
     (dotimes (ix n)
       (fwd x :dest y)))))
(compile 'doits)

|#

;; -- end of mac-sfft.lisp -- ;;

