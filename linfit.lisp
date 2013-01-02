;; linfit.lisp -- Linear regression with Peter Stetson's iteration
;;
;; DM/RAL  11/10
;; ---------------------------------------------------------------

(defpackage :linfit
  (:use #:common-lisp)
  (:export
   #:regression))

(in-package :linfit)

(defun stetson-refinement (wts devs sigma &key (alpha 2) (beta 2))
  ;; take the resids and wts and produce new wts
  (map 'vector (lambda (dev wt)
                 (if (zerop dev)
                     wt
                   (/ wt (+ 1 (expt (/ (abs dev) (* alpha sigma)) beta)))))
       devs wts))

(defun regression (xs ys wts &key (alpha 2) (beta 2) (eps 1d-3))
  ;; return y-weighted mean and slope though data centroid along with stdev of fit.
  ;; iterate until the relative improvement in overall chisq
  ;; is less than eps. Be careful not to set eps below the sqrt machine precision.
  (let* ((xmn   (vm:mean xs))
         (xsmrm (matrix:- xs xmn))
         (wts   (if (numberp wts)
                    (make-array (length ys)
                                :initial-element wts)
                  wts)))
    (um:nlet-tail iter ((swts    wts)
                        (sigprev nil)
                        (niter   0))
      (let* ((ywmn   (/ (matrix:<*> ys swts) (matrix:sum swts)))
             (ysmrm  (matrix:- ys ywmn))
             (slope  (/ (matrix:<*> ysmrm xsmrm swts)
                        (matrix:<*> xsmrm xsmrm swts)))
             (devs   (matrix:- ysmrm (matrix:* slope xsmrm)))
             (wsigma (sqrt (/ (matrix:<*> swts devs devs)
                              (matrix:sum swts)))))
        (if (and (< niter 100)
                 (plusp wsigma)
                 (or (null sigprev)
                     (> (abs (- wsigma sigprev)) (* eps sigprev))))
            (iter (stetson-refinement wts devs wsigma
                                      :alpha alpha
                                      :beta  beta)
                  wsigma
                  (1+ niter))
          ;; else
          (progn
            (format t "~%linfit:regression niter = ~A" niter)
            (values xmn ywmn slope wsigma niter) ))))))
  
