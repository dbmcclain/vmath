;; lmfit.lisp -- Levenberg-Marquardt nonlinear fitting
;;
;; DM/RAL  11/10
;; -------------------------------------------------------

(defpackage :lmfit
  (:use #:common-lisp)
  (:export
   #:fit))

(in-package :lmfit)

;; ------------------------------------------------------------
;; Levenberg-Marquardt Nonlinear LSQ Fitting

(defun fit (v errfn derfn)
  (let* ((nparm  (length v))
         (errs   (funcall errfn v))
         (ndof   (- (length errs) nparm)))
    (labels ((chisq (v)
               (let ((errs (funcall errfn v)))
                 (/ (matrix:<*> errs errs) ndof)))
             
             (iter (chisq0 v lam)
               ;; (format t "~%iter: chisq = ~g, lam = ~g" chisq0 lam)
               (let* ((malpha  (matrix:make-matrix
                                :rows (matrix:as-vector
                                       (loop repeat nparm collect
                                             (make-array nparm)) )))
                      (vbeta   (make-array nparm)))
                 (let ((err0 (funcall errfn v))
                       (vder (funcall derfn v)))
                   (loop for ix from 0 below nparm do
                         (let ((der (aref vder ix)))
                           (setf (aref vbeta ix) (matrix:<*> der err0)
                                 (matrix:aref malpha ix ix) (matrix:<*> der der))
                           (loop for jx from 0 below ix do
                                 (let ((tot (matrix:<*> der (aref vder jx))))
                                   (setf (matrix:aref malpha ix jx) tot
                                         (matrix:aref malpha jx ix) tot)) ))) )
                 
                 (labels ((solve (lam)
                            ;; (format t "~%Solve: lam = ~g" lam)
                            (let ((malphax (matrix:copy malpha)))
                              (loop for ix from 0 below nparm do
                                    (setf (matrix:aref malphax ix ix)
                                          (if (zerop (matrix:aref malpha ix ix))
                                              1e-3
                                            (* (matrix:aref malpha ix ix)
                                               (+ 1d0 lam))) ))
                              (let* ((dv (matrix:cholsl malphax vbeta))
                                     (vx (vops:vsub v dv))
                                     (chisqx (chisq vx))
                                     (dchisq (abs (- chisq0 chisqx))))
                                ;; (format t "~%ChiSq = ~g" chisqx)
                                (if (or (< dchisq (* 1d-6 chisq0))
                                        (< chisqx 1d-12))
                                    (let ((cov (ignore-errors
                                                 (matrix:sqrt
                                                  (matrix:* chisqx
                                                            (matrix:abs
                                                             (matrix:inv malpha)))))))
                                      (unless cov
                                        (print "LMFIT: Singular COV"))
                                      ;; (format t "~%Finished...")
                                      (values chisqx vx cov))
                                  ;; else
                                  (if (> chisqx chisq0)
                                      (solve (* 10d0 lam))
                                    ;; else
                                    (iter chisqx vx (/ lam 10d0))) ))) ))
                   (solve lam)))))
      (iter (chisq v) v 0.001d0))))

