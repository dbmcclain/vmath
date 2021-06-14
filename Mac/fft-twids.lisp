
(in-package :fft)

;; --------------------------------------------------------------

(defvar *twids-lock*  (mp:make-lock))
(defvar *ka-queue*    nil)

(defvar *stwids*
  (make-array 22
              :initial-element nil
              :adjustable t
              :weak t))

(defvar *dtwids*
  (make-array 22
              :initial-element nil
              :adjustable t
              :weak t))

(defstruct twids
  psetup)

(defstruct (stwids
            (:include twids)
            (:constructor stwids (psetup))
            ))
(defstruct (dtwids
            (:include twids)
            (:constructor dtwids (psetup))
            ))

(defun free-twids (obj)
  (cond
   ((stwids-p obj)
    (mp:with-lock (*twids-lock*)
      (cond ((assoc obj *ka-queue*)
             (hcl:flag-special-free-action obj))
            (t
             (nsubstitute nil obj *stwids*)
             (mp:funcall-async #'destroy-fft-setup (stwids-psetup obj)))
            )))

   ((dtwids-p obj)
    (mp:with-lock (*twids-lock*)
      (cond ((assoc obj *ka-queue*)
             (hcl:flag-special-free-action obj))
            (t
             (nsubstitute nil obj *dtwids*)
             (mp:funcall-async #'destroy-fft-setupD (dtwids-psetup obj)))
            )))
   ))
(hcl:add-special-free-action 'free-twids)

;; -------------------------------------------------------
;; Keep-Alive Queue... keep a reference on a TWIDS for 30 sec before
;; turning over to GC

(defvar *twids-timer* nil)

(defun scan-twids ()
  (mp:with-lock (*twids-lock*)
    (let ((expiration (- (get-universal-time) 30)))
      (if (um:deletef-if *ka-queue* (lambda (pair)
                                      (> expiration (cdr pair))))
          (mp:schedule-timer-relative *twids-timer* 30)
        (setf *twids-timer* nil))
      )))

(defun mark-twids (twids)
  (let ((pair (assoc twids *ka-queue*)))
    (if pair
        (setf (cdr pair) (get-universal-time))
      (um:aconsf *ka-queue* twids (get-universal-time)))
    (unless *twids-timer*
      (setf *twids-timer* (mp:make-timer #'mp:funcall-async #'scan-twids))
      (mp:schedule-timer-relative *twids-timer* 31))
    ))
                   
;; ---------------------------------------------------

(defun get-twids (nx cache create-fn)
  (let ((lg2nx (um:ceiling-log2 nx)))
    (assert (and (>= lg2nx 3)
                 (<= lg2nx 24)))
    (mp:with-lock (*twids-lock*)
      (let* ((slot  (- lg2nx 3))
             (twids (or (find-if #'identity cache :start slot)
                        (let ((twids (funcall create-fn lg2nx)))
                          (hcl:flag-special-free-action twids)
                          (setf (aref cache slot) twids))
                        )))
        (mark-twids twids)
        (twids-psetup twids)))))

(defun get-stwids (nx)
  (get-twids nx *stwids*
             (lambda (lg2nx)
               (stwids (create-fft-setup lg2nx)))
             ))

(defun get-dtwids (nx)
  (get-twids nx *dtwids*
             (lambda (lg2nx)
               (dtwids (create-fft-setupD lg2nx)))
             ))

;; ---------------------------------------------------------------------

(defun get-process-split-tmp ()
  (mp:process-private-property 'fft-split-tmp))

(defun set-process-split-tmp (tmp)
  (setf (mp:process-private-property 'fft-split-tmp) tmp))

(defsetf get-process-split-tmp  set-process-split-tmp)

(defun get-split-temp-array (nx type)
  (let ((tmp (get-process-split-tmp)))
    (unless (and (fft-buffer-p tmp)
                 (eql type (array-element-type (fft-buffer-r tmp)))
                 (>= (array-total-size (fft-buffer-r tmp)) nx))
      (setf tmp (make-fft-buffer nx type)
            (get-process-split-tmp) tmp))
    (unless (= (fft-buffer-nx tmp) nx)
      (setf (fft-buffer-nx tmp) nx
            (fft-buffer-hr tmp) (make-array (half-dim nx)
                                            :element-type type
                                            :displaced-to (fft-buffer-r tmp)
                                            :displaced-index-offset (fft-buffer-roff tmp))))
    ;; (chk-buf tmp)
    (values (fft-buffer-r tmp) (fft-buffer-roff tmp) (fft-buffer-pr tmp)
            (fft-buffer-i tmp) (fft-buffer-ioff tmp) (fft-buffer-pi tmp))
    ))

;; ---------------------------------------------------------------------
;; Per-process tmp buffer for FFT's. We never peek inside. Just need
;; an aligned split array of size equal or greater for the NX.

(defun get-process-tmp ()
  (mp:process-private-property 'fft-tmp))

(defun set-process-tmp (tmp)
  (setf (mp:process-private-property 'fft-tmp) tmp))

(defsetf get-process-tmp  set-process-tmp)

(defun get-temp-array (nx)
  (let ((tmp  (get-process-tmp)))
    (unless (and (fft-buffer-p tmp)
                 (>= (array-total-size (fft-buffer-r tmp)) nx))
      (setf tmp (make-fft-buffer nx 'single-float)
            (get-process-tmp) tmp))
    ;; (chk-buf tmp)
    (values (fft-buffer-pr tmp)
            (fft-buffer-pi tmp))
    ))

(defun get-stmp (nx)
  ;; buffer is max of 16384 bytes
  (get-temp-array (min nx 4096)))

(defun get-dtmp (nx)
  (get-stmp (* 2 nx)))

;; ---------------------------------------------------------------------

