#| DATE           : 28Dec1 
 | USER           : Dave 
 | PROCESSED FILE : D:\Projects\LispWorks\vmath\fft_intf.h
 |#

(in-package :FFT)

;;; Derived from file : "C:\\TEMP\\PCL188.h"
;;; Hand edited for module linkage DM/MCFA 12/01

(defvar *fftx-library*
  (fli:register-module :fftlib
                       :real-name
                       (translate-logical-pathname
                        #+:LISPWORKS-32BIT "PROJECTS:DYLIB;liblispfft.dylib"
                        #+:LISPWORKS-64BIT "PROJECTS:DYLIB64;liblispfft-64.dylib")))

;; ------------------------------------------------------------
(fli:define-foreign-function (getfftversionstring
                              "GetFFTVersionString"
                              :source)
    ((buf  (:reference-return (:ef-mb-string :limit 256)))
     (:constant 256 :int))
  :result-type :void
  :language :ansi-c
  :module *fftx-library*)

;; ------------------------------------------------------------

(defvar $fftw-forward 0)
(defvar $fftw-inverse 1)

;; ------------------------------------------------------------
;; Double-precision 1D FFTW routines

(fli:define-foreign-function (d2zfft "d2zfft" :source)
    ((nx :int)
     (src (:pointer :double))
     (dst (:pointer :double)))
  :result-type :int
  :language :ansi-c
  :module *fftx-library*)

(fli:define-foreign-function (z2dfft "z2dfft" :source)
    ((nx :int)
     (src (:pointer :double))
     (dst (:pointer :double)))
  :result-type :int
  :language :ansi-c
  :module *fftx-library*)

(fli:define-foreign-function (z2zfft "z2zfft" :source)
    ((nx :int)
     (src (:pointer :double))
     (dst (:pointer :double))
     (direction :int))
  :result-type :int
  :language :ansi-c
  :module *fftx-library*)

(fli:define-foreign-function (unsafe-z2zfft "unsafe_z2zfft" :source)
    ((nx :int)
     (src-r :lisp-simple-1d-array)
     (src-roff :int)
     (src-i :lisp-simple-1d-array)
     (src-ioff :int)
     (direction :int))
  :result-type :int
  :language :ansi-c
  :module *fftx-library*)

;; ------------------------------------------------------------
;; Double precision 2-D routines

(fli:define-foreign-function (d2zfft2d "d2zfft2d" :source)
    ((nx :int)
     (ny :int)
     (src (:pointer :double))
     (dst (:pointer :double)))
  :result-type :int
  :language :ansi-c
  :module *fftx-library*)

(fli:define-foreign-function (z2dfft2d "z2dfft2d" :source)
    ((nx :int)
     (ny :int)
     (src (:pointer :double))
     (dst (:pointer :double)))
  :result-type :int
  :language :ansi-c
  :module *fftx-library*)

(fli:define-foreign-function (z2zfft2d "z2zfft2d" :source)
    ((nx :int)
     (ny :int)
     (src (:pointer :double))
     (dst (:pointer :double))
     (direction :int))
  :result-type :int
  :language :ansi-c
  :module *fftx-library*)

(fli:define-foreign-function (unsafe-z2zfft2d "unsafe_z2zfft2d" :source)
    ((nx :int)
     (ny :int)
     (src-r :lisp-simple-1d-array)
     (src-roff :int)
     (src-i :lisp-simple-1d-array)
     (src-ioff :int)
     (direction :int))
  :result-type :int
  :language :ansi-c
  :module *fftx-library*)

;; ------------------------------------------------------------
;; Single-precision Altivec 1D routines
(fli:define-foreign-function (r2cfft "r2cfft" :source)
    ((nx :int)
     (src (:pointer :float))
     (dst (:pointer :float)))
  :result-type :int
  :language :ansi-c
  :module *fftx-library*)

(fli:define-foreign-function (c2rfft "c2rfft" :source)
    ((nx :int)
     (src (:pointer :float))
     (dst (:pointer :float)))
  :result-type :int
  :language :ansi-c
  :module *fftx-library*)

(fli:define-foreign-function (c2cfft "c2cfft" :source)
    ((nx :int)
     (src (:pointer :float))
     (dst (:pointer :float))
     (direction :int))
  :result-type :int
  :language :ansi-c
  :module *fftx-library*)

(fli:define-foreign-function (unsafe-c2cfft "unsafe_c2cfft" :source)
    ((nx :int)
     (src-r :lisp-simple-1d-array)
     (src-roff :int)
     (src-i :lisp-simple-1d-array)
     (src-ioff :int)
     (direction :int))
  :result-type :int
  :language :ansi-c
  :module *fftx-library*)

;; ------------------------------------------------------------
;; Single-precision Altivec 2D routines
(fli:define-foreign-function (r2cfft2d "r2cfft2d" :source)
    ((nx :int)
     (ny :int)
     (src (:pointer :float))
     (dst (:pointer :float)))
  :result-type :int
  :language :ansi-c
  :module *fftx-library*)

(fli:define-foreign-function (c2rfft2d "c2rfft2d" :source)
    ((nx :int)
     (ny :int)
     (src (:pointer :float))
     (dst (:pointer :float)))
  :result-type :int
  :language :ansi-c
  :module *fftx-library*)

(fli:define-foreign-function (c2cfft2d "c2cfft2d" :source)
    ((nx :int)
     (ny :int)
     (src (:pointer :float))
     (dst (:pointer :float))
     (direction :int))
  :result-type :int
  :language :ansi-c
  :module *fftx-library*)

(fli:define-foreign-function (unsafe-c2cfft2d "unsafe_c2cfft2d" :source)
    ((nx :int)
     (ny :int)
     (src-r :lisp-simple-1d-array)
     (src-roff :int)
     (src-i :lisp-simple-1d-array)
     (src-ioff :int)
     (direction :int))
  :result-type :int
  :language :ansi-c
  :module *fftx-library*)

;; -------------------------------------------------------

(fli:define-foreign-function (get-align16-offset "get_align16_offset" :source)
    ((buf  :lisp-simple-1d-array))
  :result-type :int
  :language :ansi-c
  :module *fftx-library*)

;; -------------------------------------------------------

(fli:define-foreign-function (siglab_sbFFT
                              "siglab_sbFFT"
                              :source)
    ((rsrc (:pointer :float))
     (cdst (:pointer :float))
     (nfft :int))
  :result-type :void
  :language :ansi-c
  :module *fftx-library*)

;; --------------------------------------------------------

(in-package :VMATH)

(defvar *siglab-library*
  (fli:register-module :siglab
                       :real-name
                       (translate-logical-pathname
                        #+:LISPWORKS-32BIT "PROJECTS:DYLIB;liblispsiglab.dylib"
                        #+:LISPWORKS-64BIT "PROJECTS:DYLIB64;liblispsiglab-64.dylib")))

(fli:define-foreign-function (disable-denormals
                              "siglab_disable_denormals"
                              :source)
    ()
  :result-type :int
  :language :ansi-c
  :module *siglab-library*)

(fli:define-foreign-function (restore-denormals
                              "siglab_restore_denormals"
                              :source)
    ((savemxcsr :int))
  :result-type :void
  :language :ansi-c
  :module *siglab-library*)

