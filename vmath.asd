
(asdf:defsystem "vmath"
  :description "vmath: a system for vectorized math"
  :version     "1.0"
  :author      "D.McClain <dbm@spectrodynamics.com>"
  :license     "Copyright (c) 2008 by SpectroDynamics, LLC. All rights reserved."
  :components  ((:file "packages")
                (:file "simple-vector-ops")
                (:file "vmath")
                (:file "matrix")
                (:file "lmfit")
                (:file "linfit")
                (:file "kalman")
                (:file "locate")
                (:file "interpolate")
                (:file "monotonic-spline")
                (:file "integrate")
                
                #+:MACOSX
                (:MODULE "Mac"
                 :COMPONENTS ((:file "mac_fft_intf")
                              (:file "mac-fft" :depends-on ("mac_fft_intf"))
                              (:file "mac-sfft")
                              (:file "mac-dfft")
                              (:file "mac-fft2d")
                              ))

                #+:win32
                (:MODULE "Win"
                 :components ((:file "win_fft_intf")
                              (:file "win-fft"  :depends-on ("win_fft_intf"))))
                
                #+:MACOSX (:file "dgesvd"))
  :serial t
  :depends-on  ("data-objects"
                #+:WIN32 "safearrays"
                #-:WIN32 "c-arrays"))

