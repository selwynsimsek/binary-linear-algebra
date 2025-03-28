;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;; Short form aliases for functions.

(in-package #:com.selwynsimsek.binary-linear-algebra)

(defmacro function-alias (function-name &rest aliases)
  (assert (symbolp function-name))
  `(progn ,@(loop for alias in aliases collect `(setf (fdefinition ',alias) (function ,function-name)))
          (export '(,function-name ,@aliases))))

(function-alias binary-identity-matrix bim)
(function-alias binary-zero-matrix b0m)
(function-alias binary-ones-matrix b1m)
(function-alias binary-scalar-vector bsv)
(function-alias binary-scalar-matrix bsm)
(function-alias binary-symplectic-matrix bsym)

(function-alias binary-zero-vector b0v)
(function-alias binary-ones-vector b1v)
(function-alias binary-one-hot-vector b1hv)
(function-alias binary-pad-vector bpv)
(function-alias binary-matrix-matrix-equal-p bmm=)
(function-alias binary-vector-vector-equal-p bvv=)
(function-alias mod2 m2)
(function-alias binary-vector-inner-product bv. bvv.)
(function-alias binary-vector-symplectic-inner-product bsv. bsvv.)
(function-alias popcount popcnt)
(function-alias binary-vector-outer-product bvo)

(function-alias random-bit rb)
(function-alias random-binary-vector rbv)
(function-alias random-binary-matrix rbm)
(function-alias random-binary-diagonal-matrix rbdm)
(function-alias random-binary-upper-triangular-unit-matrix rbutum)
(function-alias random-binary-lower-triangular-unit-matrix rbltum)
(function-alias random-binary-upper-triangular-matrix rbutm)
(function-alias random-binary-lower-triangular-matrix rbltm)
(function-alias random-binary-symplectic-matrix rbsm)
(function-alias random-binary-invertible-matrix rbim)

(function-alias binary-transpose bt)
(function-alias binary-matrix-diagonal-vector bmdv)

(function-alias binary-matrix-matrix-product bmm*)
(function-alias binary-matrix-vector-product bmv*)
(function-alias binary-scalar-vector-product bsv*)
(function-alias binary-matrix-matrix-sum bmm+)
(function-alias binary-vector-vector-sum bvv+)

(function-alias binary-vector-as-column-vector bvcv)
(function-alias binary-vector-as-row-vector bvrv)

(function-alias pluq-decomposition pluq pluqr)
(function-alias invert-binary-matrix ibm)

(function-alias binary-vector-p bv?)
(function-alias binary-matrix-p bm?)
(function-alias binary-matrix-permutation-p bmperm?)
(function-alias binary-matrix-square-p bmsq?)
(function-alias binary-matrix-lower-triangular-p bmlt?)
(function-alias binary-matrix-upper-triangular-p bmut?)
(function-alias binary-matrix-lower-triangular-unit-p bmltu?)
(function-alias binary-matrix-upper-triangular-unit-p bmutu?)
(function-alias binary-matrix-invertible-p bmi?)
(function-alias binary-matrix-symplectic-p bmsym?)

(function-alias null-space)
(function-alias row-echelon-null-space)
(function-alias upper-trapezoidal-null-space)
(function-alias trsm)
(function-alias upper-left-matrix ulm)
(function-alias solve-matrix-system)

(function-alias row-echelon-p re? rep)
(function-alias reduced-row-echelon-p rre? rrep)

(function-alias row-echelon-form ref)
(function-alias reduced-row-echelon-form rref)
