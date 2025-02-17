(uiop:define-package #:com.selwynsimsek.binary-linear-algebra
  (:use #:cl)
  (:shadowing-import-from #:zr-utils #:define-function)
  (:shadowing-import-from #:metabang-bind #:bind))

(in-package #:com.selwynsimsek.binary-linear-algebra)


(defun binary-transpose (matrix)
  "Returns the bt of a binary matrix."
  (destructuring-bind (m n) (array-dimensions matrix)
    (let ((result (bm0 n m)))
      (loop for i from 0 below m
            do (loop for j from 0 below n
                     do (setf (aref result j i) (aref matrix i j))))
      result)))
